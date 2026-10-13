
from pathlib import Path
import copy
import json
import random

import numpy as np
import pandas as pd
import torch
from torch import nn
import torch.nn.functional as F

from torch_geometric.nn import RGCNConv
from sklearn.metrics import roc_auc_score, average_precision_score


# =========================================================
# 1. Configuration and reproducibility
# =========================================================

SEED = 42
HIDDEN_DIM = 64
DROPOUT = 0.2
LEARNING_RATE = 0.005
WEIGHT_DECAY = 0.001
MAX_EPOCHS = 100
PATIENCE = 12

random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.set_num_threads(4)

DEVICE = torch.device("cpu")

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "processed"
MODEL_DIR = ROOT / "data" / "model"
SPLIT_DIR = ROOT / "data" / "splits"
OUTPUT_DIR = ROOT / "reports" / "model_evaluation"
CHECKPOINT_DIR = ROOT / "models"

for folder in [OUTPUT_DIR, CHECKPOINT_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

print("Device:", DEVICE)
print("Random seed:", SEED)


# =========================================================
# 2. Load nodes and leakage-safe graph
# =========================================================

nodes = pd.read_parquet(DATA_DIR / "nodes.parquet")
graph = pd.read_parquet(
    MODEL_DIR / "training_graph_edges.parquet"
)

nodes["node_id"] = nodes["node_id"].astype(str)
graph["source_id"] = graph["source_id"].astype(str)
graph["target_id"] = graph["target_id"].astype(str)

node_to_idx = {
    node_id: i
    for i, node_id in enumerate(nodes["node_id"])
}

node_types = nodes["node_type"].tolist()
node_type_names = ["drug", "disease", "protein", "pathway"]
node_type_to_idx = {
    name: i for i, name in enumerate(node_type_names)
}

unknown_types = set(node_types) - set(node_type_to_idx)

if unknown_types:
    raise ValueError(f"Unexpected node types: {unknown_types}")

node_type_tensor = torch.tensor(
    [node_type_to_idx[t] for t in node_types],
    dtype=torch.long,
    device=DEVICE,
)

print("Nodes:", len(nodes))
print("Node types:", nodes["node_type"].value_counts().to_dict())


# =========================================================
# 3. Build relation-aware message-passing graph
# =========================================================

# Use biological relations only. Treatment labels and
# contraindications are not message-passing edges.
SUPPORT_RELATIONS = [
    "targets",
    "participates_in",
    "associated_with",
]

relation_to_idx = {
    relation: i
    for i, relation in enumerate(SUPPORT_RELATIONS)
}

num_base_relations = len(SUPPORT_RELATIONS)
num_relations = num_base_relations * 2

src_indices = []
dst_indices = []
relation_indices = []

used_counts = {r: 0 for r in SUPPORT_RELATIONS}

for row in graph.itertuples(index=False):
    relation = row.relation

    if relation not in relation_to_idx:
        continue

    source = row.source_id
    target = row.target_id

    if source not in node_to_idx or target not in node_to_idx:
        raise ValueError(
            f"Unknown endpoint: {source} -> {target}"
        )

    src = node_to_idx[source]
    dst = node_to_idx[target]
    rel = relation_to_idx[relation]

    # Original direction.
    src_indices.append(src)
    dst_indices.append(dst)
    relation_indices.append(rel)

    # Reverse direction, with its own relation ID.
    src_indices.append(dst)
    dst_indices.append(src)
    relation_indices.append(rel + num_base_relations)

    used_counts[relation] += 1

if not src_indices:
    raise ValueError("No supporting biological edges were found.")

edge_index = torch.tensor(
    [src_indices, dst_indices],
    dtype=torch.long,
    device=DEVICE,
)

edge_type = torch.tensor(
    relation_indices,
    dtype=torch.long,
    device=DEVICE,
)

print("\nSupporting graph")
print("Edges by relation:", used_counts)
print("Directed message-passing edges:", edge_index.shape[1])
print("R-GCN relation types, including inverse:", num_relations)


# =========================================================
# 4. Load labelled drug-disease pairs
# =========================================================

def load_pairs(split_name):
    path = SPLIT_DIR / f"{split_name}_pairs.csv"

    if not path.exists():
        raise FileNotFoundError(f"Missing split: {path}")

    df = pd.read_csv(
        path,
        dtype={"drug_id": str, "disease_id": str},
    )

    required = {"drug_id", "disease_id", "label"}
    if not required.issubset(df.columns):
        raise ValueError(f"Unexpected columns in {path}")

    if df[["drug_id", "disease_id"]].isna().any().any():
        raise ValueError(f"Missing pair IDs in {path}")

    if df.duplicated(["drug_id", "disease_id"]).any():
        raise ValueError(f"Duplicate pairs in {path}")

    unknown_ids = (
        set(df["drug_id"]) | set(df["disease_id"])
    ) - set(node_to_idx)

    if unknown_ids:
        raise ValueError(
            f"Pair IDs not found in nodes: {list(unknown_ids)[:5]}"
        )

    if any(
        nodes.iloc[node_to_idx[d]]["node_type"] != "drug"
        for d in df["drug_id"]
    ):
        raise ValueError(f"Non-drug IDs found in {path}")

    if any(
        nodes.iloc[node_to_idx[d]]["node_type"] != "disease"
        for d in df["disease_id"]
    ):
        raise ValueError(f"Non-disease IDs found in {path}")

    return df


splits = {
    name: load_pairs(name)
    for name in ["train", "validation", "test"]
}


def make_pair_tensors(df):
    drugs = torch.tensor(
        [node_to_idx[x] for x in df["drug_id"]],
        dtype=torch.long,
        device=DEVICE,
    )

    diseases = torch.tensor(
        [node_to_idx[x] for x in df["disease_id"]],
        dtype=torch.long,
        device=DEVICE,
    )

    labels = torch.tensor(
        df["label"].to_numpy(dtype=np.float32),
        dtype=torch.float32,
        device=DEVICE,
    )

    return drugs, diseases, labels


pair_tensors = {
    name: make_pair_tensors(df)
    for name, df in splits.items()
}

for name, df in splits.items():
    print(
        f"{name.capitalize()} pairs: {len(df)}, "
        f"positives: {(df['label'] == 1).sum()}, "
        f"negatives: {(df['label'] == 0).sum()}"
    )


# =========================================================
# 5. Define the R-GCN and prediction decoder
# =========================================================

class ReMedyRGCN(nn.Module):
    def __init__(
        self,
        num_nodes,
        num_node_types,
        hidden_dim,
        num_relations,
        dropout,
    ):
        super().__init__()

        self.dropout = dropout

        # Learn an embedding for each entity and a shared
        # embedding for each node type.
        self.entity_embedding = nn.Embedding(
            num_nodes, hidden_dim
        )
        self.type_embedding = nn.Embedding(
            num_node_types, hidden_dim
        )

        nn.init.xavier_uniform_(
            self.entity_embedding.weight
        )
        nn.init.xavier_uniform_(
            self.type_embedding.weight
        )

        self.conv1 = RGCNConv(
            hidden_dim,
            hidden_dim,
            num_relations=num_relations,
            num_bases=3,
        )

        self.conv2 = RGCNConv(
            hidden_dim,
            hidden_dim,
            num_relations=num_relations,
            num_bases=3,
        )

        # The decoder receives drug, disease, product and
        # absolute-difference representations.
        self.decoder = nn.Sequential(
            nn.Linear(hidden_dim * 4, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def encode(self, node_type_ids, edge_index, edge_type):
        x = (
            self.entity_embedding.weight
            + self.type_embedding(node_type_ids)
        )

        x = self.conv1(x, edge_index, edge_type)
        x = F.relu(x)
        x = F.dropout(
            x, p=self.dropout, training=self.training
        )

        x = self.conv2(x, edge_index, edge_type)
        return x

    def decode(self, embeddings, drug_idx, disease_idx):
        drug_vectors = embeddings[drug_idx]
        disease_vectors = embeddings[disease_idx]

        pair_features = torch.cat(
            [
                drug_vectors,
                disease_vectors,
                drug_vectors * disease_vectors,
                torch.abs(drug_vectors - disease_vectors),
            ],
            dim=1,
        )

        return self.decoder(pair_features).squeeze(-1)


model = ReMedyRGCN(
    num_nodes=len(nodes),
    num_node_types=len(node_type_names),
    hidden_dim=HIDDEN_DIM,
    num_relations=num_relations,
    dropout=DROPOUT,
).to(DEVICE)

optimizer = torch.optim.Adam(
    model.parameters(),
    lr=LEARNING_RATE,
    weight_decay=WEIGHT_DECAY,
)

criterion = nn.BCEWithLogitsLoss()

train_drugs, train_diseases, train_labels = pair_tensors["train"]
val_drugs, val_diseases, val_labels = pair_tensors["validation"]
test_drugs, test_diseases, test_labels = pair_tensors["test"]


# =========================================================
# 6. Evaluation utilities
# =========================================================

def calculate_metrics(labels, scores):
    labels_np = labels.detach().cpu().numpy().astype(int)
    scores_np = scores.detach().cpu().numpy()

    auc = roc_auc_score(labels_np, scores_np)
    ap = average_precision_score(labels_np, scores_np)

    order = np.argsort(-scores_np)
    k = min(10, len(labels_np))
    top_k = labels_np[order[:k]]

    precision_at_k = float(top_k.mean())
    positives = int(labels_np.sum())
    recall_at_k = (
        float(top_k.sum() / positives)
        if positives > 0 else 0.0
    )

    return {
        "roc_auc": float(auc),
        "auprc_average_precision": float(ap),
        "precision_at_10": precision_at_k,
        "recall_at_10": recall_at_k,
    }


def predict_pairs(df):
    drugs, diseases, labels = make_pair_tensors(df)

    model.eval()
    with torch.no_grad():
        embeddings = model.encode(
            node_type_tensor, edge_index, edge_type
        )
        logits = model.decode(embeddings, drugs, diseases)
        probabilities = torch.sigmoid(logits)

    return labels, logits, probabilities


# =========================================================
# 7. Train using validation AUPRC for early stopping
# =========================================================

best_validation_ap = -1.0
best_epoch = 0
best_state = None
epochs_without_improvement = 0
training_log = []

print("\nStarting R-GCN training...")
print("Maximum epochs:", MAX_EPOCHS)
print("Early stopping patience:", PATIENCE)

for epoch in range(1, MAX_EPOCHS + 1):
    model.train()
    optimizer.zero_grad()

    embeddings = model.encode(
        node_type_tensor, edge_index, edge_type
    )

    train_logits = model.decode(
        embeddings, train_drugs, train_diseases
    )

    loss = criterion(train_logits, train_labels)
    loss.backward()

    torch.nn.utils.clip_grad_norm_(
        model.parameters(), max_norm=2.0
    )

    optimizer.step()

    # Evaluate on validation pairs without changing parameters.
    model.eval()

    with torch.no_grad():
        validation_embeddings = model.encode(
            node_type_tensor, edge_index, edge_type
        )

        validation_logits = model.decode(
            validation_embeddings, val_drugs, val_diseases
        )

        validation_probabilities = torch.sigmoid(
            validation_logits
        )

        val_metrics = calculate_metrics(
            val_labels, validation_probabilities
        )

    val_ap = val_metrics["auprc_average_precision"]

    training_log.append({
        "epoch": epoch,
        "training_loss": float(loss.item()),
        "validation_roc_auc": val_metrics["roc_auc"],
        "validation_auprc": val_ap,
    })

    if val_ap > best_validation_ap + 1e-5:
        best_validation_ap = val_ap
        best_epoch = epoch
        best_state = copy.deepcopy(model.state_dict())
        epochs_without_improvement = 0
    else:
        epochs_without_improvement += 1

    if epoch == 1 or epoch % 5 == 0:
        print(
            f"Epoch {epoch:03d} | "
            f"Loss {loss.item():.4f} | "
            f"Validation AP {val_ap:.4f} | "
            f"Best AP {best_validation_ap:.4f}"
        )

    if epochs_without_improvement >= PATIENCE:
        print(f"\nEarly stopping at epoch {epoch}.")
        break

if best_state is None:
    raise RuntimeError("Training failed to produce a model.")

# Restore the checkpoint selected by validation performance.
model.load_state_dict(best_state)

pd.DataFrame(training_log).to_csv(
    OUTPUT_DIR / "rgcn_training_log.csv",
    index=False,
)


# =========================================================
# 8. Final validation and held-out test evaluation
# =========================================================

val_labels_final, _, val_probabilities = predict_pairs(
    splits["validation"]
)
test_labels_final, _, test_probabilities = predict_pairs(
    splits["test"]
)

val_metrics = calculate_metrics(
    val_labels_final, val_probabilities
)
test_metrics = calculate_metrics(
    test_labels_final, test_probabilities
)

metrics_rows = [
    {
        "model": "R-GCN",
        "split": "validation",
        **val_metrics,
    },
    {
        "model": "R-GCN",
        "split": "test",
        **test_metrics,
    },
]

pd.DataFrame(metrics_rows).to_csv(
    OUTPUT_DIR / "rgcn_metrics.csv",
    index=False,
)

test_predictions = splits["test"].copy()
test_predictions["rgcn_probability"] = (
    test_probabilities.cpu().numpy()
)
test_predictions = test_predictions.sort_values(
    "rgcn_probability", ascending=False
)

test_predictions.to_csv(
    OUTPUT_DIR / "rgcn_test_predictions.csv",
    index=False,
)


# =========================================================
# 9. Save the model checkpoint and metadata
# =========================================================

checkpoint_path = CHECKPOINT_DIR / "rgcn_baseline.pt"

torch.save(
    {
        "model_state_dict": model.state_dict(),
        "node_to_idx": node_to_idx,
        "node_type_to_idx": node_type_to_idx,
        "relation_to_idx": relation_to_idx,
        "num_relations": num_relations,
        "hidden_dim": HIDDEN_DIM,
        "dropout": DROPOUT,
        "dataset_version": "ReMedy-KG v1.0",
        "best_epoch": best_epoch,
        "best_validation_auprc": best_validation_ap,
        "seed": SEED,
    },
    checkpoint_path,
)

metadata = {
    "model": "R-GCN",
    "dataset_version": "ReMedy-KG v1.0",
    "seed": SEED,
    "device": str(DEVICE),
    "hidden_dimension": HIDDEN_DIM,
    "dropout": DROPOUT,
    "learning_rate": LEARNING_RATE,
    "weight_decay": WEIGHT_DECAY,
    "best_epoch": best_epoch,
    "best_validation_auprc": best_validation_ap,
    "support_relations": SUPPORT_RELATIONS,
    "inverse_relations_added": True,
    "treatment_relations_used_in_message_passing": False,
    "contraindications_used_in_message_passing": False,
    "validation_metrics": val_metrics,
    "test_metrics": test_metrics,
}

with open(
    OUTPUT_DIR / "rgcn_run_metadata.json",
    "w",
    encoding="utf-8",
) as file:
    json.dump(metadata, file, indent=2)

# Create a comparison table with the previous baseline results.
baseline_path = OUTPUT_DIR / "baseline_metrics.csv"

if baseline_path.exists():
    baseline = pd.read_csv(baseline_path)
    baseline_test = baseline.loc[
        baseline["split"] == "test"
    ].copy()

    rgcn_test = pd.DataFrame([
        {
            "model": "R-GCN",
            "split": "test",
            **test_metrics,
        }
    ])

    comparison = pd.concat(
        [baseline_test, rgcn_test],
        ignore_index=True,
        sort=False,
    )

    comparison.to_csv(
        OUTPUT_DIR / "model_comparison.csv",
        index=False,
    )


# =========================================================
# 10. Print results
# =========================================================

print("\n" + "=" * 60)
print("REMEDY R-GCN RESULTS")
print("=" * 60)

print("Best epoch:", best_epoch)
print("Best validation AUPRC:", round(best_validation_ap, 4))

print("\nValidation metrics:")
for metric, value in val_metrics.items():
    print(f"{metric}: {value:.4f}")

print("\nHeld-out test metrics:")
for metric, value in test_metrics.items():
    print(f"{metric}: {value:.4f}")

print("\nSaved files:")
print("models/rgcn_baseline.pt")
print("reports/model_evaluation/rgcn_training_log.csv")
print("reports/model_evaluation/rgcn_metrics.csv")
print("reports/model_evaluation/rgcn_test_predictions.csv")
print("reports/model_evaluation/rgcn_run_metadata.json")

if (OUTPUT_DIR / "model_comparison.csv").exists():
    print("reports/model_evaluation/model_comparison.csv")

print("\nR-GCN training completed.")
