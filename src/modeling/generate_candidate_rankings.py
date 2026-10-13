
from pathlib import Path
from collections import defaultdict
import json

import numpy as np
import pandas as pd
import torch
from torch import nn
import torch.nn.functional as F
from torch_geometric.nn import RGCNConv


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "processed"
MODEL_DIR = ROOT / "data" / "model"
SPLIT_DIR = ROOT / "data" / "splits"
OUTPUT_DIR = ROOT / "reports" / "model_evaluation"
CHECKPOINT_PATH = ROOT / "models" / "rgcn_baseline.pt"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
DEVICE = torch.device("cpu")
BATCH_SIZE = 1024

if not CHECKPOINT_PATH.exists():
    raise FileNotFoundError(
        f"Missing checkpoint: {CHECKPOINT_PATH}"
    )


# ---------------------------------------------------------
# 1. Re-create the exact model architecture used in training
# ---------------------------------------------------------

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

        self.entity_embedding = nn.Embedding(
            num_nodes, hidden_dim
        )
        self.type_embedding = nn.Embedding(
            num_node_types, hidden_dim
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

        return self.conv2(x, edge_index, edge_type)

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


# ---------------------------------------------------------
# 2. Load the trained checkpoint and frozen KG
# ---------------------------------------------------------

checkpoint = torch.load(
    CHECKPOINT_PATH,
    map_location=DEVICE,
    weights_only=False,
)

nodes = pd.read_parquet(DATA_DIR / "nodes.parquet")
edges = pd.read_parquet(DATA_DIR / "edges.parquet")
graph = pd.read_parquet(
    MODEL_DIR / "training_graph_edges.parquet"
)

nodes["node_id"] = nodes["node_id"].astype(str)
edges["source_id"] = edges["source_id"].astype(str)
edges["target_id"] = edges["target_id"].astype(str)
graph["source_id"] = graph["source_id"].astype(str)
graph["target_id"] = graph["target_id"].astype(str)

node_to_idx = checkpoint["node_to_idx"]
node_type_to_idx = checkpoint["node_type_to_idx"]
relation_to_idx = checkpoint["relation_to_idx"]

if len(node_to_idx) != len(nodes):
    raise ValueError("Checkpoint and current node table disagree.")

if set(node_to_idx) != set(nodes["node_id"]):
    raise ValueError("Node IDs differ from the trained checkpoint.")

node_type_by_id = dict(
    zip(nodes["node_id"], nodes["node_type"])
)

ordered_node_ids = [
    node_id
    for node_id, index in sorted(
        node_to_idx.items(), key=lambda item: item[1]
    )
]

node_type_tensor = torch.tensor(
    [
        node_type_to_idx[node_type_by_id[node_id]]
        for node_id in ordered_node_ids
    ],
    dtype=torch.long,
    device=DEVICE,
)

node_labels = nodes.set_index("node_id")["label"].to_dict()

model = ReMedyRGCN(
    num_nodes=len(node_to_idx),
    num_node_types=len(node_type_to_idx),
    hidden_dim=checkpoint["hidden_dim"],
    num_relations=checkpoint["num_relations"],
    dropout=checkpoint["dropout"],
).to(DEVICE)

model.load_state_dict(checkpoint["model_state_dict"])
model.eval()


# ---------------------------------------------------------
# 3. Build the same biological message-passing graph
# ---------------------------------------------------------

support_relations = [
    "targets",
    "participates_in",
    "associated_with",
]

base_relation_count = len(relation_to_idx)

src_indices = []
dst_indices = []
relation_indices = []

for row in graph.itertuples(index=False):
    if row.relation not in relation_to_idx:
        continue

    src = node_to_idx[row.source_id]
    dst = node_to_idx[row.target_id]
    rel = relation_to_idx[row.relation]

    # Forward relation.
    src_indices.append(src)
    dst_indices.append(dst)
    relation_indices.append(rel)

    # Reverse relation with a separate relation ID.
    src_indices.append(dst)
    dst_indices.append(src)
    relation_indices.append(rel + base_relation_count)

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

with torch.no_grad():
    embeddings = model.encode(
        node_type_tensor, edge_index, edge_type
    )


# ---------------------------------------------------------
# 4. Find all drug-disease pairs not labelled in the KG
# ---------------------------------------------------------

drug_ids = sorted(
    nodes.loc[nodes["node_type"] == "drug", "node_id"].tolist()
)
disease_ids = sorted(
    nodes.loc[nodes["node_type"] == "disease", "node_id"].tolist()
)

known_relations = {
    "indication",
    "off-label use",
    "contraindication",
}

# Exclude every pair that already has any of these relations.
known_pairs = set(
    edges.loc[
        edges["relation"].isin(known_relations),
        ["source_id", "target_id"],
    ].itertuples(index=False, name=None)
)

candidate_pairs = [
    (drug, disease)
    for drug in drug_ids
    for disease in disease_ids
    if (drug, disease) not in known_pairs
]

expected_candidates = len(drug_ids) * len(disease_ids) - len(
    known_pairs
)

if len(candidate_pairs) != expected_candidates:
    raise AssertionError("Candidate-pair count check failed.")

print("Drugs:", len(drug_ids))
print("Diseases:", len(disease_ids))
print("Known treatment-use/contraindication pairs:", len(known_pairs))
print("Unobserved candidate pairs:", len(candidate_pairs))


# ---------------------------------------------------------
# 5. Score candidates in batches
# ---------------------------------------------------------

scores = []

with torch.no_grad():
    for start in range(0, len(candidate_pairs), BATCH_SIZE):
        batch = candidate_pairs[start:start + BATCH_SIZE]

        drug_tensor = torch.tensor(
            [node_to_idx[drug] for drug, _ in batch],
            dtype=torch.long,
            device=DEVICE,
        )

        disease_tensor = torch.tensor(
            [node_to_idx[disease] for _, disease in batch],
            dtype=torch.long,
            device=DEVICE,
        )

        logits = model.decode(
            embeddings, drug_tensor, disease_tensor
        )

        # Ranking score, not a clinically calibrated probability.
        batch_scores = torch.sigmoid(logits).cpu().numpy()
        scores.extend(batch_scores.tolist())


# ---------------------------------------------------------
# 6. Save the full ranking and top 100 candidates
# ---------------------------------------------------------

ranking = pd.DataFrame({
    "drug_id": [drug for drug, _ in candidate_pairs],
    "disease_id": [disease for _, disease in candidate_pairs],
    "rgcn_score": scores,
})

ranking["drug_name"] = ranking["drug_id"].map(node_labels)
ranking["disease_name"] = ranking["disease_id"].map(node_labels)

ranking = ranking.sort_values(
    ["rgcn_score", "drug_id", "disease_id"],
    ascending=[False, True, True],
).reset_index(drop=True)

ranking.insert(0, "rank", np.arange(1, len(ranking) + 1))

full_path = OUTPUT_DIR / "rgcn_unobserved_candidate_ranking.csv"
top_path = OUTPUT_DIR / "rgcn_top_100_candidates.csv"

ranking.to_csv(full_path, index=False)
ranking.head(100).to_csv(top_path, index=False)

metadata = {
    "dataset_version": "ReMedy-KG v1.0",
    "checkpoint": "models/rgcn_baseline.pt",
    "score_interpretation": (
        "Sigmoid-transformed R-GCN ranking score; "
        "not a calibrated probability or clinical efficacy estimate."
    ),
    "candidate_definition": (
        "Drug-disease pairs with no indication, off-label use, "
        "or contraindication edge in the frozen KG."
    ),
    "candidate_count": len(ranking),
    "model_selected_epoch": checkpoint["best_epoch"],
    "best_validation_auprc": checkpoint["best_validation_auprc"],
}

with open(
    OUTPUT_DIR / "rgcn_candidate_ranking_metadata.json",
    "w",
    encoding="utf-8",
) as file:
    json.dump(metadata, file, indent=2)

print("\nTop 20 unobserved candidates:")
print(
    ranking[
        ["rank", "drug_name", "disease_name", "rgcn_score"]
    ].head(20).to_string(index=False)
)

print("\nSaved:")
print(full_path.relative_to(ROOT))
print(top_path.relative_to(ROOT))
print(
    (OUTPUT_DIR / "rgcn_candidate_ranking_metadata.json")
    .relative_to(ROOT)
)
