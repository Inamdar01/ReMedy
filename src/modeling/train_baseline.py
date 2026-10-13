
from pathlib import Path
import json

import joblib
import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# =========================================================
# 1. Paths and configuration
# =========================================================

ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = ROOT / "data" / "processed"
MODEL_DIR = ROOT / "data" / "model"
SPLIT_DIR = ROOT / "data" / "splits"
OUTPUT_DIR = ROOT / "reports" / "model_evaluation"
CHECKPOINT_DIR = ROOT / "models"

for directory in [MODEL_DIR, OUTPUT_DIR, CHECKPOINT_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

SEED = 42

FEATURE_COLUMNS = [
    "shared_protein_count",
    "shared_protein_jaccard",
    "shared_pathway_count",
    "shared_pathway_jaccard",
    "drug_target_count",
    "disease_protein_count",
    "drug_pathway_count",
    "disease_pathway_count",
]


# =========================================================
# 2. Load nodes and leakage-safe training graph
# =========================================================

nodes = pd.read_parquet(DATA_DIR / "nodes.parquet")

graph_path = MODEL_DIR / "training_graph_edges.parquet"

if not graph_path.exists():
    raise FileNotFoundError(
        "Run create_link_prediction_splits.py first."
    )

graph_edges = pd.read_parquet(graph_path)

node_types = nodes.set_index("node_id")["node_type"].to_dict()

drug_ids = set(
    nodes.loc[nodes["node_type"] == "drug", "node_id"]
)
disease_ids = set(
    nodes.loc[nodes["node_type"] == "disease", "node_id"]
)

# We intentionally derive features only from supporting
# biological relationships, never from treatment labels.
drug_targets = {}
disease_proteins = {}
protein_pathways = {}

# Use sets to prevent duplicate relationships inflating counts.
from collections import defaultdict

drug_targets = defaultdict(set)
disease_proteins = defaultdict(set)
protein_pathways = defaultdict(set)

used_relation_counts = defaultdict(int)

for row in graph_edges.itertuples(index=False):
    source = row.source_id
    target = row.target_id
    relation = row.relation

    if relation not in {
        "targets",
        "participates_in",
        "associated_with",
    }:
        continue

    source_type = node_types.get(source)
    target_type = node_types.get(target)

    if relation == "targets":
        expected_types = ("drug", "protein")
    elif relation == "participates_in":
        expected_types = ("protein", "pathway")
    else:
        expected_types = ("disease", "protein")

    if (source_type, target_type) != expected_types:
        raise ValueError(
            f"Unexpected endpoint types for {relation}: "
            f"{source} ({source_type}) -> "
            f"{target} ({target_type})"
        )

    if relation == "targets":
        drug_targets[source].add(target)

    elif relation == "participates_in":
        protein_pathways[source].add(target)

    elif relation == "associated_with":
        disease_proteins[source].add(target)

    used_relation_counts[relation] += 1


# =========================================================
# 3. Derive pathway sets for drugs and diseases
# =========================================================

drug_pathways = {}

for drug in drug_ids:
    pathways = set()

    for protein in drug_targets[drug]:
        pathways.update(protein_pathways[protein])

    drug_pathways[drug] = pathways


disease_pathways = {}

for disease in disease_ids:
    pathways = set()

    for protein in disease_proteins[disease]:
        pathways.update(protein_pathways[protein])

    disease_pathways[disease] = pathways


# =========================================================
# 4. Define biological features for a drug-disease pair
# =========================================================

def jaccard_similarity(set_a, set_b):
    union = set_a | set_b

    if not union:
        return 0.0

    return len(set_a & set_b) / len(union)


def extract_features(drug, disease):
    drug_protein_set = drug_targets[drug]
    disease_protein_set = disease_proteins[disease]

    drug_pathway_set = drug_pathways[drug]
    disease_pathway_set = disease_pathways[disease]

    shared_proteins = drug_protein_set & disease_protein_set
    shared_pathways = drug_pathway_set & disease_pathway_set

    return {
        "shared_protein_count": len(shared_proteins),
        "shared_protein_jaccard": jaccard_similarity(
            drug_protein_set, disease_protein_set
        ),
        "shared_pathway_count": len(shared_pathways),
        "shared_pathway_jaccard": jaccard_similarity(
            drug_pathway_set, disease_pathway_set
        ),
        "drug_target_count": len(drug_protein_set),
        "disease_protein_count": len(disease_protein_set),
        "drug_pathway_count": len(drug_pathway_set),
        "disease_pathway_count": len(disease_pathway_set),
    }


# =========================================================
# 5. Generate feature tables from saved pair splits
# =========================================================

feature_tables = {}

for split_name in ["train", "validation", "test"]:
    input_path = SPLIT_DIR / f"{split_name}_pairs.csv"

    if not input_path.exists():
        raise FileNotFoundError(f"Missing split: {input_path}")

    pairs = pd.read_csv(input_path)

    required_columns = {"drug_id", "disease_id", "label"}

    if not required_columns.issubset(pairs.columns):
        raise ValueError(
            f"Unexpected split schema in {input_path}"
        )

    feature_rows = []

    for row in pairs.itertuples(index=False):
        feature_rows.append(
            extract_features(row.drug_id, row.disease_id)
        )

    features = pd.DataFrame(feature_rows)

    data = pd.concat(
        [pairs.reset_index(drop=True), features],
        axis=1,
    )

    feature_tables[split_name] = data

    data.to_csv(
        MODEL_DIR / f"{split_name}_pair_features.csv",
        index=False,
    )

    print(
        f"{split_name.capitalize()} pairs: {len(data)}, "
        f"positive labels: {(data['label'] == 1).sum()}, "
        f"negative labels: {(data['label'] == 0).sum()}"
    )


# =========================================================
# 6. Define evaluation metrics
# =========================================================

def evaluate_scores(y_true, scores, model_name, split_name):
    y_true = np.asarray(y_true, dtype=int)
    scores = np.asarray(scores, dtype=float)

    order = np.argsort(-scores)
    k = min(10, len(y_true))

    top_k_labels = y_true[order[:k]]
    total_positives = int(y_true.sum())

    precision_at_k = float(top_k_labels.mean())

    recall_at_k = (
        float(top_k_labels.sum() / total_positives)
        if total_positives > 0
        else 0.0
    )

    return {
        "model": model_name,
        "split": split_name,
        "roc_auc": float(roc_auc_score(y_true, scores)),
        "auprc_average_precision": float(
            average_precision_score(y_true, scores)
        ),
        "precision_at_10": precision_at_k,
        "recall_at_10": recall_at_k,
    }


# =========================================================
# 7. Create similarity and degree-based baselines
# =========================================================

metric_rows = []
prediction_tables = {}

for split_name, data in feature_tables.items():
    data = data.copy()

    data["similarity_score"] = (
        data["shared_protein_jaccard"]
        + data["shared_pathway_jaccard"]
    )

    data["degree_popularity_score"] = (
        np.log1p(data["drug_target_count"])
        + np.log1p(data["disease_protein_count"])
    )

    feature_tables[split_name] = data

    if split_name in {"validation", "test"}:
        metric_rows.append(
            evaluate_scores(
                data["label"],
                data["similarity_score"],
                "Biological similarity",
                split_name,
            )
        )

        metric_rows.append(
            evaluate_scores(
                data["label"],
                data["degree_popularity_score"],
                "Degree popularity",
                split_name,
            )
        )


# =========================================================
# 8. Select logistic-regression settings using validation
# =========================================================

train = feature_tables["train"]
validation = feature_tables["validation"]
test = feature_tables["test"]

X_train = train[FEATURE_COLUMNS]
y_train = train["label"].astype(int)

X_validation = validation[FEATURE_COLUMNS]
y_validation = validation["label"].astype(int)

best_c = None
best_validation_ap = -1.0
best_validation_model = None
best_validation_probabilities = None

for c_value in [0.01, 0.1, 1.0, 10.0]:
    candidate_model = Pipeline([
        ("scaler", StandardScaler()),
        ("classifier", LogisticRegression(
            C=c_value,
            max_iter=2000,
            random_state=SEED,
        )),
    ])

    candidate_model.fit(X_train, y_train)

    validation_probabilities = (
        candidate_model.predict_proba(X_validation)[:, 1]
    )

    validation_ap = average_precision_score(
        y_validation,
        validation_probabilities,
    )

    print(
        f"C={c_value:<5} "
        f"validation AUPRC={validation_ap:.4f}"
    )

    if validation_ap > best_validation_ap:
        best_validation_ap = validation_ap
        best_c = c_value
        best_validation_model = candidate_model
        best_validation_probabilities = validation_probabilities


metric_rows.append(
    evaluate_scores(
        y_validation,
        best_validation_probabilities,
        "Logistic regression",
        "validation",
    )
)


# =========================================================
# 9. Refit on training + validation; evaluate final test
# =========================================================

train_validation = pd.concat(
    [train, validation],
    ignore_index=True,
)

final_model = Pipeline([
    ("scaler", StandardScaler()),
    ("classifier", LogisticRegression(
        C=best_c,
        max_iter=2000,
        random_state=SEED,
    )),
])

final_model.fit(
    train_validation[FEATURE_COLUMNS],
    train_validation["label"].astype(int),
)

X_test = test[FEATURE_COLUMNS]
y_test = test["label"].astype(int)

test_probabilities = final_model.predict_proba(X_test)[:, 1]

metric_rows.append(
    evaluate_scores(
        y_test,
        test_probabilities,
        "Logistic regression",
        "test",
    )
)


# =========================================================
# 10. Save predictions, metrics and model
# =========================================================

test_predictions = test[
    ["drug_id", "disease_id", "label", "label_source"]
].copy()

test_predictions["similarity_score"] = (
    test["similarity_score"].to_numpy()
)

test_predictions["degree_popularity_score"] = (
    test["degree_popularity_score"].to_numpy()
)

test_predictions["logistic_probability"] = test_probabilities

test_predictions = test_predictions.sort_values(
    "logistic_probability",
    ascending=False,
)

test_predictions.to_csv(
    OUTPUT_DIR / "baseline_test_predictions.csv",
    index=False,
)

metrics = pd.DataFrame(metric_rows)

metrics.to_csv(
    OUTPUT_DIR / "baseline_metrics.csv",
    index=False,
)

joblib.dump(
    {
        "model": final_model,
        "feature_columns": FEATURE_COLUMNS,
        "selected_C": best_c,
        "dataset_version": "ReMedy-KG v1.0",
    },
    CHECKPOINT_DIR / "logistic_regression_baseline.joblib",
)

metadata = {
    "dataset_version": "ReMedy-KG v1.0",
    "random_seed": SEED,
    "features": FEATURE_COLUMNS,
    "logistic_regression_selected_C": best_c,
    "best_validation_auprc": best_validation_ap,
    "training_positive_count": int(y_train.sum()),
    "validation_positive_count": int(y_validation.sum()),
    "test_positive_count": int(y_test.sum()),
    "feature_relations_used": dict(used_relation_counts),
    "treatment_labels_used_as_features": False,
}

with open(
    OUTPUT_DIR / "baseline_run_metadata.json",
    "w",
    encoding="utf-8",
) as file:
    json.dump(metadata, file, indent=2)


# =========================================================
# 11. Print a summary
# =========================================================

print("\n" + "=" * 65)
print("REMEDY BASELINE RESULTS")
print("=" * 65)

print("\nRelationship types used for biological features:")
for relation, count in sorted(used_relation_counts.items()):
    print(f"{relation}: {count}")

print("\nSelected logistic-regression C:", best_c)
print("Best validation AUPRC:", round(best_validation_ap, 4))

print("\nMetrics:")
print(
    metrics.to_string(
        index=False,
        float_format=lambda value: f"{value:.4f}",
    )
)

print("\nSaved files:")
print("data/model/train_pair_features.csv")
print("data/model/validation_pair_features.csv")
print("data/model/test_pair_features.csv")
print("reports/model_evaluation/baseline_metrics.csv")
print("reports/model_evaluation/baseline_test_predictions.csv")
print("reports/model_evaluation/baseline_run_metadata.json")
print("models/logistic_regression_baseline.joblib")

print("\nBaseline development completed.")
