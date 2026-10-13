
from pathlib import Path
import json
import random
import pandas as pd

SEED = 42
TREATMENT_RELATIONS = {"indication", "off-label use"}
CONTRAINDICATION = "contraindication"

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "processed"
SPLITS_DIR = ROOT / "data" / "splits"
MODEL_DIR = ROOT / "data" / "model"

SPLITS_DIR.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# 1. Load the frozen dataset without modifying it.
# --------------------------------------------------
nodes = pd.read_parquet(DATA / "nodes.parquet")
edges = pd.read_parquet(DATA / "edges.parquet")
evidence = pd.read_parquet(DATA / "evidence.parquet")

node_types = nodes.set_index("node_id")["node_type"]

drug_ids = sorted(
    nodes.loc[nodes["node_type"] == "drug", "node_id"].unique()
)
disease_ids = sorted(
    nodes.loc[nodes["node_type"] == "disease", "node_id"].unique()
)

# --------------------------------------------------
# 2. Extract drug-disease pair sets.
# --------------------------------------------------
def extract_pairs(relation):
    selected = edges.loc[
        edges["relation"] == relation,
        ["source_id", "target_id"]
    ].drop_duplicates()

    # Every selected pair must be drug -> disease.
    for row in selected.itertuples(index=False):
        assert node_types.get(row.source_id) == "drug"
        assert node_types.get(row.target_id) == "disease"

    return set(map(tuple, selected.to_numpy()))


indication_pairs = extract_pairs("indication")
off_label_pairs = extract_pairs("off-label use")
contra_pairs = extract_pairs(CONTRAINDICATION)

treatment_pairs = indication_pairs | off_label_pairs

# Contradictory pairs are not assigned a supervised label.
conflict_pairs = treatment_pairs & contra_pairs

# Use only unambiguous known treatment associations as positives.
positive_pairs = treatment_pairs - contra_pairs

# Sample negatives only from pairs with no known treatment-use
# AND no known contraindication relationship.
all_possible_pairs = {
    (drug, disease)
    for drug in drug_ids
    for disease in disease_ids
}

known_pairs = treatment_pairs | contra_pairs
negative_pool = sorted(all_possible_pairs - known_pairs)

# --------------------------------------------------
# 3. Split positive pairs using a fixed random seed.
# --------------------------------------------------
rng = random.Random(SEED)

positive_list = sorted(positive_pairs)
rng.shuffle(positive_list)

n_positive = len(positive_list)
n_train = int(n_positive * 0.70)
n_validation = (n_positive - n_train) // 2
n_test = n_positive - n_train - n_validation

positive_splits = {
    "train": positive_list[:n_train],
    "validation": positive_list[
        n_train:n_train + n_validation
    ],
    "test": positive_list[n_train + n_validation:],
}

# --------------------------------------------------
# 4. Sample distinct, unobserved negative pairs.
#    Use 1 negative for every positive in each split.
# --------------------------------------------------
rng.shuffle(negative_pool)

required_negatives = n_train + n_validation + n_test
if len(negative_pool) < required_negatives:
    raise ValueError("Not enough candidate negatives.")

negative_splits = {
    "train": negative_pool[:n_train],
    "validation": negative_pool[
        n_train:n_train + n_validation
    ],
    "test": negative_pool[n_train + n_validation:
                          n_train + n_validation + n_test],
}

# --------------------------------------------------
# 5. Save labelled pairs for each split.
# --------------------------------------------------
def positive_relation_label(pair):
    labels = []
    if pair in indication_pairs:
        labels.append("indication")
    if pair in off_label_pairs:
        labels.append("off-label use")
    return ";".join(labels)


split_sizes = {}

for split_name in ["train", "validation", "test"]:
    records = []

    for drug, disease in positive_splits[split_name]:
        records.append({
            "drug_id": drug,
            "disease_id": disease,
            "label": 1,
            "label_source": positive_relation_label(
                (drug, disease)
            ),
            "split": split_name,
        })

    for drug, disease in negative_splits[split_name]:
        records.append({
            "drug_id": drug,
            "disease_id": disease,
            "label": 0,
            "label_source": "sampled_unobserved_pair",
            "split": split_name,
        })

    result = pd.DataFrame(records)
    result = result.sample(
        frac=1, random_state=SEED
    ).reset_index(drop=True)

    if result.duplicated(["drug_id", "disease_id"]).any():
        raise AssertionError(
            f"Duplicate drug-disease pairs in {split_name}"
        )

    output_path = SPLITS_DIR / f"{split_name}_pairs.csv"
    result.to_csv(output_path, index=False)

    split_sizes[split_name] = {
        "positive": int((result["label"] == 1).sum()),
        "negative": int((result["label"] == 0).sum()),
        "total": len(result),
    }

# --------------------------------------------------
# 6. Record contradictory pairs excluded from labels.
# --------------------------------------------------
conflict_records = []

for drug, disease in sorted(conflict_pairs):
    labels = []
    if (drug, disease) in indication_pairs:
        labels.append("indication")
    if (drug, disease) in off_label_pairs:
        labels.append("off-label use")

    conflict_records.append({
        "drug_id": drug,
        "disease_id": disease,
        "treatment_relations": ";".join(labels),
        "contraindication_present": True,
        "exclusion_reason": "conflicting treatment/contraindication labels",
    })

pd.DataFrame(
    conflict_records,
    columns=[
        "drug_id",
        "disease_id",
        "treatment_relations",
        "contraindication_present",
        "exclusion_reason",
    ],
).to_csv(
    SPLITS_DIR / "conflicting_pairs_excluded.csv",
    index=False
)

# --------------------------------------------------
# 7. Make a leakage-safe graph for model training.
#
# Keep training positive edges and supporting biology.
# Remove all indication/off-label edges for validation
# and test positives. Remove all three target-related
# edges for contradictory pairs.
# --------------------------------------------------
holdout_pairs = (
    set(positive_splits["validation"])
    | set(positive_splits["test"])
)
excluded_target_pairs = holdout_pairs | conflict_pairs

remove_mask = []

for row in edges.itertuples(index=False):
    pair = (row.source_id, row.target_id)

    remove_edge = (
        row.relation in TREATMENT_RELATIONS
        and pair in excluded_target_pairs
    )

    if (
        row.relation == CONTRAINDICATION
        and pair in conflict_pairs
    ):
        remove_edge = True

    remove_mask.append(remove_edge)

remove_mask = pd.Series(remove_mask, index=edges.index)

training_graph_edges = edges.loc[~remove_mask].copy()

# Verify that held-out treatment labels cannot be seen
# directly in the graph used for training.
leaked = training_graph_edges.loc[
    training_graph_edges["relation"].isin(TREATMENT_RELATIONS)
    & training_graph_edges.apply(
        lambda row: (
            row["source_id"], row["target_id"]
        ) in holdout_pairs,
        axis=1,
    )
]

if not leaked.empty:
    raise AssertionError("Validation/test treatment-edge leakage detected.")

training_edge_ids = set(training_graph_edges["edge_id"])
training_graph_evidence = evidence.loc[
    evidence["edge_id"].isin(training_edge_ids)
].copy()

# Confirm every retained graph edge has supporting evidence.
evidence_edge_ids = set(training_graph_evidence["edge_id"])
missing_evidence = training_edge_ids - evidence_edge_ids

if missing_evidence:
    raise AssertionError(
        f"{len(missing_evidence)} training graph edges lack evidence."
    )

training_graph_edges.to_parquet(
    MODEL_DIR / "training_graph_edges.parquet",
    index=False
)
training_graph_evidence.to_parquet(
    MODEL_DIR / "training_graph_evidence.parquet",
    index=False
)

# --------------------------------------------------
# 8. Save split metadata for reproducibility.
# --------------------------------------------------
metadata = {
    "dataset_version": "ReMedy-KG v1.0",
    "random_seed": SEED,
    "positive_definition": (
        "indication or off-label use, excluding pairs "
        "that also have a contraindication"
    ),
    "negative_definition": (
        "sampled pairs with no treatment-use and no "
        "contraindication relationship in the frozen KG; "
        "these are unobserved candidates, not proven negatives"
    ),
    "split_method": "random pair-level split",
    "positive_split_counts": {
        name: len(pairs)
        for name, pairs in positive_splits.items()
    },
    "split_sizes": split_sizes,
    "indication_unique_pairs": len(indication_pairs),
    "off_label_unique_pairs": len(off_label_pairs),
    "unique_treatment_pairs_before_conflict_exclusion": len(
        treatment_pairs
    ),
    "conflicting_pairs_excluded": len(conflict_pairs),
    "unambiguous_positive_pairs": len(positive_pairs),
    "contraindication_only_pairs": len(
        contra_pairs - treatment_pairs
    ),
    "available_unobserved_negative_pool": len(negative_pool),
    "training_graph_edge_count": len(training_graph_edges),
    "training_graph_evidence_rows": len(training_graph_evidence),
    "training_graph_edges_without_evidence": len(missing_evidence),
    "heldout_treatment_edges_leakage": len(leaked),
}

with open(
    SPLITS_DIR / "split_metadata.json",
    "w",
    encoding="utf-8",
) as file:
    json.dump(metadata, file, indent=2)

print("=" * 60)
print("LINK-PREDICTION SPLITS CREATED")
print("=" * 60)
print("Positive pairs:", len(positive_pairs))
print("Conflicting pairs excluded:", len(conflict_pairs))
print("Available negative candidates:", len(negative_pool))

for name, counts in split_sizes.items():
    print(
        f"{name.capitalize()}: "
        f"{counts['positive']} positives, "
        f"{counts['negative']} negatives, "
        f"{counts['total']} total"
    )

print("\nTraining graph edges:", len(training_graph_edges))
print("Training graph evidence rows:", len(training_graph_evidence))
print("Edges without evidence:", len(missing_evidence))
print("Held-out target-edge leakage:", len(leaked))
print("\nSaved files:")
print("data/splits/train_pairs.csv")
print("data/splits/validation_pairs.csv")
print("data/splits/test_pairs.csv")
print("data/splits/conflicting_pairs_excluded.csv")
print("data/splits/split_metadata.json")
print("data/model/training_graph_edges.parquet")
print("data/model/training_graph_evidence.parquet")
print("\nOriginal frozen data files were not modified.")
