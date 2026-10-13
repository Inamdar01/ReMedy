
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "processed"

nodes = pd.read_parquet(DATA / "nodes.parquet")
edges = pd.read_parquet(DATA / "edges.parquet")

# Map each canonical node ID to its entity type.
node_types = nodes.set_index("node_id")["node_type"]

# Keep only the three drug-disease relationship types.
relations = {"indication", "off-label use", "contraindication"}
df = edges[edges["relation"].isin(relations)].copy()

df["source_type"] = df["source_id"].map(node_types)
df["target_type"] = df["target_id"].map(node_types)

# Check that these relations connect drugs to diseases.
valid = (
    (df["source_type"] == "drug")
    & (df["target_type"] == "disease")
)

print("=" * 60)
print("TREATMENT LABEL AUDIT")
print("=" * 60)

print("\nRaw edge counts:")
print(df["relation"].value_counts().to_string())

print("\nEndpoint type audit:")
print("Selected edges:", len(df))
print("Drug -> disease:", int(valid.sum()))
print("Other endpoint combinations:", int((~valid).sum()))

if (~valid).any():
    print("\nUnexpected endpoint examples:")
    print(
        df.loc[
            ~valid,
            ["source_id", "relation", "target_id",
             "source_type", "target_type"]
        ].head(10).to_string(index=False)
    )

# Continue with correctly typed drug -> disease edges only.
df = df.loc[valid].copy()

def get_pairs(relation):
    selected = df.loc[
        df["relation"] == relation,
        ["source_id", "target_id"]
    ].drop_duplicates()

    return set(map(tuple, selected.to_numpy()))

indication = get_pairs("indication")
off_label = get_pairs("off-label use")
contraindication = get_pairs("contraindication")

# A treatment-positive pair has at least one treatment-use relation.
treatment = indication | off_label

drugs = set(nodes.loc[nodes["node_type"] == "drug", "node_id"])
diseases = set(nodes.loc[nodes["node_type"] == "disease", "node_id"])

all_possible_pairs = {
    (drug, disease)
    for drug in drugs
    for disease in diseases
}

known_treatment_or_contra = treatment | contraindication
unknown_pairs = all_possible_pairs - known_treatment_or_contra

print("\nUnique drug-disease pair counts:")
print("Indication pairs:", len(indication))
print("Off-label use pairs:", len(off_label))
print("Unique treatment-use pairs:", len(treatment))
print("Contraindication pairs:", len(contraindication))
print("Pairs labelled both treatment-use and contraindication:",
      len(treatment & contraindication))
print("Pairs labelled both indication and off-label use:",
      len(indication & off_label))
print("All possible drug-disease pairs:", len(all_possible_pairs))
print("Unlabelled, non-contraindicated pairs:", len(unknown_pairs))

print("\nAudit completed.")
