
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "reports" / "model_evaluation"

input_path = RESULTS / "candidate_rank_comparison.csv"

if not input_path.exists():
    raise FileNotFoundError(f"Missing file: {input_path}")

df = pd.read_csv(input_path)

required = {
    "drug_id",
    "drug_name",
    "disease_id",
    "disease_name",
    "rgcn_score",
    "shared_protein_count",
    "shared_pathway_count",
    "biological_similarity_score",
}

missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing columns: {sorted(missing)}")

# A candidate has direct graph-derived support when at least
# one protein or pathway is shared between the drug and disease.
df["has_shared_protein"] = df["shared_protein_count"] > 0
df["has_shared_pathway"] = df["shared_pathway_count"] > 0

df["has_graph_support"] = (
    df["has_shared_protein"] | df["has_shared_pathway"]
)

supported = df[df["has_graph_support"]].copy()
unsupported = df[~df["has_graph_support"]].copy()

supported = supported.sort_values(
    ["rgcn_score", "biological_similarity_score"],
    ascending=False,
)

unsupported = unsupported.sort_values(
    "rgcn_score",
    ascending=False,
)

# Save the full subsets and a compact top-100 shortlist.
supported.to_csv(
    RESULTS / "rgcn_candidates_with_graph_support.csv",
    index=False,
)

unsupported.to_csv(
    RESULTS / "rgcn_candidates_without_shared_support.csv",
    index=False,
)

supported.head(100).to_csv(
    RESULTS / "rgcn_top_100_supported_candidates.csv",
    index=False,
)

unsupported.head(100).to_csv(
    RESULTS / "rgcn_top_100_without_shared_support.csv",
    index=False,
)

print("=" * 65)
print("EVIDENCE-SUPPORTED CANDIDATE AUDIT")
print("=" * 65)

print("Total candidates:", len(df))
print("With shared-protein or pathway support:", len(supported))
print("Without either type of support:", len(unsupported))

print(
    "Supported candidate percentage:",
    round(100 * len(supported) / len(df), 2),
    "%",
)

print("\nTop 20 R-GCN candidates WITH graph support:")
columns = [
    "drug_name",
    "disease_name",
    "rgcn_score",
    "shared_protein_count",
    "shared_pathway_count",
    "biological_similarity_score",
]

print(supported[columns].head(20).to_string(index=False))

print("\nTop diseases in the supported top 100:")
print(supported.head(100)["disease_name"].value_counts().to_string())

print("\nSaved files:")
print("reports/model_evaluation/rgcn_candidates_with_graph_support.csv")
print("reports/model_evaluation/rgcn_candidates_without_shared_support.csv")
print("reports/model_evaluation/rgcn_top_100_supported_candidates.csv")
print("reports/model_evaluation/rgcn_top_100_without_shared_support.csv")

print("\nAudit completed.")
