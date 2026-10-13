
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "reports" / "model_evaluation"

ranking_path = RESULTS / "rgcn_unobserved_candidate_ranking.csv"
metrics_path = RESULTS / "model_comparison.csv"

if not ranking_path.exists():
    raise FileNotFoundError(
        f"Candidate ranking not found: {ranking_path}"
    )

df = pd.read_csv(ranking_path)

if df.empty:
    raise ValueError("Candidate ranking file is empty.")

required = {
    "rank", "drug_id", "disease_id", "rgcn_score"
}

if not required.issubset(df.columns):
    raise ValueError(
        f"Missing columns: {required - set(df.columns)}"
    )

scores = df["rgcn_score"].astype(float)

print("=" * 65)
print("REMEDY CANDIDATE RANKING AUDIT")
print("=" * 65)

print("\nDataset checks")
print("Candidate rows:", len(df))
print("Unique drugs:", df["drug_id"].nunique())
print("Unique diseases:", df["disease_id"].nunique())
print(
    "Duplicate drug-disease pairs:",
    df.duplicated(["drug_id", "disease_id"]).sum()
)
print("Missing scores:", scores.isna().sum())

print("\nScore distribution")
print(scores.describe(
    percentiles=[0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
).to_string())

print("\nNumber of unique scores, rounded:")
for decimals in [3, 4, 6]:
    print(
        f"{decimals} decimal places:",
        scores.round(decimals).nunique()
    )

print("\nScore separation")
print("Top 10 mean:", scores.head(10).mean())
print("Top 100 mean:", scores.head(100).mean())
print("Top 1000 mean:", scores.head(1000).mean())
print("Overall mean:", scores.mean())
print("Top score minus median:", scores.max() - scores.median())
print("Top score minus minimum:", scores.max() - scores.min())

print("\nMost common diseases in the top 100")
print(
    df.head(100)["disease_name"]
    .value_counts()
    .to_string()
)

print("\nMost common drugs in the top 100")
print(
    df.head(100)["drug_name"]
    .value_counts()
    .head(20)
    .to_string()
)

print("\nFirst 10 candidates")
columns = [
    col for col in [
        "rank", "drug_id", "drug_name",
        "disease_id", "disease_name", "rgcn_score"
    ] if col in df.columns
]
print(df[columns].head(10).to_string(index=False))

if metrics_path.exists():
    print("\nBaseline vs R-GCN metrics")
    comparison = pd.read_csv(metrics_path)
    print(comparison.to_string(index=False))

print("\nRanking audit completed.")
