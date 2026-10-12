from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
EVIDENCE_PATH = ROOT / "data" / "processed" / "evidence.parquet"
OUTPUT_DIR = ROOT / "reports" / "visualizations"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Load frozen KG evidence (read-only)
evidence = pd.read_parquet(EVIDENCE_PATH)


# Count evidence records by source
source_counts = (
    evidence["source"]
    .value_counts()
    .sort_values(ascending=False)
)


# Save counts
source_counts.to_csv(
    OUTPUT_DIR / "source_contribution_counts.csv",
    header=["evidence_count"]
)


# Create visualization
plt.figure(figsize=(9, 6))

source_counts.plot(kind="bar")

plt.title("ReMedy-KG v1.0 — Evidence Source Contribution")
plt.xlabel("Evidence Source")
plt.ylabel("Number of Evidence Records")
plt.xticks(rotation=0)
plt.tight_layout()


# Save figure
output_path = OUTPUT_DIR / "source_contribution.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()


print("=" * 60)
print("ReMedy-KG v1.0 — SOURCE CONTRIBUTION")
print("=" * 60)
print()
print(source_counts.to_string())
print()
print(f"Total evidence records: {len(evidence):,}")
print(f"Saved: {output_path}")
print()
print("Frozen Parquet files were read only; no dataset files were modified.")
