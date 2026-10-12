from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
EDGES_PATH = ROOT / "data" / "processed" / "edges.parquet"
OUTPUT_DIR = ROOT / "reports" / "visualizations"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Load frozen KG edges (read-only)
edges = pd.read_parquet(EDGES_PATH)


# Count relations
relation_counts = (
    edges["relation"]
    .value_counts()
    .sort_values(ascending=False)
)


# Save counts
relation_counts.to_csv(
    OUTPUT_DIR / "relation_distribution_counts.csv",
    header=["edge_count"]
)


# Create visualization
plt.figure(figsize=(11, 6))

relation_counts.plot(kind="bar")

plt.title("ReMedy-KG v1.0 — Relation Distribution")
plt.xlabel("Relation")
plt.ylabel("Number of Edges")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()


# Save figure
output_path = OUTPUT_DIR / "relation_distribution.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()


print("=" * 60)
print("ReMedy-KG v1.0 — RELATION DISTRIBUTION")
print("=" * 60)
print()
print(relation_counts.to_string())
print()
print(f"Total edges: {len(edges):,}")
print(f"Saved: {output_path}")
print()
print("Frozen Parquet files were read only; no dataset files were modified.")
