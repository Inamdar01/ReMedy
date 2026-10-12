from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
EDGES_PATH = ROOT / "data" / "processed" / "edges.parquet"
OUTPUT_DIR = ROOT / "reports" / "visualizations"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Load frozen KG edges — read only
edges = pd.read_parquet(EDGES_PATH)


# Count total degree for every node.
# Each edge contributes one connection to both its source and target.
source_degree = edges["source_id"].value_counts()
target_degree = edges["target_id"].value_counts()

degree = source_degree.add(target_degree, fill_value=0)
degree = degree.astype(int).sort_values(ascending=False)


# Save complete degree statistics
degree_df = degree.rename("degree").reset_index()
degree_df.columns = ["node_id", "degree"]

degree_df.to_csv(
    OUTPUT_DIR / "degree_distribution_counts.csv",
    index=False,
)


# Basic statistics
print("=" * 60)
print("ReMedy-KG v1.0 — DEGREE DISTRIBUTION")
print("=" * 60)
print()

print(f"Nodes with at least one edge: {len(degree):,}")
print(f"Maximum degree: {degree.max():,}")
print(f"Minimum degree: {degree.min():,}")
print(f"Median degree: {degree.median():,.0f}")
print(f"Mean degree: {degree.mean():,.2f}")

print()
print("Top 20 most connected nodes:")
print(degree_df.head(20).to_string(index=False))


# Plot degree distribution
plt.figure(figsize=(10, 6))

plt.hist(
    degree,
    bins=50,
)

plt.title("ReMedy-KG v1.0 — Node Degree Distribution")
plt.xlabel("Node Degree")
plt.ylabel("Number of Nodes")
plt.tight_layout()


output_path = OUTPUT_DIR / "degree_distribution.png"

plt.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()


print()
print(f"Saved: {output_path}")
print(f"Saved: {OUTPUT_DIR / 'degree_distribution_counts.csv'}")
print()
print("Frozen Parquet files were read only; no dataset files were modified.")
