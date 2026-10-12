from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# Project paths
ROOT = Path(__file__).resolve().parents[2]
NODES_PATH = ROOT / "data" / "processed" / "nodes.parquet"
OUTPUT_DIR = ROOT / "reports" / "visualizations"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Load frozen KG nodes (read-only)
nodes = pd.read_parquet(NODES_PATH)


# Count nodes by type
node_counts = (
    nodes["node_type"]
    .value_counts()
    .sort_values(ascending=False)
)


# Save the underlying counts for reproducibility
node_counts.to_csv(
    OUTPUT_DIR / "kg_composition_counts.csv",
    header=["node_count"]
)


# Create visualization
plt.figure(figsize=(10, 6))

node_counts.plot(kind="bar")

plt.title("ReMedy-KG v1.0 — Node Type Composition")
plt.xlabel("Node Type")
plt.ylabel("Number of Nodes")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()


# Save figure
output_path = OUTPUT_DIR / "kg_composition.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()


print("=" * 60)
print("ReMedy-KG v1.0 — KG COMPOSITION")
print("=" * 60)
print()
print(node_counts.to_string())
print()
print(f"Total nodes: {len(nodes):,}")
print(f"Saved: {output_path}")
print()
print("Frozen Parquet files were read only; no dataset files were modified.")
