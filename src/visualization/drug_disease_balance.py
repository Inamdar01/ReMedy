from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
EDGES_PATH = ROOT / "data" / "processed" / "edges.parquet"
OUTPUT_DIR = ROOT / "reports" / "visualizations"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Load frozen KG edges — read only
edges = pd.read_parquet(EDGES_PATH)


# Drug-Disease relations only
drug_disease_relations = edges[
    edges["relation"].isin(
        ["indication", "off-label use", "contraindication"]
    )
].copy()


# Count each Drug-Disease relation
label_counts = (
    drug_disease_relations["relation"]
    .value_counts()
    .reindex(
        ["indication", "off-label use", "contraindication"],
        fill_value=0,
    )
)


# Save counts
label_counts.to_csv(
    OUTPUT_DIR / "drug_disease_label_balance_counts.csv",
    header=["edge_count"],
)


# Create visualization
plt.figure(figsize=(9, 6))

label_counts.plot(kind="bar")

plt.title("ReMedy-KG v1.0 — Drug–Disease Relation Distribution")
plt.xlabel("Drug–Disease Relation")
plt.ylabel("Number of Edges")
plt.xticks(rotation=0)
plt.tight_layout()


output_path = OUTPUT_DIR / "drug_disease_label_balance.png"

plt.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()


print("=" * 60)
print("ReMedy-KG v1.0 — DRUG–DISEASE LABEL BALANCE")
print("=" * 60)
print()

print(label_counts.to_string())

print()
print(f"Total Drug–Disease relation edges: {len(drug_disease_relations):,}")

treatment_count = (
    label_counts["indication"] +
    label_counts["off-label use"]
)

print(f"Treatment-use edges: {treatment_count:,}")
print(
    f"Contraindication edges: "
    f"{label_counts['contraindication']:,}"
)

print()
print(f"Saved: {output_path}")
print(
    f"Saved: "
    f"{OUTPUT_DIR / 'drug_disease_label_balance_counts.csv'}"
)
print()
print("Frozen Parquet files were read only; no dataset files were modified.")
