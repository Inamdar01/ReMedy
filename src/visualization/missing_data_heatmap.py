from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "reports" / "visualizations"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Load frozen KG files — read only
nodes = pd.read_parquet(ROOT / "data" / "processed" / "nodes.parquet")
edges = pd.read_parquet(ROOT / "data" / "processed" / "edges.parquet")
evidence = pd.read_parquet(ROOT / "data" / "processed" / "evidence.parquet")


datasets = {
    "Nodes": nodes,
    "Edges": edges,
    "Evidence": evidence,
}


# Calculate missing-value percentage for every field
rows = []

for dataset_name, df in datasets.items():
    for column in df.columns:
        series = df[column]

        # Actual null values
        missing = series.isna().sum()

        # Empty strings are also treated as missing
        if pd.api.types.is_object_dtype(series):
            missing += (series.astype("string").str.strip() == "").sum()

        total = len(df)
        percentage = (missing / total) * 100 if total else 0

        rows.append(
            {
                "dataset": dataset_name,
                "field": column,
                "missing_count": int(missing),
                "total_count": int(total),
                "missing_percentage": percentage,
            }
        )


missing_df = pd.DataFrame(rows)


# Save detailed missing-data statistics
missing_df.to_csv(
    OUTPUT_DIR / "missing_data_statistics.csv",
    index=False,
)


# Create matrix for heatmap
heatmap_data = missing_df.pivot(
    index="dataset",
    columns="field",
    values="missing_percentage",
)


# Create heatmap
fig, ax = plt.subplots(figsize=(16, 5))

image = ax.imshow(
    heatmap_data.values,
    aspect="auto",
    interpolation="nearest",
)

ax.set_title("ReMedy-KG v1.0 — Missing Data Heatmap")
ax.set_xlabel("Field")
ax.set_ylabel("Dataset")

ax.set_xticks(range(len(heatmap_data.columns)))
ax.set_xticklabels(
    heatmap_data.columns,
    rotation=45,
    ha="right",
)

ax.set_yticks(range(len(heatmap_data.index)))
ax.set_yticklabels(heatmap_data.index)

# Add percentages to cells
for row in range(heatmap_data.shape[0]):
    for col in range(heatmap_data.shape[1]):
        value = heatmap_data.iloc[row, col]
        ax.text(
            col,
            row,
            f"{value:.1f}%",
            ha="center",
            va="center",
            fontsize=8,
        )

colorbar = fig.colorbar(image, ax=ax)
colorbar.set_label("Missing values (%)")

plt.tight_layout()


output_path = OUTPUT_DIR / "missing_data_heatmap.png"

plt.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()


print("=" * 60)
print("ReMedy-KG v1.0 — MISSING DATA ANALYSIS")
print("=" * 60)
print()

print(
    missing_df[
        ["dataset", "field", "missing_count", "total_count", "missing_percentage"]
    ].to_string(index=False)
)

print()
print(f"Saved: {output_path}")
print(
    f"Saved: {OUTPUT_DIR / 'missing_data_statistics.csv'}"
)
print()
print("Frozen Parquet files were read only; no dataset files were modified.")
