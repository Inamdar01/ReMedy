
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "processed"

files = {
    "nodes": DATA / "nodes.parquet",
    "edges": DATA / "edges.parquet",
    "evidence": DATA / "evidence.parquet",
}

tables = {}

for name, path in files.items():
    print("\n" + "=" * 60)
    print(f"FILE: {name.upper()}")

    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    df = pd.read_parquet(path)
    tables[name] = df

    print("Shape:", df.shape)
    print("Columns:", df.columns.tolist())
    print("\nFirst 3 rows:")
    print(df.head(3).to_string(index=False))

nodes = tables["nodes"]
edges = tables["edges"]

for col in ["node_type", "type", "entity_type"]:
    if col in nodes.columns:
        print(f"\nNODE COUNTS BY {col}:")
        print(nodes[col].value_counts(dropna=False))
        break

for col in ["relation", "relation_type"]:
    if col in edges.columns:
        print(f"\nRELATION COUNTS BY {col}:")
        print(edges[col].value_counts(dropna=False))
        break

print("\nInspection completed successfully.")
