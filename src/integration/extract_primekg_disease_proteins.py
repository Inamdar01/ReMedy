from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

PRIMEKG_FILE = ROOT / "data/raw/primekg/primekg.csv"
OUT_FILE = ROOT / "data/interim/primekg_disease_protein_relations.csv"
REPORT_FILE = ROOT / "reports/primekg_disease_protein_report.md"

SELECTED_DISEASES = [
    "epilepsy",
    "anxiety disorder",
    "stroke disorder",
    "Parkinson disease",
    "Alzheimer disease",
    "migraine disorder",
    "bipolar disorder",
    "schizophrenia",
    "dementia (disease)",
    "major depressive disorder",
    "attention deficit-hyperactivity disorder",
    "multiple sclerosis",
]


def main():
    print("=" * 70)
    print("REMEDY - PRIMEKG DISEASE → PROTEIN EXTRACTION")
    print("=" * 70)

    print("\nLoading PrimeKG...")

    columns = [
        "relation",
        "x_id",
        "x_type",
        "x_name",
        "x_source",
        "y_id",
        "y_type",
        "y_name",
        "y_source",
    ]

    df = pd.read_csv(
        PRIMEKG_FILE,
        usecols=columns,
        low_memory=False
    )

    print(f"PrimeKG rows loaded: {len(df):,}")

    # ---------------------------------------------------------
    # Keep only disease_protein relationships
    # ---------------------------------------------------------
    df = df[df["relation"] == "disease_protein"].copy()

    print(f"Disease-protein rows: {len(df):,}")

    # ---------------------------------------------------------
    # Disease can appear on either side.
    # Normalize into:
    # disease -> protein
    # ---------------------------------------------------------

    rows = []

    # Case 1:
    # x = gene/protein
    # y = disease
    mask = (
        (df["x_type"] == "gene/protein")
        & (df["y_type"] == "disease")
        & (df["y_name"].isin(SELECTED_DISEASES))
    )

    a = df.loc[
        mask,
        [
            "x_id",
            "x_name",
            "x_source",
            "y_id",
            "y_name",
            "y_source",
        ]
    ].copy()

    a = a.rename(
        columns={
            "x_id": "protein_entrez_id",
            "x_name": "protein_gene_symbol",
            "x_source": "protein_source",
            "y_id": "disease_primekg_id",
            "y_name": "disease_name",
            "y_source": "disease_source",
        }
    )

    rows.append(a)

    # Case 2:
    # x = disease
    # y = gene/protein
    mask = (
        (df["x_type"] == "disease")
        & (df["y_type"] == "gene/protein")
        & (df["x_name"].isin(SELECTED_DISEASES))
    )

    b = df.loc[
        mask,
        [
            "x_id",
            "x_name",
            "x_source",
            "y_id",
            "y_name",
            "y_source",
        ]
    ].copy()

    b = b.rename(
        columns={
            "x_id": "disease_primekg_id",
            "x_name": "disease_name",
            "x_source": "disease_source",
            "y_id": "protein_entrez_id",
            "y_name": "protein_gene_symbol",
            "y_source": "protein_source",
        }
    )

    rows.append(b)

    result = pd.concat(rows, ignore_index=True)

    # ---------------------------------------------------------
    # Normalize IDs
    # ---------------------------------------------------------
    result["disease_primekg_id"] = (
        result["disease_primekg_id"]
        .astype(str)
        .str.strip()
    )

    result["protein_entrez_id"] = (
        result["protein_entrez_id"]
        .astype(str)
        .str.strip()
    )

    result["disease_name"] = (
        result["disease_name"]
        .astype(str)
        .str.strip()
    )

    result["protein_gene_symbol"] = (
        result["protein_gene_symbol"]
        .astype(str)
        .str.strip()
    )

    # ---------------------------------------------------------
    # Remove exact duplicates
    # ---------------------------------------------------------
    before = len(result)

    result = result.drop_duplicates(
        subset=[
            "disease_primekg_id",
            "protein_entrez_id",
            "disease_source",
            "protein_source",
        ]
    )

    removed = before - len(result)

    # ---------------------------------------------------------
    # Sort
    # ---------------------------------------------------------
    result = result.sort_values(
        [
            "disease_name",
            "disease_source",
            "protein_gene_symbol",
        ]
    )

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT_FILE, index=False)

    # ---------------------------------------------------------
    # Statistics
    # ---------------------------------------------------------
    disease_counts = (
        result.groupby(
            ["disease_name", "disease_source"],
            dropna=False
        )
        .agg(
            protein_count=("protein_entrez_id", "nunique"),
            relation_rows=("protein_entrez_id", "size"),
        )
        .reset_index()
        .sort_values("disease_name")
    )

    unique_diseases = result["disease_name"].nunique()
    unique_proteins = result["protein_entrez_id"].nunique()

    # ---------------------------------------------------------
    # Report
    # ---------------------------------------------------------
    report_lines = [
        "# PrimeKG Disease-Protein Extraction Report",
        "",
        "## Input",
        "",
        "- Source: `data/raw/primekg/primekg.csv`",
        "- Relation: `disease_protein`",
        "",
        "## Selected diseases",
        "",
    ]

    for disease in SELECTED_DISEASES:
        report_lines.append(f"- {disease}")

    report_lines.extend(
        [
            "",
            "## Results",
            "",
            f"- Disease-protein rows before deduplication: {before:,}",
            f"- Exact duplicates removed: {removed:,}",
            f"- Final disease-protein rows: {len(result):,}",
            f"- Selected diseases represented: {unique_diseases}",
            f"- Unique protein identifiers: {unique_proteins:,}",
            "",
            "## Coverage",
            "",
            disease_counts.to_markdown(index=False),
            "",
            "## Identifier note",
            "",
            "PrimeKG protein identifiers are retained as supplied by PrimeKG.",
            "Disease identifiers are also retained as supplied and are not",
            "treated as canonical MONDO identifiers until an explicit mapping",
            "step is completed.",
            "",
            "## Output",
            "",
            "`data/interim/primekg_disease_protein_relations.csv`",
        ]
    )

    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)
    REPORT_FILE.write_text(
        "\n".join(report_lines),
        encoding="utf-8"
    )

    # ---------------------------------------------------------
    # Terminal output
    # ---------------------------------------------------------
    print()
    print("=" * 70)
    print("DISEASE-PROTEIN EXTRACTION COMPLETE")
    print("=" * 70)
    print(f"Selected diseases represented: {unique_diseases}")
    print(f"Unique proteins: {unique_proteins:,}")
    print(f"Final relation rows: {len(result):,}")
    print(f"Duplicates removed: {removed:,}")

    print("\nDisease coverage:")
    print(disease_counts.to_string(index=False))

    print()
    print(f"Saved: {OUT_FILE}")
    print(f"Saved: {REPORT_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    main()