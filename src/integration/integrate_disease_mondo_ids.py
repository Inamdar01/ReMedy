from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

DISEASE_PROTEIN_FILE = (
    ROOT / "data/interim/primekg_disease_protein_uniprot.csv"
)

MONDO_MAPPING_FILE = (
    ROOT / "data/interim/selected_diseases_mondo_mapping.csv"
)

OUT_FILE = (
    ROOT / "data/interim/disease_uniprot_relations_canonical.csv"
)

UNRESOLVED_FILE = (
    ROOT / "reports/unresolved_disease_mondo_mapping.csv"
)

REPORT_FILE = (
    ROOT / "reports/disease_mondo_integration_report.md"
)


def normalize_name(value):
    value = str(value).strip().lower()

    # PrimeKG uses "dementia (disease)"
    # while MONDO uses "dementia".
    if value == "dementia (disease)":
        value = "dementia"

    return value


def main():

    print("=" * 70)
    print("REMEDY - DISEASE → CANONICAL MONDO INTEGRATION")
    print("=" * 70)

    # ---------------------------------------------------------
    # 1. Load disease-protein-UniProt layer
    # ---------------------------------------------------------
    print("\nLoading disease-protein-UniProt relations...")

    disease_protein = pd.read_csv(
        DISEASE_PROTEIN_FILE,
        dtype=str
    ).fillna("")

    print(
        f"Input disease-UniProt rows: "
        f"{len(disease_protein):,}"
    )

    # ---------------------------------------------------------
    # 2. Load MONDO mappings
    # ---------------------------------------------------------
    print("Loading MONDO mappings...")

    mondo = pd.read_csv(
        MONDO_MAPPING_FILE,
        dtype=str
    ).fillna("")

    # Create normalized join key
    disease_protein["disease_join_name"] = (
        disease_protein["disease_name"]
        .apply(normalize_name)
    )

    mondo["disease_join_name"] = (
        mondo["selected_disease"]
        .apply(normalize_name)
    )

    # ---------------------------------------------------------
    # 3. Check MONDO mapping uniqueness
    # ---------------------------------------------------------
    duplicate_names = (
        mondo.groupby("disease_join_name")["mondo_id"]
        .nunique()
    )

    ambiguous_names = duplicate_names[
        duplicate_names > 1
    ]

    print(
        f"Ambiguous disease names in MONDO mapping: "
        f"{len(ambiguous_names):,}"
    )

    # Keep one row per disease mapping
    mondo_unique = mondo.drop_duplicates(
        subset=["disease_join_name"]
    )[
        [
            "disease_join_name",
            "selected_disease",
            "mondo_id",
            "mondo_name",
        ]
    ]

    # ---------------------------------------------------------
    # 4. Join
    # ---------------------------------------------------------
    result = disease_protein.merge(
        mondo_unique,
        on="disease_join_name",
        how="left"
    )

    # ---------------------------------------------------------
    # 5. Find unresolved diseases
    # ---------------------------------------------------------
    unresolved = (
        result[
            result["mondo_id"].isna()
            | (result["mondo_id"] == "")
        ][
            [
                "disease_primekg_id",
                "disease_name",
                "disease_source",
            ]
        ]
        .drop_duplicates()
    )

    # ---------------------------------------------------------
    # 6. Keep only canonicalized rows
    # ---------------------------------------------------------
    canonical = result[
        result["mondo_id"].notna()
        & (result["mondo_id"] != "")
    ].copy()

    # ---------------------------------------------------------
    # 7. Select clean columns
    # ---------------------------------------------------------
    canonical = canonical[
        [
            "mondo_id",
            "mondo_name",
            "selected_disease",
            "disease_primekg_id",
            "disease_name",
            "disease_source",
            "protein_entrez_id",
            "protein_gene_symbol",
            "protein_source",
            "hgnc_id",
            "hgnc_symbol",
            "hgnc_name",
            "uniprot_accession",
        ]
    ].drop_duplicates()

    # ---------------------------------------------------------
    # 8. Save
    # ---------------------------------------------------------
    OUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    canonical.to_csv(
        OUT_FILE,
        index=False
    )

    UNRESOLVED_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    unresolved.to_csv(
        UNRESOLVED_FILE,
        index=False
    )

    # ---------------------------------------------------------
    # 9. Statistics
    # ---------------------------------------------------------
    unique_diseases = canonical["mondo_id"].nunique()
    unique_proteins = canonical["uniprot_accession"].nunique()
    rows = len(canonical)

    # ---------------------------------------------------------
    # 10. Disease coverage
    # ---------------------------------------------------------
    coverage = (
        canonical.groupby(
            ["mondo_id", "mondo_name"]
        )
        .agg(
            unique_uniprot=(
                "uniprot_accession",
                "nunique"
            ),
            relation_rows=(
                "uniprot_accession",
                "size"
            ),
        )
        .reset_index()
        .sort_values("mondo_name")
    )

    # ---------------------------------------------------------
    # 11. Report
    # ---------------------------------------------------------
    report = f"""# Disease → Canonical MONDO Integration Report

## Inputs

- Disease-UniProt layer:
  `data/interim/primekg_disease_protein_uniprot.csv`
- Selected MONDO mappings:
  `data/interim/selected_diseases_mondo_mapping.csv`

## Results

- Canonical diseases represented: {unique_diseases}
- Unique UniProt proteins: {unique_proteins:,}
- Canonical disease-UniProt rows: {rows:,}
- Unresolved disease mappings: {len(unresolved):,}

## Disease coverage

{coverage.to_markdown(index=False)}

## Identifier policy

PrimeKG disease identifiers are preserved as provenance fields.

The canonical disease identifier used by ReMedy is the MONDO ID.

No disease identifier is inferred from the disease name beyond
the explicit reviewed mapping table.

## Outputs

- `data/interim/disease_uniprot_relations_canonical.csv`
- `reports/unresolved_disease_mondo_mapping.csv`
"""

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8"
    )

    # ---------------------------------------------------------
    # 12. Terminal output
    # ---------------------------------------------------------
    print()
    print("=" * 70)
    print("DISEASE → MONDO INTEGRATION COMPLETE")
    print("=" * 70)
    print(f"Canonical diseases: {unique_diseases}")
    print(f"Unique UniProt proteins: {unique_proteins:,}")
    print(f"Canonical disease-UniProt rows: {rows:,}")
    print(f"Unresolved disease mappings: {len(unresolved):,}")

    print("\nDisease coverage:")
    print(
        coverage.to_string(index=False)
    )

    print()
    print(f"Saved: {OUT_FILE}")
    print(f"Saved: {UNRESOLVED_FILE}")
    print(f"Saved: {REPORT_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    main()
    