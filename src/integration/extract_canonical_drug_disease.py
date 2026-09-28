from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

PRIMEKG_FILE = ROOT / "data/raw/primekg/primekg.csv"

DRUG_MAPPING_FILE = (
    ROOT / "data/raw/chembl/drugbank_to_chembl_validated.csv"
)

MONDO_MAPPING_FILE = (
    ROOT / "data/interim/selected_diseases_mondo_mapping.csv"
)

OUT_FILE = (
    ROOT / "data/interim/drug_disease_relations_canonical.csv"
)

UNRESOLVED_FILE = (
    ROOT / "reports/unresolved_drug_disease_edges.csv"
)

REPORT_FILE = (
    ROOT / "reports/drug_disease_canonical_report.md"
)


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


def normalize_disease_name(value):
    value = str(value).strip().lower()

    # PrimeKG uses this label while MONDO uses "dementia"
    if value == "dementia (disease)":
        return "dementia"

    return value


def main():

    print("=" * 70)
    print("REMEDY - CANONICAL DRUG → DISEASE EDGE EXTRACTION")
    print("=" * 70)

    # ---------------------------------------------------------
    # 1. Load validated DrugBank → ChEMBL mapping
    # ---------------------------------------------------------
    print("\nLoading DrugBank → ChEMBL mapping...")

    drug_mapping = pd.read_csv(
        DRUG_MAPPING_FILE,
        dtype=str
    ).fillna("")

    drug_mapping = drug_mapping[
        [
            "drugbank_id",
            "drug_name",
            "chembl_id",
            "chembl_pref_name",
            "mapping_status",
            "mapping_reason",
        ]
    ].drop_duplicates()

    print(
        f"Validated DrugBank rows: {len(drug_mapping):,}"
    )

    # ---------------------------------------------------------
    # 2. Load disease → MONDO mapping
    # ---------------------------------------------------------
    print("Loading disease → MONDO mapping...")

    mondo = pd.read_csv(
        MONDO_MAPPING_FILE,
        dtype=str
    ).fillna("")

    mondo["disease_join_name"] = (
        mondo["selected_disease"]
        .apply(normalize_disease_name)
    )

    mondo = mondo[
        [
            "disease_join_name",
            "selected_disease",
            "mondo_id",
            "mondo_name",
        ]
    ].drop_duplicates()

    print(
        f"Canonical disease mappings: "
        f"{mondo['mondo_id'].nunique():,}"
    )

    # ---------------------------------------------------------
    # 3. Load only required PrimeKG columns
    # ---------------------------------------------------------
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

    print(
        f"PrimeKG rows loaded: {len(df):,}"
    )

    # ---------------------------------------------------------
    # 4. Keep only selected drug-disease relationships
    # ---------------------------------------------------------
    selected = set(SELECTED_DISEASES)

    forward = (
        (df["x_type"] == "drug")
        & (df["y_type"] == "disease")
        & (df["y_name"].isin(selected))
    )

    reverse = (
        (df["x_type"] == "disease")
        & (df["y_type"] == "drug")
        & (df["x_name"].isin(selected))
    )

    edges = df[forward | reverse].copy()

    print(
        f"Selected PrimeKG drug-disease rows: "
        f"{len(edges):,}"
    )

    # ---------------------------------------------------------
    # 5. Normalize orientation to DrugBank → Disease
    # ---------------------------------------------------------
    edges["drugbank_id"] = edges["x_id"].where(
        edges["x_type"] == "drug",
        edges["y_id"]
    )

    edges["drug_name"] = edges["x_name"].where(
        edges["x_type"] == "drug",
        edges["y_name"]
    )

    edges["primekg_disease_id"] = edges["y_id"].where(
        edges["y_type"] == "disease",
        edges["x_id"]
    )

    edges["disease_name"] = edges["y_name"].where(
        edges["y_type"] == "disease",
        edges["x_name"]
    )

    edges["disease_source"] = edges["y_source"].where(
        edges["y_type"] == "disease",
        edges["x_source"]
    )

    # ---------------------------------------------------------
    # 6. Normalize disease join key
    # ---------------------------------------------------------
    edges["disease_join_name"] = (
        edges["disease_name"]
        .apply(normalize_disease_name)
    )

    # ---------------------------------------------------------
    # 7. Remove duplicate orientations
    # ---------------------------------------------------------
    before_dedup = len(edges)

    edges = edges.drop_duplicates(
        subset=[
            "drugbank_id",
            "primekg_disease_id",
            "disease_name",
            "relation",
        ]
    )

    duplicate_rows_removed = (
        before_dedup - len(edges)
    )

    print(
        f"Duplicate orientation rows removed: "
        f"{duplicate_rows_removed:,}"
    )

    # ---------------------------------------------------------
    # 8. Join DrugBank → ChEMBL
    # ---------------------------------------------------------
    result = edges.merge(
        drug_mapping,
        on="drugbank_id",
        how="left",
        suffixes=("", "_mapping")
    )

    # ---------------------------------------------------------
    # 9. Join disease → MONDO
    # ---------------------------------------------------------
    result = result.merge(
        mondo,
        on="disease_join_name",
        how="left"
    )

    # ---------------------------------------------------------
    # 10. Identify unresolved rows
    # ---------------------------------------------------------
    unresolved = result[
        result["chembl_id"].isna()
        | (result["chembl_id"] == "")
        | result["mondo_id"].isna()
        | (result["mondo_id"] == "")
    ].copy()

    unresolved = unresolved[
        [
            "drugbank_id",
            "drug_name",
            "relation",
            "primekg_disease_id",
            "disease_name",
            "chembl_id",
            "mondo_id",
        ]
    ].drop_duplicates()

    # ---------------------------------------------------------
    # 11. Keep only completely resolved edges
    # ---------------------------------------------------------
    result = result[
        result["chembl_id"].notna()
        & (result["chembl_id"] != "")
        & result["mondo_id"].notna()
        & (result["mondo_id"] != "")
    ].copy()

    # ---------------------------------------------------------
    # 12. Canonical edge table
    # ---------------------------------------------------------
    result = result[
        [
            "chembl_id",
            "chembl_pref_name",
            "drugbank_id",
            "drug_name",
            "mondo_id",
            "mondo_name",
            "disease_name",
            "primekg_disease_id",
            "relation",
            "disease_source",
            "mapping_status",
            "mapping_reason",
        ]
    ].copy()

    # ---------------------------------------------------------
    # 12A. Collapse duplicate canonical edges
    #
    # Different DrugBank records can map to the same ChEMBL
    # drug. The KG should contain one canonical edge:
    #
    # ChEMBL Drug ──relation──> MONDO Disease
    #
    # while retaining all DrugBank records as provenance.
    # ---------------------------------------------------------
    canonical_before = len(result)

    def join_unique(values):
        values = sorted({
            str(v).strip()
            for v in values
            if str(v).strip()
        })
        return " | ".join(values)

    result = (
        result
        .groupby(
            [
                "chembl_id",
                "mondo_id",
                "relation",
            ],
            as_index=False
        )
        .agg(
            chembl_pref_name=(
                "chembl_pref_name",
                join_unique
            ),
            drugbank_id=(
                "drugbank_id",
                join_unique
            ),
            drug_name=(
                "drug_name",
                join_unique
            ),
            mondo_name=(
                "mondo_name",
                join_unique
            ),
            disease_name=(
                "disease_name",
                join_unique
            ),
            primekg_disease_id=(
                "primekg_disease_id",
                join_unique
            ),
            disease_source=(
                "disease_source",
                join_unique
            ),
            mapping_status=(
                "mapping_status",
                join_unique
            ),
            mapping_reason=(
                "mapping_reason",
                join_unique
            ),
        )
    )

    canonical_duplicates_removed = (
        canonical_before - len(result)
    )

    result = result.sort_values(
        [
            "mondo_name",
            "relation",
            "chembl_id",
        ]
    ).reset_index(drop=True)

    # ---------------------------------------------------------
    # 13. Save outputs
    # ---------------------------------------------------------
    OUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    result.to_csv(
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
    # 14. Statistics
    # ---------------------------------------------------------
    relation_counts = (
        result["relation"]
        .value_counts()
        .rename_axis("relation")
        .reset_index(name="rows")
    )

    disease_counts = (
        result.groupby(
            ["mondo_id", "mondo_name"]
        )
        .agg(
            unique_drugs=(
                "chembl_id",
                "nunique"
            ),
            rows=(
                "chembl_id",
                "size"
            ),
        )
        .reset_index()
        .sort_values("mondo_name")
    )

    unique_drugs = result["chembl_id"].nunique()
    unique_diseases = result["mondo_id"].nunique()

    # ---------------------------------------------------------
    # 15. Report
    # ---------------------------------------------------------
    report = f"""# Canonical Drug → Disease Edge Report

## Inputs

- PrimeKG:
  `data/raw/primekg/primekg.csv`
- DrugBank → ChEMBL:
  `data/raw/chembl/drugbank_to_chembl_validated.csv`
- Disease → MONDO:
  `data/interim/selected_diseases_mondo_mapping.csv`

## Results

- PrimeKG selected rows before deduplication: {len(edges) + duplicate_rows_removed:,}
- Duplicate orientation rows removed: {duplicate_rows_removed:,}
    - Canonical duplicate edges collapsed: {canonical_duplicates_removed:,}
- Canonical resolved edges: {len(result):,}
- Unique ChEMBL drugs: {unique_drugs:,}
- Unique MONDO diseases: {unique_diseases:,}
- Unresolved edges: {len(unresolved):,}

## Relation counts

{relation_counts.to_markdown(index=False)}

## Disease coverage

{disease_counts.to_markdown(index=False)}

## Relation semantics

PrimeKG relations are preserved separately:

- `indication`
- `off-label use`
- `contraindication`

They are not collapsed into a single treatment relation.

## Identifier policy

- Drug node identifier: ChEMBL ID
- Disease node identifier: MONDO ID
- PrimeKG DrugBank and disease identifiers are retained as provenance.
- No unmapped identifier is guessed.

## Outputs

- `data/interim/drug_disease_relations_canonical.csv`
- `reports/unresolved_drug_disease_edges.csv`
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
    # 16. Terminal output
    # ---------------------------------------------------------
    print()
    print("=" * 70)
    print("CANONICAL DRUG → DISEASE EXTRACTION COMPLETE")
    print("=" * 70)

    print(
        f"Canonical resolved edges: "
        f"{len(result):,}"
    )

    print(
        f"Unique ChEMBL drugs: "
        f"{unique_drugs:,}"
    )

    print(
        f"Unique MONDO diseases: "
        f"{unique_diseases:,}"
    )

    print(
        f"Unresolved edges: "
        f"{len(unresolved):,}"
    )

    print("\nRelation counts:")
    print(
        relation_counts.to_string(index=False)
    )

    print("\nDisease coverage:")
    print(
        disease_counts.to_string(index=False)
    )

    print()
    print(f"Saved: {OUT_FILE}")
    print(f"Saved: {UNRESOLVED_FILE}")
    print(f"Saved: {REPORT_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    main()