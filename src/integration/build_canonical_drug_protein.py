from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

RELATIONS_FILE = (
    ROOT / "data/interim/chembl_drug_protein_relations.csv"
)

MAPPING_FILE = (
    ROOT / "data/interim/chembl_target_uniprot_mapping_extended.csv"
)

OUT_FILE = (
    ROOT / "data/interim/drug_protein_relations_canonical.csv"
)

NONHUMAN_FILE = (
    ROOT / "reports/nonhuman_drug_protein_relations.csv"
)

UNRESOLVED_FILE = (
    ROOT / "reports/unresolved_drug_protein_edges.csv"
)

REPORT_FILE = (
    ROOT / "reports/drug_protein_canonical_report.md"
)


def main():

    print("=" * 70)
    print("REMEDY - CANONICAL DRUG → UNIPROT PROTEIN RELATIONS")
    print("=" * 70)

    # ---------------------------------------------------------
    # 1. Load ChEMBL drug-target summary
    # ---------------------------------------------------------
    relations = pd.read_csv(
        RELATIONS_FILE,
        dtype=str
    ).fillna("")

    print(
        f"Drug-target rows: {len(relations):,}"
    )

    # Convert counts back to numeric.
    for col in [
        "evidence_count",
        "human_evidence_count",
        "non_human_evidence_count",
        "unspecified_evidence_count",
    ]:
        relations[col] = pd.to_numeric(
            relations[col],
            errors="coerce"
        ).fillna(0).astype(int)

    # ---------------------------------------------------------
    # 2. Load extended target → UniProt mapping
    # ---------------------------------------------------------
    mapping = pd.read_csv(
        MAPPING_FILE,
        dtype=str
    ).fillna("")

    print(
        f"Target-UniProt mapping rows: "
        f"{len(mapping):,}"
    )

    # ---------------------------------------------------------
    # 3. Separate mapped human / non-human / unresolved
    # ---------------------------------------------------------
    human_mapping = mapping[
        (mapping["human_target"].str.lower() == "true")
        &
        (mapping["mapping_status"] == "mapped")
        &
        (mapping["uniprot_accession"] != "")
    ].copy()

    nonhuman_mapping = mapping[
        (mapping["human_target"].str.lower() != "true")
        &
        (mapping["mapping_status"] == "mapped")
        &
        (mapping["uniprot_accession"] != "")
    ].copy()

    unresolved_mapping = mapping[
        (
            (mapping["mapping_status"] != "mapped")
            |
            (mapping["uniprot_accession"] == "")
        )
    ].copy()

    print(
        f"Human target mapping rows: "
        f"{len(human_mapping):,}"
    )

    print(
        f"Non-human target mapping rows: "
        f"{len(nonhuman_mapping):,}"
    )

    # ---------------------------------------------------------
    # 4. Join all ChEMBL drug-target relations
    #    with target → UniProt mapping
    # ---------------------------------------------------------
    joined = relations.merge(
        mapping[
            [
                "target_id",
                "target_name",
                "target_organism",
                "uniprot_accession",
                "mapping_status",
                "human_target",
            ]
        ],
        on="target_id",
        how="left",
        suffixes=("", "_mapping")
    )

    # ---------------------------------------------------------
    # 5. Identify unresolved drug-target edges
    # ---------------------------------------------------------
    unresolved = joined[
        (
            joined["uniprot_accession"].isna()
            | (joined["uniprot_accession"] == "")
            | (joined["mapping_status"] != "mapped")
        )
    ].copy()

    unresolved = unresolved[
        [
            "chembl_id",
            "drugbank_id",
            "drug_name",
            "target_id",
            "target_name",
            "evidence_count",
            "human_evidence_count",
            "non_human_evidence_count",
            "unspecified_evidence_count",
            "target_organism",
            "mapping_status",
        ]
    ].drop_duplicates()

    # ---------------------------------------------------------
    # 6. Identify non-human target edges
    # ---------------------------------------------------------
    nonhuman = joined[
        (
            joined["mapping_status"] == "mapped"
        )
        &
        (
            joined["human_target"].str.lower() != "true"
        )
        &
        (
            joined["uniprot_accession"] != ""
        )
    ].copy()

    nonhuman = nonhuman[
        [
            "chembl_id",
            "drugbank_id",
            "drug_name",
            "target_id",
            "target_name",
            "target_organism",
            "uniprot_accession",
            "evidence_count",
            "human_evidence_count",
            "non_human_evidence_count",
            "unspecified_evidence_count",
        ]
    ].drop_duplicates()

    # ---------------------------------------------------------
    # 7. Keep human target relationships for core graph
    # ---------------------------------------------------------
    canonical = joined[
        (
            joined["mapping_status"] == "mapped"
        )
        &
        (
            joined["human_target"].str.lower() == "true"
        )
        &
        (
            joined["uniprot_accession"] != ""
        )
    ].copy()

    # ---------------------------------------------------------
    # 8. Canonical relation table
    # ---------------------------------------------------------
    canonical = canonical[
        [
            "chembl_id",
            "drugbank_id",
            "drug_name",
            "target_id",
            "target_name",
            "uniprot_accession",
            "target_organism",
            "evidence_count",
            "human_evidence_count",
            "non_human_evidence_count",
            "unspecified_evidence_count",
        ]
    ].drop_duplicates()

    canonical = canonical.sort_values(
        [
            "chembl_id",
            "uniprot_accession",
            "target_id",
        ]
    ).reset_index(drop=True)

    # ---------------------------------------------------------
    # 9. Save
    # ---------------------------------------------------------
    OUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    canonical.to_csv(
        OUT_FILE,
        index=False
    )

    NONHUMAN_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    nonhuman.to_csv(
        NONHUMAN_FILE,
        index=False
    )

    unresolved.to_csv(
        UNRESOLVED_FILE,
        index=False
    )

    # ---------------------------------------------------------
    # 10. Statistics
    # ---------------------------------------------------------
    unique_drugs = canonical["chembl_id"].nunique()

    unique_targets = canonical["target_id"].nunique()

    unique_uniprot = canonical[
        "uniprot_accession"
    ].nunique()

    canonical_rows = len(canonical)

    nonhuman_rows = len(nonhuman)

    unresolved_rows = len(unresolved)

    human_evidence = canonical[
        "human_evidence_count"
    ].sum()

    total_evidence = canonical[
        "evidence_count"
    ].sum()

    # ---------------------------------------------------------
    # 11. Report
    # ---------------------------------------------------------
    report = f"""# Canonical Drug → Protein Report

## Inputs

- ChEMBL drug-target summary:
  `data/interim/chembl_drug_protein_relations.csv`
- Extended target-UniProt mapping:
  `data/interim/chembl_target_uniprot_mapping_extended.csv`

## Core human graph

Only ChEMBL targets explicitly classified as
`Homo sapiens` are included in the canonical core
Drug → UniProt layer.

- Canonical rows: {canonical_rows:,}
- Unique ChEMBL drugs: {unique_drugs:,}
- Unique ChEMBL targets: {unique_targets:,}
- Unique UniProt accessions: {unique_uniprot:,}
- Human assay evidence represented: {human_evidence:,}
- Total evidence represented on these target edges: {total_evidence:,}

## Non-human evidence

Non-human target relationships are preserved separately
rather than used to create human biological links.

- Non-human relation rows: {nonhuman_rows:,}

## Unresolved target mappings

- Unresolved drug-target rows: {unresolved_rows:,}

## Important interpretation rule

The current source table supports a generic ChEMBL
`Drug → target` relation.

We do NOT infer `inhibits` or `activates` solely from the
activity type because doing so would require an explicit
direction/interaction interpretation.

Activity values remain in the ChEMBL evidence layer.

## Outputs

- `data/interim/drug_protein_relations_canonical.csv`
- `reports/nonhuman_drug_protein_relations.csv`
- `reports/unresolved_drug_protein_edges.csv`
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
    print("CANONICAL DRUG → PROTEIN COMPLETE")
    print("=" * 70)

    print(
        f"Canonical human rows: "
        f"{canonical_rows:,}"
    )

    print(
        f"Unique ChEMBL drugs: "
        f"{unique_drugs:,}"
    )

    print(
        f"Unique ChEMBL targets: "
        f"{unique_targets:,}"
    )

    print(
        f"Unique UniProt proteins: "
        f"{unique_uniprot:,}"
    )

    print(
        f"Non-human rows preserved: "
        f"{nonhuman_rows:,}"
    )

    print(
        f"Unresolved rows: "
        f"{unresolved_rows:,}"
    )

    print()
    print(f"Saved: {OUT_FILE}")
    print(f"Saved: {NONHUMAN_FILE}")
    print(f"Saved: {UNRESOLVED_FILE}")
    print(f"Saved: {REPORT_FILE}")

    print("=" * 70)


if __name__ == "__main__":
    main()