from pathlib import Path
import re

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]

PRIMEKG_FILE = (
    ROOT / "data/raw/primekg/primekg.csv"
)

DRUG_MAPPING_FILE = (
    ROOT / "data/raw/chembl/drugbank_to_chembl_validated.csv"
)

SELECTED_DRUGS_FILE = (
    ROOT / "data/interim/drug_disease_relations_canonical.csv"
)

HGNC_FILE = (
    ROOT / "data/raw/hgnc/hgnc_complete_set.txt"
)

OUTPUT_FILE = (
    ROOT / "data/interim/primekg_drug_protein_relations_canonical.csv"
)

REPORT_FILE = (
    ROOT / "reports/primekg_drug_protein_canonical_report.md"
)


def clean_string(value):
    if pd.isna(value):
        return ""

    return str(value).strip()


def normalize_entrez(value):
    value = clean_string(value)

    if not value:
        return ""

    # Keep numeric Entrez IDs only.
    value = re.sub(r"\.0$", "", value)

    return value


def split_uniprot(value):
    """
    HGNC may contain multiple UniProt accessions in one field.
    Keep every valid accession rather than arbitrarily selecting one.
    """

    value = clean_string(value)

    if not value:
        return []

    parts = re.split(
        r"[;,| ]+",
        value
    )

    return sorted(
        {
            p.strip()
            for p in parts
            if p.strip()
        }
    )


def join_unique(values):
    values = sorted(
        {
            clean_string(v)
            for v in values
            if clean_string(v)
        }
    )

    return " | ".join(values)


def main():

    print("=" * 70)
    print("REMEDY - PRIMEKG DRUG → PROTEIN CANONICAL LAYER")
    print("=" * 70)

    # ---------------------------------------------------------
    # 1. Load selected ChEMBL drugs
    # ---------------------------------------------------------

    selected = pd.read_csv(
        SELECTED_DRUGS_FILE,
        dtype=str,
    ).fillna("")

    selected_chembl = set(
        selected["chembl_id"]
        .astype(str)
        .str.strip()
    )

    selected_chembl.discard("")

    print(
        f"\nSelected ChEMBL drugs in current scope: "
        f"{len(selected_chembl):,}"
    )

    # ---------------------------------------------------------
    # 2. Load validated DrugBank → ChEMBL mapping
    # ---------------------------------------------------------

    drug_mapping = pd.read_csv(
        DRUG_MAPPING_FILE,
        dtype=str,
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

    # Only DrugBank records corresponding to our selected
    # ChEMBL drug scope.
    drug_mapping = drug_mapping[
        drug_mapping["chembl_id"]
        .isin(selected_chembl)
    ].copy()

    print(
        f"Relevant DrugBank mappings: "
        f"{len(drug_mapping):,}"
    )

    # ---------------------------------------------------------
    # 3. Load HGNC Entrez → UniProt mapping
    # ---------------------------------------------------------

    print(
        "\nLoading HGNC..."
    )

    hgnc = pd.read_csv(
        HGNC_FILE,
        sep="\t",
        dtype=str,
        low_memory=False,
    ).fillna("")

    required_hgnc = [
        "hgnc_id",
        "symbol",
        "name",
        "entrez_id",
        "uniprot_ids",
    ]

    missing_hgnc = [
        c
        for c in required_hgnc
        if c not in hgnc.columns
    ]

    if missing_hgnc:

        raise ValueError(
            "HGNC is missing required columns: "
            f"{missing_hgnc}"
        )

    hgnc = hgnc[
        required_hgnc
    ].copy()

    hgnc["entrez_id"] = (
        hgnc["entrez_id"]
        .apply(normalize_entrez)
    )

    hgnc = hgnc[
        hgnc["entrez_id"] != ""
    ].copy()

    # ---------------------------------------------------------
    # Expand all UniProt IDs
    # ---------------------------------------------------------

    hgnc["uniprot_accessions"] = (
        hgnc["uniprot_ids"]
        .apply(split_uniprot)
    )

    hgnc_expanded = (
        hgnc[
            [
                "entrez_id",
                "hgnc_id",
                "symbol",
                "name",
                "uniprot_accessions",
            ]
        ]
        .explode(
            "uniprot_accessions"
        )
        .rename(
            columns={
                "uniprot_accessions":
                    "uniprot_accession"
            }
        )
    )

    hgnc_expanded["uniprot_accession"] = (
        hgnc_expanded[
            "uniprot_accession"
        ]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    hgnc_expanded = hgnc_expanded[
        hgnc_expanded[
            "uniprot_accession"
        ] != ""
    ].copy()

    hgnc_expanded = (
        hgnc_expanded
        .drop_duplicates(
            subset=[
                "entrez_id",
                "uniprot_accession",
            ]
        )
    )

    print(
        f"HGNC Entrez → UniProt rows: "
        f"{len(hgnc_expanded):,}"
    )

    # ---------------------------------------------------------
    # 4. Load PrimeKG drug-protein rows
    # ---------------------------------------------------------

    print(
        "\nLoading PrimeKG..."
    )

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
        low_memory=False,
    )

    # ---------------------------------------------------------
    # Keep DrugBank → NCBI orientation only
    # ---------------------------------------------------------

    selected_mask = (
        (df["x_type"] == "drug")
        &
        (df["x_source"] == "DrugBank")
        &
        (df["y_type"] == "gene/protein")
        &
        (df["y_source"] == "NCBI")
        &
        (df["relation"] == "drug_protein")
    )

    edges = df[
        selected_mask
    ].copy()

    print(
        f"PrimeKG DrugBank → NCBI rows: "
        f"{len(edges):,}"
    )

    # ---------------------------------------------------------
    # 5. Normalize source IDs
    # ---------------------------------------------------------

    edges["drugbank_id"] = (
        edges["x_id"]
        .astype(str)
        .str.strip()
    )

    edges["drug_name_primekg"] = (
        edges["x_name"]
        .astype(str)
        .str.strip()
    )

    edges["entrez_id"] = (
        edges["y_id"]
        .astype(str)
        .apply(normalize_entrez)
    )

    edges["protein_gene_symbol_primekg"] = (
        edges["y_name"]
        .astype(str)
        .str.strip()
    )

    # ---------------------------------------------------------
    # 6. Restrict to selected drugs
    # ---------------------------------------------------------

    edges = edges[
        edges["drugbank_id"]
        .isin(
            set(
                drug_mapping["drugbank_id"]
            )
        )
    ].copy()

    print(
        f"Rows within selected drug scope: "
        f"{len(edges):,}"
    )

    # ---------------------------------------------------------
    # 7. Remove exact duplicate source rows
    # ---------------------------------------------------------

    before_dedup = len(edges)

    edges = edges.drop_duplicates(
        subset=[
            "drugbank_id",
            "entrez_id",
            "relation",
        ]
    ).copy()

    print(
        f"Exact duplicate source rows removed: "
        f"{before_dedup - len(edges):,}"
    )

    # ---------------------------------------------------------
    # 8. Map DrugBank → ChEMBL
    # ---------------------------------------------------------

    result = edges.merge(
        drug_mapping,
        on="drugbank_id",
        how="left",
        suffixes=(
            "",
            "_mapping",
        ),
    )

    # ---------------------------------------------------------
    # 9. Keep only successful ChEMBL mappings
    # ---------------------------------------------------------

    unresolved_drugs = result[
        result["chembl_id"]
        .isna()
        |
        (
            result["chembl_id"]
            .astype(str)
            .str.strip()
            == ""
        )
    ].copy()

    print(
        f"Rows with unresolved DrugBank → ChEMBL: "
        f"{len(unresolved_drugs):,}"
    )

    result = result[
        result["chembl_id"]
        .notna()
        &
        (
            result["chembl_id"]
            .astype(str)
            .str.strip()
            != ""
        )
    ].copy()

    # ---------------------------------------------------------
    # 10. Map Entrez → UniProt through HGNC
    # ---------------------------------------------------------

    result = result.merge(
        hgnc_expanded[
            [
                "entrez_id",
                "hgnc_id",
                "symbol",
                "name",
                "uniprot_accession",
            ]
        ],
        on="entrez_id",
        how="left",
    )

    # ---------------------------------------------------------
    # 11. Identify unresolved protein mappings
    # ---------------------------------------------------------

    unresolved_proteins = result[
        result["uniprot_accession"]
        .isna()
        |
        (
            result["uniprot_accession"]
            .astype(str)
            .str.strip()
            == ""
        )
    ].copy()

    print(
        f"Rows with unresolved Entrez → UniProt: "
        f"{len(unresolved_proteins):,}"
    )

    result = result[
        result["uniprot_accession"]
        .notna()
        &
        (
            result["uniprot_accession"]
            .astype(str)
            .str.strip()
            != ""
        )
    ].copy()

    # ---------------------------------------------------------
    # 12. Create canonical Drug → Protein relation
    # ---------------------------------------------------------

    canonical = (
        result[
            [
                "chembl_id",
                "chembl_pref_name",
                "drugbank_id",
                "drug_name",
                "entrez_id",
                "protein_gene_symbol_primekg",
                "hgnc_id",
                "symbol",
                "name",
                "uniprot_accession",
                "x_source",
                "y_source",
                "mapping_status",
                "mapping_reason",
            ]
        ]
        .rename(
            columns={
                "symbol": "hgnc_symbol",
                "name": "hgnc_name",
                "x_source": "drug_source",
                "y_source": "protein_source",
            }
        )
        .copy()
    )

    canonical["relation"] = (
        "targets"
    )

    # ---------------------------------------------------------
    # Canonical deduplication
    #
    # Multiple DrugBank / Entrez records can represent the
    # same ChEMBL → UniProt canonical edge.
    #
    # Keep one edge and preserve source identifiers.
    # ---------------------------------------------------------

    canonical = (
        canonical
        .groupby(
            [
                "chembl_id",
                "uniprot_accession",
            ],
            as_index=False,
        )
        .agg(
            chembl_pref_name=(
                "chembl_pref_name",
                join_unique,
            ),
            drugbank_id=(
                "drugbank_id",
                join_unique,
            ),
            drug_name=(
                "drug_name",
                join_unique,
            ),
            entrez_id=(
                "entrez_id",
                join_unique,
            ),
            protein_gene_symbol=(
                "protein_gene_symbol_primekg",
                join_unique,
            ),
            hgnc_id=(
                "hgnc_id",
                join_unique,
            ),
            hgnc_symbol=(
                "hgnc_symbol",
                join_unique,
            ),
            hgnc_name=(
                "hgnc_name",
                join_unique,
            ),
            drug_source=(
                "drug_source",
                join_unique,
            ),
            protein_source=(
                "protein_source",
                join_unique,
            ),
            mapping_status=(
                "mapping_status",
                join_unique,
            ),
            mapping_reason=(
                "mapping_reason",
                join_unique,
            ),
            relation=(
                "relation",
                "first",
            ),
        )
    )

    # ---------------------------------------------------------
    # 13. Statistics
    # ---------------------------------------------------------

    unique_drugs = canonical[
        "chembl_id"
    ].nunique()

    unique_proteins = canonical[
        "uniprot_accession"
    ].nunique()

    unique_edges = len(
        canonical
    )

    mapped_drugbank = canonical[
        "drugbank_id"
    ].nunique()

    mapped_entrez = canonical[
        "entrez_id"
    ].nunique()

    # ---------------------------------------------------------
    # 14. Save canonical layer
    # ---------------------------------------------------------

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    canonical = canonical[
        [
            "chembl_id",
            "chembl_pref_name",
            "drugbank_id",
            "drug_name",
            "uniprot_accession",
            "entrez_id",
            "protein_gene_symbol",
            "hgnc_id",
            "hgnc_symbol",
            "hgnc_name",
            "relation",
            "drug_source",
            "protein_source",
            "mapping_status",
            "mapping_reason",
        ]
    ].sort_values(
        [
            "chembl_id",
            "uniprot_accession",
        ]
    ).reset_index(
        drop=True
    )

    canonical.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    # ---------------------------------------------------------
    # 15. Report
    # ---------------------------------------------------------

    selected_drugbank_count = (
        drug_mapping[
            "drugbank_id"
        ].nunique()
    )

    drugbank_with_primekg = (
        edges[
            "drugbank_id"
        ].nunique()
    )

    entrez_with_primekg = (
        edges[
            "entrez_id"
        ].nunique()
    )

    report = f"""# PrimeKG Drug → Protein Canonical Report

## Source relation

PrimeKG relation:

`drug -- drug_protein --> gene/protein`

Only the DrugBank → NCBI orientation was retained.

## Scope

- Selected ChEMBL drugs: {len(selected_chembl):,}
- Relevant DrugBank mappings: {selected_drugbank_count:,}

## PrimeKG extraction

- DrugBank → NCBI rows: {len(edges):,}
- DrugBank IDs represented: {drugbank_with_primekg:,}
- Entrez IDs represented: {entrez_with_primekg:,}

## Canonical mapping

DrugBank → ChEMBL → Entrez → HGNC → UniProt

- Canonical Drug → Protein edges: {unique_edges:,}
- Unique ChEMBL drugs: {unique_drugs:,}
- Unique UniProt proteins: {unique_proteins:,}
- DrugBank IDs represented in canonical layer: {mapped_drugbank:,}
- Entrez IDs represented in canonical layer: {mapped_entrez:,}

## Relation

`drug -- targets --> protein`

## Provenance

PrimeKG source identifiers are retained:

- DrugBank ID
- Entrez ID
- HGNC ID
- PrimeKG source fields

These can later be attached to the canonical Drug → Protein
edge as PrimeKG evidence.

## Important policy

The relation is mapped to the generic `targets` edge.

No `inhibits` or `activates` relation is inferred from the
PrimeKG `drug_protein` relation.

## Outputs

`data/interim/primekg_drug_protein_relations_canonical.csv`

`reports/primekg_drug_protein_canonical_report.md`

Unresolved mappings are not silently invented or reassigned.
"""

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    # ---------------------------------------------------------
    # 16. Terminal output
    # ---------------------------------------------------------

    print()
    print("=" * 70)
    print("PRIMEKG DRUG → PROTEIN CANONICAL COMPLETE")
    print("=" * 70)

    print(
        f"Canonical Drug → Protein edges: "
        f"{unique_edges:,}"
    )

    print(
        f"Unique ChEMBL drugs: "
        f"{unique_drugs:,}"
    )

    print(
        f"Unique UniProt proteins: "
        f"{unique_proteins:,}"
    )

    print(
        f"DrugBank IDs represented: "
        f"{mapped_drugbank:,}"
    )

    print(
        f"Entrez IDs represented: "
        f"{mapped_entrez:,}"
    )

    print()
    print(
        f"Saved: {OUTPUT_FILE}"
    )

    print(
        f"Saved: {REPORT_FILE}"
    )

    print("=" * 70)


if __name__ == "__main__":
    main()