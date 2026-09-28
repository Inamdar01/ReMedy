from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

DISEASE_PROTEIN_FILE = (
    ROOT / "data/interim/primekg_disease_protein_relations.csv"
)

HGNC_FILE = (
    ROOT / "data/raw/hgnc/hgnc_complete_set.txt"
)

MAPPING_OUT = (
    ROOT / "data/interim/primekg_protein_uniprot_mapping.csv"
)

HARMONIZED_OUT = (
    ROOT / "data/interim/primekg_disease_protein_uniprot.csv"
)

UNRESOLVED_OUT = (
    ROOT / "reports/unresolved_primekg_proteins.csv"
)

REPORT_OUT = (
    ROOT / "reports/primekg_protein_uniprot_mapping_report.md"
)


def clean_entrez(value):
    """Normalize Entrez IDs such as 1234, 1234.0, or whitespace."""
    if pd.isna(value):
        return ""

    value = str(value).strip()

    if value.endswith(".0"):
        value = value[:-2]

    return value


def split_uniprot(value):
    """
    HGNC may contain multiple UniProt accessions.
    Preserve all supplied accessions instead of choosing one arbitrarily.
    """
    if pd.isna(value):
        return []

    value = str(value).strip()

    if not value or value.lower() == "nan":
        return []

    parts = re.split(r"[|,;\s]+", value)

    return sorted(
        {
            p.strip()
            for p in parts
            if p.strip()
            and p.strip().lower() != "nan"
        }
    )


def main():

    print("=" * 70)
    print("REMEDY - PRIMEKG ENTrez → HGNC → UNIPROT MAPPING")
    print("=" * 70)

    # ---------------------------------------------------------
    # 1. Load PrimeKG disease-protein relationships
    # ---------------------------------------------------------
    print("\nLoading PrimeKG disease-protein relations...")

    disease_protein = pd.read_csv(
        DISEASE_PROTEIN_FILE,
        dtype=str
    )

    disease_protein = disease_protein.fillna("")

    disease_protein["protein_entrez_id"] = (
        disease_protein["protein_entrez_id"]
        .apply(clean_entrez)
    )

    print(
        f"Disease-protein relation rows: "
        f"{len(disease_protein):,}"
    )

    unique_entrez = (
        disease_protein["protein_entrez_id"]
        .replace("", pd.NA)
        .dropna()
        .nunique()
    )

    print(f"Unique PrimeKG Entrez IDs: {unique_entrez:,}")

    # ---------------------------------------------------------
    # 2. Load HGNC
    # ---------------------------------------------------------
    print("\nLoading HGNC complete set...")

    hgnc_columns = [
        "hgnc_id",
        "symbol",
        "name",
        "locus_group",
        "locus_type",
        "status",
        "entrez_id",
        "uniprot_ids",
    ]

    hgnc = pd.read_csv(
        HGNC_FILE,
        sep="\t",
        usecols=hgnc_columns,
        dtype=str
    )

    hgnc = hgnc.fillna("")

    hgnc["entrez_id"] = (
        hgnc["entrez_id"]
        .apply(clean_entrez)
    )

    # ---------------------------------------------------------
    # 3. Restrict to approved protein-coding genes
    # ---------------------------------------------------------
    hgnc_protein = hgnc[
        (hgnc["status"].str.strip().str.lower() == "approved")
        &
        (
            (hgnc["locus_group"].str.strip().str.lower()
             == "protein-coding gene")
            |
            (hgnc["locus_type"].str.strip().str.lower()
             == "gene with protein product")
        )
        &
        (hgnc["entrez_id"] != "")
    ].copy()

    print(
        f"Approved protein-coding HGNC records: "
        f"{len(hgnc_protein):,}"
    )

    # ---------------------------------------------------------
    # 4. Keep only PrimeKG Entrez IDs
    # ---------------------------------------------------------
    primekg_ids = set(
        disease_protein["protein_entrez_id"]
        .loc[
            disease_protein["protein_entrez_id"] != ""
        ]
        .unique()
    )

    mapping = hgnc_protein[
        hgnc_protein["entrez_id"].isin(primekg_ids)
    ].copy()

    # ---------------------------------------------------------
    # 5. Detect ambiguous Entrez → HGNC mappings
    # ---------------------------------------------------------
    entrez_hgnc_counts = (
        mapping.groupby("entrez_id")["hgnc_id"]
        .nunique()
    )

    ambiguous_entrez = set(
        entrez_hgnc_counts[
            entrez_hgnc_counts > 1
        ].index
    )

    print(
        f"Entrez IDs with multiple HGNC mappings: "
        f"{len(ambiguous_entrez):,}"
    )

    # ---------------------------------------------------------
    # 6. Expand UniProt accessions
    # ---------------------------------------------------------
    mapping["uniprot_accession"] = (
        mapping["uniprot_ids"]
        .apply(split_uniprot)
    )

    mapping = mapping.explode(
        "uniprot_accession",
        ignore_index=True
    )

    mapping["uniprot_accession"] = (
        mapping["uniprot_accession"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    # ---------------------------------------------------------
    # 7. Assign mapping status
    # ---------------------------------------------------------
    mapping["mapping_status"] = "mapped"

    mapping.loc[
        mapping["uniprot_accession"] == "",
        "mapping_status"
    ] = "hgnc_no_uniprot"

    if ambiguous_entrez:
        mapping.loc[
            mapping["entrez_id"].isin(ambiguous_entrez),
            "mapping_status"
        ] = "ambiguous_hgnc"

    # ---------------------------------------------------------
    # 8. Create clean mapping table
    # ---------------------------------------------------------
    mapping_out = mapping[
        [
            "entrez_id",
            "hgnc_id",
            "symbol",
            "name",
            "status",
            "locus_group",
            "locus_type",
            "uniprot_accession",
            "mapping_status",
        ]
    ].copy()

    mapping_out = mapping_out.rename(
        columns={
            "entrez_id": "protein_entrez_id",
            "symbol": "hgnc_symbol",
            "name": "hgnc_name",
        }
    )

    mapping_out = mapping_out.drop_duplicates()

    # ---------------------------------------------------------
    # 9. Determine completely unresolved PrimeKG proteins
    # ---------------------------------------------------------
    mapped_entrez = set(
        mapping_out.loc[
            (mapping_out["mapping_status"] == "mapped")
            &
            (mapping_out["uniprot_accession"] != ""),
            "protein_entrez_id"
        ]
    )

    unresolved_ids = sorted(
        primekg_ids - mapped_entrez
    )

    unresolved = disease_protein[
        disease_protein["protein_entrez_id"].isin(
            unresolved_ids
        )
    ][
        [
            "protein_entrez_id",
            "protein_gene_symbol",
        ]
    ].drop_duplicates()

    # ---------------------------------------------------------
    # 10. Join disease relations to UniProt mapping
    # ---------------------------------------------------------
    harmonized = disease_protein.merge(
        mapping_out[
            [
                "protein_entrez_id",
                "hgnc_id",
                "hgnc_symbol",
                "hgnc_name",
                "uniprot_accession",
                "mapping_status",
            ]
        ],
        on="protein_entrez_id",
        how="left"
    )

    # Keep only successfully mapped UniProt relationships
    harmonized = harmonized[
        (harmonized["uniprot_accession"].notna())
        &
        (harmonized["uniprot_accession"] != "")
        &
        (harmonized["mapping_status"] == "mapped")
    ].copy()

    harmonized = harmonized.drop_duplicates()

    # ---------------------------------------------------------
    # 11. Save outputs
    # ---------------------------------------------------------
    MAPPING_OUT.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    mapping_out.to_csv(
        MAPPING_OUT,
        index=False
    )

    harmonized.to_csv(
        HARMONIZED_OUT,
        index=False
    )

    UNRESOLVED_OUT.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    unresolved.to_csv(
        UNRESOLVED_OUT,
        index=False
    )

    # ---------------------------------------------------------
    # 12. Statistics
    # ---------------------------------------------------------
    mapped_unique_entrez = (
        harmonized["protein_entrez_id"]
        .nunique()
    )

    unique_hgnc = (
        harmonized["hgnc_id"]
        .nunique()
    )

    unique_uniprot = (
        harmonized["uniprot_accession"]
        .nunique()
    )

    unmapped_count = len(unresolved_ids)

    mapping_rate = (
        mapped_unique_entrez / unique_entrez * 100
        if unique_entrez
        else 0
    )

    # ---------------------------------------------------------
    # 13. Disease coverage after harmonization
    # ---------------------------------------------------------
    disease_coverage = (
        harmonized.groupby("disease_name")
        .agg(
            unique_entrez=(
                "protein_entrez_id",
                "nunique"
            ),
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
        .sort_values("disease_name")
    )

    # ---------------------------------------------------------
    # 14. Report
    # ---------------------------------------------------------
    report = f"""# PrimeKG Protein → HGNC → UniProt Mapping Report

## Input

- PrimeKG disease-protein relations:
  `data/interim/primekg_disease_protein_relations.csv`
- HGNC complete set:
  `data/raw/hgnc/hgnc_complete_set.txt`

## Mapping strategy

PrimeKG protein identifiers are treated as NCBI/Entrez Gene IDs.

The mapping uses exact Entrez ID matching against HGNC.

Only approved protein-coding HGNC records are retained.

All UniProt accessions supplied by HGNC are preserved.

No identifier is guessed from gene name alone.

## Results

- PrimeKG unique Entrez IDs: {unique_entrez:,}
- Entrez IDs mapped to UniProt: {mapped_unique_entrez:,}
- HGNC records represented: {unique_hgnc:,}
- Unique UniProt accessions: {unique_uniprot:,}
- Unresolved Entrez IDs: {unmapped_count:,}
- Entrez → UniProt coverage: {mapping_rate:.2f}%

## Ambiguity

- Entrez IDs with multiple HGNC mappings:
  {len(ambiguous_entrez):,}

## Harmonized disease-protein layer

- Original disease-protein rows:
  {len(disease_protein):,}
- Harmonized disease-UniProt rows:
  {len(harmonized):,}

## Disease coverage

{disease_coverage.to_markdown(index=False)}

## Outputs

- `data/interim/primekg_protein_uniprot_mapping.csv`
- `data/interim/primekg_disease_protein_uniprot.csv`
- `reports/unresolved_primekg_proteins.csv`

## Important limitation

The original PrimeKG relation represents a combined
gene/protein concept. This mapping provides identifier
harmonization to HGNC and UniProt; it does not independently
change the biological meaning or strength of the PrimeKG
association.
"""

    REPORT_OUT.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    REPORT_OUT.write_text(
        report,
        encoding="utf-8"
    )

    # ---------------------------------------------------------
    # 15. Terminal output
    # ---------------------------------------------------------
    print()
    print("=" * 70)
    print("PRIMEKG PROTEIN HARMONIZATION COMPLETE")
    print("=" * 70)

    print(
        f"PrimeKG unique Entrez IDs: "
        f"{unique_entrez:,}"
    )

    print(
        f"Mapped Entrez IDs: "
        f"{mapped_unique_entrez:,}"
    )

    print(
        f"Mapping coverage: "
        f"{mapping_rate:.2f}%"
    )

    print(
        f"Unique HGNC genes: "
        f"{unique_hgnc:,}"
    )

    print(
        f"Unique UniProt accessions: "
        f"{unique_uniprot:,}"
    )

    print(
        f"Unresolved Entrez IDs: "
        f"{unmapped_count:,}"
    )

    print(
        f"Ambiguous Entrez → HGNC IDs: "
        f"{len(ambiguous_entrez):,}"
    )

    print(
        f"Harmonized disease-UniProt rows: "
        f"{len(harmonized):,}"
    )

    print("\nDisease coverage:")
    print(
        disease_coverage.to_string(
            index=False
        )
    )

    print()
    print(f"Saved: {MAPPING_OUT}")
    print(f"Saved: {HARMONIZED_OUT}")
    print(f"Saved: {UNRESOLVED_OUT}")
    print(f"Saved: {REPORT_OUT}")

    print("=" * 70)


if __name__ == "__main__":
    main()