from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = (
    ROOT / "data/interim/chembl_reactome_protein_pathways.csv"
)

OUTPUT_FILE = (
    ROOT / "data/interim/protein_pathway_relations_canonical.csv"
)

REPORT_FILE = (
    ROOT / "reports/protein_pathway_canonical_report.md"
)


def main():

    print("=" * 70)
    print("REMEDY - CANONICAL UNIPROT → REACTOME PATHWAY RELATIONS")
    print("=" * 70)

    df = pd.read_csv(
        INPUT_FILE,
        dtype=str
    ).fillna("")

    print(
        f"Input mapping rows: {len(df):,}"
    )

    # ---------------------------------------------------------
    # Keep only valid UniProt + pathway records
    # ---------------------------------------------------------
    df = df[
        (df["uniprot_accession"] != "")
        &
        (df["reactome_pathway_id"] != "")
    ].copy()

    # ---------------------------------------------------------
    # Canonical protein → pathway relation
    #
    # Keep evidence_code because multiple Reactome records may
    # support the same protein-pathway relationship.
    # ---------------------------------------------------------
    canonical = df[
        [
            "uniprot_accession",
            "reactome_pathway_id",
            "pathway_name",
            "evidence_code",
        ]
    ].drop_duplicates()

    # ---------------------------------------------------------
    # Aggregate evidence per protein-pathway edge
    # ---------------------------------------------------------
    edges = (
        canonical.groupby(
            [
                "uniprot_accession",
                "reactome_pathway_id",
                "pathway_name",
            ],
            dropna=False
        )
        .agg(
            evidence_count=(
                "evidence_code",
                "nunique"
            ),
            evidence_codes=(
                "evidence_code",
                lambda s: "|".join(
                    sorted(
                        {
                            str(x).strip()
                            for x in s
                            if str(x).strip()
                        }
                    )
                )
            ),
        )
        .reset_index()
    )

    edges["relation"] = "participates_in"

    edges["source"] = "Reactome"

    edges = edges[
        [
            "uniprot_accession",
            "relation",
            "reactome_pathway_id",
            "pathway_name",
            "evidence_count",
            "evidence_codes",
            "source",
        ]
    ]

    edges = edges.sort_values(
        [
            "uniprot_accession",
            "reactome_pathway_id",
        ]
    ).reset_index(drop=True)

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------
    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    edges.to_csv(
        OUTPUT_FILE,
        index=False
    )

    # ---------------------------------------------------------
    # Statistics
    # ---------------------------------------------------------
    unique_proteins = edges[
        "uniprot_accession"
    ].nunique()

    unique_pathways = edges[
        "reactome_pathway_id"
    ].nunique()

    total_edges = len(edges)

    # ---------------------------------------------------------
    # Report
    # ---------------------------------------------------------
    report = f"""# Canonical Protein → Pathway Report

## Input

`data/interim/chembl_reactome_protein_pathways.csv`

## Results

- Canonical UniProt → Reactome edges: {total_edges:,}
- Unique UniProt proteins: {unique_proteins:,}
- Unique Reactome pathways: {unique_pathways:,}

## Relation

Canonical relation:

`protein -- participates_in --> pathway`

Protein identifier:
`UniProt accession`

Pathway identifier:
`Reactome pathway ID`

Source:
`Reactome`

## Evidence

Reactome evidence codes are preserved and aggregated
per protein-pathway edge.

## Output

`data/interim/protein_pathway_relations_canonical.csv`
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
    # Terminal output
    # ---------------------------------------------------------
    print()
    print("=" * 70)
    print("CANONICAL PROTEIN → PATHWAY COMPLETE")
    print("=" * 70)

    print(
        f"Canonical edges: {total_edges:,}"
    )

    print(
        f"Unique UniProt proteins: "
        f"{unique_proteins:,}"
    )

    print(
        f"Unique Reactome pathways: "
        f"{unique_pathways:,}"
    )

    print()
    print(f"Saved: {OUTPUT_FILE}")
    print(f"Saved: {REPORT_FILE}")

    print("=" * 70)


if __name__ == "__main__":
    main()