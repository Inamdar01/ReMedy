from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

PATHWAY_FILE = ROOT / "data/raw/reactome/ReactomePathways.txt"
RELATION_FILE = ROOT / "data/raw/reactome/ReactomePathwaysRelation.txt"
PROTEIN_PATHWAY_FILE = ROOT / "data/interim/chembl_reactome_protein_pathways.csv"

OUT_FILE = ROOT / "data/interim/reactome_pathway_hierarchy.csv"
REPORT_FILE = ROOT / "reports/reactome_hierarchy_report.md"


def main():
    print("=" * 70)
    print("REMEDY - REACTOME PATHWAY HIERARCHY INTEGRATION")
    print("=" * 70)

    # ---------------------------------------------------------
    # 1. Load pathways already connected to our proteins
    # ---------------------------------------------------------
    protein_pathways = pd.read_csv(PROTEIN_PATHWAY_FILE)

    print(f"Protein-pathway rows: {len(protein_pathways):,}")

    pathway_ids = set(
        protein_pathways["reactome_pathway_id"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    print(f"Unique pathways from protein mapping: {len(pathway_ids):,}")

    # ---------------------------------------------------------
    # 2. Load official Reactome pathway definitions
    # ---------------------------------------------------------
    pathways = pd.read_csv(
        PATHWAY_FILE,
        sep="\t",
        header=None,
        names=["pathway_id", "pathway_name", "species"],
        dtype=str
    )

    pathways = pathways.fillna("")

    # Keep only the pathways that appear in our protein mapping
    pathways_used = pathways[
        pathways["pathway_id"].isin(pathway_ids)
    ].copy()

    print(f"Pathways found in Reactome definitions: {len(pathways_used):,}")

    # ---------------------------------------------------------
    # 3. Load pathway hierarchy
    # ---------------------------------------------------------
    relations = pd.read_csv(
        RELATION_FILE,
        sep="\t",
        header=None,
        names=["parent_pathway_id", "child_pathway_id"],
        dtype=str
    )

    relations = relations.fillna("")

    # Only keep relationships relevant to our current pathways.
    relevant_relations = relations[
        relations["parent_pathway_id"].isin(pathway_ids)
        | relations["child_pathway_id"].isin(pathway_ids)
    ].copy()

    # ---------------------------------------------------------
    # 4. Add pathway names
    # ---------------------------------------------------------
    name_map = pathways.drop_duplicates("pathway_id").set_index("pathway_id")

    relevant_relations["parent_pathway_name"] = (
        relevant_relations["parent_pathway_id"]
        .map(name_map["pathway_name"])
    )

    relevant_relations["child_pathway_name"] = (
        relevant_relations["child_pathway_id"]
        .map(name_map["pathway_name"])
    )

    # ---------------------------------------------------------
    # 5. Clean and deduplicate
    # ---------------------------------------------------------
    relevant_relations = relevant_relations[
        [
            "parent_pathway_id",
            "parent_pathway_name",
            "child_pathway_id",
            "child_pathway_name",
        ]
    ].drop_duplicates()

    relevant_relations = relevant_relations.sort_values(
        ["parent_pathway_id", "child_pathway_id"]
    )

    # ---------------------------------------------------------
    # 6. Save
    # ---------------------------------------------------------
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    relevant_relations.to_csv(OUT_FILE, index=False)

    # ---------------------------------------------------------
    # 7. Report
    # ---------------------------------------------------------
    total_relations = len(relations)
    relevant_count = len(relevant_relations)

    report = f"""# Reactome Pathway Hierarchy Integration Report

## Input

- Protein-pathway mapping: `data/interim/chembl_reactome_protein_pathways.csv`
- Reactome pathway definitions: `data/raw/reactome/ReactomePathways.txt`
- Reactome hierarchy: `data/raw/reactome/ReactomePathwaysRelation.txt`

## Results

- Unique pathways connected to ReMedy proteins: {len(pathway_ids):,}
- Pathways found in Reactome definitions: {len(pathways_used):,}
- Total Reactome hierarchy relations: {total_relations:,}
- Relevant hierarchy relations: {relevant_count:,}
- Unique parent pathways: {relevant_relations["parent_pathway_id"].nunique():,}
- Unique child pathways: {relevant_relations["child_pathway_id"].nunique():,}

## Output

`data/interim/reactome_pathway_hierarchy.csv`

The hierarchy is retained as explicit parent-to-child pathway relationships.
No disease relationships are inferred from the pathway hierarchy.
"""

    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)
    REPORT_FILE.write_text(report, encoding="utf-8")

    print()
    print("=" * 70)
    print("REACTOME HIERARCHY INTEGRATION COMPLETE")
    print("=" * 70)
    print(f"Relevant hierarchy relations: {relevant_count:,}")
    print(f"Unique parent pathways: {relevant_relations['parent_pathway_id'].nunique():,}")
    print(f"Unique child pathways: {relevant_relations['child_pathway_id'].nunique():,}")
    print()
    print(f"Saved: {OUT_FILE}")
    print(f"Saved: {REPORT_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    main()