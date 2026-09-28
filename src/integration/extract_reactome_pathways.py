import pandas as pd
from pathlib import Path

# ============================================================
# FILES
# ============================================================

MAPPING_FILE = (
    "data/interim/chembl_target_uniprot_mapping_extended.csv"
)

REACTOME_MAPPING_FILE = (
    "data/raw/reactome/UniProt2Reactome_All_Levels.txt"
)

REACTOME_PATHWAYS_FILE = (
    "data/raw/reactome/ReactomePathways.txt"
)

OUTPUT_FILE = (
    "data/interim/chembl_reactome_protein_pathways.csv"
)

REPORT_FILE = (
    "reports/reactome_pathway_mapping_report.md"
)

Path("data/interim").mkdir(
    parents=True,
    exist_ok=True
)

Path("reports").mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# LOAD CHemBL → UNIPROT
# ============================================================

print("=" * 70)
print("REMEDY - CHEMBL / UNIPROT → REACTOME PATHWAY MAPPING")
print("=" * 70)

mapping = pd.read_csv(
    MAPPING_FILE,
    low_memory=False
)

mapping = mapping[
    mapping["uniprot_accession"].notna()
].copy()

mapping["uniprot_accession"] = (
    mapping["uniprot_accession"]
    .astype(str)
    .str.strip()
)

# Human targets only for the primary Reactome layer
human_mapping = mapping[
    mapping["human_target"] == True
].copy()

print(
    f"ChEMBL mappings: {len(mapping)}"
)

print(
    f"Human target mappings: {len(human_mapping)}"
)

print(
    f"Unique human UniProt accessions: "
    f"{human_mapping['uniprot_accession'].nunique()}"
)


# ============================================================
# LOAD REACTOME UNIPROT MAPPING
# ============================================================

print()
print(
    "Loading Reactome UniProt mapping..."
)

reactome_map = pd.read_csv(
    REACTOME_MAPPING_FILE,
    sep="\t",
    header=None,
    names=[
        "uniprot_accession",
        "reactome_pathway_id",
        "pathway_url",
        "pathway_name",
        "evidence_code",
        "species"
    ],
    dtype=str
)

print(
    f"Reactome mapping rows: "
    f"{len(reactome_map)}"
)

# ============================================================
# KEEP HUMAN REACTOME ENTRIES
# ============================================================

reactome_human = reactome_map[
    reactome_map["species"]
    == "Homo sapiens"
].copy()

print(
    f"Human Reactome rows: "
    f"{len(reactome_human)}"
)

# ============================================================
# NORMALIZE ACCESSIONS
# ============================================================

reactome_human["uniprot_accession"] = (
    reactome_human["uniprot_accession"]
    .astype(str)
    .str.strip()
)


# ============================================================
# JOIN
# ============================================================

print()
print("Joining UniProt accessions...")

joined = human_mapping.merge(
    reactome_human,
    on="uniprot_accession",
    how="left"
)

# Only keep actual pathway matches
matched = joined[
    joined["reactome_pathway_id"].notna()
].copy()

print(
    f"Matched mapping rows: "
    f"{len(matched)}"
)

print(
    f"Unique ChEMBL targets matched: "
    f"{matched['target_id'].nunique()}"
)

print(
    f"Unique UniProt proteins matched: "
    f"{matched['uniprot_accession'].nunique()}"
)

print(
    f"Unique Reactome pathways: "
    f"{matched['reactome_pathway_id'].nunique()}"
)


# ============================================================
# REMOVE DUPLICATES
# ============================================================

matched = matched.drop_duplicates(
    subset=[
        "target_id",
        "uniprot_accession",
        "reactome_pathway_id"
    ]
)


# ============================================================
# SAVE
# ============================================================

matched.to_csv(
    OUTPUT_FILE,
    index=False
)


# ============================================================
# UNMATCHED PROTEINS
# ============================================================

all_human_targets = set(
    human_mapping["target_id"]
)

matched_targets = set(
    matched["target_id"]
)

unmatched_targets = (
    all_human_targets
    - matched_targets
)

unmatched_df = human_mapping[
    human_mapping["target_id"]
    .isin(unmatched_targets)
][
    [
        "target_id",
        "target_name",
        "uniprot_accession",
        "target_organism"
    ]
].drop_duplicates()

UNMATCHED_FILE = (
    "reports/"
    "reactome_unmatched_human_targets.csv"
)

unmatched_df.to_csv(
    UNMATCHED_FILE,
    index=False
)


# ============================================================
# REPORT
# ============================================================

report = f"""# Reactome Protein → Pathway Mapping Report

## Input

- ChEMBL → UniProt mapping:
  `{MAPPING_FILE}`
- Reactome bulk mapping:
  `{REACTOME_MAPPING_FILE}`

## Results

- ChEMBL target mappings: {len(mapping)}
- Human target mappings: {len(human_mapping)}
- Unique human UniProt accessions:
  {human_mapping['uniprot_accession'].nunique()}

- Human target-pathway mapping rows:
  {len(matched)}

- Unique ChEMBL targets matched:
  {matched['target_id'].nunique()}

- Unique UniProt proteins matched:
  {matched['uniprot_accession'].nunique()}

- Unique Reactome pathways:
  {matched['reactome_pathway_id'].nunique()}

- Human ChEMBL targets without a Reactome pathway:
  {len(unmatched_targets)}

## Filtering rule

Only Reactome entries with:

`species = Homo sapiens`

were used for the primary ReMedy human pathway layer.

## Output

`{OUTPUT_FILE}`

Unmatched targets:

`{UNMATCHED_FILE}`
"""

Path(REPORT_FILE).write_text(
    report,
    encoding="utf-8"
)


# ============================================================
# FINAL
# ============================================================

print()
print("=" * 70)
print("REACTOME MAPPING COMPLETE")
print("=" * 70)

print(
    f"Protein-pathway rows: "
    f"{len(matched)}"
)

print(
    f"Unique targets: "
    f"{matched['target_id'].nunique()}"
)

print(
    f"Unique UniProt: "
    f"{matched['uniprot_accession'].nunique()}"
)

print(
    f"Unique pathways: "
    f"{matched['reactome_pathway_id'].nunique()}"
)

print(
    f"Unmatched human targets: "
    f"{len(unmatched_targets)}"
)

print()
print(f"Saved: {OUTPUT_FILE}")
print(f"Saved: {REPORT_FILE}")
print(f"Saved: {UNMATCHED_FILE}")

print("=" * 70)