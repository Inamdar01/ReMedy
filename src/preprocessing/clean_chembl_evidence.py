import pandas as pd
from pathlib import Path

INPUT = "data/interim/chembl_target_evidence_combined.csv"

OUTPUT_DIR = Path("data/interim")
OUTPUT = OUTPUT_DIR / "chembl_evidence_clean.csv"
RELATIONS = OUTPUT_DIR / "chembl_drug_protein_relations.csv"
REPORT = Path("reports/chembl_cleaning_report.md")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
REPORT.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    INPUT,
    low_memory=False
)

print("=" * 70)
print("REMEDY - CHEMBL EVIDENCE CLEANING")
print("=" * 70)

print(f"Input rows: {len(df)}")

# ------------------------------------------------------------
# 1. Keep only SINGLE PROTEIN evidence
# ------------------------------------------------------------

df = df[
    df["target_type"].astype(str).str.upper()
    == "SINGLE PROTEIN"
].copy()

# ------------------------------------------------------------
# 2. Require target and measured value
# ------------------------------------------------------------

df = df[
    df["target_id"].notna()
    & df["standard_value"].notna()
].copy()

# ------------------------------------------------------------
# 3. Normalize IDs / strings
# ------------------------------------------------------------

for col in [
    "drugbank_id",
    "drug_name",
    "chembl_id",
    "activity_id",
    "target_id",
    "target_name",
    "standard_type",
    "standard_units",
    "standard_relation",
    "assay_id",
    "assay_type",
    "assay_organism",
    "target_organism",
    "document_id",
    "publication",
    "doi",
    "pmid",
    "source"
]:
    if col in df.columns:
        df[col] = df[col].where(
            df[col].notna(),
            None
        )

# ------------------------------------------------------------
# 4. Standard experimental organism classification
# ------------------------------------------------------------
# IMPORTANT:
# Use assay_organism only.
# Do NOT use target_organism to classify the experiment.

def classify_assay_organism(value):

    if pd.isna(value) or str(value).strip() == "":
        return "organism_unspecified"

    value = str(value).strip()

    if value.lower() == "homo sapiens":
        return "human"

    return "non_human"


df["assay_organism_class"] = (
    df["assay_organism"]
    .apply(classify_assay_organism)
)

# ------------------------------------------------------------
# 5. Preserve target organism separately
# ------------------------------------------------------------

df["target_organism"] = df[
    "target_organism"
].where(
    df["target_organism"].notna(),
    None
)

# ------------------------------------------------------------
# 6. Flag suspicious values
# ------------------------------------------------------------

df["value_quality"] = df[
    "value_quality"
].fillna("usable")

# ------------------------------------------------------------
# 7. Remove exact duplicate activities
# ------------------------------------------------------------

before = len(df)

df = df.drop_duplicates(
    subset=["activity_id"]
)

duplicate_rows_removed = (
    before - len(df)
)

# ------------------------------------------------------------
# 8. Create cleaner drug-protein relation table
# ------------------------------------------------------------

relation_df = (
    df.groupby(
        [
            "chembl_id",
            "drugbank_id",
            "drug_name",
            "target_id",
            "target_name"
        ],
        dropna=False
    )
    .agg(
        evidence_count=("activity_id", "nunique"),
        human_evidence_count=(
            "assay_organism_class",
            lambda s: (s == "human").sum()
        ),
        non_human_evidence_count=(
            "assay_organism_class",
            lambda s: (s == "non_human").sum()
        ),
        unspecified_evidence_count=(
            "assay_organism_class",
            lambda s: (
                s == "organism_unspecified"
            ).sum()
        )
    )
    .reset_index()
)

# ------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------

df.to_csv(
    OUTPUT,
    index=False
)

relation_df.to_csv(
    RELATIONS,
    index=False
)

# ------------------------------------------------------------
# 10. Report
# ------------------------------------------------------------

human = (
    df["assay_organism_class"]
    == "human"
).sum()

non_human = (
    df["assay_organism_class"]
    == "non_human"
).sum()

unknown = (
    df["assay_organism_class"]
    == "organism_unspecified"
).sum()

flagged = (
    df["value_quality"]
    == "flagged"
).sum()

usable = (
    df["value_quality"]
    == "usable"
).sum()

unique_drugs = df["chembl_id"].nunique()
unique_targets = df["target_id"].nunique()

report = f"""# ChEMBL Evidence Cleaning Report

## Input

Source:
`{INPUT}`

## Final cleaned evidence

- Evidence records: {len(df)}
- Unique drugs with evidence: {unique_drugs}
- Unique protein targets: {unique_targets}
- Exact duplicate activity rows removed: {duplicate_rows_removed}

## Experimental organism classification

Classification uses `assay_organism`.

- Human: {human}
- Non-human: {non_human}
- Organism unspecified: {unknown}

`target_organism` is retained separately and is not used to classify the experimental assay.

## Value quality

- Usable: {usable}
- Flagged: {flagged}

## Activity types

"""

for activity_type, count in (
    df["standard_type"]
    .value_counts()
    .items()
):
    report += f"- {activity_type}: {count}\n"

report += """
## Important interpretation rule

ChEMBL activity values are preserved with their original
activity type, value, units, and relation. Different activity
types are not collapsed into one numerical score.

Flagged values are retained and identified through the
`value_quality` field.
"""

REPORT.write_text(
    report,
    encoding="utf-8"
)

print()
print("=" * 70)
print("CLEANING COMPLETE")
print("=" * 70)

print(f"Clean evidence rows : {len(df)}")
print(f"Unique drugs        : {unique_drugs}")
print(f"Unique targets      : {unique_targets}")
print(f"Human assays       : {human}")
print(f"Non-human assays   : {non_human}")
print(f"Unspecified assays : {unknown}")
print(f"Usable values      : {usable}")
print(f"Flagged values     : {flagged}")

print()
print(f"Saved: {OUTPUT}")
print(f"Saved: {RELATIONS}")
print(f"Saved: {REPORT}")
print("=" * 70)