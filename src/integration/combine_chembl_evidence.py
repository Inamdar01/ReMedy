from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

ORIGINAL_FILE = (
    ROOT / "data/raw/chembl/evidence/chembl_target_evidence_final.csv"
)

SIX_FILE = (
    ROOT / "data/raw/chembl/evidence/chembl_target_evidence_six.csv"
)

OUT_FILE = (
    ROOT / "data/interim/chembl_target_evidence_combined.csv"
)

REPORT_FILE = (
    ROOT / "reports/chembl_evidence_combination_report.md"
)


def main():

    print("=" * 70)
    print("REMEDY - COMBINE CHEMBL EVIDENCE EXTRACTS")
    print("=" * 70)

    # ---------------------------------------------------------
    # Load both extracts
    # ---------------------------------------------------------
    print("\nLoading original evidence...")
    original = pd.read_csv(
        ORIGINAL_FILE,
        dtype=str
    ).fillna("")

    print(
        f"Original evidence rows: "
        f"{len(original):,}"
    )

    print("\nLoading targeted six-drug evidence...")
    six = pd.read_csv(
        SIX_FILE,
        dtype=str
    ).fillna("")

    print(
        f"Targeted evidence rows: "
        f"{len(six):,}"
    )

    # ---------------------------------------------------------
    # Schema check
    # ---------------------------------------------------------
    if list(original.columns) != list(six.columns):
        missing = sorted(
            set(original.columns) - set(six.columns)
        )

        extra = sorted(
            set(six.columns) - set(original.columns)
        )

        raise ValueError(
            f"Schema mismatch.\n"
            f"Missing from six: {missing}\n"
            f"Extra in six: {extra}"
        )

    print("\nSchema: identical")

    # ---------------------------------------------------------
    # Check activity IDs before combining
    # ---------------------------------------------------------
    original_activity_ids = set(
        original["activity_id"]
        .loc[original["activity_id"] != ""]
    )

    six_activity_ids = set(
        six["activity_id"]
        .loc[six["activity_id"] != ""]
    )

    overlap = (
        original_activity_ids
        & six_activity_ids
    )

    print(
        f"Activity IDs overlapping between extracts: "
        f"{len(overlap):,}"
    )

    if overlap:
        print("\nOverlapping activity IDs:")
        print("\n".join(sorted(overlap)[:20]))

        raise ValueError(
            "Activity IDs overlap. Do not combine until resolved."
        )

    # ---------------------------------------------------------
    # Combine
    # ---------------------------------------------------------
    combined = pd.concat(
        [original, six],
        ignore_index=True
    )

    # Only remove completely identical rows.
    before = len(combined)

    combined = combined.drop_duplicates()

    exact_duplicates_removed = (
        before - len(combined)
    )

    # ---------------------------------------------------------
    # Validate unique activity IDs
    # ---------------------------------------------------------
    duplicate_activity_ids = (
        combined["activity_id"]
        .loc[combined["activity_id"] != ""]
        .duplicated()
        .sum()
    )

    if duplicate_activity_ids:
        raise ValueError(
            f"Duplicate activity IDs after combination: "
            f"{duplicate_activity_ids}"
        )

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------
    OUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    combined.to_csv(
        OUT_FILE,
        index=False
    )

    # ---------------------------------------------------------
    # Statistics
    # ---------------------------------------------------------
    unique_drugs = combined["chembl_id"].nunique()

    unique_targets = combined["target_id"].nunique()

    human = (
        combined["organism_class"] == "human"
    ).sum()

    non_human = (
        combined["organism_class"] == "non_human"
    ).sum()

    unspecified = (
        combined["organism_class"]
        == "organism_unspecified"
    ).sum()

    # ---------------------------------------------------------
    # Drug coverage
    # ---------------------------------------------------------
    drug_counts = (
        combined.groupby(
            ["chembl_id", "drug_name"]
        )
        .size()
        .sort_values(
            ascending=False
        )
    )

    # ---------------------------------------------------------
    # Report
    # ---------------------------------------------------------
    report = f"""# ChEMBL Evidence Combination Report

## Inputs

- Original:
  `data/raw/chembl/evidence/chembl_target_evidence_final.csv`
- Targeted six-drug extract:
  `data/raw/chembl/evidence/chembl_target_evidence_six.csv`

## Results

- Original evidence rows: {len(original):,}
- Targeted evidence rows: {len(six):,}
- Exact duplicate rows removed: {exact_duplicates_removed:,}
- Combined evidence rows: {len(combined):,}
- Unique ChEMBL drugs with evidence: {unique_drugs:,}
- Unique ChEMBL targets: {unique_targets:,}

## Assay organism classification

- Human: {human:,}
- Non-human: {non_human:,}
- Unspecified: {unspecified:,}

## Activity IDs

- Overlapping activity IDs before combination: {len(overlap):,}
- Duplicate activity IDs after combination: {duplicate_activity_ids:,}

## Evidence by drug

{drug_counts.to_string()}

## Output

`data/interim/chembl_target_evidence_combined.csv`
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
    print("CHEMBL EVIDENCE COMBINATION COMPLETE")
    print("=" * 70)

    print(
        f"Combined evidence rows: "
        f"{len(combined):,}"
    )

    print(
        f"Unique drugs with evidence: "
        f"{unique_drugs:,}"
    )

    print(
        f"Unique targets: "
        f"{unique_targets:,}"
    )

    print(
        f"Human evidence: "
        f"{human:,}"
    )

    print(
        f"Non-human evidence: "
        f"{non_human:,}"
    )

    print(
        f"Unspecified evidence: "
        f"{unspecified:,}"
    )

    print(
        f"Exact duplicates removed: "
        f"{exact_duplicates_removed:,}"
    )

    print()
    print(f"Saved: {OUT_FILE}")
    print(f"Saved: {REPORT_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    main()
    