
from pathlib import Path
import pandas as pd

# Project paths
ROOT = Path(__file__).resolve().parents[2]
INPUT_PATH = (
    ROOT
    / "reports"
    / "model_evaluation"
    / "top20_candidate_direct_evidence_report.csv"
)
OUTPUT_PATH = (
    ROOT
    / "reports"
    / "model_evaluation"
    / "drug_target_evidence_summary.csv"
)


def join_unique(values):
    """Join unique, non-empty values into a readable string."""
    cleaned = values.dropna().astype(str).str.strip()
    cleaned = cleaned[cleaned != ""]
    return ", ".join(sorted(set(cleaned)))


def main():
    if not INPUT_PATH.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_PATH}"
        )

    data = pd.read_csv(INPUT_PATH)

    required_columns = [
        "drug_name",
        "target_name",
        "target_organism",
        "disease_name",
        "activity_id",
        "activity_type",
        "activity_value",
        "publication",
    ]

    missing = set(required_columns) - set(data.columns)

    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}"
        )

    # Group by drug, target, and target organism.
    grouped = data.groupby(
        ["drug_name", "target_name", "target_organism"],
        dropna=False,
    )

    summary = grouped.agg(
        candidate_diseases=("disease_name", join_unique),
        assay_records=("activity_id", "nunique"),
        assay_types=("activity_type", join_unique),
        activity_values=("activity_value", join_unique),
        publications=("publication", join_unique),
    ).reset_index()

    # Sort the report for readability.
    summary = summary.sort_values(
        ["drug_name", "target_name"]
    )

    # Save the report.
    OUTPUT_PATH.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    summary.to_csv(OUTPUT_PATH, index=False)

    # Print results.
    print("=" * 70)
    print("DRUG-TARGET EVIDENCE SUMMARY")
    print("=" * 70)

    print("Input records:", len(data))
    print("Unique drug-target groups:", len(summary))
    print("Unique drugs:", summary["drug_name"].nunique())
    print("Unique targets:", summary["target_name"].nunique())

    print("\nSummary:")
    print(summary.to_string(index=False))

    print("\nSaved report:")
    print(OUTPUT_PATH)

    print(
        "\nNote: assay evidence supports drug-target relationships; "
        "it does not independently establish treatment effectiveness "
        "for the listed diseases."
    )


if __name__ == "__main__":
    main()
