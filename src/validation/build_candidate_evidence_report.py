
from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INPUT_PATH = (
    ROOT
    / "reports"
    / "model_evaluation"
    / "top20_candidate_edge_evidence.csv"
)
OUTPUT_DIR = ROOT / "reports" / "model_evaluation"
OUTPUT_PATH = OUTPUT_DIR / "top20_candidate_direct_evidence_report.csv"


def extract_context_value(context, key):
    """Extract a key=value field from the context string."""
    if pd.isna(context):
        return ""

    pattern = rf"(?:^|;\s*){re.escape(key)}=([^;]*)"
    match = re.search(pattern, str(context))

    return match.group(1).strip() if match else ""


def main():
    if not INPUT_PATH.exists():
        raise FileNotFoundError(
            f"Input evidence file not found: {INPUT_PATH}"
        )

    data = pd.read_csv(INPUT_PATH)

    required = {
        "rank",
        "drug_id",
        "drug_name",
        "disease_id",
        "disease_name",
        "rgcn_score",
        "edge_id",
        "relation",
        "source_name",
        "target_name",
        "evidence_source",
        "publication",
        "assay",
        "organism",
        "activity_value",
        "context",
    }

    missing = required - set(data.columns)
    if missing:
        raise ValueError(
            f"Missing columns: {sorted(missing)}"
        )

    # This report focuses on direct drug-target assay evidence.
    direct = data[
        (data["relation"] == "targets")
        & (data["evidence_source"] == "ChEMBL")
    ].copy()

    if direct.empty:
        raise ValueError(
            "No ChEMBL drug-target records were found."
        )

    # Extract detailed metadata preserved in the context field.
    direct["activity_id"] = direct["context"].apply(
        lambda x: extract_context_value(x, "activity_id")
    )
    direct["target_chembl_id"] = direct["context"].apply(
        lambda x: extract_context_value(x, "target_id")
    )
    direct["target_description"] = direct["context"].apply(
        lambda x: extract_context_value(x, "target_name")
    )
    direct["target_type"] = direct["context"].apply(
        lambda x: extract_context_value(x, "target_type")
    )
    direct["activity_type"] = direct["context"].apply(
        lambda x: extract_context_value(x, "activity_type_raw")
    )
    direct["assay_type"] = direct["context"].apply(
        lambda x: extract_context_value(x, "assay_type")
    )
    direct["assay_description"] = direct["context"].apply(
        lambda x: extract_context_value(x, "assay_description")
    )
    direct["target_organism"] = direct["context"].apply(
        lambda x: extract_context_value(x, "target_organism")
    )
    direct["organism_class"] = direct["context"].apply(
        lambda x: extract_context_value(x, "organism_class")
    )
    direct["pchembl_value"] = direct["context"].apply(
        lambda x: extract_context_value(x, "pchembl_value")
    )
    direct["value_quality"] = direct["context"].apply(
        lambda x: extract_context_value(x, "value_quality")
    )

    # Keep the assay-level fields most useful for explanations.
    columns = [
        "rank",
        "drug_id",
        "drug_name",
        "disease_id",
        "disease_name",
        "rgcn_score",
        "source_name",
        "target_name",
        "target_chembl_id",
        "target_description",
        "target_type",
        "activity_id",
        "publication",
        "assay",
        "assay_type",
        "assay_description",
        "activity_type",
        "activity_value",
        "pchembl_value",
        "target_organism",
        "organism_class",
        "value_quality",
        "edge_id",
    ]

    report = direct[columns].copy()

    # Remove repeated copies of the same candidate/assay record.
    report = report.drop_duplicates(
        subset=[
            "drug_id",
            "disease_id",
            "activity_id",
            "assay",
            "activity_value",
        ]
    )

    # Put the highest-ranked candidates first.
    report = report.sort_values(
        ["rank", "drug_name", "target_name", "activity_id"]
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    report.to_csv(OUTPUT_PATH, index=False)

    print("=" * 70)
    print("REMEDY DIRECT DRUG-TARGET EVIDENCE REPORT")
    print("=" * 70)
    print("Input ChEMBL rows:", len(direct))
    print("Report rows after deduplication:", len(report))
    print("Unique candidates:", report[
        ["drug_id", "disease_id"]
    ].drop_duplicates().shape[0])
    print("Unique drugs:", report["drug_id"].nunique())
    print("Unique targets:", report["target_name"].nunique())
    print("Unique assay IDs:", report["assay"].nunique())
    print("Unique activity IDs:", report["activity_id"].nunique())
    print(
        "Rows with human target organism:",
        int(
            report["target_organism"]
            .str.contains("Homo sapiens", case=False, na=False)
            .sum()
        ),
    )

    print("\nSample evidence records:")
    print(
        report[
            [
                "rank",
                "drug_name",
                "disease_name",
                "target_name",
                "activity_type",
                "activity_value",
                "target_organism",
                "publication",
            ]
        ].head(15).to_string(index=False)
    )

    print("\nSaved report:")
    print(OUTPUT_PATH)

    print(
        "\nInterpretation note: these are drug-target assay records "
        "associated with ranked candidates. They do not independently "
        "prove drug efficacy for the candidate diseases."
    )


if __name__ == "__main__":
    main()
