
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = ROOT / "reports" / "model_evaluation"

REPORT_PATH = (
    RESULTS_DIR / "top20_explainable_candidate_report.csv"
)
DIRECT_PATH = (
    RESULTS_DIR / "top20_candidate_direct_evidence_report.csv"
)

OUTPUT_PATH = RESULTS_DIR / "explainable_report_audit.txt"


def split_values(value):
    """Split semicolon-separated values safely."""
    if pd.isna(value) or not str(value).strip():
        return []
    return [
        item.strip()
        for item in str(value).split(";")
        if item.strip()
    ]


def main():
    for path in [REPORT_PATH, DIRECT_PATH]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required file not found: {path}"
            )

    report = pd.read_csv(REPORT_PATH)
    direct = pd.read_csv(DIRECT_PATH)

    required_report_columns = [
        "rank",
        "drug_id",
        "disease_id",
        "drug_name",
        "disease_name",
        "rgcn_score",
        "shared_protein_count",
        "shared_pathway_count",
        "shared_protein_ids",
        "shared_pathway_ids",
        "direct_assay_record_count",
    ]

    missing = set(required_report_columns) - set(report.columns)
    if missing:
        raise ValueError(
            f"Missing report columns: {sorted(missing)}"
        )

    findings = []

    def check(title, passed, details=""):
        status = "PASS" if passed else "REVIEW"
        findings.append(f"[{status}] {title}")
        if details:
            findings.append(f"       {details}")

    # --------------------------------------------------
    # 1. Candidate identity and duplication
    # --------------------------------------------------

    pair_duplicates = report.duplicated(
        subset=["drug_id", "disease_id"]
    ).sum()

    check(
        "Candidate drug-disease pairs are unique",
        pair_duplicates == 0,
        f"Duplicate pairs: {pair_duplicates}",
    )

    check(
        "Report contains 20 candidates",
        len(report) == 20,
        f"Rows found: {len(report)}",
    )

    # --------------------------------------------------
    # 2. Ranking and score checks
    # --------------------------------------------------

    scores = pd.to_numeric(
        report["rgcn_score"], errors="coerce"
    )

    check(
        "All R-GCN scores are numeric",
        scores.notna().all(),
        f"Missing/non-numeric scores: {scores.isna().sum()}",
    )

    if scores.notna().all():
        descending = scores.is_monotonic_decreasing
        check(
            "Candidates are sorted by descending score",
            descending,
            f"Score range: {scores.min():.6f} to {scores.max():.6f}",
        )

        check(
            "Scores are within the sigmoid output range",
            scores.between(0, 1).all(),
            "This range check does not establish calibration.",
        )

    # --------------------------------------------------
    # 3. Shared protein and pathway count consistency
    # --------------------------------------------------

    protein_mismatches = []
    pathway_mismatches = []

    for row in report.to_dict(orient="records"):
        protein_ids = split_values(row["shared_protein_ids"])
        pathway_ids = split_values(row["shared_pathway_ids"])

        if len(protein_ids) != int(row["shared_protein_count"]):
            protein_mismatches.append(
                f"{row['drug_name']} / {row['disease_name']}"
            )

        if len(pathway_ids) != int(row["shared_pathway_count"]):
            pathway_mismatches.append(
                f"{row['drug_name']} / {row['disease_name']}"
            )

    check(
        "Shared protein IDs match reported counts",
        len(protein_mismatches) == 0,
        f"Mismatches: {len(protein_mismatches)}",
    )

    check(
        "Shared pathway IDs match reported counts",
        len(pathway_mismatches) == 0,
        f"Mismatches: {len(pathway_mismatches)}",
    )

    # --------------------------------------------------
    # 4. Assay count checks against the direct evidence file
    # --------------------------------------------------

    direct = direct.copy()

    for column in ["drug_id", "disease_id"]:
        direct[column] = direct[column].astype(str).str.strip()
        report[column] = report[column].astype(str).str.strip()

    expected_counts = (
        direct.groupby(["drug_id", "disease_id"])
        .size()
        .to_dict()
    )

    assay_mismatches = []

    for row in report.to_dict(orient="records"):
        key = (row["drug_id"], row["disease_id"])
        expected = expected_counts.get(key, 0)
        actual = int(row["direct_assay_record_count"])

        if actual != expected:
            assay_mismatches.append(
                f"{row['drug_name']} / {row['disease_name']}: "
                f"report={actual}, source={expected}"
            )

    check(
        "Candidate assay counts match the direct evidence file",
        len(assay_mismatches) == 0,
        f"Mismatches: {len(assay_mismatches)}",
    )

    # --------------------------------------------------
    # 5. Identify assays reused across different diseases
    # --------------------------------------------------

    required_direct_columns = {
        "drug_id",
        "disease_id",
        "activity_id",
    }

    missing_direct = required_direct_columns - set(direct.columns)

    if missing_direct:
        raise ValueError(
            f"Missing direct evidence columns: {sorted(missing_direct)}"
        )

    valid_activity = direct[
        direct["activity_id"].notna()
        & direct["activity_id"].astype(str).str.strip().ne("")
    ].copy()

    reused = (
        valid_activity.groupby(["drug_id", "activity_id"])
        ["disease_id"]
        .nunique()
        .reset_index(name="disease_count")
    )

    reused = reused[reused["disease_count"] > 1]

    check(
        "Repeated assays across diseases are identified",
        True,
        f"Drug/activity combinations reused across diseases: {len(reused)}. "
        "This can be expected; it is not automatically an error.",
    )

    # --------------------------------------------------
    # 6. Explainability coverage
    # --------------------------------------------------

    assay_count = pd.to_numeric(
        report["direct_assay_record_count"], errors="coerce"
    ).fillna(0)

    protein_count = pd.to_numeric(
        report["shared_protein_count"], errors="coerce"
    ).fillna(0)

    pathway_count = pd.to_numeric(
        report["shared_pathway_count"], errors="coerce"
    ).fillna(0)

    check(
        "Candidates without ChEMBL records are clearly countable",
        True,
        f"No ChEMBL records in this export: {(assay_count == 0).sum()}",
    )

    check(
        "Candidates with no shared proteins are identified",
        True,
        f"No shared proteins: {(protein_count == 0).sum()}",
    )

    check(
        "Candidates with no shared pathways are identified",
        True,
        f"No shared pathways: {(pathway_count == 0).sum()}",
    )

    # --------------------------------------------------
    # 7. Write audit results
    # --------------------------------------------------

    lines = [
        "=" * 70,
        "REMEDY EXPLAINABLE CANDIDATE REPORT AUDIT",
        "=" * 70,
        f"Report rows: {len(report)}",
        f"Unique drugs: {report['drug_id'].nunique()}",
        f"Unique diseases: {report['disease_id'].nunique()}",
        "",
        *findings,
        "",
        "ASSAY COUNT MISMATCH DETAILS",
        "-" * 70,
        *(assay_mismatches or ["None"]),
        "",
        "PROTEIN COUNT MISMATCH DETAILS",
        "-" * 70,
        *(protein_mismatches or ["None"]),
        "",
        "PATHWAY COUNT MISMATCH DETAILS",
        "-" * 70,
        *(pathway_mismatches or ["None"]),
        "",
        "INTERPRETATION",
        "-" * 70,
        "This audit checks report consistency, not biological validity.",
        "A passing audit does not prove candidate efficacy or model accuracy.",
    ]

    output_text = "\n".join(lines)
    OUTPUT_PATH.write_text(output_text, encoding="utf-8")

    print(output_text)
    print("\nAudit saved to:")
    print(OUTPUT_PATH)


if __name__ == "__main__":
    main()
