
from pathlib import Path
import pandas as pd

# ============================================================
# REMEDY: EXPLAINABLE AI CANDIDATE REPORT
# Uses existing outputs; does not retrain the model or alter
# the frozen knowledge graph.
# ============================================================

ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = ROOT / "reports" / "model_evaluation"

TRACE_PATH = RESULTS_DIR / "top20_candidate_evidence_traces.csv"
DIRECT_PATH = RESULTS_DIR / "top20_candidate_direct_evidence_report.csv"

CSV_OUTPUT = RESULTS_DIR / "top20_explainable_candidate_report.csv"
MD_OUTPUT = RESULTS_DIR / "top20_explainable_candidate_report.md"


def unique_values(values):
    """Return sorted, unique, non-empty values as text."""
    cleaned = values.dropna().astype(str).str.strip()
    cleaned = cleaned[cleaned != ""]
    cleaned = cleaned[~cleaned.str.lower().isin(
        ["unknown", "none", "nan", "null", "n/a", "na"]
    )]
    return "; ".join(sorted(set(cleaned)))


def load_inputs():
    """Load and validate the existing evidence files."""
    for path in [TRACE_PATH, DIRECT_PATH]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required input file not found: {path}"
            )

    traces = pd.read_csv(TRACE_PATH)
    direct = pd.read_csv(DIRECT_PATH)

    trace_required = {
        "rank",
        "drug_id",
        "drug_name",
        "disease_id",
        "disease_name",
        "rgcn_score",
        "shared_protein_count",
        "shared_pathway_count",
        "shared_protein_names",
        "shared_pathway_names",
    }

    direct_required = {
        "drug_id",
        "disease_id",
        "target_name",
        "activity_id",
        "activity_type",
        "activity_value",
        "assay",
        "publication",
        "assay_description",
        "target_organism",
    }

    missing_trace = trace_required - set(traces.columns)
    missing_direct = direct_required - set(direct.columns)

    if missing_trace:
        raise ValueError(
            f"Missing trace columns: {sorted(missing_trace)}"
        )

    if missing_direct:
        raise ValueError(
            f"Missing direct-evidence columns: {sorted(missing_direct)}"
        )

    # Normalize IDs before joining tables.
    for column in ["drug_id", "disease_id"]:
        traces[column] = traces[column].astype(str).str.strip()
        direct[column] = direct[column].astype(str).str.strip()

    return traces, direct


def summarize_direct_evidence(direct):
    """Summarize available ChEMBL evidence per drug-disease pair."""
    direct = direct.copy()

    # Build readable assay descriptions.
    def describe_assay(row):
        target = str(row["target_name"])
        activity = str(row["activity_value"])
        assay = str(row["assay"])
        document = str(row["publication"])

        parts = [f"{target}: {activity}"]

        if assay and assay.lower() not in ["nan", "unknown"]:
            parts.append(f"assay={assay}")

        if document and document.lower() not in [
            "nan", "unknown", "none"
        ]:
            parts.append(f"document={document}")

        return " | ".join(parts)

    direct["assay_detail"] = direct.apply(
        describe_assay, axis=1
    )

    grouped_rows = []

    for (drug_id, disease_id), group in direct.groupby(
        ["drug_id", "disease_id"], sort=False
    ):
        # Deduplicate assay details for the summary.
        assay_details = sorted(set(group["assay_detail"]))

        grouped_rows.append({
            "drug_id": drug_id,
            "disease_id": disease_id,
            "direct_assay_record_count": len(group),
            "direct_target_count": group["target_name"].nunique(),
            "direct_target_names": unique_values(
                group["target_name"]
            ),
            "direct_assay_ids": unique_values(group["assay"]),
            "direct_activity_types": unique_values(
                group["activity_type"]
            ),
            "direct_assay_details": "; ".join(assay_details),
            "direct_document_ids": unique_values(
                group["publication"]
            ),
            "direct_target_organisms": unique_values(
                group["target_organism"]
            ),
            "direct_assay_descriptions": unique_values(
                group["assay_description"]
            ),
        })

    return pd.DataFrame(grouped_rows)


def make_explanation(row):
    """Create a cautious, readable explanation for a candidate."""
    protein_count = int(row["shared_protein_count"])
    pathway_count = int(row["shared_pathway_count"])
    assay_count = int(row["direct_assay_record_count"])

    parts = []

    parts.append(
        f"The R-GCN ranked this drug-disease pair at position "
        f"{int(row['rank'])} with a score of "
        f"{float(row['rgcn_score']):.6f}. This score is an "
        f"uncalibrated ranking score, not a probability."
    )

    if protein_count:
        parts.append(
            f"The graph contains {protein_count} shared protein "
            f"connection(s) in the drug and disease neighborhoods."
        )
    else:
        parts.append(
            "No shared protein connections were found by the "
            "current neighborhood-tracing method."
        )

    if pathway_count:
        parts.append(
            f"The two graph neighborhoods share {pathway_count} "
            f"pathway connection(s). These are indirect graph "
            f"connections, not proof of a causal mechanism."
        )
    else:
        parts.append(
            "No shared pathway connections were found by the "
            "current neighborhood-tracing method."
        )

    if assay_count:
        parts.append(
            f"The current evidence export contains {assay_count} "
            f"ChEMBL drug-target assay record(s) associated with "
            f"this drug-disease candidate. These support drug-target "
            f"relationships, not disease-treatment efficacy."
        )
    else:
        parts.append(
            "No ChEMBL drug-target assay records were found for "
            "this pair in the current top-20 evidence export. "
            "This does not mean that no such evidence exists elsewhere."
        )

    parts.append(
        "This candidate is a research hypothesis requiring "
        "independent validation, not a clinical recommendation."
    )

    return " ".join(parts)


def build_markdown(report):
    """Create a readable Markdown report."""
    lines = []

    lines.append("# ReMedy — Explainable AI Candidate Report")
    lines.append("")
    lines.append("## Purpose")
    lines.append("")
    lines.append(
        "This report summarizes the top 20 candidates from the "
        "existing R-GCN ranking and connects each candidate to "
        "available knowledge-graph and ChEMBL drug-target evidence."
    )
    lines.append("")
    lines.append("## Important interpretation limits")
    lines.append("")
    lines.append(
        "- R-GCN scores are ranking scores, not calibrated probabilities."
    )
    lines.append(
        "- Shared proteins and pathways indicate graph connections, "
        "not proven causal mechanisms."
    )
    lines.append(
        "- Drug-target assay evidence does not establish efficacy "
        "for the candidate disease."
    )
    lines.append(
        "- Missing evidence in this export is not proof that no "
        "evidence exists elsewhere."
    )
    lines.append(
        "- All candidates are research hypotheses, not clinical advice."
    )
    lines.append("")
    lines.append("## Candidate overview")
    lines.append("")
    lines.append(
        "| Rank | Drug | Disease | R-GCN score | "
        "Shared proteins | Shared pathways | ChEMBL assay records |"
    )
    lines.append(
        "|---:|---|---|---:|---:|---:|---:|"
    )

    for row in report.itertuples(index=False):
        lines.append(
            f"| {row.rank} | {row.drug_name} | {row.disease_name} "
            f"| {row.rgcn_score:.6f} | {row.shared_protein_count} "
            f"| {row.shared_pathway_count} "
            f"| {row.direct_assay_record_count} |"
        )

    lines.append("")
    lines.append("## Candidate explanations")
    lines.append("")

    for row in report.itertuples(index=False):
        lines.append(
            f"### {row.rank}. {row.drug_name} — {row.disease_name}"
        )
        lines.append("")
        lines.append(row.explanation)
        lines.append("")

        if row.shared_protein_names:
            lines.append(
                f"**Shared proteins:** {row.shared_protein_names}"
            )
            lines.append("")

        if row.shared_pathway_names:
            lines.append(
                f"**Shared pathways:** {row.shared_pathway_names}"
            )
            lines.append("")

        if row.direct_target_names:
            lines.append(
                f"**Targets with ChEMBL assay records:** "
                f"{row.direct_target_names}"
            )
            lines.append("")

        if row.direct_activity_types:
            lines.append(
                f"**Reported activity types:** "
                f"{row.direct_activity_types}"
            )
            lines.append("")

        if row.direct_assay_details:
            lines.append("**Available assay details:**")
            lines.append("")

            # Limit Markdown detail to keep the report readable.
            details = row.direct_assay_details.split("; ")
            for detail in details[:5]:
                lines.append(f"- {detail}")

            if len(details) > 5:
                lines.append(
                    f"- Plus {len(details) - 5} additional assay record(s); "
                    "see the CSV report."
                )

            lines.append("")

        lines.append("---")
        lines.append("")

    return "\n".join(lines)


def main():
    print("Loading existing ReMedy candidate evidence...")

    traces, direct = load_inputs()
    direct_summary = summarize_direct_evidence(direct)

    print("Candidate traces:", len(traces))
    print("Drug-disease pairs with ChEMBL records:", len(direct_summary))

    # Join using canonical drug and disease IDs.
    report = traces.merge(
        direct_summary,
        on=["drug_id", "disease_id"],
        how="left",
        validate="one_to_one",
    )

    # Missing assay data is represented explicitly.
    numeric_columns = [
        "direct_assay_record_count",
        "direct_target_count",
    ]

    for column in numeric_columns:
        report[column] = report[column].fillna(0).astype(int)

    text_columns = [
        "direct_target_names",
        "direct_assay_ids",
        "direct_activity_types",
        "direct_assay_details",
        "direct_document_ids",
        "direct_target_organisms",
        "direct_assay_descriptions",
    ]

    for column in text_columns:
        report[column] = report[column].fillna("")

    report["explanation"] = report.apply(
        make_explanation, axis=1
    )

    report = report.sort_values("rank").reset_index(drop=True)

    # Save the structured report.
    report.to_csv(CSV_OUTPUT, index=False)

    # Save the readable report.
    markdown = build_markdown(report)
    MD_OUTPUT.write_text(markdown, encoding="utf-8")

    # Print a concise audit.
    print("\n" + "=" * 70)
    print("REMEDY EXPLAINABLE CANDIDATE REPORT")
    print("=" * 70)

    print("Candidates reported:", len(report))
    print(
        "Candidates with ChEMBL assay records:",
        int((report["direct_assay_record_count"] > 0).sum()),
    )
    print(
        "Candidates without ChEMBL records in this export:",
        int((report["direct_assay_record_count"] == 0).sum()),
    )
    print(
        "Candidates with shared protein connections:",
        int((report["shared_protein_count"] > 0).sum()),
    )
    print(
        "Candidates with shared pathway connections:",
        int((report["shared_pathway_count"] > 0).sum()),
    )

    print("\nTop 10 candidate overview:")
    print(
        report[
            [
                "rank",
                "drug_name",
                "disease_name",
                "rgcn_score",
                "shared_protein_count",
                "shared_pathway_count",
                "direct_assay_record_count",
            ]
        ].head(10).to_string(index=False)
    )

    print("\nFiles created:")
    print(CSV_OUTPUT)
    print(MD_OUTPUT)

    print(
        "\nCompleted. No model training or frozen graph files "
        "were modified."
    )


if __name__ == "__main__":
    main()
