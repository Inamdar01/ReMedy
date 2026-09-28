from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

EDGES_FILE = ROOT / "data/processed/edges.parquet"
EVIDENCE_FILE = ROOT / "data/processed/evidence.parquet"

REPORT_FILE = ROOT / "reports/H3_source_overlap_check.md"


def main():

    print("=" * 70)
    print("REMEDY - H3 MULTI-SOURCE EVIDENCE CHECK")
    print("=" * 70)

    edges = pd.read_parquet(
        EDGES_FILE
    )

    evidence = pd.read_parquet(
        EVIDENCE_FILE
    )

    print(
        f"\nEdges: {len(edges):,}"
    )

    print(
        f"Evidence rows: {len(evidence):,}"
    )

    # ---------------------------------------------------------
    # Count distinct evidence sources per edge
    # ---------------------------------------------------------

    source_counts = (
        evidence
        .groupby("edge_id")["source"]
        .nunique()
        .rename("source_count")
        .reset_index()
    )

    edge_source = edges.merge(
        source_counts,
        on="edge_id",
        how="left",
    )

    edge_source["source_count"] = (
        edge_source["source_count"]
        .fillna(0)
        .astype(int)
    )

    # ---------------------------------------------------------
    # Multi-source definition
    # ---------------------------------------------------------

    edge_source["evidence_class"] = (
        edge_source["source_count"]
        .apply(
            lambda x:
                "multi-source"
                if x >= 2
                else "single-source"
        )
    )

    class_counts = (
        edge_source["evidence_class"]
        .value_counts()
    )

    relation_source_summary = (
        edge_source
        .groupby(
            [
                "relation",
                "evidence_class",
            ]
        )
        .size()
        .reset_index(
            name="edge_count"
        )
    )

    # ---------------------------------------------------------
    # Source combinations
    # ---------------------------------------------------------

    source_combinations = (
        evidence
        .groupby("edge_id")["source"]
        .apply(
            lambda s:
                " | ".join(
                    sorted(
                        set(
                            str(x).strip()
                            for x in s
                            if str(x).strip()
                        )
                    )
                )
        )
        .rename("source_combination")
        .reset_index()
    )

    combination_counts = (
        source_combinations[
            "source_combination"
        ]
        .value_counts()
        .rename_axis(
            "source_combination"
        )
        .reset_index(
            name="edge_count"
        )
    )

    # ---------------------------------------------------------
    # Terminal output
    # ---------------------------------------------------------

    multi_source = int(
        class_counts.get(
            "multi-source",
            0
        )
    )

    single_source = int(
        class_counts.get(
            "single-source",
            0
        )
    )

    print()
    print("=" * 70)
    print("H3 SOURCE OVERLAP CHECK")
    print("=" * 70)

    print(
        f"Single-source edges: "
        f"{single_source:,}"
    )

    print(
        f"Multi-source edges: "
        f"{multi_source:,}"
    )

    if len(edges) > 0:

        print(
            f"Multi-source coverage: "
            f"{multi_source / len(edges) * 100:.4f}%"
        )

    print("\nSource combinations:")

    print(
        combination_counts.to_string(
            index=False
        )
    )

    print("\nBy relation:")

    print(
        relation_source_summary.to_string(
            index=False
        )
    )

    # ---------------------------------------------------------
    # Report
    # ---------------------------------------------------------

    report = f"""# H3 Multi-Source Evidence Check

## Purpose

H3 asks whether multi-source relations are more reliable than
single-source relations.

## Current final KG

- Total edges: {len(edges):,}
- Total evidence records: {len(evidence):,}
- Single-source edges: {single_source:,}
- Multi-source edges: {multi_source:,}
- Multi-source coverage: {
        (multi_source / len(edges) * 100)
        if len(edges)
        else 0
    :.4f}%

## Source combinations

{combination_counts.to_markdown(index=False)}

## Relation/source distribution

{relation_source_summary.to_markdown(index=False)}

## Decision

"""

    if multi_source == 0:

        report += """No canonical edge currently has evidence from
two or more distinct source databases.

Therefore H3 cannot be statistically tested using the current
final evidence table without adding another independent source
for at least some existing canonical relations.

An H3 result should not be fabricated from single-source edges.
"""

        decision = (
            "H3 NOT CURRENTLY TESTABLE"
        )

    else:

        report += """Multi-source edges are present, so a formal
comparison of multi-source versus single-source evidence can
proceed.
"""

        decision = (
            "H3 CAN PROCEED"
        )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8"
    )

    print()
    print(
        f"Decision: {decision}"
    )

    print(
        f"Saved: {REPORT_FILE}"
    )

    print("=" * 70)


if __name__ == "__main__":
    main()