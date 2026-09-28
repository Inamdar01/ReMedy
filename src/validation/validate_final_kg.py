from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

NODES_FILE = ROOT / "data/processed/nodes.parquet"
EDGES_FILE = ROOT / "data/processed/edges.parquet"
EVIDENCE_FILE = ROOT / "data/processed/evidence.parquet"

REPORT_FILE = ROOT / "reports/final_kg_validation_report.md"


EXPECTED_NODE_COLUMNS = [
    "node_id",
    "node_type",
    "label",
    "source_database",
]

EXPECTED_EDGE_COLUMNS = [
    "edge_id",
    "source_id",
    "relation",
    "target_id",
    "direction",
    "evidence_count",
]

EXPECTED_EVIDENCE_COLUMNS = [
    "edge_id",
    "source",
    "publication",
    "assay",
    "organism",
    "cell_line",
    "activity_value",
    "context",
]


def check_columns(df, expected, name):
    actual = list(df.columns)

    missing = [
        col for col in expected
        if col not in actual
    ]

    extra = [
        col for col in actual
        if col not in expected
    ]

    return {
        "name": name,
        "missing": missing,
        "extra": extra,
        "valid": not missing,
    }


def main():

    print("=" * 70)
    print("REMEDY - FINAL KNOWLEDGE GRAPH VALIDATION")
    print("=" * 70)

    # ---------------------------------------------------------
    # 1. Load Parquet files
    # ---------------------------------------------------------

    print("\nLoading final Parquet files...")

    nodes = pd.read_parquet(NODES_FILE)
    edges = pd.read_parquet(EDGES_FILE)
    evidence = pd.read_parquet(EVIDENCE_FILE)

    print(
        f"Nodes: {len(nodes):,}"
    )

    print(
        f"Edges: {len(edges):,}"
    )

    print(
        f"Evidence: {len(evidence):,}"
    )

    failures = []

    # ---------------------------------------------------------
    # 2. Schema checks
    # ---------------------------------------------------------

    schema_results = [
        check_columns(
            nodes,
            EXPECTED_NODE_COLUMNS,
            "nodes",
        ),
        check_columns(
            edges,
            EXPECTED_EDGE_COLUMNS,
            "edges",
        ),
        check_columns(
            evidence,
            EXPECTED_EVIDENCE_COLUMNS,
            "evidence",
        ),
    ]

    for result in schema_results:

        if result["missing"]:

            failures.append(
                f"{result['name']} missing columns: "
                f"{result['missing']}"
            )

        if result["extra"]:

            print(
                f"Note: {result['name']} has extra columns: "
                f"{result['extra']}"
            )

    # ---------------------------------------------------------
    # 3. Node ID uniqueness
    # ---------------------------------------------------------

    duplicate_node_ids = int(
        nodes["node_id"].duplicated().sum()
    )

    print(
        f"\nDuplicate node IDs: "
        f"{duplicate_node_ids:,}"
    )

    if duplicate_node_ids:
        failures.append(
            f"Duplicate node IDs: {duplicate_node_ids:,}"
        )

    # ---------------------------------------------------------
    # 4. Edge ID uniqueness
    # ---------------------------------------------------------

    duplicate_edge_ids = int(
        edges["edge_id"].duplicated().sum()
    )

    print(
        f"Duplicate edge IDs: "
        f"{duplicate_edge_ids:,}"
    )

    if duplicate_edge_ids:
        failures.append(
            f"Duplicate edge IDs: {duplicate_edge_ids:,}"
        )

    # ---------------------------------------------------------
    # 5. Canonical edge uniqueness
    # ---------------------------------------------------------

    duplicate_edges = int(
        edges.duplicated(
            subset=[
                "source_id",
                "relation",
                "target_id",
            ]
        ).sum()
    )

    print(
        f"Duplicate canonical edges: "
        f"{duplicate_edges:,}"
    )

    if duplicate_edges:
        failures.append(
            f"Duplicate canonical edges: "
            f"{duplicate_edges:,}"
        )

    # ---------------------------------------------------------
    # 6. Missing node values
    # ---------------------------------------------------------

    node_missing = {}

    for col in EXPECTED_NODE_COLUMNS:

        count = int(
            nodes[col]
            .isna()
            .sum()
        )

        blank = int(
            nodes[col]
            .astype(str)
            .str.strip()
            .eq("")
            .sum()
        )

        node_missing[col] = count + blank

    # ---------------------------------------------------------
    # 7. Missing edge values
    # ---------------------------------------------------------

    edge_missing = {}

    for col in EXPECTED_EDGE_COLUMNS:

        count = int(
            edges[col]
            .isna()
            .sum()
        )

        blank = int(
            edges[col]
            .astype(str)
            .str.strip()
            .eq("")
            .sum()
        )

        edge_missing[col] = count + blank

    # ---------------------------------------------------------
    # 8. Missing evidence values
    # ---------------------------------------------------------

    evidence_missing = {}

    for col in EXPECTED_EVIDENCE_COLUMNS:

        count = int(
            evidence[col]
            .isna()
            .sum()
        )

        blank = int(
            evidence[col]
            .astype(str)
            .str.strip()
            .eq("")
            .sum()
        )

        evidence_missing[col] = count + blank

    # ---------------------------------------------------------
    # 9. Edge endpoints must exist as nodes
    # ---------------------------------------------------------

    node_ids = set(
        nodes["node_id"]
        .astype(str)
    )

    orphan_sources = (
        set(
            edges["source_id"]
            .astype(str)
        )
        - node_ids
    )

    orphan_targets = (
        set(
            edges["target_id"]
            .astype(str)
        )
        - node_ids
    )

    print(
        f"\nOrphan edge source IDs: "
        f"{len(orphan_sources):,}"
    )

    print(
        f"Orphan edge target IDs: "
        f"{len(orphan_targets):,}"
    )

    if orphan_sources:
        failures.append(
            f"Orphan edge source IDs: "
            f"{len(orphan_sources):,}"
        )

    if orphan_targets:
        failures.append(
            f"Orphan edge target IDs: "
            f"{len(orphan_targets):,}"
        )

    # ---------------------------------------------------------
    # 10. Evidence edge IDs must exist
    # ---------------------------------------------------------

    edge_ids = set(
        edges["edge_id"]
        .astype(str)
    )

    evidence_edge_ids = set(
        evidence["edge_id"]
        .astype(str)
    )

    orphan_evidence = (
        evidence_edge_ids
        - edge_ids
    )

    edges_without_evidence = (
        edge_ids
        - evidence_edge_ids
    )

    print(
        f"Orphan evidence edge IDs: "
        f"{len(orphan_evidence):,}"
    )

    print(
        f"Edges without evidence records: "
        f"{len(edges_without_evidence):,}"
    )

    if orphan_evidence:
        failures.append(
            f"Orphan evidence edge IDs: "
            f"{len(orphan_evidence):,}"
        )

    if edges_without_evidence:
        failures.append(
            f"Edges without evidence: "
            f"{len(edges_without_evidence):,}"
        )

    # ---------------------------------------------------------
    # 11. Evidence count must match evidence.parquet
    # ---------------------------------------------------------

    actual_evidence_counts = (
        evidence
        .groupby("edge_id")
        .size()
        .rename("actual_evidence_count")
    )

    check_edges = edges[
        [
            "edge_id",
            "evidence_count",
        ]
    ].copy()

    check_edges["actual_evidence_count"] = (
        check_edges["edge_id"]
        .map(actual_evidence_counts)
        .fillna(0)
        .astype(int)
    )

    evidence_count_mismatches = int(
        (
            check_edges["evidence_count"].astype(int)
            !=
            check_edges["actual_evidence_count"]
        ).sum()
    )

    print(
        f"Evidence-count mismatches: "
        f"{evidence_count_mismatches:,}"
    )

    if evidence_count_mismatches:
        failures.append(
            f"Evidence-count mismatches: "
            f"{evidence_count_mismatches:,}"
        )

    # ---------------------------------------------------------
    # 12. Validate relation endpoint types
    # ---------------------------------------------------------

    node_type_lookup = dict(
        zip(
            nodes["node_id"].astype(str),
            nodes["node_type"].astype(str),
        )
    )

    relation_rules = {
        "indication": ("drug", "disease"),
        "off-label use": ("drug", "disease"),
        "contraindication": ("drug", "disease"),
        "targets": ("drug", "protein"),
        "associated_with": ("disease", "protein"),
        "participates_in": ("protein", "pathway"),
    }

    invalid_relation_endpoints = []

    for _, row in edges.iterrows():

        relation = str(
            row["relation"]
        )

        expected = relation_rules.get(
            relation
        )

        if expected is None:

            invalid_relation_endpoints.append({
                "edge_id": row["edge_id"],
                "relation": relation,
                "reason": "unknown relation",
            })

            continue

        actual_source_type = node_type_lookup.get(
            str(row["source_id"])
        )

        actual_target_type = node_type_lookup.get(
            str(row["target_id"])
        )

        if (
            actual_source_type != expected[0]
            or
            actual_target_type != expected[1]
        ):

            invalid_relation_endpoints.append({
                "edge_id": row["edge_id"],
                "relation": relation,
                "expected": f"{expected[0]}->{expected[1]}",
                "actual": (
                    f"{actual_source_type}"
                    f"->{actual_target_type}"
                ),
            })

    print(
        f"Invalid relation endpoint types: "
        f"{len(invalid_relation_endpoints):,}"
    )

    if invalid_relation_endpoints:

        failures.append(
            "Invalid relation endpoint types: "
            f"{len(invalid_relation_endpoints):,}"
        )

    # ---------------------------------------------------------
    # 13. Node-type distribution
    # ---------------------------------------------------------

    node_counts = (
        nodes["node_type"]
        .value_counts()
        .sort_index()
    )

    # ---------------------------------------------------------
    # 14. Relation distribution
    # ---------------------------------------------------------

    relation_counts = (
        edges["relation"]
        .value_counts()
        .sort_index()
    )

    # ---------------------------------------------------------
    # 15. Evidence source distribution
    # ---------------------------------------------------------

    evidence_sources = (
        evidence["source"]
        .value_counts()
        .sort_index()
    )

    # ---------------------------------------------------------
    # 16. Evidence coverage
    # ---------------------------------------------------------

    total_edges = len(edges)

    covered_edges = (
        total_edges
        - len(edges_without_evidence)
    )

    evidence_coverage_pct = (
        covered_edges
        / total_edges
        * 100
        if total_edges
        else 0
    )

    # ---------------------------------------------------------
    # 17. Duplicate evidence records
    # ---------------------------------------------------------

    duplicate_evidence = int(
        evidence.duplicated().sum()
    )

    print(
        f"Duplicate evidence rows: "
        f"{duplicate_evidence:,}"
    )

    # Duplicate evidence is reported but not automatically
    # treated as a graph failure because two records may
    # legitimately contain different provenance/context.
    #
    # Exact duplicate rows, however, should normally be zero.

    if duplicate_evidence:
        failures.append(
            f"Duplicate evidence rows: "
            f"{duplicate_evidence:,}"
        )

    # ---------------------------------------------------------
    # 18. Write report
    # ---------------------------------------------------------

    report = []

    report.append(
        "# ReMedy Final KG Validation Report\n\n"
    )

    report.append(
        "Validation performed directly against the final "
        "Parquet KG files.\n\n"
    )

    report.append(
        "## 1. Files\n\n"
    )

    report.append(
        "- `data/processed/nodes.parquet`\n"
    )

    report.append(
        "- `data/processed/edges.parquet`\n"
    )

    report.append(
        "- `data/processed/evidence.parquet`\n"
    )

    report.append(
        "\n## 2. Counts\n\n"
    )

    report.append(
        f"- Nodes: {len(nodes):,}\n"
    )

    report.append(
        f"- Edges: {len(edges):,}\n"
    )

    report.append(
        f"- Evidence records: {len(evidence):,}\n"
    )

    report.append(
        f"- Edges with evidence: {covered_edges:,}\n"
    )

    report.append(
        f"- Evidence coverage: "
        f"{evidence_coverage_pct:.2f}%\n"
    )

    report.append(
        "\n## 3. Node types\n\n"
    )

    report.append(
        "| Node type | Count |\n"
    )

    report.append(
        "|---|---:|\n"
    )

    for node_type, count in node_counts.items():

        report.append(
            f"| {node_type} | {count:,} |\n"
        )

    report.append(
        "\n## 4. Relations\n\n"
    )

    report.append(
        "| Relation | Count |\n"
    )

    report.append(
        "|---|---:|\n"
    )

    for relation, count in relation_counts.items():

        report.append(
            f"| {relation} | {count:,} |\n"
        )

    report.append(
        "\n## 5. Evidence sources\n\n"
    )

    report.append(
        "| Source | Records |\n"
    )

    report.append(
        "|---|---:|\n"
    )

    for source, count in evidence_sources.items():

        report.append(
            f"| {source} | {count:,} |\n"
        )

    report.append(
        "\n## 6. Integrity checks\n\n"
    )

    checks = [
        (
            "Duplicate node IDs",
            duplicate_node_ids,
        ),
        (
            "Duplicate edge IDs",
            duplicate_edge_ids,
        ),
        (
            "Duplicate canonical edges",
            duplicate_edges,
        ),
        (
            "Orphan edge source IDs",
            len(orphan_sources),
        ),
        (
            "Orphan edge target IDs",
            len(orphan_targets),
        ),
        (
            "Orphan evidence edge IDs",
            len(orphan_evidence),
        ),
        (
            "Edges without evidence",
            len(edges_without_evidence),
        ),
        (
            "Evidence-count mismatches",
            evidence_count_mismatches,
        ),
        (
            "Invalid relation endpoint types",
            len(invalid_relation_endpoints),
        ),
        (
            "Duplicate evidence rows",
            duplicate_evidence,
        ),
    ]

    report.append(
        "| Check | Count |\n"
    )

    report.append(
        "|---|---:|\n"
    )

    for name, value in checks:

        report.append(
            f"| {name} | {value:,} |\n"
        )

    report.append(
        "\n## 7. Missing-value counts\n\n"
    )

    report.append(
        "### Nodes\n\n"
    )

    report.append(
        "| Column | Missing/blank |\n"
    )

    report.append(
        "|---|---:|\n"
    )

    for col, value in node_missing.items():

        report.append(
            f"| {col} | {value:,} |\n"
        )

    report.append(
        "\n### Edges\n\n"
    )

    report.append(
        "| Column | Missing/blank |\n"
    )

    report.append(
        "|---|---:|\n"
    )

    for col, value in edge_missing.items():

        report.append(
            f"| {col} | {value:,} |\n"
        )

    report.append(
        "\n### Evidence\n\n"
    )

    report.append(
        "| Column | Missing/blank |\n"
    )

    report.append(
        "|---|---:|\n"
    )

    for col, value in evidence_missing.items():

        report.append(
            f"| {col} | {value:,} |\n"
        )

    report.append(
        "\n## 8. Validation status\n\n"
    )

    if failures:

        report.append(
            "**FAILED**\n\n"
        )

        report.append(
            "Issues detected:\n\n"
        )

        for failure in failures:

            report.append(
                f"- {failure}\n"
            )

    else:

        report.append(
            "**PASSED**\n\n"
        )

        report.append(
            "All structural integrity checks passed.\n"
        )

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    REPORT_FILE.write_text(
        "".join(report),
        encoding="utf-8",
    )

    # ---------------------------------------------------------
    # Terminal summary
    # ---------------------------------------------------------

    print()
    print("=" * 70)

    if failures:
        print("FINAL KG VALIDATION: ISSUES FOUND")
    else:
        print("FINAL KG VALIDATION: PASSED")

    print("=" * 70)

    print(
        f"Nodes: {len(nodes):,}"
    )

    print(
        f"Edges: {len(edges):,}"
    )

    print(
        f"Evidence records: {len(evidence):,}"
    )

    print(
        f"Evidence coverage: "
        f"{evidence_coverage_pct:.2f}%"
    )

    if failures:

        print("\nIssues:")

        for failure in failures:
            print(
                f"- {failure}"
            )

    print()
    print(
        f"Saved: {REPORT_FILE}"
    )

    print("=" * 70)

    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()