from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

FILES = {
    "drug_disease": ROOT / "data/interim/drug_disease_relations_canonical.csv",
    "drug_protein": ROOT / "data/interim/drug_protein_relations_canonical.csv",
    "disease_protein": ROOT / "data/interim/disease_uniprot_relations_canonical.csv",
    "protein_pathway": ROOT / "data/interim/protein_pathway_relations_canonical.csv",
}

OUTPUT_REPORT = ROOT / "reports/data_quality_report.md"
OUTPUT_STATS = ROOT / "reports/graph_statistics.csv"


def find_column(df, candidates):
    """Return the first matching column from a list of candidates."""
    for col in candidates:
        if col in df.columns:
            return col
    return None


def load_file(path):
    if not path.exists():
        print(f"WARNING: Missing file: {path}")
        return None

    df = pd.read_csv(path, dtype=str).fillna("")
    return df


def main():

    print("=" * 70)
    print("REMEDY - DATA QUALITY CHECK")
    print("=" * 70)

    tables = {}

    for name, path in FILES.items():
        print(f"\nLoading: {name}")
        df = load_file(path)

        if df is not None:
            tables[name] = df
            print(f"Rows: {len(df):,}")
            print(f"Columns: {len(df.columns)}")

    # ---------------------------------------------------------
    # Column-level missing values
    # ---------------------------------------------------------
    missing_rows = []

    for name, df in tables.items():

        for col in df.columns:

            missing_count = (
                df[col]
                .astype(str)
                .str.strip()
                .eq("")
                .sum()
            )

            missing_pct = (
                missing_count / len(df) * 100
                if len(df) > 0 else 0
            )

            missing_rows.append({
                "table": name,
                "column": col,
                "rows": len(df),
                "missing_count": int(missing_count),
                "missing_percentage": round(missing_pct, 2)
            })

    missing_df = pd.DataFrame(missing_rows)

    # ---------------------------------------------------------
    # Define canonical endpoint columns
    # ---------------------------------------------------------
    schemas = {
        "drug_disease": {
            "source_node": [
                "chembl_id",
                "drug_chembl_id",
            ],
            "target_node": [
                "mondo_id",
                "disease_mondo_id",
            ],
            "relation": ["relation"],
            "source": ["source", "provenance_source"],
        },

        "drug_protein": {
            "source_node": [
                "chembl_id",
                "drug_chembl_id",
            ],
            "target_node": [
                "uniprot_accession",
                "uniprot_id",
            ],
            "relation": ["relation"],
            "source": ["source", "provenance_source"],
        },

        "disease_protein": {
            "source_node": [
                "mondo_id",
                "disease_mondo_id",
            ],
            "target_node": [
                "uniprot_accession",
                "uniprot_id",
            ],
            "relation": ["relation"],
            "source": ["source", "provenance_source"],
        },

        "protein_pathway": {
            "source_node": [
                "uniprot_accession",
                "uniprot_id",
            ],
            "target_node": [
                "reactome_pathway_id",
                "pathway_id",
            ],
            "relation": ["relation"],
            "source": ["source", "provenance_source"],
        },
    }

    stats = []
    relation_rows = []

    all_nodes = {
        "drug": set(),
        "disease": set(),
        "protein": set(),
        "pathway": set(),
        "gene": set(),
    }

    duplicate_rows = []

    # ---------------------------------------------------------
    # Analyze each relation table
    # ---------------------------------------------------------
    for name, df in tables.items():

        schema = schemas[name]

        src_col = find_column(
            df,
            schema["source_node"]
        )

        dst_col = find_column(
            df,
            schema["target_node"]
        )

        rel_col = find_column(
            df,
            schema["relation"]
        )

        source_col = find_column(
            df,
            schema["source"]
        )

        # Missing endpoint information
        endpoint_missing = 0

        if src_col and dst_col:
            endpoint_missing = (
                (df[src_col].str.strip() == "")
                |
                (df[dst_col].str.strip() == "")
            ).sum()

        # Duplicate edge count
        duplicate_count = 0

        if src_col and dst_col and rel_col:

            edge_df = df[
                [src_col, rel_col, dst_col]
            ].copy()

            duplicate_count = int(
                edge_df.duplicated().sum()
            )

        duplicate_rows.append({
            "table": name,
            "duplicate_edge_count": duplicate_count,
        })

        # Provenance coverage
        #
        # Each canonical relation has a source-specific field:
        # drug_disease    -> disease_source
        # drug_protein    -> evidence_count
        # disease_protein -> disease_source + protein_source
        # protein_pathway -> source
        #
        # We measure whether a row has identifiable provenance
        # rather than requiring a literal "source" column.

        if name == "drug_disease":

            provenance_count = (
                df["disease_source"]
                .astype(str)
                .str.strip()
                .ne("")
                .sum()
            )

        elif name == "drug_protein":

            provenance_count = (
                pd.to_numeric(
                    df["evidence_count"],
                    errors="coerce"
                )
                .fillna(0)
                .gt(0)
                .sum()
            )

        elif name == "disease_protein":

            provenance_count = (
                df["disease_source"]
                .astype(str)
                .str.strip()
                .ne("")
                |
                df["protein_source"]
                .astype(str)
                .str.strip()
                .ne("")
            ).sum()

        elif name == "protein_pathway":

            provenance_count = (
                df["source"]
                .astype(str)
                .str.strip()
                .ne("")
            ).sum()

        else:
            provenance_count = 0

        provenance_pct = (
            provenance_count / len(df) * 100
            if len(df) else 0
        )

        # Store relationship statistics
        if rel_col:
            counts = (
                df[rel_col]
                .replace("", pd.NA)
                .dropna()
                .value_counts()
            )

            for relation, count in counts.items():

                relation_rows.append({
                    "table": name,
                    "relation": relation,
                    "edge_count": int(count),
                })

        # -----------------------------------------------------
        # Node collection
        # -----------------------------------------------------
        if src_col and dst_col:

            if name == "drug_disease":
                all_nodes["drug"].update(
                    df[src_col][
                        df[src_col].str.strip() != ""
                    ]
                )

                all_nodes["disease"].update(
                    df[dst_col][
                        df[dst_col].str.strip() != ""
                    ]
                )

            elif name == "drug_protein":
                all_nodes["drug"].update(
                    df[src_col][
                        df[src_col].str.strip() != ""
                    ]
                )

                all_nodes["protein"].update(
                    df[dst_col][
                        df[dst_col].str.strip() != ""
                    ]
                )

            elif name == "disease_protein":
                all_nodes["disease"].update(
                    df[src_col][
                        df[src_col].str.strip() != ""
                    ]
                )

                all_nodes["protein"].update(
                    df[dst_col][
                        df[dst_col].str.strip() != ""
                    ]
                )

            elif name == "protein_pathway":
                all_nodes["protein"].update(
                    df[src_col][
                        df[src_col].str.strip() != ""
                    ]
                )

                all_nodes["pathway"].update(
                    df[dst_col][
                        df[dst_col].str.strip() != ""
                    ]
                )

        stats.append({
            "table": name,
            "rows": len(df),
            "source_column": source_col or "",
            "source_coverage_pct": round(
                provenance_pct,
                2
            ),
            "endpoint_missing_rows": int(
                endpoint_missing
            ),
            "duplicate_edge_count": int(
                duplicate_count
            ),
        })

    # ---------------------------------------------------------
    # Total graph statistics
    # ---------------------------------------------------------
    total_edges = sum(
        len(df)
        for df in tables.values()
    )

    total_nodes = sum(
        len(values)
        for values in all_nodes.values()
    )

    # ---------------------------------------------------------
    # Drug-disease label imbalance
    # ---------------------------------------------------------
    label_counts = {}

    if "drug_disease" in tables:

        df = tables["drug_disease"]

        rel_col = find_column(
            df,
            ["relation"]
        )

        if rel_col:
            label_counts = (
                df[rel_col]
                .replace("", pd.NA)
                .dropna()
                .value_counts()
                .to_dict()
            )

    # ---------------------------------------------------------
    # High-degree nodes
    # ---------------------------------------------------------
    high_degree_sections = {}

    for name, df in tables.items():

        schema = schemas[name]

        src_col = find_column(
            df,
            schema["source_node"]
        )

        dst_col = find_column(
            df,
            schema["target_node"]
        )

        if not src_col or not dst_col:
            continue

        combined = pd.concat(
            [
                df[src_col],
                df[dst_col],
            ],
            ignore_index=True,
        )

        combined = combined[
            combined.astype(str).str.strip() != ""
        ]

        degree_counts = (
            combined.astype(str)
            .value_counts()
            .head(10)
        )
        high_degree_sections[name] = degree_counts

    # ---------------------------------------------------------
    # DrugBank → ChEMBL mapping statistics
    # ---------------------------------------------------------
    mapping_file = (
        ROOT
        / "data/raw/chembl/drugbank_to_chembl_validated.csv"
    )

    mapping_count = None

    if mapping_file.exists():

        mapping_df = pd.read_csv(
            mapping_file,
            dtype=str
        ).fillna("")

        db_col = find_column(
            mapping_df,
            ["drugbank_id"]
        )

        ch_col = find_column(
            mapping_df,
            ["chembl_id"]
        )

        if db_col and ch_col:
            mapping_count = mapping_df[
                (mapping_df[db_col].str.strip() != "")
                &
                (mapping_df[ch_col].str.strip() != "")
            ][db_col].nunique()

    # ---------------------------------------------------------
    # Save graph statistics CSV
    # ---------------------------------------------------------
    graph_stats = [
        {
            "metric": "total_nodes",
            "value": total_nodes,
        },
        {
            "metric": "total_edges",
            "value": total_edges,
        },
        {
            "metric": "unique_drugs",
            "value": len(all_nodes["drug"]),
        },
        {
            "metric": "unique_diseases",
            "value": len(all_nodes["disease"]),
        },
        {
            "metric": "unique_genes",
            "value": len(all_nodes["gene"]),
        },
        {
            "metric": "unique_proteins",
            "value": len(all_nodes["protein"]),
        },
        {
            "metric": "unique_pathways",
            "value": len(all_nodes["pathway"]),
        },
    ]

    if mapping_count is not None:
        graph_stats.append({
            "metric": "validated_drugbank_chembl_mappings",
            "value": mapping_count,
        })

    graph_stats.extend(
        relation_rows
    )

    graph_stats_df = pd.DataFrame(
        graph_stats
    )

    OUTPUT_STATS.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    graph_stats_df.to_csv(
        OUTPUT_STATS,
        index=False
    )

    # ---------------------------------------------------------
    # Write Markdown report
    # ---------------------------------------------------------
    report = []

    report.append("# ReMedy Data Quality Report\n")
    report.append(
        "This report was generated from the current canonical "
        "relation layers.\n"
    )

    report.append("## 1. Overall graph statistics\n")

    report.append(
        f"- Total unique nodes: {total_nodes:,}\n"
    )
    report.append(
        f"- Total relation rows: {total_edges:,}\n"
    )
    report.append(
        f"- Drugs: {len(all_nodes['drug']):,}\n"
    )
    report.append(
        f"- Diseases: {len(all_nodes['disease']):,}\n"
    )
    report.append(
        f"- Genes: {len(all_nodes['gene']):,} "
        "(no separate canonical gene relation layer yet)\n"
    )
    report.append(
        f"- Proteins: {len(all_nodes['protein']):,}\n"
    )
    report.append(
        f"- Pathways: {len(all_nodes['pathway']):,}\n"
    )

    report.append("\n## 2. Relation-table quality\n")
    report.append(
        "| Table | Rows | Missing endpoints | Duplicate edges | Provenance coverage |\n"
    )
    report.append(
        "|---|---:|---:|---:|---:|\n"
    )

    for row in stats:

        report.append(
            f"| {row['table']} | "
            f"{row['rows']:,} | "
            f"{row['endpoint_missing_rows']:,} | "
            f"{row['duplicate_edge_count']:,} | "
            f"{row['source_coverage_pct']:.2f}% |\n"
        )

    report.append("\n## 3. Relation distribution\n")

    report.append(
        "| Table | Relation | Edge count |\n"
    )
    report.append(
        "|---|---|---:|\n"
    )

    for row in relation_rows:

        report.append(
            f"| {row['table']} | "
            f"{row['relation']} | "
            f"{row['edge_count']:,} |\n"
        )

    report.append("\n## 4. Missing-value percentage by column\n")

    report.append(
        "| Table | Column | Missing | Missing % |\n"
    )
    report.append(
        "|---|---|---:|---:|\n"
    )

    for _, row in missing_df.iterrows():

        report.append(
            f"| {row['table']} | "
            f"{row['column']} | "
            f"{row['missing_count']:,} | "
            f"{row['missing_percentage']:.2f}% |\n"
        )

    report.append("\n## 5. Drug–disease label distribution\n")

    if label_counts:

        total_labels = sum(
            label_counts.values()
        )

        report.append(
            "| Relation | Count | Percentage |\n"
        )
        report.append(
            "|---|---:|---:|\n"
        )

        for label, count in label_counts.items():

            pct = (
                count / total_labels * 100
                if total_labels
                else 0
            )

            report.append(
                f"| {label} | "
                f"{count:,} | "
                f"{pct:.2f}% |\n"
            )

    else:
        report.append(
            "Drug–disease relation labels could not be detected.\n"
        )

    report.append("\n## 6. High-degree nodes\n")

    for table_name, degree_counts in high_degree_sections.items():

        report.append(
            f"\n### {table_name}\n\n"
        )

        report.append(
            "| Node | Degree |\n"
        )
        report.append(
            "|---|---:|\n"
        )

        for node, degree in degree_counts.items():

            report.append(
                f"| {node} | {degree:,} |\n"
            )

    report.append("\n## 7. Identifier mapping\n")

    if mapping_count is not None:

        report.append(
            f"- Validated DrugBank → ChEMBL mappings: "
            f"{mapping_count:,}\n"
        )

    report.append(
        "- PrimeKG disease-protein → UniProt coverage: "
        "96.50% (from the completed mapping stage)\n"
    )

    report.append(
        "- ChEMBL target → UniProt unresolved mappings: "
        "0 (from the completed mapping stage)\n"
    )

    report.append(
        "- Reactome unmatched human ChEMBL targets: "
        "165\n"
    )

    report.append("\n## 8. Notes\n")

    report.append(
        "- Raw source files were not modified.\n"
    )
    report.append(
        "- Canonical identifiers are used in the current relation layers.\n"
    )
    report.append(
        "- Drug–disease relations remain separated as "
        "indication, off-label use, and contraindication.\n"
    )
    report.append(
        "- Non-human ChEMBL evidence is preserved separately "
        "rather than being treated as human evidence.\n"
    )
    report.append(
        "- A separate canonical gene relation layer has not yet "
        "been created, so genes are reported as 0 rather than inferred "
        "from protein identifiers.\n"
    )

    OUTPUT_REPORT.write_text(
        "".join(report),
        encoding="utf-8"
    )

    # ---------------------------------------------------------
    # Terminal summary
    # ---------------------------------------------------------
    print()
    print("=" * 70)
    print("DATA QUALITY CHECK COMPLETE")
    print("=" * 70)

    print(f"Total unique nodes: {total_nodes:,}")
    print(f"Total relation rows: {total_edges:,}")
    print(f"Drugs: {len(all_nodes['drug']):,}")
    print(f"Diseases: {len(all_nodes['disease']):,}")
    print(f"Genes: {len(all_nodes['gene']):,}")
    print(f"Proteins: {len(all_nodes['protein']):,}")
    print(f"Pathways: {len(all_nodes['pathway']):,}")

    print()
    print(f"Saved: {OUTPUT_REPORT}")
    print(f"Saved: {OUTPUT_STATS}")

    print("=" * 70)


if __name__ == "__main__":
    main()