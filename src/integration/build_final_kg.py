from pathlib import Path
import hashlib
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

# ---------------------------------------------------------
# INPUT FILES
# ---------------------------------------------------------

DRUG_DISEASE_FILE = (
    ROOT / "data/interim/drug_disease_relations_canonical.csv"
)

DRUG_PROTEIN_FILE = (
    ROOT / "data/interim/drug_protein_relations_canonical.csv"
)

DISEASE_PROTEIN_FILE = (
    ROOT / "data/interim/disease_uniprot_relations_canonical.csv"
)

PROTEIN_PATHWAY_FILE = (
    ROOT / "data/interim/protein_pathway_relations_canonical.csv"
)

CHEMBL_EVIDENCE_FILE = (
    ROOT / "data/interim/chembl_evidence_clean.csv"
)

CHEMBL_MAPPING_FILE = (
    ROOT / "data/interim/chembl_target_uniprot_mapping_extended.csv"
)

# ---------------------------------------------------------
# OUTPUT FILES
# ---------------------------------------------------------

PROCESSED_DIR = ROOT / "data/processed"

NODES_FILE = (
    PROCESSED_DIR / "nodes.parquet"
)

EDGES_FILE = (
    PROCESSED_DIR / "edges.parquet"
)

EVIDENCE_FILE = (
    PROCESSED_DIR / "evidence.parquet"
)

REPORT_FILE = (
    ROOT / "reports/final_kg_build_report.md"
)


# ---------------------------------------------------------
# HELPERS
# ---------------------------------------------------------

def clean_string(value):
    """Normalize a value to a clean string."""
    if pd.isna(value):
        return ""

    return str(value).strip()


def first_nonempty(*values):
    """Return the first non-empty value."""
    for value in values:
        value = clean_string(value)

        if value:
            return value

    return ""


def join_unique(values):
    """
    Join unique non-empty values using '|'.
    """
    cleaned = sorted(
        {
            clean_string(value)
            for value in values
            if clean_string(value)
        }
    )

    return " | ".join(cleaned)


def make_edge_id(source_id, relation, target_id):
    """
    Generate a stable edge ID from canonical endpoints.
    """
    key = (
        f"{source_id}|{relation}|{target_id}"
    )

    digest = hashlib.sha1(
        key.encode("utf-8")
    ).hexdigest()[:16]

    return f"E_{digest}"


# ---------------------------------------------------------
# NODE REGISTRY
# ---------------------------------------------------------

def add_node_records(
    records,
    ids,
    node_type,
    labels,
    source_database,
    priority=1,
):
    """
    Add candidate node records to the registry.
    """

    for node_id, label in zip(ids, labels):

        node_id = clean_string(node_id)
        label = clean_string(label)

        if not node_id:
            continue

        records.append(
            {
                "node_id": node_id,
                "node_type": node_type,
                "label": label,
                "source_database": source_database,
                "priority": priority,
            }
        )


def build_nodes(
    drug_disease,
    drug_protein,
    disease_protein,
    protein_pathway,
):

    records = []

    # -----------------------------------------------------
    # Drugs
    # -----------------------------------------------------

    add_node_records(
        records,
        drug_disease["chembl_id"],
        "drug",
        drug_disease["chembl_pref_name"],
        "ChEMBL|PrimeKG",
        priority=1,
    )

    add_node_records(
        records,
        drug_protein["chembl_id"],
        "drug",
        drug_protein["drug_name"],
        "ChEMBL",
        priority=2,
    )

    # -----------------------------------------------------
    # Diseases
    # -----------------------------------------------------

    add_node_records(
        records,
        drug_disease["mondo_id"],
        "disease",
        drug_disease["mondo_name"],
        "MONDO|PrimeKG",
        priority=1,
    )

    add_node_records(
        records,
        disease_protein["mondo_id"],
        "disease",
        disease_protein["mondo_name"],
        "MONDO|PrimeKG",
        priority=1,
    )

    # -----------------------------------------------------
    # Proteins
    # -----------------------------------------------------

    add_node_records(
        records,
        disease_protein["uniprot_accession"],
        "protein",
        disease_protein["hgnc_symbol"],
        "UniProt|HGNC|PrimeKG",
        priority=1,
    )

    add_node_records(
        records,
        disease_protein["uniprot_accession"],
        "protein",
        disease_protein["hgnc_name"],
        "UniProt|HGNC|PrimeKG",
        priority=2,
    )

    add_node_records(
        records,
        drug_protein["uniprot_accession"],
        "protein",
        drug_protein["target_name"],
        "UniProt|ChEMBL",
        priority=2,
    )

    add_node_records(
        records,
        drug_protein["uniprot_accession"],
        "protein",
        drug_protein["uniprot_accession"],
        "UniProt|ChEMBL",
        priority=3,
    )

    add_node_records(
        records,
        protein_pathway["uniprot_accession"],
        "protein",
        protein_pathway["uniprot_accession"],
        "UniProt|Reactome",
        priority=3,
    )

    # -----------------------------------------------------
    # Pathways
    # -----------------------------------------------------

    add_node_records(
        records,
        protein_pathway["reactome_pathway_id"],
        "pathway",
        protein_pathway["pathway_name"],
        "Reactome",
        priority=1,
    )

    nodes_raw = pd.DataFrame(records)

    if nodes_raw.empty:
        raise ValueError(
            "No nodes were generated."
        )

    # -----------------------------------------------------
    # Choose best label by priority
    # -----------------------------------------------------

    nodes_raw = nodes_raw.sort_values(
        [
            "node_id",
            "priority",
            "label",
        ]
    )

    labels = (
        nodes_raw[
            nodes_raw["label"].str.strip() != ""
        ]
        .groupby("node_id")["label"]
        .first()
    )

    sources = (
        nodes_raw
        .groupby("node_id")[
            "source_database"
        ]
        .apply(join_unique)
    )

    types = (
        nodes_raw
        .groupby("node_id")[
            "node_type"
        ]
        .first()
    )

    nodes = pd.DataFrame(
        {
            "node_id": types.index,
            "node_type": types.values,
            "label": [
                labels.get(
                    node_id,
                    node_id,
                )
                for node_id in types.index
            ],
            "source_database": [
                sources.get(
                    node_id,
                    "unknown",
                )
                for node_id in types.index
            ],
        }
    )

    # Ensure exact schema
    nodes = nodes[
        [
            "node_id",
            "node_type",
            "label",
            "source_database",
        ]
    ]

    return nodes


# ---------------------------------------------------------
# EDGE BUILDER
# ---------------------------------------------------------

def build_edges(
    drug_disease,
    drug_protein,
    disease_protein,
    protein_pathway,
):

    edge_frames = []

    # -----------------------------------------------------
    # Drug → Disease
    # -----------------------------------------------------

    dd = pd.DataFrame(
        {
            "source_id": drug_disease["chembl_id"],
            "relation": drug_disease["relation"],
            "target_id": drug_disease["mondo_id"],
            "direction": "forward",
        }
    )

    dd["_key"] = (
        dd["source_id"]
        + "|"
        + dd["relation"]
        + "|"
        + dd["target_id"]
    )

    dd["_evidence_type"] = "drug_disease"

    edge_frames.append(dd)

    # -----------------------------------------------------
    # Drug → Protein
    # -----------------------------------------------------

    dp = pd.DataFrame(
        {
            "source_id": drug_protein["chembl_id"],
            "relation": "targets",
            "target_id": drug_protein["uniprot_accession"],
            "direction": "forward",
        }
    )

    dp["_key"] = (
        dp["source_id"]
        + "|"
        + dp["relation"]
        + "|"
        + dp["target_id"]
    )

    dp["_evidence_type"] = "drug_protein"

    edge_frames.append(dp)

    # -----------------------------------------------------
    # Disease → Protein
    #
    # Generic association relation is used because the
    # source table represents disease-associated proteins.
    # We do not infer a more specific mechanism.
    # -----------------------------------------------------

    dpp = pd.DataFrame(
        {
            "source_id": disease_protein["mondo_id"],
            "relation": "associated_with",
            "target_id": disease_protein[
                "uniprot_accession"
            ],
            "direction": "forward",
        }
    )

    dpp["_key"] = (
        dpp["source_id"]
        + "|"
        + dpp["relation"]
        + "|"
        + dpp["target_id"]
    )

    dpp["_evidence_type"] = "disease_protein"

    edge_frames.append(dpp)

    # -----------------------------------------------------
    # Protein → Pathway
    # -----------------------------------------------------

    pp = pd.DataFrame(
        {
            "source_id": protein_pathway[
                "uniprot_accession"
            ],
            "relation": "participates_in",
            "target_id": protein_pathway[
                "reactome_pathway_id"
            ],
            "direction": "forward",
        }
    )

    pp["_key"] = (
        pp["source_id"]
        + "|"
        + pp["relation"]
        + "|"
        + pp["target_id"]
    )

    pp["_evidence_type"] = "protein_pathway"

    edge_frames.append(pp)

    all_edges = pd.concat(
        edge_frames,
        ignore_index=True,
    )

    # -----------------------------------------------------
    # Remove any accidental duplicate canonical edges
    # -----------------------------------------------------

    all_edges = (
        all_edges
        .drop_duplicates(
            subset=[
                "source_id",
                "relation",
                "target_id",
            ]
        )
        .copy()
    )

    # Stable edge IDs
    all_edges["edge_id"] = [
        make_edge_id(
            source_id,
            relation,
            target_id,
        )
        for source_id, relation, target_id
        in zip(
            all_edges["source_id"],
            all_edges["relation"],
            all_edges["target_id"],
        )
    ]

    return all_edges


# ---------------------------------------------------------
# EVIDENCE BUILDER
# ---------------------------------------------------------

def build_evidence(
    edges,
    drug_disease,
    disease_protein,
    protein_pathway,
    chembl_evidence,
    chembl_mapping,
):

    evidence_frames = []

    edge_lookup = (
        edges[
            [
                "edge_id",
                "source_id",
                "relation",
                "target_id",
            ]
        ]
        .copy()
    )

    edge_lookup["key"] = (
        edge_lookup["source_id"]
        + "|"
        + edge_lookup["relation"]
        + "|"
        + edge_lookup["target_id"]
    )

    key_to_edge_id = dict(
        zip(
            edge_lookup["key"],
            edge_lookup["edge_id"],
        )
    )

    # =====================================================
    # 1. Drug → Disease evidence
    # =====================================================

    dd_rows = []

    for _, row in drug_disease.iterrows():

        edge_key = (
            f"{row['chembl_id']}|"
            f"{row['relation']}|"
            f"{row['mondo_id']}"
        )

        edge_id = key_to_edge_id.get(
            edge_key
        )

        if not edge_id:
            continue

        context_parts = [
            f"drugbank_id={row['drugbank_id']}",
            f"primekg_disease_id={row['primekg_disease_id']}",
            f"disease_source={row['disease_source']}",
            f"mapping_status={row['mapping_status']}",
            f"mapping_reason={row['mapping_reason']}",
        ]

        dd_rows.append(
            {
                "edge_id": edge_id,
                "source": "PrimeKG",
                "publication": "unknown",
                "assay": "unknown",
                "organism": "unknown",
                "cell_line": "unknown",
                "activity_value": "unknown",
                "context": "; ".join(
                    context_parts
                ),
            }
        )

    if dd_rows:
        evidence_frames.append(
            pd.DataFrame(dd_rows)
        )

    # =====================================================
    # 2. Disease → Protein evidence
    # =====================================================

    dpr_rows = []

    for _, row in disease_protein.iterrows():

        edge_key = (
            f"{row['mondo_id']}|"
            f"associated_with|"
            f"{row['uniprot_accession']}"
        )

        edge_id = key_to_edge_id.get(
            edge_key
        )

        if not edge_id:
            continue

        context_parts = [
            f"disease_primekg_id={row['disease_primekg_id']}",
            f"disease_source={row['disease_source']}",
            f"protein_entrez_id={row['protein_entrez_id']}",
            f"protein_gene_symbol={row['protein_gene_symbol']}",
            f"protein_source={row['protein_source']}",
            f"hgnc_id={row['hgnc_id']}",
            f"hgnc_symbol={row['hgnc_symbol']}",
        ]

        dpr_rows.append(
            {
                "edge_id": edge_id,
                "source": "PrimeKG",
                "publication": "unknown",
                "assay": "unknown",
                "organism": "unknown",
                "cell_line": "unknown",
                "activity_value": "unknown",
                "context": "; ".join(
                    context_parts
                ),
            }
        )

    if dpr_rows:
        evidence_frames.append(
            pd.DataFrame(dpr_rows)
        )

    # =====================================================
    # 3. Protein → Pathway evidence
    # =====================================================

    pp_rows = []

    for _, row in protein_pathway.iterrows():

        edge_key = (
            f"{row['uniprot_accession']}|"
            f"participates_in|"
            f"{row['reactome_pathway_id']}"
        )

        edge_id = key_to_edge_id.get(
            edge_key
        )

        if not edge_id:
            continue

        context_parts = [
            f"pathway_name={row['pathway_name']}",
            f"evidence_count={row['evidence_count']}",
            f"evidence_codes={row['evidence_codes']}",
        ]

        pp_rows.append(
            {
                "edge_id": edge_id,
                "source": "Reactome",
                "publication": "unknown",
                "assay": "unknown",
                "organism": "Homo sapiens",
                "cell_line": "unknown",
                "activity_value": "unknown",
                "context": "; ".join(
                    context_parts
                ),
            }
        )

    if pp_rows:
        evidence_frames.append(
            pd.DataFrame(pp_rows)
        )

    # =====================================================
    # 4. Drug → Protein ChEMBL evidence
    # =====================================================

    print(
        "\nMapping ChEMBL evidence to canonical "
        "Drug → Protein edges..."
    )

    mapping = chembl_mapping.copy()

    mapping["human_target"] = (
        mapping["human_target"]
        .astype(str)
        .str.lower()
        .eq("true")
    )

    mapping = mapping[
        mapping["human_target"]
        &
        mapping["uniprot_accession"]
        .astype(str)
        .str.strip()
        .ne("")
    ][
        [
            "target_id",
            "uniprot_accession",
        ]
    ].drop_duplicates()

    chembl_ev = chembl_evidence.copy()

    chembl_ev["target_id"] = (
        chembl_ev["target_id"]
        .astype(str)
        .str.strip()
    )

    mapped_ev = chembl_ev.merge(
        mapping,
        on="target_id",
        how="inner",
    )

    print(
        f"ChEMBL evidence rows before mapping: "
        f"{len(chembl_ev):,}"
    )

    print(
        f"ChEMBL evidence rows after target mapping: "
        f"{len(mapped_ev):,}"
    )

    # Avoid exact repeated evidence for the same
    # activity → protein mapping.
    mapped_ev = mapped_ev.drop_duplicates(
        subset=[
            "activity_id",
            "uniprot_accession",
        ]
    )

    # -----------------------------------------------------
    # Match only edges that actually exist in the
    # canonical human Drug → Protein layer.
    # -----------------------------------------------------

    mapped_ev["edge_key"] = (
        mapped_ev["chembl_id"]
        + "|targets|"
        + mapped_ev["uniprot_accession"]
    )

    mapped_ev = mapped_ev[
        mapped_ev["edge_key"].isin(
            key_to_edge_id
        )
    ].copy()

    print(
        f"ChEMBL evidence rows linked to final edges: "
        f"{len(mapped_ev):,}"
    )

    chembl_rows = []

    for _, row in mapped_ev.iterrows():

        edge_id = key_to_edge_id.get(
            row["edge_key"]
        )

        if not edge_id:
            continue

        publication = first_nonempty(
            row["publication"],
            (
                f"PMID:{row['pmid']}"
                if clean_string(row["pmid"])
                else ""
            ),
            (
                f"DOI:{row['doi']}"
                if clean_string(row["doi"])
                else ""
            ),
            (
                f"Document:{row['document_id']}"
                if clean_string(
                    row["document_id"]
                )
                else ""
            ),
            "unknown",
        )

        assay = first_nonempty(
            (
                f"Assay:{row['assay_id']}"
                if clean_string(
                    row["assay_id"]
                )
                else ""
            ),
            "unknown",
        )

        organism = first_nonempty(
            row["assay_organism"],
            "unknown",
        )

        cell_line = first_nonempty(
            row["cell_line"],
            "unknown",
        )

        activity_value = "unknown"

        if clean_string(
            row["standard_value"]
        ):

            activity_parts = [
                clean_string(
                    row["standard_type"]
                ),
                clean_string(
                    row["standard_relation"]
                ),
                clean_string(
                    row["standard_value"]
                ),
                clean_string(
                    row["standard_units"]
                ),
            ]

            activity_value = " ".join(
                x
                for x in activity_parts
                if x
            )

        context_parts = [
            f"activity_id={row['activity_id']}",
            f"target_id={row['target_id']}",
            f"target_name={row['target_name']}",
            f"target_type={row['target_type']}",
            f"activity_type_raw={row['activity_type_raw']}",
            f"assay_type={row['assay_type']}",
            f"assay_description={row['assay_description']}",
            f"organism_class={row['organism_class']}",
            f"assay_organism_class={row['assay_organism_class']}",
            f"pchembl_value={row['pchembl_value']}",
            f"data_validity_comment={row['data_validity_comment']}",
            f"value_quality={row['value_quality']}",
            f"target_organism={row['target_organism']}",
            f"organism_source={row['organism_source']}",
        ]

        chembl_rows.append(
            {
                "edge_id": edge_id,
                "source": "ChEMBL",
                "publication": publication,
                "assay": assay,
                "organism": organism,
                "cell_line": cell_line,
                "activity_value": activity_value,
                "context": "; ".join(
                    context_parts
                ),
            }
        )

    if chembl_rows:
        evidence_frames.append(
            pd.DataFrame(chembl_rows)
        )

    # =====================================================
    # Combine evidence
    # =====================================================

    if evidence_frames:

        evidence = pd.concat(
            evidence_frames,
            ignore_index=True,
        )

    else:

        evidence = pd.DataFrame(
            columns=[
                "edge_id",
                "source",
                "publication",
                "assay",
                "organism",
                "cell_line",
                "activity_value",
                "context",
            ]
        )

    # -----------------------------------------------------
    # Final evidence cleanup
    # -----------------------------------------------------

    for column in [
        "source",
        "publication",
        "assay",
        "organism",
        "cell_line",
        "activity_value",
        "context",
    ]:

        evidence[column] = (
            evidence[column]
            .fillna("")
            .astype(str)
            .str.strip()
        )

        evidence.loc[
            evidence[column] == "",
            column,
        ] = "unknown"

    evidence = evidence[
        [
            "edge_id",
            "source",
            "publication",
            "assay",
            "organism",
            "cell_line",
            "activity_value",
            "context",
        ]
    ]

    return evidence


# ---------------------------------------------------------
# MAIN
# ---------------------------------------------------------

def main():

    print("=" * 70)
    print("REMEDY - BUILD FINAL DERIVED KNOWLEDGE GRAPH")
    print("=" * 70)

    # -----------------------------------------------------
    # Load canonical relation tables
    # -----------------------------------------------------

    print("\nLoading canonical relation tables...")

    drug_disease = pd.read_csv(
        DRUG_DISEASE_FILE,
        dtype=str,
    ).fillna("")

    drug_protein = pd.read_csv(
        DRUG_PROTEIN_FILE,
        dtype=str,
    ).fillna("")

    disease_protein = pd.read_csv(
        DISEASE_PROTEIN_FILE,
        dtype=str,
    ).fillna("")

    protein_pathway = pd.read_csv(
        PROTEIN_PATHWAY_FILE,
        dtype=str,
    ).fillna("")

    chembl_evidence = pd.read_csv(
        CHEMBL_EVIDENCE_FILE,
        dtype=str,
    ).fillna("")

    chembl_mapping = pd.read_csv(
        CHEMBL_MAPPING_FILE,
        dtype=str,
    ).fillna("")

    print(
        f"Drug-Disease: {len(drug_disease):,}"
    )

    print(
        f"Drug-Protein: {len(drug_protein):,}"
    )

    print(
        f"Disease-Protein: {len(disease_protein):,}"
    )

    print(
        f"Protein-Pathway: {len(protein_pathway):,}"
    )

    print(
        f"ChEMBL evidence: {len(chembl_evidence):,}"
    )

    # -----------------------------------------------------
    # Build nodes
    # -----------------------------------------------------

    print("\nBuilding nodes...")

    nodes = build_nodes(
        drug_disease,
        drug_protein,
        disease_protein,
        protein_pathway,
    )

    print(
        f"Unique nodes: {len(nodes):,}"
    )

    # -----------------------------------------------------
    # Build edges
    # -----------------------------------------------------

    print("\nBuilding edges...")

    edges = build_edges(
        drug_disease,
        drug_protein,
        disease_protein,
        protein_pathway,
    )

    print(
        f"Unique canonical edges: {len(edges):,}"
    )

    # -----------------------------------------------------
    # Build evidence
    # -----------------------------------------------------

    print("\nBuilding evidence...")

    evidence = build_evidence(
        edges,
        drug_disease,
        disease_protein,
        protein_pathway,
        chembl_evidence,
        chembl_mapping,
    )

    print(
        f"Evidence rows: {len(evidence):,}"
    )

    # -----------------------------------------------------
    # Calculate edge evidence counts
    # -----------------------------------------------------

    evidence_counts = (
        evidence
        .groupby("edge_id")
        .size()
        .rename("evidence_count")
    )

    edges["evidence_count"] = (
        edges["edge_id"]
        .map(evidence_counts)
        .fillna(0)
        .astype(int)
    )

    edges = edges[
        [
            "edge_id",
            "source_id",
            "relation",
            "target_id",
            "direction",
            "evidence_count",
        ]
    ]

    # -----------------------------------------------------
    # Final deterministic sorting
    # -----------------------------------------------------

    nodes = nodes.sort_values(
        [
            "node_type",
            "node_id",
        ]
    ).reset_index(drop=True)

    edges = edges.sort_values(
        "edge_id"
    ).reset_index(drop=True)

    evidence = evidence.sort_values(
        [
            "edge_id",
            "source",
            "activity_value",
        ]
    ).reset_index(drop=True)

    # -----------------------------------------------------
    # Ensure output directory
    # -----------------------------------------------------

    PROCESSED_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    # -----------------------------------------------------
    # Save Parquet files
    # -----------------------------------------------------

    nodes.to_parquet(
        NODES_FILE,
        index=False,
    )

    edges.to_parquet(
        EDGES_FILE,
        index=False,
    )

    evidence.to_parquet(
        EVIDENCE_FILE,
        index=False,
    )

    # -----------------------------------------------------
    # Statistics
    # -----------------------------------------------------

    node_counts = (
        nodes["node_type"]
        .value_counts()
        .sort_index()
    )

    relation_counts = (
        edges["relation"]
        .value_counts()
        .sort_index()
    )

    evidence_source_counts = (
        evidence["source"]
        .value_counts()
        .sort_index()
    )

    edges_with_evidence = int(
        (edges["evidence_count"] > 0)
        .sum()
    )

    edges_without_evidence = int(
        (edges["evidence_count"] == 0)
        .sum()
    )

    # -----------------------------------------------------
    # Report
    # -----------------------------------------------------

    report = f"""# ReMedy Final Knowledge Graph Build Report

## Final files

- `data/processed/nodes.parquet`
- `data/processed/edges.parquet`
- `data/processed/evidence.parquet`

## Nodes

- Total nodes: {len(nodes):,}
- Drugs: {node_counts.get("drug", 0):,}
- Diseases: {node_counts.get("disease", 0):,}
- Proteins: {node_counts.get("protein", 0):,}
- Pathways: {node_counts.get("pathway", 0):,}
- Genes: {node_counts.get("gene", 0):,}

## Edges

- Total canonical edges: {len(edges):,}

### Relation distribution

{relation_counts.to_frame("edge_count").to_markdown()}

## Evidence

- Total evidence rows: {len(evidence):,}
- Edges with >=1 evidence row: {edges_with_evidence:,}
- Edges with 0 evidence rows: {edges_without_evidence:,}

### Evidence source distribution

{evidence_source_counts.to_frame("evidence_rows").to_markdown()}

## Canonical identifier policy

- Drug: ChEMBL ID
- Disease: MONDO ID
- Protein: UniProt accession
- Pathway: Reactome ID

## Relation semantics

- Drug → Disease relations remain separate:
  - indication
  - off-label use
  - contraindication
- Drug → Protein uses the generic `targets` relation.
- Disease → Protein uses the generic `associated_with` relation.
- Protein → Pathway uses `participates_in`.

Mechanistic relations such as `inhibits` or `activates`
were not inferred from activity types.

## Evidence policy

- ChEMBL assay-level evidence is linked to canonical
  human ChEMBL-target → UniProt edges.
- Assay organism is preserved as recorded.
- Missing context is represented as `unknown`.
- PrimeKG provenance is retained in the evidence context.
- Reactome evidence codes are retained in the evidence context.
- Non-human assay evidence is not relabeled as human evidence.

## Important scope note

No separate gene relation layer is constructed because the
current canonical disease-protein source combines gene/protein
information. Gene edges are not inferred from protein records.

## Integrity

The final edge table uses stable deterministic edge IDs
generated from:

`source_id | relation | target_id`

The final evidence table references these edge IDs.
"""

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    # -----------------------------------------------------
    # Terminal summary
    # -----------------------------------------------------

    print()
    print("=" * 70)
    print("FINAL RE MEDY KG BUILD COMPLETE")
    print("=" * 70)

    print(
        f"Nodes: {len(nodes):,}"
    )

    print(
        f"Edges: {len(edges):,}"
    )

    print(
        f"Evidence rows: {len(evidence):,}"
    )

    print(
        f"Edges with evidence: "
        f"{edges_with_evidence:,}"
    )

    print(
        f"Edges without evidence: "
        f"{edges_without_evidence:,}"
    )

    print()
    print(f"Saved: {NODES_FILE}")
    print(f"Saved: {EDGES_FILE}")
    print(f"Saved: {EVIDENCE_FILE}")
    print(f"Saved: {REPORT_FILE}")

    print("=" * 70)


if __name__ == "__main__":
    main()