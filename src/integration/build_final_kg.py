from pathlib import Path
import hashlib

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]

# =========================================================
# INPUT FILES
# =========================================================

DRUG_DISEASE_FILE = (
    ROOT / "data/interim/drug_disease_relations_canonical.csv"
)

DRUG_PROTEIN_FILE = (
    ROOT / "data/interim/drug_protein_relations_canonical.csv"
)

PRIMEKG_DRUG_PROTEIN_FILE = (
    ROOT / "data/interim/primekg_drug_protein_relations_canonical.csv"
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

# =========================================================
# OUTPUT FILES
# =========================================================

PROCESSED_DIR = ROOT / "data/processed"

NODES_FILE = PROCESSED_DIR / "nodes.parquet"
EDGES_FILE = PROCESSED_DIR / "edges.parquet"
EVIDENCE_FILE = PROCESSED_DIR / "evidence.parquet"

REPORT_FILE = (
    ROOT / "reports/final_kg_build_report.md"
)


# =========================================================
# HELPERS
# =========================================================

def clean_string(value):
    if pd.isna(value):
        return ""

    return str(value).strip()


def first_nonempty(*values):

    for value in values:

        value = clean_string(value)

        if value:
            return value

    return ""


def join_unique(values):

    cleaned = sorted(
        {
            clean_string(value)
            for value in values
            if clean_string(value)
        }
    )

    return " | ".join(cleaned)


def make_edge_id(
    source_id,
    relation,
    target_id,
):

    key = (
        f"{source_id}|"
        f"{relation}|"
        f"{target_id}"
    )

    digest = hashlib.sha1(
        key.encode("utf-8")
    ).hexdigest()[:16]

    return f"E_{digest}"


# =========================================================
# NODE BUILDER
# =========================================================

def add_nodes(
    records,
    node_ids,
    node_type,
    labels,
    source_database,
    priority=1,
):

    for node_id, label in zip(
        node_ids,
        labels,
    ):

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
    primekg_drug_protein,
    disease_protein,
    protein_pathway,
):

    records = []

    # -----------------------------------------------------
    # Drugs
    # -----------------------------------------------------

    add_nodes(
        records,
        drug_disease["chembl_id"],
        "drug",
        drug_disease["chembl_pref_name"],
        "ChEMBL|PrimeKG",
        1,
    )

    add_nodes(
        records,
        drug_protein["chembl_id"],
        "drug",
        drug_protein["drug_name"],
        "ChEMBL",
        2,
    )

    add_nodes(
        records,
        primekg_drug_protein["chembl_id"],
        "drug",
        primekg_drug_protein["chembl_pref_name"],
        "ChEMBL|PrimeKG",
        1,
    )

    # -----------------------------------------------------
    # Diseases
    # -----------------------------------------------------

    add_nodes(
        records,
        drug_disease["mondo_id"],
        "disease",
        drug_disease["mondo_name"],
        "MONDO|PrimeKG",
        1,
    )

    add_nodes(
        records,
        disease_protein["mondo_id"],
        "disease",
        disease_protein["mondo_name"],
        "MONDO|PrimeKG",
        1,
    )

    # -----------------------------------------------------
    # Proteins
    # -----------------------------------------------------

    add_nodes(
        records,
        disease_protein["uniprot_accession"],
        "protein",
        disease_protein["hgnc_symbol"],
        "UniProt|HGNC|PrimeKG",
        1,
    )

    add_nodes(
        records,
        disease_protein["uniprot_accession"],
        "protein",
        disease_protein["hgnc_name"],
        "UniProt|HGNC|PrimeKG",
        2,
    )

    add_nodes(
        records,
        drug_protein["uniprot_accession"],
        "protein",
        drug_protein["target_name"],
        "UniProt|ChEMBL",
        2,
    )

    add_nodes(
        records,
        primekg_drug_protein["uniprot_accession"],
        "protein",
        primekg_drug_protein["hgnc_symbol"],
        "UniProt|HGNC|PrimeKG",
        1,
    )

    add_nodes(
        records,
        protein_pathway["uniprot_accession"],
        "protein",
        protein_pathway["pathway_name"],
        "UniProt|Reactome",
        3,
    )

    # -----------------------------------------------------
    # Pathways
    # -----------------------------------------------------

    add_nodes(
        records,
        protein_pathway["reactome_pathway_id"],
        "pathway",
        protein_pathway["pathway_name"],
        "Reactome",
        1,
    )

    nodes_raw = pd.DataFrame(records)

    if nodes_raw.empty:
        raise ValueError(
            "No nodes generated."
        )

    nodes_raw = nodes_raw.sort_values(
        [
            "node_id",
            "priority",
            "label",
        ]
    )

    # One canonical node per ID.
    # Node IDs are globally unique across node types because
    # the identifier namespaces are distinct.
    nodes = (
        nodes_raw
        .groupby(
            "node_id",
            as_index=False,
        )
        .agg(
            node_type=(
                "node_type",
                "first",
            ),
            label=(
                "label",
                lambda s: next(
                    (
                        x
                        for x in s
                        if clean_string(x)
                    ),
                    "",
                ),
            ),
            source_database=(
                "source_database",
                join_unique,
            ),
        )
    )

    nodes = nodes[
        [
            "node_id",
            "node_type",
            "label",
            "source_database",
        ]
    ]

    return nodes


# =========================================================
# EDGE BUILDER
# =========================================================

def build_edges(
    drug_disease,
    drug_protein,
    primekg_drug_protein,
    disease_protein,
    protein_pathway,
):

    frames = []

    # -----------------------------------------------------
    # Drug → Disease
    # -----------------------------------------------------

    dd = pd.DataFrame(
        {
            "source_id": drug_disease[
                "chembl_id"
            ],
            "relation": drug_disease[
                "relation"
            ],
            "target_id": drug_disease[
                "mondo_id"
            ],
        }
    )

    frames.append(dd)

    # -----------------------------------------------------
    # Drug → Protein from ChEMBL
    # -----------------------------------------------------

    dp_chembl = pd.DataFrame(
        {
            "source_id": drug_protein[
                "chembl_id"
            ],
            "relation": "targets",
            "target_id": drug_protein[
                "uniprot_accession"
            ],
        }
    )

    frames.append(dp_chembl)

    # -----------------------------------------------------
    # Drug → Protein from PrimeKG
    # -----------------------------------------------------

    dp_primekg = pd.DataFrame(
        {
            "source_id": primekg_drug_protein[
                "chembl_id"
            ],
            "relation": "targets",
            "target_id": primekg_drug_protein[
                "uniprot_accession"
            ],
        }
    )

    frames.append(dp_primekg)

    # -----------------------------------------------------
    # Disease → Protein
    # -----------------------------------------------------

    disease_prot = pd.DataFrame(
        {
            "source_id": disease_protein[
                "mondo_id"
            ],
            "relation": "associated_with",
            "target_id": disease_protein[
                "uniprot_accession"
            ],
        }
    )

    frames.append(disease_prot)

    # -----------------------------------------------------
    # Protein → Pathway
    # -----------------------------------------------------

    prot_path = pd.DataFrame(
        {
            "source_id": protein_pathway[
                "uniprot_accession"
            ],
            "relation": "participates_in",
            "target_id": protein_pathway[
                "reactome_pathway_id"
            ],
        }
    )

    frames.append(prot_path)

    all_edges = pd.concat(
        frames,
        ignore_index=True,
    )

    # -----------------------------------------------------
    # Remove duplicate canonical edges.
    #
    # This is crucial for ChEMBL + PrimeKG overlap:
    #
    # Drug → Protein
    #      one canonical edge
    #      multiple evidence sources
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

    all_edges["direction"] = (
        "forward"
    )

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


# =========================================================
# EVIDENCE BUILDER
# =========================================================

def build_evidence(
    edges,
    drug_disease,
    primekg_drug_protein,
    disease_protein,
    protein_pathway,
    chembl_evidence,
    chembl_mapping,
):

    evidence_frames = []

    edge_keys = (
        edges["source_id"]
        + "|"
        + edges["relation"]
        + "|"
        + edges["target_id"]
    )

    edge_lookup = dict(
        zip(
            edge_keys,
            edges["edge_id"],
        )
    )

    # =====================================================
    # 1. Drug → Disease : PrimeKG
    # =====================================================

    rows = []

    for _, row in drug_disease.iterrows():

        key = (
            f"{row['chembl_id']}|"
            f"{row['relation']}|"
            f"{row['mondo_id']}"
        )

        edge_id = edge_lookup.get(
            key
        )

        if not edge_id:
            continue

        rows.append(
            {
                "edge_id": edge_id,
                "source": "PrimeKG",
                "publication": "unknown",
                "assay": "unknown",
                "organism": "unknown",
                "cell_line": "unknown",
                "activity_value": "unknown",
                "context": "; ".join(
                    [
                        f"drugbank_id={row['drugbank_id']}",
                        f"primekg_disease_id={row['primekg_disease_id']}",
                        f"disease_source={row['disease_source']}",
                        f"mapping_status={row['mapping_status']}",
                        f"mapping_reason={row['mapping_reason']}",
                    ]
                ),
            }
        )

    if rows:
        evidence_frames.append(
            pd.DataFrame(rows)
        )

    # =====================================================
    # 2. Drug → Protein : PrimeKG
    # =====================================================

    rows = []

    for _, row in primekg_drug_protein.iterrows():

        key = (
            f"{row['chembl_id']}|"
            f"targets|"
            f"{row['uniprot_accession']}"
        )

        edge_id = edge_lookup.get(
            key
        )

        if not edge_id:
            continue

        rows.append(
            {
                "edge_id": edge_id,
                "source": "PrimeKG",
                "publication": "unknown",
                "assay": "unknown",
                "organism": "unknown",
                "cell_line": "unknown",
                "activity_value": "unknown",
                "context": "; ".join(
                    [
                        f"drugbank_id={row['drugbank_id']}",
                        f"entrez_id={row['entrez_id']}",
                        f"protein_gene_symbol={row['protein_gene_symbol']}",
                        f"hgnc_id={row['hgnc_id']}",
                        f"hgnc_symbol={row['hgnc_symbol']}",
                        f"hgnc_name={row['hgnc_name']}",
                        f"drug_source={row['drug_source']}",
                        f"protein_source={row['protein_source']}",
                    ]
                ),
            }
        )

    if rows:
        evidence_frames.append(
            pd.DataFrame(rows)
        )

    # =====================================================
    # 3. Disease → Protein : PrimeKG
    # =====================================================

    rows = []

    for _, row in disease_protein.iterrows():

        key = (
            f"{row['mondo_id']}|"
            f"associated_with|"
            f"{row['uniprot_accession']}"
        )

        edge_id = edge_lookup.get(
            key
        )

        if not edge_id:
            continue

        rows.append(
            {
                "edge_id": edge_id,
                "source": "PrimeKG",
                "publication": "unknown",
                "assay": "unknown",
                "organism": "unknown",
                "cell_line": "unknown",
                "activity_value": "unknown",
                "context": "; ".join(
                    [
                        f"disease_primekg_id={row['disease_primekg_id']}",
                        f"disease_source={row['disease_source']}",
                        f"protein_entrez_id={row['protein_entrez_id']}",
                        f"protein_gene_symbol={row['protein_gene_symbol']}",
                        f"protein_source={row['protein_source']}",
                        f"hgnc_id={row['hgnc_id']}",
                        f"hgnc_symbol={row['hgnc_symbol']}",
                    ]
                ),
            }
        )

    if rows:
        evidence_frames.append(
            pd.DataFrame(rows)
        )

    # =====================================================
    # 4. Protein → Pathway : Reactome
    # =====================================================

    rows = []

    for _, row in protein_pathway.iterrows():

        key = (
            f"{row['uniprot_accession']}|"
            f"participates_in|"
            f"{row['reactome_pathway_id']}"
        )

        edge_id = edge_lookup.get(
            key
        )

        if not edge_id:
            continue

        rows.append(
            {
                "edge_id": edge_id,
                "source": "Reactome",
                "publication": "unknown",
                "assay": "unknown",
                "organism": "Homo sapiens",
                "cell_line": "unknown",
                "activity_value": "unknown",
                "context": "; ".join(
                    [
                        f"pathway_name={row['pathway_name']}",
                        f"evidence_count={row['evidence_count']}",
                        f"evidence_codes={row['evidence_codes']}",
                    ]
                ),
            }
        )

    if rows:
        evidence_frames.append(
            pd.DataFrame(rows)
        )

    # =====================================================
    # 5. Drug → Protein : ChEMBL assay evidence
    # =====================================================

    print(
        "\nMapping ChEMBL evidence..."
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

    mapped_ev = chembl_evidence.merge(
        mapping,
        on="target_id",
        how="inner",
    )

    print(
        f"ChEMBL evidence rows before mapping: "
        f"{len(chembl_evidence):,}"
    )

    print(
        f"ChEMBL evidence rows after target mapping: "
        f"{len(mapped_ev):,}"
    )

    mapped_ev = mapped_ev.drop_duplicates(
        subset=[
            "activity_id",
            "uniprot_accession",
        ]
    )

    mapped_ev["edge_key"] = (
        mapped_ev["chembl_id"]
        + "|targets|"
        + mapped_ev["uniprot_accession"]
    )

    mapped_ev = mapped_ev[
        mapped_ev["edge_key"].isin(
            set(edge_lookup.keys())
        )
    ].copy()

    print(
        f"ChEMBL evidence linked to final edges: "
        f"{len(mapped_ev):,}"
    )

    rows = []

    for _, row in mapped_ev.iterrows():

        edge_id = edge_lookup.get(
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
                if clean_string(row["document_id"])
                else ""
            ),
            "unknown",
        )

        assay = first_nonempty(
            (
                f"Assay:{row['assay_id']}"
                if clean_string(row["assay_id"])
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

        activity_value = first_nonempty(
            " ".join(
                x
                for x in [
                    clean_string(row["standard_type"]),
                    clean_string(row["standard_relation"]),
                    clean_string(row["standard_value"]),
                    clean_string(row["standard_units"]),
                ]
                if x
            ),
            "unknown",
        )

        context = "; ".join(
            [
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
        )

        rows.append(
            {
                "edge_id": edge_id,
                "source": "ChEMBL",
                "publication": publication,
                "assay": assay,
                "organism": organism,
                "cell_line": cell_line,
                "activity_value": activity_value,
                "context": context,
            }
        )

    if rows:
        evidence_frames.append(
            pd.DataFrame(rows)
        )

    # =====================================================
    # Combine all evidence
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


# =========================================================
# MAIN
# =========================================================

def main():

    print("=" * 70)
    print("REMEDY - BUILD FINAL DERIVED KNOWLEDGE GRAPH")
    print("=" * 70)

    # -----------------------------------------------------
    # Load
    # -----------------------------------------------------

    print(
        "\nLoading canonical relation tables..."
    )

    drug_disease = pd.read_csv(
        DRUG_DISEASE_FILE,
        dtype=str,
    ).fillna("")

    drug_protein = pd.read_csv(
        DRUG_PROTEIN_FILE,
        dtype=str,
    ).fillna("")

    primekg_drug_protein = pd.read_csv(
        PRIMEKG_DRUG_PROTEIN_FILE,
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
        f"Drug-Disease: "
        f"{len(drug_disease):,}"
    )

    print(
        f"ChEMBL Drug-Protein: "
        f"{len(drug_protein):,}"
    )

    print(
        f"PrimeKG Drug-Protein: "
        f"{len(primekg_drug_protein):,}"
    )

    print(
        f"Disease-Protein: "
        f"{len(disease_protein):,}"
    )

    print(
        f"Protein-Pathway: "
        f"{len(protein_pathway):,}"
    )

    print(
        f"ChEMBL evidence: "
        f"{len(chembl_evidence):,}"
    )

    # -----------------------------------------------------
    # Nodes
    # -----------------------------------------------------

    print(
        "\nBuilding nodes..."
    )

    nodes = build_nodes(
        drug_disease,
        drug_protein,
        primekg_drug_protein,
        disease_protein,
        protein_pathway,
    )

    print(
        f"Unique nodes: "
        f"{len(nodes):,}"
    )

    # -----------------------------------------------------
    # Edges
    # -----------------------------------------------------

    print(
        "\nBuilding edges..."
    )

    edges = build_edges(
        drug_disease,
        drug_protein,
        primekg_drug_protein,
        disease_protein,
        protein_pathway,
    )

    print(
        f"Unique canonical edges: "
        f"{len(edges):,}"
    )

    # -----------------------------------------------------
    # Evidence
    # -----------------------------------------------------

    print(
        "\nBuilding evidence..."
    )

    evidence = build_evidence(
        edges,
        drug_disease,
        primekg_drug_protein,
        disease_protein,
        protein_pathway,
        chembl_evidence,
        chembl_mapping,
    )

    print(
        f"Evidence rows: "
        f"{len(evidence):,}"
    )

    # -----------------------------------------------------
    # Edge evidence counts
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
    # Deterministic sorting
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
    # Save
    # -----------------------------------------------------

    PROCESSED_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

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

    evidence_sources = (
        evidence["source"]
        .value_counts()
        .sort_index()
    )

    # Distinct source count per edge.
    edge_source_counts = (
        evidence
        .groupby("edge_id")["source"]
        .nunique()
    )

    multi_source_edges = int(
        (
            edge_source_counts
            >= 2
        ).sum()
    )

    single_source_edges = int(
        (
            edge_source_counts
            == 1
        ).sum()
    )

    edges_with_evidence = int(
        (
            edges["evidence_count"]
            > 0
        ).sum()
    )

    # Specifically Drug → Protein multi-source coverage.
    dp_edge_ids = set(
        edges.loc[
            edges["relation"] == "targets",
            "edge_id",
        ]
    )

    dp_source_counts = (
        evidence[
            evidence["edge_id"]
            .isin(dp_edge_ids)
        ]
        .groupby("edge_id")["source"]
        .nunique()
    )

    dp_multi_source = int(
        (
            dp_source_counts
            >= 2
        ).sum()
    )

    dp_single_source = int(
        (
            dp_source_counts
            == 1
        ).sum()
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
- Edges with >=1 evidence record: {edges_with_evidence:,}
- Edges without evidence: {len(edges) - edges_with_evidence:,}

### Evidence source distribution

{evidence_sources.to_frame("evidence_rows").to_markdown()}

## Multi-source coverage

- Single-source edges: {single_source_edges:,}
- Multi-source edges: {multi_source_edges:,}

### Drug → Protein specifically

- Drug → Protein edges: {len(dp_edge_ids):,}
- Single-source Drug → Protein edges: {dp_single_source:,}
- Multi-source Drug → Protein edges: {dp_multi_source:,}
- Multi-source Drug → Protein coverage: {
    (dp_multi_source / len(dp_edge_ids) * 100)
    if len(dp_edge_ids)
    else 0
:.2f}%

## Canonical identifiers

- Drug: ChEMBL ID
- Disease: MONDO ID
- Protein: UniProt accession
- Pathway: Reactome ID

## Relation semantics

- Drug → Disease:
  - indication
  - off-label use
  - contraindication
- Drug → Protein:
  - targets
- Disease → Protein:
  - associated_with
- Protein → Pathway:
  - participates_in

Mechanistic `inhibits` or `activates` relations are not inferred
from activity values.

## Multi-source Drug → Protein policy

ChEMBL and PrimeKG records that refer to the same canonical
ChEMBL → UniProt edge are represented as one edge with
multiple evidence sources.

Different source records are retained in `evidence.parquet`.

## Evidence policy

- ChEMBL assay-level evidence is preserved.
- PrimeKG provenance is preserved.
- Reactome provenance and evidence codes are preserved.
- Assay organism is preserved as recorded.
- Missing context is represented as `unknown`.
- Non-human assay evidence is not relabeled as human evidence.

## Scope note

No separate gene relation layer is constructed because the
current canonical disease-protein layer combines gene/protein
information. Gene edges are not inferred.

## Integrity

Edge IDs are deterministic hashes of:

`source_id | relation | target_id`

The evidence table references those edge IDs.
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
    # Terminal
    # -----------------------------------------------------

    print()
    print("=" * 70)
    print("FINAL REMEDY KG BUILD COMPLETE")
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
        f"Multi-source edges: "
        f"{multi_source_edges:,}"
    )

    print(
        f"Multi-source Drug → Protein edges: "
        f"{dp_multi_source:,}"
    )

    print()
    print(
        f"Saved: {NODES_FILE}"
    )

    print(
        f"Saved: {EDGES_FILE}"
    )

    print(
        f"Saved: {EVIDENCE_FILE}"
    )

    print(
        f"Saved: {REPORT_FILE}"
    )

    print("=" * 70)


if __name__ == "__main__":
    main()