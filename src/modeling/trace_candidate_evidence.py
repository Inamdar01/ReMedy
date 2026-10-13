
from pathlib import Path
import pandas as pd

# --------------------------------------------------
# 1. PATHS
# --------------------------------------------------

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "processed"
RESULTS_DIR = ROOT / "reports" / "model_evaluation"

RANKING_PATH = RESULTS_DIR / "rgcn_top_100_supported_candidates.csv"

NODES_PATH = DATA_DIR / "nodes.parquet"
EDGES_PATH = DATA_DIR / "edges.parquet"
EVIDENCE_PATH = DATA_DIR / "evidence.parquet"


# --------------------------------------------------
# 2. LOAD DATA
# --------------------------------------------------

print("Loading ReMedy knowledge graph files...")

for path in [NODES_PATH, EDGES_PATH, EVIDENCE_PATH, RANKING_PATH]:
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")

nodes = pd.read_parquet(NODES_PATH)
edges = pd.read_parquet(EDGES_PATH)
evidence = pd.read_parquet(EVIDENCE_PATH)
ranking = pd.read_csv(RANKING_PATH)

print(f"Nodes loaded: {len(nodes):,}")
print(f"Edges loaded: {len(edges):,}")
print(f"Evidence records loaded: {len(evidence):,}")
print(f"Candidate rows loaded: {len(ranking):,}")

# --------------------------------------------------
# 3. VALIDATE COLUMNS AND CREATE RANK IF MISSING
# --------------------------------------------------

required_columns = {"drug_id", "disease_id", "rgcn_score"}

missing_columns = required_columns - set(ranking.columns)

if missing_columns:
    raise ValueError(
        f"Missing columns in candidate CSV: {sorted(missing_columns)}\n"
        f"Available columns: {list(ranking.columns)}"
    )

# Preserve existing ranks; otherwise assign ranks by current row order.
if "rank" not in ranking.columns:
    ranking.insert(0, "rank", range(1, len(ranking) + 1))

ranking["rank"] = pd.to_numeric(ranking["rank"], errors="raise")
ranking["rgcn_score"] = pd.to_numeric(
    ranking["rgcn_score"], errors="raise"
)

# Keep the current ranking order and select the top 20.
selected = ranking.head(20).copy()

if selected.empty:
    raise ValueError("The candidate ranking file contains no candidates.")

# --------------------------------------------------
# 4. CREATE NODE LABEL LOOKUP
# --------------------------------------------------

node_labels = nodes.set_index("node_id")["label"].to_dict()


def readable(node_id):
    """Return a readable label for a node ID."""
    return node_labels.get(node_id, "Unknown")


# --------------------------------------------------
# 5. INDEX EVIDENCE BY EDGE ID
# --------------------------------------------------

evidence_by_edge = {
    edge_id: group
    for edge_id, group in evidence.groupby("edge_id", sort=False)
}


def edge_evidence(edge_row):
    """Return up to five provenance records for an edge."""
    records = evidence_by_edge.get(edge_row.edge_id)

    if records is None or records.empty:
        return []

    output = []

    for record in records.head(5).itertuples(index=False):
        output.append({
            "evidence_source": getattr(record, "source", None),
            "publication": getattr(record, "publication", None),
            "assay": getattr(record, "assay", None),
            "organism": getattr(record, "organism", None),
            "cell_line": getattr(record, "cell_line", None),
            "activity_value": getattr(record, "activity_value", None),
            "context": getattr(record, "context", None),
        })

    return output


# --------------------------------------------------
# 6. SELECT KNOWLEDGE GRAPH RELATIONS
# --------------------------------------------------

targets = edges[edges["relation"] == "targets"]
pathways = edges[edges["relation"] == "participates_in"]
disease_associations = edges[
    edges["relation"] == "associated_with"
]

# Drug -> target proteins
drug_targets = {
    drug_id: set(group["target_id"])
    for drug_id, group in targets.groupby("source_id")
}

# Protein -> pathways
protein_pathways = {
    protein_id: set(group["target_id"])
    for protein_id, group in pathways.groupby("source_id")
}

# Disease -> associated proteins
disease_proteins = {
    disease_id: set(group["target_id"])
    for disease_id, group in disease_associations.groupby("source_id")
}

# Lookup for canonical edges
edge_lookup = {
    (row.source_id, row.relation, row.target_id): row
    for row in edges.itertuples(index=False)
}

# --------------------------------------------------
# 7. TRACE TOP 20 CANDIDATES
# --------------------------------------------------

trace_rows = []
evidence_rows = []

print("\nTracing evidence for the top 20 candidates...")

for candidate in selected.to_dict(orient="records"):

    # Dictionary access avoids the Pandas .rank attribute conflict.
    rank = candidate["rank"]
    drug = candidate["drug_id"]
    disease = candidate["disease_id"]
    score = candidate["rgcn_score"]

    drug_protein_set = drug_targets.get(drug, set())
    disease_protein_set = disease_proteins.get(disease, set())

    # Proteins directly linked to both drug and disease
    shared_proteins = drug_protein_set & disease_protein_set

    # Pathways linked to the drug's target proteins
    drug_pathway_set = set()

    for protein in drug_protein_set:
        drug_pathway_set.update(
            protein_pathways.get(protein, set())
        )

    # Pathways linked to the disease-associated proteins
    disease_pathway_set = set()

    for protein in disease_protein_set:
        disease_pathway_set.update(
            protein_pathways.get(protein, set())
        )

    # Common pathways in the two graph neighborhoods
    shared_pathways = drug_pathway_set & disease_pathway_set

    # Save one summary row per candidate
    trace_rows.append({
        "rank": rank,
        "drug_id": drug,
        "drug_name": readable(drug),
        "disease_id": disease,
        "disease_name": readable(disease),
        "rgcn_score": score,
        "drug_target_protein_count": len(drug_protein_set),
        "disease_associated_protein_count": len(disease_protein_set),
        "shared_protein_count": len(shared_proteins),
        "shared_pathway_count": len(shared_pathways),
        "shared_protein_ids": ";".join(sorted(shared_proteins)),
        "shared_protein_names": ";".join(
            readable(protein)
            for protein in sorted(shared_proteins)
        ),
        "shared_pathway_ids": ";".join(sorted(shared_pathways)),
        "shared_pathway_names": ";".join(
            readable(pathway)
            for pathway in sorted(shared_pathways)
        ),
        "interpretation": (
            "Shared protein and/or pathway links found in KG"
            if shared_proteins or shared_pathways
            else "No shared protein or pathway links found"
        ),
    })

    # ----------------------------------------------
    # A. Drug -> target protein evidence
    # ----------------------------------------------

    for protein in sorted(drug_protein_set):

        edge_row = edge_lookup.get(
            (drug, "targets", protein)
        )

        if edge_row is None:
            continue

        for record in edge_evidence(edge_row):
            evidence_rows.append({
                "rank": rank,
                "drug_id": drug,
                "drug_name": readable(drug),
                "disease_id": disease,
                "disease_name": readable(disease),
                "rgcn_score": score,
                "edge_id": edge_row.edge_id,
                "source_id": drug,
                "relation": "targets",
                "target_id": protein,
                "source_name": readable(drug),
                "target_name": readable(protein),
                **record,
            })

    # ----------------------------------------------
    # B. Disease -> associated protein evidence
    # ----------------------------------------------

    for protein in sorted(disease_protein_set):

        edge_row = edge_lookup.get(
            (disease, "associated_with", protein)
        )

        if edge_row is None:
            continue

        for record in edge_evidence(edge_row):
            evidence_rows.append({
                "rank": rank,
                "drug_id": drug,
                "drug_name": readable(drug),
                "disease_id": disease,
                "disease_name": readable(disease),
                "rgcn_score": score,
                "edge_id": edge_row.edge_id,
                "source_id": disease,
                "relation": "associated_with",
                "target_id": protein,
                "source_name": readable(disease),
                "target_name": readable(protein),
                **record,
            })

    # ----------------------------------------------
    # C. Protein -> shared pathway evidence
    # ----------------------------------------------

    relevant_proteins = drug_protein_set | disease_protein_set

    for protein in sorted(relevant_proteins):

        protein_pathway_set = protein_pathways.get(
            protein, set()
        )

        for pathway in sorted(protein_pathway_set):

            if pathway not in shared_pathways:
                continue

            edge_row = edge_lookup.get(
                (protein, "participates_in", pathway)
            )

            if edge_row is None:
                continue

            for record in edge_evidence(edge_row):
                evidence_rows.append({
                    "rank": rank,
                    "drug_id": drug,
                    "drug_name": readable(drug),
                    "disease_id": disease,
                    "disease_name": readable(disease),
                    "rgcn_score": score,
                    "edge_id": edge_row.edge_id,
                    "source_id": protein,
                    "relation": "participates_in",
                    "target_id": pathway,
                    "source_name": readable(protein),
                    "target_name": readable(pathway),
                    **record,
                })


# --------------------------------------------------
# 8. SAVE OUTPUT FILES
# --------------------------------------------------

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

trace_df = pd.DataFrame(trace_rows)
evidence_df = pd.DataFrame(evidence_rows)

trace_path = RESULTS_DIR / "top20_candidate_evidence_traces.csv"
evidence_path = RESULTS_DIR / "top20_candidate_edge_evidence.csv"

trace_df.to_csv(trace_path, index=False)
evidence_df.to_csv(evidence_path, index=False)


# --------------------------------------------------
# 9. PRINT RESULTS
# --------------------------------------------------

print("\n" + "=" * 75)
print("REMedy: TOP-20 CANDIDATE EVIDENCE TRACE")
print("=" * 75)

print("\nCandidate summaries:")

summary_columns = [
    "rank",
    "drug_name",
    "disease_name",
    "rgcn_score",
    "shared_protein_count",
    "shared_pathway_count",
    "interpretation",
]

print(
    trace_df[summary_columns].to_string(index=False)
)

print("\nTotal traceable evidence records:", len(evidence_df))

if not evidence_df.empty:

    print("\nEvidence records by relationship:")
    print(
        evidence_df["relation"]
        .value_counts()
        .to_string()
    )

    print("\nEvidence records by source database:")
    print(
        evidence_df["evidence_source"]
        .fillna("Unknown")
        .value_counts()
        .to_string()
    )

else:
    print(
        "No evidence records were returned. "
        "Review edge IDs and evidence-table contents."
    )

print("\nOutput files saved successfully:")
print(trace_path)
print(evidence_path)

print("\nEvidence tracing completed.")

print(
    "\nScientific caution: R-GCN scores are ranking scores, "
    "not calibrated probabilities. Shared proteins and pathways "
    "show graph connections, not proof of a causal mechanism "
    "or clinical effectiveness."
)
