from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx


ROOT = Path(__file__).resolve().parents[2]

NODES_PATH = ROOT / "data" / "processed" / "nodes.parquet"
EDGES_PATH = ROOT / "data" / "processed" / "edges.parquet"
OUTPUT_DIR = ROOT / "reports" / "visualizations"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ================================================================
# Load frozen ReMedy-KG v1.0 — READ ONLY
# ================================================================

nodes = pd.read_parquet(NODES_PATH)
edges = pd.read_parquet(EDGES_PATH)

node_info = nodes.set_index("node_id").to_dict("index")


# ================================================================
# Extract required relation types
# ================================================================

drug_protein = edges[
    edges["relation"] == "targets"
][
    ["source_id", "target_id", "evidence_count"]
].rename(
    columns={
        "source_id": "drug_id",
        "target_id": "protein_id",
        "evidence_count": "drug_protein_evidence",
    }
)

protein_pathway = edges[
    edges["relation"] == "participates_in"
][
    ["source_id", "target_id", "evidence_count"]
].rename(
    columns={
        "source_id": "protein_id",
        "target_id": "pathway_id",
        "evidence_count": "protein_pathway_evidence",
    }
)

disease_protein = edges[
    edges["relation"] == "associated_with"
][
    ["source_id", "target_id", "evidence_count"]
].rename(
    columns={
        "source_id": "disease_id",
        "target_id": "protein_id",
        "evidence_count": "disease_protein_evidence",
    }
)

drug_disease = edges[
    edges["relation"].isin(["indication", "off-label use"])
][
    ["source_id", "relation", "target_id", "evidence_count"]
].rename(
    columns={
        "source_id": "drug_id",
        "target_id": "disease_id",
        "evidence_count": "drug_disease_evidence",
    }
)


# ================================================================
# Find the exact 4-edge connected pattern
#
# Drug -> Protein -> Pathway
# Drug -> Disease -> Protein
# ================================================================

# First ensure the same protein connects to a pathway
mechanism = drug_protein.merge(
    protein_pathway,
    on="protein_id",
    how="inner",
)

# Add the disease-protein relationship
mechanism = mechanism.merge(
    disease_protein,
    on="protein_id",
    how="inner",
)

# Require that the SAME drug has a treatment-use relationship
# with the SAME disease
connected = mechanism.merge(
    drug_disease,
    on=["drug_id", "disease_id"],
    how="inner",
)


if connected.empty:
    raise RuntimeError(
        "No Drug -> Protein -> Pathway + "
        "Drug -> Disease -> Protein example exists "
        "in ReMedy-KG v1.0."
    )


# ================================================================
# Prefer an indication example and stronger combined evidence
# ================================================================

connected["relation_priority"] = connected["relation"].map(
    {
        "indication": 0,
        "off-label use": 1,
    }
)

connected["combined_evidence"] = (
    connected["drug_protein_evidence"]
    + connected["protein_pathway_evidence"]
    + connected["disease_protein_evidence"]
    + connected["drug_disease_evidence"]
)

connected = connected.sort_values(
    [
        "relation_priority",
        "combined_evidence",
    ],
    ascending=[True, False],
)


example = connected.iloc[0]


drug_id = example["drug_id"]
protein_id = example["protein_id"]
pathway_id = example["pathway_id"]
disease_id = example["disease_id"]
treatment_relation = example["relation"]


# ================================================================
# Get labels
# ================================================================

def get_label(node_id):
    return node_info.get(node_id, {}).get("label", node_id)


drug_label = get_label(drug_id)
protein_label = get_label(protein_id)
pathway_label = get_label(pathway_id)
disease_label = get_label(disease_id)


# ================================================================
# Save exact selected relationships
# ================================================================

example_rows = [
    {
        "source_id": drug_id,
        "source_label": drug_label,
        "source_type": node_info[drug_id]["node_type"],
        "relation": "targets",
        "target_id": protein_id,
        "target_label": protein_label,
        "target_type": node_info[protein_id]["node_type"],
        "evidence_count": example["drug_protein_evidence"],
    },
    {
        "source_id": protein_id,
        "source_label": protein_label,
        "source_type": node_info[protein_id]["node_type"],
        "relation": "participates_in",
        "target_id": pathway_id,
        "target_label": pathway_label,
        "target_type": node_info[pathway_id]["node_type"],
        "evidence_count": example["protein_pathway_evidence"],
    },
    {
        "source_id": drug_id,
        "source_label": drug_label,
        "source_type": node_info[drug_id]["node_type"],
        "relation": treatment_relation,
        "target_id": disease_id,
        "target_label": disease_label,
        "target_type": node_info[disease_id]["node_type"],
        "evidence_count": example["drug_disease_evidence"],
    },
    {
        "source_id": disease_id,
        "source_label": disease_label,
        "source_type": node_info[disease_id]["node_type"],
        "relation": "associated_with",
        "target_id": protein_id,
        "target_label": protein_label,
        "target_type": node_info[protein_id]["node_type"],
        "evidence_count": example["disease_protein_evidence"],
    },
]

example_df = pd.DataFrame(example_rows)

csv_path = (
    OUTPUT_DIR
    / "connected_drug_protein_pathway_disease_example.csv"
)

example_df.to_csv(
    csv_path,
    index=False,
)


# ================================================================
# Build graph
# ================================================================

G = nx.DiGraph()

G.add_edge(
    drug_id,
    protein_id,
    relation="targets",
)

G.add_edge(
    protein_id,
    pathway_id,
    relation="participates_in",
)

G.add_edge(
    drug_id,
    disease_id,
    relation=treatment_relation,
)

G.add_edge(
    disease_id,
    protein_id,
    relation="associated_with",
)


# ================================================================
# Positions designed to resemble the intended structure
# ================================================================

pos = {
    drug_id: (0, 0),
    disease_id: (0, -3),
    protein_id: (4, 0),
    pathway_id: (8, 0),
}


# ================================================================
# Draw visualization
# ================================================================

plt.figure(figsize=(14, 8))

nx.draw_networkx_nodes(
    G,
    pos,
    node_size=3000,
)

nx.draw_networkx_labels(
    G,
    pos,
    labels={
        drug_id: f"{drug_label}\nDrug",
        protein_id: f"{protein_label}\nProtein",
        pathway_id: f"{pathway_label}\nPathway",
        disease_id: f"{disease_label}\nDisease",
    },
    font_size=9,
)

nx.draw_networkx_edges(
    G,
    pos,
    arrows=True,
    arrowsize=22,
    width=2,
    connectionstyle="arc3,rad=0.08",
)

edge_labels = {
    (drug_id, protein_id): "targets",
    (protein_id, pathway_id): "participates_in",
    (drug_id, disease_id): treatment_relation,
    (disease_id, protein_id): "associated_with",
}

nx.draw_networkx_edge_labels(
    G,
    pos,
    edge_labels=edge_labels,
    font_size=9,
    label_pos=0.5,
)

plt.title(
    "ReMedy-KG v1.0 — Connected Drug–Protein–Pathway–Disease Example"
)

plt.axis("off")
plt.tight_layout()


png_path = (
    OUTPUT_DIR
    / "connected_drug_protein_pathway_disease_example.png"
)

plt.savefig(
    png_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()


# ================================================================
# Print result
# ================================================================

print("=" * 72)
print(
    "ReMedy-KG v1.0 — "
    "CONNECTED DRUG–PROTEIN–PATHWAY–DISEASE EXAMPLE"
)
print("=" * 72)
print()

print("Selected entities:")
print(f"  Drug:     {drug_id} | {drug_label}")
print(f"  Protein:  {protein_id} | {protein_label}")
print(f"  Pathway:  {pathway_id} | {pathway_label}")
print(f"  Disease:  {disease_id} | {disease_label}")
print()

print("Relationships:")
print(
    f"  {drug_id} "
    f"--targets--> "
    f"{protein_id}"
)
print(
    f"  {protein_id} "
    f"--participates_in--> "
    f"{pathway_id}"
)
print(
    f"  {drug_id} "
    f"--{treatment_relation}--> "
    f"{disease_id}"
)
print(
    f"  {disease_id} "
    f"--associated_with--> "
    f"{protein_id}"
)
print()

print("Evidence counts:")
print(
    f"  Drug -> Protein:       "
    f"{example['drug_protein_evidence']}"
)
print(
    f"  Protein -> Pathway:    "
    f"{example['protein_pathway_evidence']}"
)
print(
    f"  Drug -> Disease:       "
    f"{example['drug_disease_evidence']}"
)
print(
    f"  Disease -> Protein:    "
    f"{example['disease_protein_evidence']}"
)
print()

print(
    f"Total matching 4-edge examples found: "
    f"{len(connected):,}"
)
print()

print("Interpretation:")
print(
    "This visualization shows one connected subgraph in which "
    "a drug targets a protein, that protein participates in a "
    "pathway, the same drug has a treatment-use relationship "
    "with the disease, and the same protein is associated with "
    "that disease."
)
print()

print(
    "All four relationships already exist in ReMedy-KG v1.0. "
    "No new biological relationship was inferred or added."
)
print()

print(f"Saved: {png_path}")
print(f"Saved: {csv_path}")
print()

print(
    "Frozen Parquet files were read only; "
    "no dataset files were modified."
)
