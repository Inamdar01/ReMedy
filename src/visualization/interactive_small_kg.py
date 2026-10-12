from pathlib import Path
import json

import pandas as pd
from pyvis.network import Network


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
# Identify the same connected 4-edge pattern used in Visualization 7
#
# Drug -> Protein -> Pathway
# Drug -> Disease -> Protein
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


# Build connected candidates
connected = drug_protein.merge(
    protein_pathway,
    on="protein_id",
    how="inner",
)

connected = connected.merge(
    disease_protein,
    on="protein_id",
    how="inner",
)

connected = connected.merge(
    drug_disease,
    on=["drug_id", "disease_id"],
    how="inner",
)


if connected.empty:
    raise RuntimeError(
        "No connected 4-edge Drug-Protein-Pathway-Disease "
        "example exists in ReMedy-KG v1.0."
    )


# Prefer indication over off-label use,
# then prefer stronger combined evidence.
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
    ["relation_priority", "combined_evidence"],
    ascending=[True, False],
)

example = connected.iloc[0]

drug_id = example["drug_id"]
protein_id = example["protein_id"]
pathway_id = example["pathway_id"]
disease_id = example["disease_id"]
treatment_relation = example["relation"]


# ================================================================
# Helper for labels
# ================================================================

def get_label(node_id):
    return node_info.get(node_id, {}).get("label", node_id)


drug_label = get_label(drug_id)
protein_label = get_label(protein_id)
pathway_label = get_label(pathway_id)
disease_label = get_label(disease_id)


# ================================================================
# Create interactive network
# ================================================================

net = Network(
    height="750px",
    width="100%",
    directed=True,
    heading="ReMedy-KG v1.0 — Interactive Evidence Neighborhood",
    cdn_resources="in_line",
)

net.set_options(
    json.dumps(
        {
            "nodes": {
                "font": {
                    "size": 18,
                },
                "borderWidth": 2,
            },
            "edges": {
                "arrows": {
                    "to": {
                        "enabled": True,
                        "scaleFactor": 0.8,
                    }
                },
                "font": {
                    "size": 14,
                    "align": "middle",
                },
                "smooth": {
                    "type": "curvedCW",
                    "roundness": 0.15,
                },
            },
            "physics": {
                "enabled": True,
                "stabilization": {
                    "iterations": 200,
                },
            },
            "interaction": {
                "hover": True,
                "navigationButtons": True,
                "zoomView": True,
                "dragView": True,
            },
        }
    )
)


# ================================================================
# Add nodes
# ================================================================

selected_nodes = [
    (drug_id, drug_label, "Drug"),
    (protein_id, protein_label, "Protein"),
    (pathway_id, pathway_label, "Pathway"),
    (disease_id, disease_label, "Disease"),
]


# Fixed positions make the initial layout easier to understand.
positions = {
    drug_id: (-350, 0),
    protein_id: (50, 0),
    pathway_id: (450, 0),
    disease_id: (-350, 250),
}


for node_id, label, node_type in selected_nodes:

    title = (
        f"<b>{node_type}</b><br>"
        f"ID: {node_id}<br>"
        f"Label: {label}"
    )

    # Distinct node shapes for entity types.
    shape = {
        "Drug": "ellipse",
        "Protein": "box",
        "Pathway": "diamond",
        "Disease": "hexagon",
    }[node_type]

    net.add_node(
        node_id,
        label=f"{label}\n({node_id})",
        title=title,
        shape=shape,
        x=positions[node_id][0],
        y=positions[node_id][1],
        physics=False,
    )


# ================================================================
# Add four source-backed edges
# ================================================================

net.add_edge(
    drug_id,
    protein_id,
    label="targets",
    title=(
        f"targets<br>"
        f"Evidence records: "
        f"{int(example['drug_protein_evidence'])}"
    ),
)

net.add_edge(
    protein_id,
    pathway_id,
    label="participates_in",
    title=(
        f"participates_in<br>"
        f"Evidence records: "
        f"{int(example['protein_pathway_evidence'])}"
    ),
)

net.add_edge(
    drug_id,
    disease_id,
    label=treatment_relation,
    title=(
        f"{treatment_relation}<br>"
        f"Evidence records: "
        f"{int(example['drug_disease_evidence'])}"
    ),
)

net.add_edge(
    disease_id,
    protein_id,
    label="associated_with",
    title=(
        f"associated_with<br>"
        f"Evidence records: "
        f"{int(example['disease_protein_evidence'])}"
    ),
)


# ================================================================
# Add explanatory text below the graph
# ================================================================

description = (
    "<div style='padding:15px; font-family:Arial;'>"
    "<h2>What this graph shows</h2>"
    "<p>"
    "This is a connected subgraph from ReMedy-KG v1.0. "
    "The drug targets a protein, the protein participates in "
    "a pathway, the same drug has a treatment-use relationship "
    "with the disease, and the same protein is associated with "
    "the disease."
    "</p>"
    "<p>"
    "<b>Important:</b> All four relationships already exist in "
    "ReMedy-KG v1.0. No new biological relationship was inferred "
    "or added by this visualization."
    "</p>"
    "</div>"
)


# ================================================================
# Save interactive HTML
# ================================================================

html_path = OUTPUT_DIR / "interactive_small_kg.html"

# Generate the HTML in memory first.
# We then write it explicitly as UTF-8 so Windows cp1252
# encoding cannot corrupt Unicode characters.
html = net.generate_html(
    notebook=False,
)

# Insert explanation before closing body tag.
html = html.replace(
    "</body>",
    description + "</body>",
)

html_path.write_text(
    html,
    encoding="utf-8",
)


# ================================================================
# Print result
# ================================================================

print("=" * 72)
print("ReMedy-KG v1.0 — INTERACTIVE SMALL KG")
print("=" * 72)
print()

print("Selected entities:")
print(f"  Drug:     {drug_id} | {drug_label}")
print(f"  Protein:  {protein_id} | {protein_label}")
print(f"  Pathway:  {pathway_id} | {pathway_label}")
print(f"  Disease:  {disease_id} | {disease_label}")
print()

print("Relationships:")
print(f"  {drug_id} --targets--> {protein_id}")
print(f"  {protein_id} --participates_in--> {pathway_id}")
print(f"  {drug_id} --{treatment_relation}--> {disease_id}")
print(f"  {disease_id} --associated_with--> {protein_id}")
print()

print(
    f"Total matching 4-edge examples found: "
    f"{len(connected):,}"
)
print()

print(f"Saved interactive graph: {html_path}")
print()

print(
    "Frozen Parquet files were read only; "
    "no dataset files were modified."
)
