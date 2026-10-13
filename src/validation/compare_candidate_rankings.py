
from pathlib import Path
from collections import defaultdict

import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = ROOT / "data" / "processed"
MODEL_DIR = ROOT / "data" / "model"
RESULTS_DIR = ROOT / "reports" / "model_evaluation"

ranking_path = (
    RESULTS_DIR / "rgcn_unobserved_candidate_ranking.csv"
)

ranking = pd.read_csv(ranking_path)
nodes = pd.read_parquet(DATA_DIR / "nodes.parquet")
graph = pd.read_parquet(
    MODEL_DIR / "training_graph_edges.parquet"
)

node_types = nodes.set_index("node_id")["node_type"].to_dict()
node_names = nodes.set_index("node_id")["label"].to_dict()

drug_targets = defaultdict(set)
disease_proteins = defaultdict(set)
protein_pathways = defaultdict(set)

# Build features only from the supporting biological graph.
for row in graph.itertuples(index=False):
    source = row.source_id
    target = row.target_id
    relation = row.relation

    if relation == "targets":
        if (
            node_types.get(source) == "drug"
            and node_types.get(target) == "protein"
        ):
            drug_targets[source].add(target)

    elif relation == "associated_with":
        if (
            node_types.get(source) == "disease"
            and node_types.get(target) == "protein"
        ):
            disease_proteins[source].add(target)

    elif relation == "participates_in":
        if (
            node_types.get(source) == "protein"
            and node_types.get(target) == "pathway"
        ):
            protein_pathways[source].add(target)

# Derive pathways reachable from each drug's targets
# and each disease's associated proteins.
drug_pathways = defaultdict(set)
disease_pathways = defaultdict(set)

for drug, proteins in drug_targets.items():
    for protein in proteins:
        drug_pathways[drug].update(protein_pathways[protein])

for disease, proteins in disease_proteins.items():
    for protein in proteins:
        disease_pathways[disease].update(protein_pathways[protein])


def jaccard(set_a, set_b):
    union = set_a | set_b
    if not union:
        return 0.0
    return len(set_a & set_b) / len(union)


# Compute graph features for each scored candidate.
rows = []

for row in ranking.itertuples(index=False):
    drug = row.drug_id
    disease = row.disease_id

    drug_proteins = drug_targets[drug]
    disease_protein_set = disease_proteins[disease]

    drug_pathway_set = drug_pathways[drug]
    disease_pathway_set = disease_pathways[disease]

    shared_proteins = drug_proteins & disease_protein_set
    shared_pathways = drug_pathway_set & disease_pathway_set

    protein_jaccard = jaccard(
        drug_proteins, disease_protein_set
    )
    pathway_jaccard = jaccard(
        drug_pathway_set, disease_pathway_set
    )

    rows.append({
        "drug_id": drug,
        "disease_id": disease,
        "rgcn_score": row.rgcn_score,
        "shared_protein_count": len(shared_proteins),
        "shared_pathway_count": len(shared_pathways),
        "shared_protein_jaccard": protein_jaccard,
        "shared_pathway_jaccard": pathway_jaccard,
        "biological_similarity_score": (
            protein_jaccard + pathway_jaccard
        ),
        "drug_target_count": len(drug_proteins),
        "disease_protein_count": len(disease_protein_set),
    })

comparison = pd.DataFrame(rows)

# Rank both approaches. Ties in biological similarity are
# kept tied in the statistical rank calculation.
comparison["rgcn_rank"] = comparison["rgcn_score"].rank(
    method="average", ascending=False
)

comparison["similarity_rank"] = (
    comparison["biological_similarity_score"].rank(
        method="average", ascending=False
    )
)

# Add human-readable labels.
comparison["drug_name"] = comparison["drug_id"].map(node_names)
comparison["disease_name"] = comparison["disease_id"].map(node_names)

comparison = comparison.sort_values(
    ["rgcn_score", "drug_id", "disease_id"],
    ascending=[False, True, True],
)

output_path = RESULTS_DIR / "candidate_rank_comparison.csv"
comparison.to_csv(output_path, index=False)

# Compare the two ranking approaches.
correlation = spearmanr(
    comparison["rgcn_score"],
    comparison["biological_similarity_score"],
)

top_rgcn = set(
    map(
        tuple,
        comparison.head(100)[["drug_id", "disease_id"]].to_numpy(),
    )
)

top_similarity_df = comparison.sort_values(
    [
        "biological_similarity_score",
        "shared_protein_count",
        "shared_pathway_count",
        "drug_id",
        "disease_id",
    ],
    ascending=[False, False, False, True, True],
).head(100)

top_similarity = set(
    map(
        tuple,
        top_similarity_df[["drug_id", "disease_id"]].to_numpy(),
    )
)

overlap = top_rgcn & top_similarity
union = top_rgcn | top_similarity

print("=" * 65)
print("R-GCN VS BIOLOGICAL SIMILARITY")
print("=" * 65)

print("\nSpearman rank correlation:", round(correlation.statistic, 4))
print("Correlation p-value:", f"{correlation.pvalue:.6g}")
print("Top-100 overlap:", len(overlap), "of 100")
print(
    "Top-100 Jaccard overlap:",
    round(len(overlap) / len(union), 4),
)

print("\nR-GCN score distribution:")
print(comparison["rgcn_score"].describe().to_string())

print("\nBiological similarity score distribution:")
print(
    comparison["biological_similarity_score"]
    .describe()
    .to_string()
)

print("\nTop 20 R-GCN candidates with supporting features:")
print(
    comparison[
        [
            "drug_name",
            "disease_name",
            "rgcn_score",
            "shared_protein_count",
            "shared_pathway_count",
            "biological_similarity_score",
        ]
    ].head(20).to_string(index=False)
)

print("\nMost frequent diseases in top 100 by R-GCN:")
print(
    comparison.head(100)["disease_name"]
    .value_counts()
    .to_string()
)

print("\nMost frequent diseases in top 100 by similarity:")
print(
    top_similarity_df["disease_name"]
    .value_counts()
    .to_string()
)

print("\nSaved:", output_path)
print("\nComparison completed.")
