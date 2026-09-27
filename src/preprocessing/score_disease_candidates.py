import pandas as pd

PATH = "data/raw/primekg/primekg.csv"

print("Loading PrimeKG...")
df = pd.read_csv(PATH, low_memory=False)

TARGET_DISEASES = [
    "epilepsy",
    "anxiety disorder",
    "stroke disorder",
    "Parkinson disease",
    "Alzheimer disease",
    "migraine disorder",
    "bipolar disorder",
    "schizophrenia",
    "dementia (disease)",
    "major depressive disorder",
    "attention deficit-hyperactivity disorder",
    "multiple sclerosis",
]

# ---------------------------------------------------------
# Find relationships involving the target diseases
# ---------------------------------------------------------

results = []

for disease_name in TARGET_DISEASES:

    disease_mask = (
        (
            (df["x_type"] == "disease") &
            (df["x_name"].str.lower() == disease_name.lower())
        )
        |
        (
            (df["y_type"] == "disease") &
            (df["y_name"].str.lower() == disease_name.lower())
        )
    )

    disease_df = df[disease_mask].copy()

    # Unique drugs
    drugs = set()

    # PrimeKG combines genes and proteins under gene/protein
    proteins_genes = set()

    # Other biological entities
    biological_entities = set()

    for _, row in disease_df.iterrows():

        if row["x_type"] == "drug":
            drugs.add(row["x_id"])

        if row["y_type"] == "drug":
            drugs.add(row["y_id"])

        if row["x_type"] == "gene/protein":
            proteins_genes.add(row["x_id"])

        if row["y_type"] == "gene/protein":
            proteins_genes.add(row["y_id"])

        if row["x_type"] in [
            "gene/protein",
            "biological_process",
            "molecular_function",
            "cellular_component",
            "pathway",
            "anatomy",
            "effect/phenotype"
        ]:
            biological_entities.add(row["x_id"])

        if row["y_type"] in [
            "gene/protein",
            "biological_process",
            "molecular_function",
            "cellular_component",
            "pathway",
            "anatomy",
            "effect/phenotype"
        ]:
            biological_entities.add(row["y_id"])

    results.append({
        "disease": disease_name,
        "total_relationships": len(disease_df),
        "unique_drugs": len(drugs),
        "unique_proteins_genes": len(proteins_genes),
        "total_biological_entities": len(biological_entities)
    })

result_df = pd.DataFrame(results)

result_df["connectivity_score"] = (
    result_df["unique_drugs"]
    + result_df["unique_proteins_genes"]
)

result_df = result_df.sort_values(
    "connectivity_score",
    ascending=False
)

print("\n========== DISEASE CONNECTIVITY ==========\n")
print(result_df.to_string(index=False))

result_df.to_csv(
    "reports/disease_candidate_scores.csv",
    index=False
)

print("\nSaved:")
print("reports/disease_candidate_scores.csv")