import pandas as pd

PATH = "data/raw/primekg/primekg.csv"

print("Loading PrimeKG...")
df = pd.read_csv(PATH, low_memory=False)

# Candidate diseases.
# Exact disease names found in PrimeKG.
TARGET_DISEASES = [
    "epilepsy",
    "anxiety disorder",
    "Parkinson disease",
    "Alzheimer disease",
    "bipolar disorder",
    "schizophrenia",
    "major depressive disorder",
    "multiple sclerosis",
    "migraine disorder",
    "dementia (disease)",
    "attention deficit-hyperactivity disorder",
    "stroke disorder",
]

# Keep only rows where one side is disease
# and the other side is drug.
disease_drug = df[
    (
        ((df["x_type"] == "disease") & (df["y_type"] == "drug")) |
        ((df["x_type"] == "drug") & (df["y_type"] == "disease"))
    )
].copy()

print("\nTotal disease-drug relationships:", len(disease_drug))

results = []

for disease_name in TARGET_DISEASES:

    # Disease can appear on either side.
    mask = (
        (
            (disease_drug["x_type"] == "disease") &
            (disease_drug["x_name"].str.lower() == disease_name.lower())
        )
        |
        (
            (disease_drug["y_type"] == "disease") &
            (disease_drug["y_name"].str.lower() == disease_name.lower())
        )
    )

    matches = disease_drug[mask]

    drugs = set()

    for _, row in matches.iterrows():

        if row["x_type"] == "drug":
            drugs.add(row["x_id"])

        if row["y_type"] == "drug":
            drugs.add(row["y_id"])

    results.append({
        "disease": disease_name,
        "primekg_name": disease_name,
        "drug_relationships": len(matches),
        "unique_drugs": len(drugs)
    })

result_df = pd.DataFrame(results)

result_df = result_df.sort_values(
    "unique_drugs",
    ascending=False
)

print("\n========== DRUG COVERAGE ==========\n")
print(result_df.to_string(index=False))

result_df.to_csv(
    "reports/disease_drug_coverage.csv",
    index=False
)

print("\nSaved:")
print("reports/disease_drug_coverage.csv")