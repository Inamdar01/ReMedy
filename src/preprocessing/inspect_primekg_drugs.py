import pandas as pd

PATH = "data/raw/primekg/primekg.csv"

TARGET_DISEASES = [
    "schizophrenia",
    "anxiety disorder",
    "bipolar disorder",
    "epilepsy",
    "major depressive disorder",
    "stroke disorder",
    "Parkinson disease",
    "Alzheimer disease",
    "migraine disorder",
    "multiple sclerosis",
]

print("Loading PrimeKG...")
df = pd.read_csv(PATH, low_memory=False)

print("\nFinding drug-disease relationships...")

mask = (
    (
        (df["x_type"] == "disease") &
        (df["y_type"] == "drug")
    )
    |
    (
        (df["x_type"] == "drug") &
        (df["y_type"] == "disease")
    )
)

dd = df[mask].copy()

all_drugs = []

for disease in TARGET_DISEASES:

    disease_mask = (
        (
            (dd["x_type"] == "disease") &
            (dd["x_name"].str.lower() == disease.lower())
        )
        |
        (
            (dd["y_type"] == "disease") &
            (dd["y_name"].str.lower() == disease.lower())
        )
    )

    matches = dd[disease_mask]

    for _, row in matches.iterrows():

        if row["x_type"] == "drug":
            all_drugs.append({
                "disease": disease,
                "drug_id": row["x_id"],
                "drug_name": row["x_name"],
                "drug_source": row["x_source"]
            })

        elif row["y_type"] == "drug":
            all_drugs.append({
                "disease": disease,
                "drug_id": row["y_id"],
                "drug_name": row["y_name"],
                "drug_source": row["y_source"]
            })

drugs = pd.DataFrame(all_drugs)

drugs = drugs.drop_duplicates()

print("\nTotal disease-drug records:", len(drugs))
print("Unique drugs:", drugs["drug_id"].nunique())

print("\n========== DRUG SOURCE COUNTS ==========\n")
print(
    drugs[["drug_id", "drug_name", "drug_source"]]
    .drop_duplicates()
    ["drug_source"]
    .value_counts()
)

print("\n========== SAMPLE DRUG IDs ==========\n")
print(
    drugs[
        ["drug_id", "drug_name", "drug_source"]
    ]
    .drop_duplicates()
    .head(50)
    .to_string(index=False)
)

drugs.to_csv(
    "reports/primekg_scope_drugs.csv",
    index=False
)

print("\nSaved:")
print("reports/primekg_scope_drugs.csv")