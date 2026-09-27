import pandas as pd

PATH = "data/raw/primekg/primekg.csv"

TARGET_DISEASES = [
    "Alzheimer",
    "Parkinson",
    "Huntington",
    "amyotrophic lateral sclerosis",
    "multiple sclerosis",
    "epilepsy",
    "schizophrenia",
    "major depressive",
    "bipolar",
    "autism"
]

print("Loading PrimeKG...")
df = pd.read_csv(PATH, low_memory=False)

# Collect all unique disease names
x_diseases = df.loc[df["x_type"] == "disease", ["x_id", "x_name", "x_source"]]
y_diseases = df.loc[df["y_type"] == "disease", ["y_id", "y_name", "y_source"]]

x_diseases.columns = ["id", "name", "source"]
y_diseases.columns = ["id", "name", "source"]

diseases = pd.concat([x_diseases, y_diseases]).drop_duplicates()

print("\nTotal unique disease records:", len(diseases))

print("\n========== TARGET DISEASE SEARCH ==========")

for target in TARGET_DISEASES:

    matches = diseases[
        diseases["name"]
        .str.contains(target, case=False, na=False, regex=False)
    ]

    print(f"\n--- {target} ---")

    if len(matches) == 0:
        print("NO MATCH FOUND")
    else:
        print(matches.to_string(index=False))

print("\nSearch complete.")