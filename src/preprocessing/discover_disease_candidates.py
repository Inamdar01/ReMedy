import pandas as pd

PATH = "data/raw/primekg/primekg.csv"

SEARCH_TERMS = [
    "migraine",
    "dementia",
    "depression",
    "anxiety",
    "attention",
    "hyperactivity",
    "stroke",
    "cerebrovascular",
    "brain"
]

print("Loading PrimeKG...")
df = pd.read_csv(PATH, low_memory=False)

# Get unique disease records from both sides
x = df.loc[
    df["x_type"] == "disease",
    ["x_id", "x_name", "x_source"]
].copy()

x.columns = ["id", "name", "source"]

y = df.loc[
    df["y_type"] == "disease",
    ["y_id", "y_name", "y_source"]
].copy()

y.columns = ["id", "name", "source"]

diseases = pd.concat([x, y]).drop_duplicates()

print("\nTotal unique disease records:", len(diseases))

for term in SEARCH_TERMS:

    matches = diseases[
        diseases["name"]
        .str.contains(
            term,
            case=False,
            na=False,
            regex=False
        )
    ]

    print("\n" + "=" * 80)
    print("SEARCH:", term)
    print("=" * 80)

    if len(matches) == 0:
        print("No matches found.")
    else:
        print(
            matches[
                ["id", "name", "source"]
            ].to_string(index=False)
        )

print("\nSearch complete.")