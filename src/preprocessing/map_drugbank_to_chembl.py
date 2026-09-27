import pandas as pd
import requests
import time
from pathlib import Path

INPUT = "reports/primekg_scope_drugs.csv"
OUTPUT = "data/raw/chembl/drugbank_to_chembl.csv"

API_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json"

Path("data/raw/chembl").mkdir(parents=True, exist_ok=True)

print("Loading PrimeKG drug list...")

df = pd.read_csv(INPUT)

drugs = (
    df[["drug_id", "drug_name"]]
    .drop_duplicates()
    .sort_values("drug_id")
    .reset_index(drop=True)
)

print("Unique drugs to map:", len(drugs))

results = []

session = requests.Session()

for i, row in drugs.iterrows():

    drugbank_id = row["drug_id"]
    drug_name = str(row["drug_name"]).strip()

    print(
        f"[{i + 1}/{len(drugs)}] "
        f"{drugbank_id} -> {drug_name}"
    )

    try:

        response = session.get(
            API_URL,
            params={
                "q": drug_name,
                "limit": 20
            },
            timeout=30
        )

        response.raise_for_status()

        data = response.json()

        molecules = data.get("molecules") or []

        # Keep only records that actually have a ChEMBL ID
        valid_molecules = []

        for molecule in molecules:

            chembl_id = molecule.get("molecule_chembl_id")

            if chembl_id:
                valid_molecules.append(molecule)

        # Remove duplicate ChEMBL IDs
        unique = {}

        for molecule in valid_molecules:
            unique[
                molecule["molecule_chembl_id"]
            ] = molecule

        molecules = list(unique.values())

        # -------------------------------------------------
        # CASE 1: No ChEMBL result
        # -------------------------------------------------

        if len(molecules) == 0:

            results.append({
                "drugbank_id": drugbank_id,
                "drug_name": drug_name,
                "chembl_id": "",
                "chembl_pref_name": "",
                "status": "not_found",
                "result_count": 0
            })

        # -------------------------------------------------
        # CASE 2: Exactly one ChEMBL result
        # -------------------------------------------------

        elif len(molecules) == 1:

            molecule = molecules[0]

            pref_name = molecule.get("pref_name")

            if pref_name is None:
                pref_name = ""

            results.append({
                "drugbank_id": drugbank_id,
                "drug_name": drug_name,
                "chembl_id": molecule.get(
                    "molecule_chembl_id", ""
                ),
                "chembl_pref_name": pref_name,
                "status": "mapped",
                "result_count": 1
            })

        # -------------------------------------------------
        # CASE 3: Multiple results
        # -------------------------------------------------

        else:

            exact_matches = []

            for molecule in molecules:

                pref_name = molecule.get("pref_name")

                if pref_name is None:
                    pref_name = ""

                if (
                    pref_name.lower().strip()
                    == drug_name.lower().strip()
                ):
                    exact_matches.append(molecule)

            # One exact preferred-name match
            if len(exact_matches) == 1:

                molecule = exact_matches[0]

                pref_name = molecule.get("pref_name")

                if pref_name is None:
                    pref_name = ""

                results.append({
                    "drugbank_id": drugbank_id,
                    "drug_name": drug_name,
                    "chembl_id": molecule.get(
                        "molecule_chembl_id", ""
                    ),
                    "chembl_pref_name": pref_name,
                    "status": "mapped_exact_name",
                    "result_count": len(molecules)
                })

            # Multiple possible matches
            else:

                results.append({
                    "drugbank_id": drugbank_id,
                    "drug_name": drug_name,
                    "chembl_id": "",
                    "chembl_pref_name": "",
                    "status": "ambiguous",
                    "result_count": len(molecules)
                })

        # Small delay between requests
        time.sleep(0.2)

    except requests.exceptions.RequestException as e:

        print("  REQUEST ERROR:", e)

        results.append({
            "drugbank_id": drugbank_id,
            "drug_name": drug_name,
            "chembl_id": "",
            "chembl_pref_name": "",
            "status": "request_error",
            "result_count": 0
        })

    except Exception as e:

        print("  ERROR:", e)

        results.append({
            "drugbank_id": drugbank_id,
            "drug_name": drug_name,
            "chembl_id": "",
            "chembl_pref_name": "",
            "status": "error",
            "result_count": 0
        })


result_df = pd.DataFrame(results)

result_df.to_csv(
    OUTPUT,
    index=False
)

print("\n========== MAPPING SUMMARY ==========\n")

print(
    result_df["status"]
    .value_counts()
    .to_string()
)

print("\nTotal records:", len(result_df))

print(
    "Successfully mapped:",
    result_df["chembl_id"].ne("").sum()
)

print("\nSaved:")
print(OUTPUT)