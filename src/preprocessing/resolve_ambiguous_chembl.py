import pandas as pd
import requests
import time
from pathlib import Path

INPUT = "data/raw/chembl/drugbank_to_chembl.csv"
OUTPUT = "data/raw/chembl/chembl_resolution_review.csv"

API_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json"

Path("data/raw/chembl").mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

problem_df = df[
    df["status"].isin(["ambiguous", "not_found"])
].copy()

print("Problematic drugs:", len(problem_df))
print()

session = requests.Session()

results = []

for i, row in problem_df.reset_index(drop=True).iterrows():

    drugbank_id = row["drugbank_id"]
    drug_name = str(row["drug_name"]).strip()

    print(
        f"[{i + 1}/{len(problem_df)}] "
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

        # Remove duplicate ChEMBL IDs
        unique = {}

        for molecule in molecules:

            chembl_id = molecule.get("molecule_chembl_id")

            if chembl_id:
                unique[chembl_id] = molecule

        molecules = list(unique.values())

        print("  Candidates:", len(molecules))

        for molecule in molecules:

            chembl_id = molecule.get(
                "molecule_chembl_id",
                ""
            )

            pref_name = molecule.get("pref_name")

            if pref_name is None:
                pref_name = ""

            molecule_type = molecule.get(
                "molecule_type",
                ""
            )

            max_phase = molecule.get(
                "max_phase"
            )

            first_approval = molecule.get(
                "first_approval"
            )

            results.append({
                "drugbank_id": drugbank_id,
                "drug_name": drug_name,
                "chembl_id": chembl_id,
                "chembl_pref_name": pref_name,
                "molecule_type": molecule_type,
                "max_phase": max_phase,
                "first_approval": first_approval,
                "original_status": row["status"]
            })

        if len(molecules) == 0:

            results.append({
                "drugbank_id": drugbank_id,
                "drug_name": drug_name,
                "chembl_id": "",
                "chembl_pref_name": "",
                "molecule_type": "",
                "max_phase": "",
                "first_approval": "",
                "original_status": row["status"]
            })

        time.sleep(0.3)

    except Exception as e:

        print("  ERROR:", e)

        results.append({
            "drugbank_id": drugbank_id,
            "drug_name": drug_name,
            "chembl_id": "",
            "chembl_pref_name": "",
            "molecule_type": "",
            "max_phase": "",
            "first_approval": "",
            "original_status": "request_error"
        })

review_df = pd.DataFrame(results)

review_df.to_csv(
    OUTPUT,
    index=False
)

print()
print("================================")
print("Resolution review complete")
print("================================")
print()

print("Rows:", len(review_df))

print()
print("Saved:")
print(OUTPUT)