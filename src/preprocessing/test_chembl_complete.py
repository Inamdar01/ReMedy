import requests
import pandas as pd
import time

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

# Load validated mappings
mapping = pd.read_csv(
    "data/raw/chembl/drugbank_to_chembl_validated.csv"
)

# Only use validated mappings
mapping = mapping[mapping["mapping_status"] == "validated"].copy()

# Test ONE drug
row = mapping[mapping["chembl_id"] == "CHEMBL1201199"].iloc[0]

drugbank_id = row["drugbank_id"]
drug_name = row["drug_name"]
chembl_id = row["chembl_id"]

print("=" * 70)
print("TESTING COMPLETE CHEMBL EVIDENCE PIPELINE")
print("=" * 70)

print(f"DrugBank ID : {drugbank_id}")
print(f"Drug Name   : {drug_name}")
print(f"ChEMBL ID   : {chembl_id}")
print()


# ---------------------------------------------------------
# 1. GET ACTIVITIES
# ---------------------------------------------------------

print("1. Fetching activities...")

url = f"{BASE_URL}/activity.json"

params = {
    "molecule_chembl_id": chembl_id,
    "limit": 5
}

response = requests.get(url, params=params, timeout=30)
response.raise_for_status()

activities = response.json().get("activities", [])

print(f"Activities found: {len(activities)}")
print()


# ---------------------------------------------------------
# 2. INSPECT ACTIVITY
# ---------------------------------------------------------

for i, activity in enumerate(activities, start=1):

    print("-" * 70)
    print(f"ACTIVITY {i}")
    print("-" * 70)

    activity_id = activity.get("activity_id")
    target_id = activity.get("target_chembl_id")
    assay_id = activity.get("assay_chembl_id")
    document_id = activity.get("document_chembl_id")

    print("Activity ID       :", activity_id)
    print("Target ID         :", target_id)
    print("Assay ID          :", assay_id)
    print("Document ID       :", document_id)

    print("Standard Type     :", activity.get("standard_type"))
    print("Standard Value    :", activity.get("standard_value"))
    print("Standard Units    :", activity.get("standard_units"))

    print("Activity Comment  :", activity.get("activity_comment"))

    print("Assay Description :", activity.get("assay_description"))

    print("Assay Organism    :", activity.get("assay_organism"))
    print("Assay Cell Type   :", activity.get("assay_cell_type"))

    print()


    # -----------------------------------------------------
    # 3. GET TARGET INFORMATION
    # -----------------------------------------------------

    if target_id:

        print("Fetching target information...")

        target_url = f"{BASE_URL}/target/{target_id}.json"

        target_response = requests.get(
            target_url,
            timeout=30
        )

        if target_response.status_code == 200:

            target = target_response.json()

            print("Target Pref Name  :", target.get("pref_name"))
            print("Target Type       :", target.get("target_type"))

            components = target.get(
                "target_components",
                []
            )

            print("Target Components :", len(components))

            for component in components[:5]:

                print(
                    "  Component:",
                    component.get("component_id"),
                    component.get("component_type")
                )

        else:

            print(
                "Could not retrieve target:",
                target_response.status_code
            )

    print()


    # -----------------------------------------------------
    # 4. GET ASSAY INFORMATION
    # -----------------------------------------------------

    if assay_id:

        print("Fetching assay information...")

        assay_url = f"{BASE_URL}/assay/{assay_id}.json"

        assay_response = requests.get(
            assay_url,
            timeout=30
        )

        if assay_response.status_code == 200:

            assay = assay_response.json()

            print("Assay Type        :", assay.get("assay_type"))
            print("Assay Description :", assay.get("description"))
            print("Assay Organism    :", assay.get("assay_organism"))
            print("Cell Type         :", assay.get("assay_cell_type"))
            print("Tissue            :", assay.get("assay_tissue"))
            print("Variant           :", assay.get("assay_variant"))

        else:

            print(
                "Could not retrieve assay:",
                assay_response.status_code
            )

    print()


    # -----------------------------------------------------
    # 5. GET PUBLICATION INFORMATION
    # -----------------------------------------------------

    if document_id:

        print("Fetching publication information...")

        document_url = (
            f"{BASE_URL}/document/{document_id}.json"
        )

        document_response = requests.get(
            document_url,
            timeout=30
        )

        if document_response.status_code == 200:

            document = document_response.json()

            print("Title             :", document.get("title"))
            print("PubMed ID         :", document.get("pubmed_id"))
            print("DOI               :", document.get("doi"))
            print("Journal           :", document.get("journal"))
            print("Year              :", document.get("year"))

        else:

            print(
                "Could not retrieve publication:",
                document_response.status_code
            )

    print()

    time.sleep(0.2)

print("=" * 70)
print("TEST COMPLETE")
print("=" * 70)