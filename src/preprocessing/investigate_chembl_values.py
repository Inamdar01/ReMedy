import requests
import json

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

CHEMBL_ID = "CHEMBL1201199"

print("=" * 70)
print("INVESTIGATING CHEMBL ACTIVITY VALUES")
print("=" * 70)
print("Drug:", CHEMBL_ID)
print()


# ---------------------------------------------------------
# 1. Get the specific activity
# ---------------------------------------------------------

activity_id = 893936

url = f"{BASE_URL}/activity/{activity_id}.json"

response = requests.get(url, timeout=30)
response.raise_for_status()

activity = response.json()

print("ACTIVITY RECORD")
print("-" * 70)

fields = [
    "activity_id",
    "molecule_chembl_id",
    "assay_chembl_id",
    "target_chembl_id",
    "standard_type",
    "standard_value",
    "standard_units",
    "standard_relation",
    "standard_text",
    "standard_upper_value",
    "standard_lower_value",
    "standard_flag",
    "value_flag",
    "pchembl_value",
    "data_validity_comment",
    "data_validity_description",
    "activity_comment",
    "record_id",
    "document_chembl_id"
]

for field in fields:
    print(f"{field:28}: {activity.get(field)}")


# ---------------------------------------------------------
# 2. Show the raw JSON
# ---------------------------------------------------------

print()
print("=" * 70)
print("RAW ACTIVITY JSON")
print("=" * 70)

print(json.dumps(activity, indent=2))


# ---------------------------------------------------------
# 3. Get assay information
# ---------------------------------------------------------

assay_id = activity.get("assay_chembl_id")

print()
print("=" * 70)
print("ASSAY INFORMATION")
print("=" * 70)

if assay_id:

    url = f"{BASE_URL}/assay/{assay_id}.json"

    response = requests.get(url, timeout=30)
    response.raise_for_status()

    assay = response.json()

    print("Assay ID:", assay_id)
    print("Assay Type:", assay.get("assay_type"))
    print("Description:", assay.get("description"))
    print("Organism:", assay.get("assay_organism"))
    print("Cell Type:", assay.get("assay_cell_type"))

    print()
    print("RAW ASSAY JSON")
    print("-" * 70)
    print(json.dumps(assay, indent=2))


print()
print("=" * 70)
print("INVESTIGATION COMPLETE")
print("=" * 70)