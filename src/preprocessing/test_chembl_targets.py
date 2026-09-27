import pandas as pd
import requests
import time
from pathlib import Path

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

INPUT_FILE = "data/raw/chembl/drugbank_to_chembl_validated.csv"
OUTPUT_FILE = "data/interim/chembl_target_evidence_test_v2.csv"

# Activity types we consider useful for drug-target evidence
KEEP_TYPES = {
    "Ki",
    "Kd",
    "IC50",
    "EC50",
    "AC50",
    "pKi",
    "pKd",
    "Activity",
    "Potency",
    "Inhibition",
    "IC20",
    "Relative IC50"
}

# Types that are not direct drug-target potency/activity evidence
EXCLUDE_TYPES = {
    "Km",
    "Km app",
    "Vmax",
    "Kcat",
    "kon",
    "koff",
    "Ka",
    "KA",
    "Log K'",
    "k",
    "Relative Vmax",
    "Inactivation"
}

session = requests.Session()

target_cache = {}
assay_cache = {}
document_cache = {}


def get_json(url, params=None, retries=3):
    for attempt in range(retries):
        try:
            r = session.get(
                url,
                params=params,
                timeout=30
            )
            r.raise_for_status()
            return r.json()

        except Exception as e:
            if attempt == retries - 1:
                print(f"    Request failed: {e}")
                return None

            time.sleep(2 ** attempt)

    return None


def get_activities(chembl_id, limit=100):
    url = f"{BASE_URL}/activity.json"

    params = {
        "molecule_chembl_id": chembl_id,
        "limit": limit,
        "offset": 0
    }

    data = get_json(url, params)

    if not data:
        return []

    return data.get("activities", [])


def get_target(target_id):
    if not target_id:
        return None

    if target_id in target_cache:
        return target_cache[target_id]

    url = f"{BASE_URL}/target/{target_id}.json"

    data = get_json(url)

    target_cache[target_id] = data

    return data


def get_assay(assay_id):
    if not assay_id:
        return None

    if assay_id in assay_cache:
        return assay_cache[assay_id]

    url = f"{BASE_URL}/assay/{assay_id}.json"

    data = get_json(url)

    assay_cache[assay_id] = data

    return data


def get_document(document_id):
    if not document_id:
        return None

    if document_id in document_cache:
        return document_cache[document_id]

    url = f"{BASE_URL}/document/{document_id}.json"

    data = get_json(url)

    document_cache[document_id] = data

    return data


# ---------------------------------------------------------
# LOAD DRUG MAPPINGS
# ---------------------------------------------------------

df = pd.read_csv(INPUT_FILE)

df = df[df["mapping_status"] == "validated"].copy()

# Unique ChEMBL molecules
drugs = (
    df[["drugbank_id", "drug_name", "chembl_id"]]
    .drop_duplicates(subset=["chembl_id"])
    .head(10)
)

print("=" * 70)
print("REMEDY - OPTIMIZED CHEMBL TARGET EVIDENCE TEST")
print("=" * 70)

print(f"Testing drugs: {len(drugs)}")

results = []


# ---------------------------------------------------------
# PROCESS DRUGS
# ---------------------------------------------------------

for _, drug in drugs.iterrows():

    drugbank_id = drug["drugbank_id"]
    drug_name = drug["drug_name"]
    chembl_id = drug["chembl_id"]

    print("\n" + "-" * 70)
    print(f"{drug_name} ({chembl_id})")
    print("-" * 70)

    activities = get_activities(
        chembl_id,
        limit=100
    )

    print(f"Activities retrieved: {len(activities)}")

    relevant = 0

    for activity in activities:

        standard_type = activity.get("standard_type")
        standard_value = activity.get("standard_value")
        standard_units = activity.get("standard_units")

        # -------------------------------------------------
        # FILTER ACTIVITY TYPE
        # -------------------------------------------------

        if standard_type in EXCLUDE_TYPES:
            continue

        if standard_type not in KEEP_TYPES:
            continue

        # Require a measurable value
        if standard_value is None:
            continue

        target_id = activity.get("target_chembl_id")

        if not target_id:
            continue

        target = get_target(target_id)

        if not target:
            continue

        target_type = target.get("target_type")

        # Only single protein targets
        if target_type != "SINGLE PROTEIN":
            continue

        target_name = target.get("pref_name")

        # -------------------------------------------------
        # ASSAY INFORMATION
        # -------------------------------------------------

        assay_id = activity.get("assay_chembl_id")

        assay = get_assay(assay_id)

        assay_type = None
        assay_description = None
        assay_organism = None

        if assay:

            assay_type = assay.get("assay_type")

            assay_description = assay.get(
                "description"
            )

            assay_organism = assay.get(
                "assay_organism"
            )

        # -------------------------------------------------
        # DOCUMENT / PUBLICATION
        # -------------------------------------------------

        document_id = activity.get(
            "document_chembl_id"
        )

        document = get_document(document_id)

        publication = None
        doi = None
        pmid = None

        if document:

            publication = document.get(
                "title"
            )

            doi = document.get(
                "doi"
            )

            pmid = document.get(
                "pubmed_id"
            )

        # -------------------------------------------------
        # STORE
        # -------------------------------------------------

        results.append({

            "drugbank_id": drugbank_id,

            "drug_name": drug_name,

            "chembl_id": chembl_id,

            "activity_id":
                activity.get("activity_id"),

            "target_id":
                target_id,

            "target_name":
                target_name,

            "target_type":
                target_type,

            "standard_type":
                standard_type,

            "standard_value":
                standard_value,

            "standard_units":
                standard_units,

            "standard_relation":
                activity.get(
                    "standard_relation"
                ),

            "pchembl_value":
                activity.get(
                    "pchembl_value"
                ),

            "assay_id":
                assay_id,

            "assay_type":
                assay_type,

            "assay_description":
                assay_description,

            "assay_organism":
                assay_organism,

            "document_id":
                document_id,

            "publication":
                publication,

            "doi":
                doi,

            "pmid":
                pmid,

            "data_validity_comment":
                activity.get(
                    "data_validity_comment"
                ),

            "data_validity_description":
                activity.get(
                    "data_validity_description"
                )
        })

        relevant += 1

    print(
        f"Relevant protein activities: {relevant}"
    )


# ---------------------------------------------------------
# REMOVE DUPLICATES
# ---------------------------------------------------------

result_df = pd.DataFrame(results)

if len(result_df) > 0:

    result_df = result_df.drop_duplicates(
        subset=["activity_id"]
    )

    # -----------------------------------------------------
    # SUMMARY
    # -----------------------------------------------------

    print("\n" + "=" * 70)
    print("TEST COMPLETE")
    print("=" * 70)

    print(
        f"Evidence records : {len(result_df)}"
    )

    print(
        f"Unique drugs     : "
        f"{result_df['chembl_id'].nunique()}"
    )

    print(
        f"Unique targets   : "
        f"{result_df['target_id'].nunique()}"
    )

    print("\nActivity types:")

    print(
        result_df["standard_type"]
        .value_counts()
    )

    print("\nOrganisms:")

    print(
        result_df["assay_organism"]
        .value_counts(dropna=False)
        .head(20)
    )

else:

    print("\nNo qualifying evidence found.")


# ---------------------------------------------------------
# SAVE
# ---------------------------------------------------------

Path(OUTPUT_FILE).parent.mkdir(
    parents=True,
    exist_ok=True
)

result_df.to_csv(
    OUTPUT_FILE,
    index=False
)

print("\nSaved to:")
print(OUTPUT_FILE)

print("=" * 70)