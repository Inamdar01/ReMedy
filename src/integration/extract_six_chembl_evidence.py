from pathlib import Path
import time
import requests
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

DRUG_MAPPING_FILE = (
    ROOT / "data/raw/chembl/drugbank_to_chembl_validated.csv"
)

OUT_FILE = (
    ROOT / "data/raw/chembl/evidence/chembl_target_evidence_six.csv"
)

PROGRESS_FILE = (
    ROOT / "data/raw/chembl/evidence/chembl_six_progress.csv"
)

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

# These are the activity types already used by the ReMedy
# ChEMBL extraction.
ALLOWED_STANDARD_TYPES = {
    "AC50",
    "IC50",
    "Ki",
    "Kd",
    "Potency",
    "Inhibition",
    "Activity",
    "EC50",
    "Relative IC50",
    "pKi",
    "IC20",
}


SIX_DRUGS = {
    "CHEMBL135": "DB00783",
    "CHEMBL3707389": "DB13396",
    "CHEMBL1094982": "DB08949",
    "CHEMBL248594": "DB12749",
    "CHEMBL1382": "DB01036",
    "CHEMBL1177": "DB01230",
    # NOTE:
    # The other two newly added drugs are also mapped, but these
    # six are the ones we are extracting here from the unresolved
    # drug-disease evidence set.
}


def get_json(session, endpoint, params=None, retries=5):
    url = f"{BASE_URL}/{endpoint}"

    for attempt in range(1, retries + 1):
        try:
            response = session.get(
                url,
                params=params,
                timeout=60
            )

            if response.status_code == 200:
                return response.json()

            if response.status_code in {429, 500, 502, 503, 504}:
                wait = 2 ** attempt
                print(
                    f"  HTTP {response.status_code}; "
                    f"retrying in {wait}s..."
                )
                time.sleep(wait)
                continue

            print(
                f"  HTTP {response.status_code} for {url}"
            )
            return None

        except requests.RequestException as exc:
            wait = 2 ** attempt

            print(
                f"  Request error: {exc}; "
                f"retrying in {wait}s..."
            )

            time.sleep(wait)

    return None


def classify_organism(organism):
    if not organism or str(organism).strip() == "":
        return "organism_unspecified"

    value = str(organism).strip().lower()

    if "homo sapiens" in value:
        return "human"

    return "non_human"


def fetch_all_activities(session, chembl_id):
    activities = []
    offset = 0
    limit = 500

    while True:

        data = get_json(
            session,
            "activity.json",
            params={
                "molecule_chembl_id": chembl_id,
                "limit": limit,
                "offset": offset,
            }
        )

        if not data:
            break

        batch = data.get("activities", [])

        if not batch:
            break

        activities.extend(batch)

        page_meta = data.get("page_meta", {})

        total = page_meta.get("total")
        next_page = page_meta.get("next")

        print(
            f"    Retrieved {len(activities):,}"
            + (f" / {total:,}" if total is not None else "")
        )

        if next_page:
            offset += limit
        else:
            break

        time.sleep(0.1)

    return activities


def fetch_target(session, target_id, cache):
    if not target_id:
        return {}

    if target_id in cache:
        return cache[target_id]

    data = get_json(
        session,
        f"target/{target_id}.json"
    )

    if not data:
        cache[target_id] = {}
        return {}

    cache[target_id] = data
    return data


def fetch_assay(session, assay_id, cache):
    if not assay_id:
        return {}

    if assay_id in cache:
        return cache[assay_id]

    data = get_json(
        session,
        f"assay/{assay_id}.json"
    )

    if not data:
        cache[assay_id] = {}
        return {}

    cache[assay_id] = data
    return data


def fetch_document(session, document_id, cache):
    if not document_id:
        return {}

    if document_id in cache:
        return cache[document_id]

    data = get_json(
        session,
        f"document/{document_id}.json"
    )

    if not data:
        cache[document_id] = {}
        return {}

    cache[document_id] = data
    return data


def main():

    print("=" * 70)
    print("REMEDY - TARGETED CHEMBL EVIDENCE EXTRACTION")
    print("=" * 70)

    mapping = pd.read_csv(
        DRUG_MAPPING_FILE,
        dtype=str
    ).fillna("")

    six_mapping = mapping[
        mapping["chembl_id"].isin(SIX_DRUGS.keys())
    ].copy()

    # Verify every requested ChEMBL ID exists.
    found = set(six_mapping["chembl_id"])

    missing = set(SIX_DRUGS.keys()) - found

    if missing:
        raise ValueError(
            f"Missing ChEMBL IDs in validated mapping: "
            f"{sorted(missing)}"
        )

    print(
        f"Drugs selected: {six_mapping['chembl_id'].nunique()}"
    )

    for _, row in six_mapping.sort_values("chembl_id").iterrows():
        print(
            f"  {row['chembl_id']} | "
            f"{row['drugbank_id']} | "
            f"{row['drug_name']}"
        )

    session = requests.Session()
    session.headers.update({
        "User-Agent": "ReMedy-Research-Pipeline/1.0"
    })

    target_cache = {}
    assay_cache = {}
    document_cache = {}

    completed = set()

    if PROGRESS_FILE.exists():
        progress = pd.read_csv(
            PROGRESS_FILE,
            dtype=str
        ).fillna("")

        if "chembl_id" in progress.columns:
            completed = set(
                progress.loc[
                    progress["status"] == "completed",
                    "chembl_id"
                ]
            )

    all_rows = []

    # ---------------------------------------------------------
    # Extract each drug
    # ---------------------------------------------------------
    for _, drug in six_mapping.iterrows():

        chembl_id = drug["chembl_id"]
        drugbank_id = drug["drugbank_id"]
        drug_name = drug["drug_name"]

        if chembl_id in completed:
            print(
                f"\nSkipping completed: "
                f"{chembl_id}"
            )
            continue

        print()
        print("-" * 70)
        print(
            f"Drug: {drug_name} | "
            f"{chembl_id}"
        )
        print("-" * 70)

        activities = fetch_all_activities(
            session,
            chembl_id
        )

        print(
            f"    Raw activities: "
            f"{len(activities):,}"
        )

        kept = []

        for activity in activities:

            standard_type = (
                activity.get("standard_type")
                or ""
            ).strip()

            standard_value = (
                activity.get("standard_value")
            )

            target_id = (
                activity.get("target_chembl_id")
                or ""
            ).strip()

            # -------------------------------------------------
            # Same core filtering rules as existing extraction
            # -------------------------------------------------
            if standard_type not in ALLOWED_STANDARD_TYPES:
                continue

            if standard_value in [None, ""]:
                continue

            # The ReMedy evidence layer requires SINGLE PROTEIN
            # target relationships.
            target = fetch_target(
                session,
                target_id,
                target_cache
            )

            target_type = (
                target.get("target_type")
                or ""
            ).strip()

            if target_type != "SINGLE PROTEIN":
                continue

            target_name = (
                target.get("pref_name")
                or activity.get("target_pref_name")
                or ""
            )

            target_organism = (
                target.get("organism")
                or ""
            )

            assay_id = (
                activity.get("assay_chembl_id")
                or ""
            )

            assay = fetch_assay(
                session,
                assay_id,
                assay_cache
            )

            assay_type = (
                assay.get("assay_type")
                or ""
            )

            assay_description = (
                assay.get("description")
                or activity.get("assay_description")
                or ""
            )

            # Prefer the assay organism because the experimental
            # context must not be inferred from target organism.
            assay_organism = (
                activity.get("assay_organism")
                or assay.get("assay_organism")
                or ""
            )

            organism_class = classify_organism(
                assay_organism
            )

            cell_line = ""

            document_id = (
                activity.get("document_chembl_id")
                or ""
            )

            document = fetch_document(
                session,
                document_id,
                document_cache
            )

            publication = (
                document.get("title")
                or ""
            )

            doi = (
                document.get("doi")
                or ""
            )

            pmid = (
                document.get("pubmed_id")
                or ""
            )

            # Preserve ChEMBL validity information.
            validity_comment = (
                activity.get("data_validity_comment")
                or ""
            )

            validity_description = (
                activity.get("data_validity_description")
                or ""
            )

            value_quality = (
                "flagged"
                if validity_comment
                else "usable"
            )

            kept.append({
                "drugbank_id": drugbank_id,
                "drug_name": drug_name,
                "chembl_id": chembl_id,
                "activity_id": activity.get(
                    "activity_id",
                    ""
                ),
                "target_id": target_id,
                "target_name": target_name,
                "target_type": target_type,
                "standard_type": standard_type,
                "standard_value": standard_value,
                "standard_units": activity.get(
                    "standard_units",
                    ""
                ),
                "standard_relation": activity.get(
                    "standard_relation",
                    ""
                ),
                "pchembl_value": activity.get(
                    "pchembl_value",
                    ""
                ),
                "activity_type_raw": activity.get(
                    "type",
                    ""
                ),
                "assay_id": assay_id,
                "assay_type": assay_type,
                "assay_description": assay_description,
                "assay_organism": assay_organism,
                "organism_class": organism_class,
                "cell_line": cell_line,
                "document_id": document_id,
                "publication": publication,
                "doi": doi,
                "pmid": pmid,
                "data_validity_comment": validity_comment,
                "data_validity_description": validity_description,
                "value_quality": value_quality,
                "source": "ChEMBL",
                "target_organism": target_organism,
                "organism": target_organism,
                "organism_source": (
                    "target"
                    if target_organism
                    else ""
                ),
            })

        print(
            f"    Qualifying SINGLE PROTEIN evidence: "
            f"{len(kept):,}"
        )

        all_rows.extend(kept)

        # Save progress immediately.
        progress_rows = []

        for chembl in SIX_DRUGS.keys():
            if chembl in completed:
                progress_rows.append({
                    "chembl_id": chembl,
                    "status": "completed"
                })
            elif chembl == chembl_id:
                progress_rows.append({
                    "chembl_id": chembl,
                    "status": "completed"
                })
            else:
                progress_rows.append({
                    "chembl_id": chembl,
                    "status": "pending"
                })

        pd.DataFrame(
            progress_rows
        ).to_csv(
            PROGRESS_FILE,
            index=False
        )

        time.sleep(0.5)

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------
    result = pd.DataFrame(all_rows)

    if result.empty:
        print("\nNo qualifying evidence found.")
    else:
        result = result.drop_duplicates(
            subset=[
                "chembl_id",
                "activity_id"
            ]
        )

    OUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    result.to_csv(
        OUT_FILE,
        index=False
    )

    # ---------------------------------------------------------
    # Statistics
    # ---------------------------------------------------------
    print()
    print("=" * 70)
    print("TARGETED CHEMBL EXTRACTION COMPLETE")
    print("=" * 70)

    print(
        f"Evidence rows: {len(result):,}"
    )

    if not result.empty:

        print(
            f"Unique drugs: "
            f"{result['chembl_id'].nunique():,}"
        )

        print(
            f"Unique targets: "
            f"{result['target_id'].nunique():,}"
        )

        print(
            f"Human assay evidence: "
            f"{(result['organism_class'] == 'human').sum():,}"
        )

        print(
            f"Non-human assay evidence: "
            f"{(result['organism_class'] == 'non_human').sum():,}"
        )

        print(
            f"Unspecified assay organism: "
            f"{(result['organism_class'] == 'organism_unspecified').sum():,}"
        )

        print("\nActivity types:")
        print(
            result["standard_type"]
            .value_counts()
            .to_string()
        )

        print("\nEvidence by drug:")
        print(
            result.groupby(
                ["chembl_id", "drug_name"]
            )
            .size()
            .sort_values(ascending=False)
            .to_string()
        )

    print()
    print(f"Saved: {OUT_FILE}")
    print(f"Saved: {PROGRESS_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    main()