import pandas as pd
import requests
import time
from pathlib import Path
from datetime import datetime

# ============================================================
# CONFIG
# ============================================================

INPUT_FILE = "data/raw/chembl/drugbank_to_chembl_validated.csv"

OUTPUT_DIR = Path("data/raw/chembl/evidence")

OUTPUT_FILE = OUTPUT_DIR / "chembl_target_evidence_final.csv"
PROGRESS_FILE = OUTPUT_DIR / "chembl_extraction_progress.csv"
FAILED_FILE = OUTPUT_DIR / "chembl_failed_drugs_final.csv"
TARGET_CACHE_FILE = OUTPUT_DIR / "chembl_target_metadata_cache.csv"
TARGET_FAILURE_FILE = OUTPUT_DIR / "chembl_unresolved_targets.csv"

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

PAGE_SIZE = 500
MAX_RETRIES = 3
REQUEST_DELAY = 0.10

# ============================================================
# ACTIVITY TYPES TO KEEP
# ============================================================

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

# ============================================================
# SETUP
# ============================================================

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

session = requests.Session()

session.headers.update({
    "User-Agent": "ReMedy/1.0"
})

# ============================================================
# TARGET CACHE
# ============================================================

target_cache = {}

if TARGET_CACHE_FILE.exists():

    try:

        cache_df = pd.read_csv(
            TARGET_CACHE_FILE
        )

        for _, row in cache_df.iterrows():

            target_id = str(
                row["target_id"]
            )

            target_cache[target_id] = {
                "target_id": target_id,
                "target_name": row.get("target_name"),
                "target_type": row.get("target_type"),
                "target_organism": row.get("target_organism")
            }

        print(
            f"Loaded cached targets: "
            f"{len(target_cache)}"
        )

    except Exception as e:

        print(
            f"Could not load target cache: {e}"
        )


def save_target_cache():

    if not target_cache:
        return

    cache_df = pd.DataFrame(
        target_cache.values()
    )

    cache_df.to_csv(
        TARGET_CACHE_FILE,
        index=False
    )


# ============================================================
# API REQUEST
# ============================================================

def get_json(url, params=None):

    for attempt in range(MAX_RETRIES):

        try:

            response = session.get(
                url,
                params=params,
                timeout=30
            )

            response.raise_for_status()

            time.sleep(
                REQUEST_DELAY
            )

            return response.json()

        except Exception as e:

            if attempt < MAX_RETRIES - 1:

                wait = 2 ** attempt

                print(
                    f"      Retry "
                    f"{attempt + 1}/{MAX_RETRIES - 1} "
                    f"in {wait}s..."
                )

                time.sleep(wait)

            else:

                print(
                    f"      API failed: {e}"
                )

                return None

    return None


# ============================================================
# TARGET METADATA
# ============================================================

def get_target_info(target_id):

    if not target_id:
        return None

    target_id = str(
        target_id
    )

    # -----------------------------------------
    # CACHE HIT
    # -----------------------------------------

    if target_id in target_cache:

        return target_cache[target_id]

    # -----------------------------------------
    # API REQUEST
    # -----------------------------------------

    url = (
        f"{BASE_URL}/target/"
        f"{target_id}.json"
    )

    data = get_json(url)

    if not data:

        return None

    info = {
        "target_id":
            target_id,

        "target_name":
            data.get("pref_name"),

        "target_type":
            data.get("target_type"),

        "target_organism":
            data.get("organism")
    }

    target_cache[target_id] = info

    return info


# ============================================================
# GET ALL ACTIVITIES
# ============================================================

def get_all_activities(chembl_id):

    activities = []

    offset = 0

    while True:

        url = (
            f"{BASE_URL}/activity.json"
        )

        params = {
            "molecule_chembl_id":
                chembl_id,

            "limit":
                PAGE_SIZE,

            "offset":
                offset
        }

        data = get_json(
            url,
            params
        )

        # API failure
        if data is None:

            return None

        page = data.get(
            "activities",
            []
        )

        if not page:

            break

        activities.extend(
            page
        )

        print(
            f"      Activities retrieved: "
            f"{len(activities)}",
            end="\r"
        )

        if len(page) < PAGE_SIZE:

            break

        offset += PAGE_SIZE

    print()

    return activities


# ============================================================
# ORGANISM CLASSIFICATION
# ============================================================

def classify_organism(organism):

    if pd.isna(organism) or not organism:

        return "organism_unspecified"

    organism = str(
        organism
    ).strip()

    if organism.lower() == "homo sapiens":

        return "human"

    return "non_human"


# ============================================================
# VALUE QUALITY
# ============================================================

def classify_value(activity):

    value = activity.get(
        "standard_value"
    )

    if value is None:

        return "missing"

    comment = activity.get(
        "data_validity_comment"
    )

    description = activity.get(
        "data_validity_description"
    )

    text = (
        str(comment or "") +
        " " +
        str(description or "")
    ).lower()

    warning_words = [
        "outside typical range",
        "unusually large",
        "unusually small",
        "may not be accurate"
    ]

    for word in warning_words:

        if word in text:

            return "flagged"

    return "usable"


# ============================================================
# LOAD DRUGS
# ============================================================

print("=" * 75)
print("REMEDY - CHEMBL EXTRACTION WITH TARGET CACHE")
print("=" * 75)

mapping = pd.read_csv(
    INPUT_FILE
)

mapping = mapping[
    mapping["mapping_status"] == "validated"
].copy()

drugs = (
    mapping[
        [
            "drugbank_id",
            "drug_name",
            "chembl_id"
        ]
    ]
    .drop_duplicates(
        subset=["chembl_id"]
    )
    .reset_index(drop=True)
)

print(
    f"Unique validated drugs: "
    f"{len(drugs)}"
)


# ============================================================
# LOAD PROGRESS
# ============================================================

completed = set()

if PROGRESS_FILE.exists():

    progress_df = pd.read_csv(
        PROGRESS_FILE
    )

    if (
        "chembl_id" in progress_df.columns
        and "status" in progress_df.columns
    ):

        completed = set(
            progress_df[
                progress_df["status"] == "completed"
            ]["chembl_id"]
            .astype(str)
        )

print(
    f"Already completed: "
    f"{len(completed)}"
)

print(
    f"Remaining: "
    f"{len(drugs) - len(completed)}"
)


# ============================================================
# LOAD EXISTING EVIDENCE
# ============================================================

if OUTPUT_FILE.exists():

    final_df = pd.read_csv(
        OUTPUT_FILE
    )

else:

    final_df = pd.DataFrame()

print(
    f"Existing evidence rows: "
    f"{len(final_df)}"
)


# ============================================================
# FAILED TARGETS
# ============================================================

unresolved_targets = []


# ============================================================
# PROCESS DRUGS
# ============================================================

total = len(drugs)

for index, drug in drugs.iterrows():

    drugbank_id = str(
        drug["drugbank_id"]
    )

    drug_name = str(
        drug["drug_name"]
    )

    chembl_id = str(
        drug["chembl_id"]
    )

    # --------------------------------------------------------
    # RESUME
    # --------------------------------------------------------

    if chembl_id in completed:

        continue

    print()
    print("-" * 75)

    print(
        f"[{index + 1}/{total}] "
        f"{drug_name} ({chembl_id})"
    )

    start_time = time.time()

    try:

        # ----------------------------------------------------
        # GET ACTIVITIES
        # ----------------------------------------------------

        activities = get_all_activities(
            chembl_id
        )

        if activities is None:

            print(
                "      Activity API failed."
            )

            continue

        print(
            f"      Total activities retrieved: "
            f"{len(activities)}"
        )

        # ----------------------------------------------------
        # FIRST FILTER
        # ----------------------------------------------------

        filtered_activities = []

        for activity in activities:

            standard_type = activity.get(
                "standard_type"
            )

            standard_value = activity.get(
                "standard_value"
            )

            if standard_type not in KEEP_TYPES:
                continue

            if standard_value is None:
                continue

            target_id = activity.get(
                "target_chembl_id"
            )

            if not target_id:
                continue

            filtered_activities.append(
                activity
            )

        print(
            f"      Activities after type/value "
            f"filter: {len(filtered_activities)}"
        )

        # ----------------------------------------------------
        # UNIQUE TARGET IDS
        # ----------------------------------------------------

        target_ids = sorted(
            {
                str(
                    a.get("target_chembl_id")
                )
                for a in filtered_activities
                if a.get("target_chembl_id")
            }
        )

        print(
            f"      Unique targets to inspect: "
            f"{len(target_ids)}"
        )

        # ----------------------------------------------------
        # GET TARGET METADATA ONCE PER TARGET
        # ----------------------------------------------------

        target_metadata = {}

        new_targets = 0
        failed_targets = 0

        for target_id in target_ids:

            info = get_target_info(
                target_id
            )

            if info is None:

                failed_targets += 1

                unresolved_targets.append({

                    "chembl_id":
                        chembl_id,

                    "drug_name":
                        drug_name,

                    "target_id":
                        target_id,

                    "reason":
                        "target_api_failed",

                    "timestamp":
                        datetime.now().isoformat()
                })

                continue

            target_metadata[
                target_id
            ] = info

            new_targets += 1

            # Save cache periodically
            if new_targets % 25 == 0:

                save_target_cache()

        save_target_cache()

        # ----------------------------------------------------
        # BUILD FINAL EVIDENCE
        # ----------------------------------------------------

        records = []

        skipped_non_single = 0

        for activity in filtered_activities:

            target_id = str(
                activity.get(
                    "target_chembl_id"
                )
            )

            info = target_metadata.get(
                target_id
            )

            if info is None:

                continue

            target_type = info.get(
                "target_type"
            )

            if target_type != "SINGLE PROTEIN":

                skipped_non_single += 1

                continue

            target_name = info.get(
                "target_name"
            )

            target_organism = info.get(
                "target_organism"
            )

            # ------------------------------------------------
            # ORGANISM
            # ------------------------------------------------

            assay_organism = activity.get(
                "assay_organism"
            )

            if assay_organism:

                organism = (
                    assay_organism
                )

                organism_source = (
                    "assay_organism"
                )

            elif target_organism:

                organism = (
                    target_organism
                )

                organism_source = (
                    "target_organism"
                )

            else:

                organism = None

                organism_source = (
                    "unknown"
                )

            organism_class = classify_organism(
                organism
            )

            # ------------------------------------------------
            # RECORD
            # ------------------------------------------------

            records.append({

                "drugbank_id":
                    drugbank_id,

                "drug_name":
                    drug_name,

                "chembl_id":
                    chembl_id,

                "activity_id":
                    activity.get(
                        "activity_id"
                    ),

                "target_id":
                    target_id,

                "target_name":
                    target_name,

                "target_type":
                    target_type,

                "standard_type":
                    activity.get(
                        "standard_type"
                    ),

                "standard_value":
                    activity.get(
                        "standard_value"
                    ),

                "standard_units":
                    activity.get(
                        "standard_units"
                    ),

                "standard_relation":
                    activity.get(
                        "standard_relation"
                    ),

                "pchembl_value":
                    activity.get(
                        "pchembl_value"
                    ),

                "activity_type_raw":
                    activity.get(
                        "type"
                    ),

                "assay_id":
                    activity.get(
                        "assay_chembl_id"
                    ),

                "assay_type":
                    activity.get(
                        "assay_type"
                    ),

                "assay_description":
                    activity.get(
                        "assay_description"
                    ),

                "assay_organism":
                    assay_organism,

                "target_organism":
                    target_organism,

                "organism":
                    organism,

                "organism_source":
                    organism_source,

                "organism_class":
                    organism_class,

                "cell_line":
                    activity.get(
                        "assay_cell_type"
                    ),

                "document_id":
                    activity.get(
                        "document_chembl_id"
                    ),

                "publication":
                    None,

                "doi":
                    None,

                "pmid":
                    None,

                "data_validity_comment":
                    activity.get(
                        "data_validity_comment"
                    ),

                "data_validity_description":
                    activity.get(
                        "data_validity_description"
                    ),

                "value_quality":
                    classify_value(
                        activity
                    ),

                "source":
                    "ChEMBL"
            })

        # ----------------------------------------------------
        # SAVE EVIDENCE
        # ----------------------------------------------------

        if records:

            batch_df = pd.DataFrame(
                records
            )

            batch_df = batch_df.drop_duplicates(
                subset=["activity_id"]
            )

            if final_df.empty:

                final_df = batch_df

            else:

                final_df = pd.concat(
                    [
                        final_df,
                        batch_df
                    ],
                    ignore_index=True
                )

            final_df = final_df.drop_duplicates(
                subset=["activity_id"]
            )

            final_df.to_csv(
                OUTPUT_FILE,
                index=False
            )

        # ----------------------------------------------------
        # STATS
        # ----------------------------------------------------

        protein_count = len(records)

        human_count = sum(
            r["organism_class"] == "human"
            for r in records
        )

        non_human_count = sum(
            r["organism_class"] == "non_human"
            for r in records
        )

        unknown_count = sum(
            r["organism_class"]
            == "organism_unspecified"
            for r in records
        )

        print(
            f"      Protein evidence: "
            f"{protein_count}"
        )

        print(
            f"      Human: "
            f"{human_count}"
        )

        print(
            f"      Non-human: "
            f"{non_human_count}"
        )

        print(
            f"      Organism unspecified: "
            f"{unknown_count}"
        )

        print(
            f"      Target lookup failures: "
            f"{failed_targets}"
        )

        print(
            f"      Non-single-protein: "
            f"{skipped_non_single}"
        )

        # ----------------------------------------------------
        # SAVE PROGRESS
        # ----------------------------------------------------

        duration = round(
            time.time() -
            start_time,
            2
        )

        progress_row = pd.DataFrame([{

            "chembl_id":
                chembl_id,

            "status":
                "completed",

            "activity_count":
                len(activities),

            "drug_name":
                drug_name,

            "activities_retrieved":
                len(activities),

            "protein_evidence":
                protein_count,

            "target_lookup_failures":
                failed_targets,

            "duration_seconds":
                duration,

            "timestamp":
                datetime.now().isoformat()
        }])

        if PROGRESS_FILE.exists():

            old_progress = pd.read_csv(
                PROGRESS_FILE
            )

            old_progress = old_progress[
                old_progress["chembl_id"].astype(str)
                != chembl_id
            ]

            progress_out = pd.concat(
                [
                    old_progress,
                    progress_row
                ],
                ignore_index=True
            )

        else:

            progress_out = progress_row

        progress_out.to_csv(
            PROGRESS_FILE,
            index=False
        )

        completed.add(
            chembl_id
        )

        print(
            f"      Time: {duration}s"
        )

    except KeyboardInterrupt:

        print()
        print(
            "Extraction interrupted."
        )
        print(
            "Previously completed drugs remain saved."
        )

        save_target_cache()

        raise

    except Exception as e:

        print(
            f"      Unexpected error: {e}"
        )

        continue


# ============================================================
# SAVE UNRESOLVED TARGETS
# ============================================================

if unresolved_targets:

    unresolved_df = pd.DataFrame(
        unresolved_targets
    )

    if TARGET_FAILURE_FILE.exists():

        old = pd.read_csv(
            TARGET_FAILURE_FILE
        )

        unresolved_df = pd.concat(
            [
                old,
                unresolved_df
            ],
            ignore_index=True
        )

        unresolved_df = unresolved_df.drop_duplicates(
            subset=[
                "chembl_id",
                "target_id"
            ]
        )

    unresolved_df.to_csv(
        TARGET_FAILURE_FILE,
        index=False
    )


# ============================================================
# FINAL STATUS
# ============================================================

print()
print("=" * 75)
print("CHEMBL EXTRACTION STATUS")
print("=" * 75)

if PROGRESS_FILE.exists():

    progress = pd.read_csv(
        PROGRESS_FILE
    )

    print(
        f"Completed drugs: "
        f"{(progress.status == 'completed').sum()}"
    )

if OUTPUT_FILE.exists():

    final_df = pd.read_csv(
        OUTPUT_FILE
    )

    print(
        f"Evidence rows: "
        f"{len(final_df)}"
    )

    print(
        f"Unique drugs: "
        f"{final_df['chembl_id'].nunique()}"
    )

    print(
        f"Unique targets: "
        f"{final_df['target_id'].nunique()}"
    )

    print()
    print("Organism:")

    print(
        final_df[
            "organism_class"
        ].value_counts(
            dropna=False
        )
    )

    print()
    print("Activity types:")

    print(
        final_df[
            "standard_type"
        ].value_counts()
    )

    print()
    print("Value quality:")

    print(
        final_df[
            "value_quality"
        ].value_counts(
            dropna=False
        )
    )

print()
print(
    f"Evidence file: {OUTPUT_FILE}"
)

print(
    f"Progress file: {PROGRESS_FILE}"
)

print(
    f"Target cache: {TARGET_CACHE_FILE}"
)

print("=" * 75)