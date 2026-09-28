import pandas as pd
import requests
import time
from pathlib import Path
from datetime import datetime

# ============================================================
# CONFIG
# ============================================================

INPUT_FILE = (
    "data/interim/chembl_evidence_clean.csv"
)

OUTPUT_FILE = (
    "data/interim/"
    "chembl_target_uniprot_mapping.csv"
)

PROGRESS_FILE = (
    "data/interim/"
    "chembl_uniprot_mapping_progress.csv"
)

REPORT_FILE = (
    "reports/"
    "chembl_uniprot_mapping_report.md"
)

BASE_URL = (
    "https://www.ebi.ac.uk/"
    "chembl/api/data"
)

BATCH_SIZE = 50
MAX_RETRIES = 3
REQUEST_TIMEOUT = 30

Path("data/interim").mkdir(
    parents=True,
    exist_ok=True
)

Path("reports").mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# SESSION
# ============================================================

session = requests.Session()

session.headers.update({
    "User-Agent": "ReMedy/1.0"
})


# ============================================================
# API REQUEST
# ============================================================

def get_targets_batch(target_ids):

    params = {
        "target_chembl_id__in":
            ",".join(target_ids),

        "limit":
            1000,

        "offset":
            0,

        "only":
            "target_chembl_id,"
            "pref_name,"
            "target_type,"
            "organism,"
            "target_components"
    }

    for attempt in range(MAX_RETRIES):

        try:

            response = session.get(
                f"{BASE_URL}/target.json",
                params=params,
                timeout=REQUEST_TIMEOUT
            )

            response.raise_for_status()

            return response.json()

        except Exception as e:

            print(
                f"      API attempt "
                f"{attempt + 1}/{MAX_RETRIES} failed: "
                f"{e}"
            )

            if attempt < MAX_RETRIES - 1:

                time.sleep(
                    2 ** attempt
                )

    return None


# ============================================================
# LOAD TARGET IDS
# ============================================================

print("=" * 70)
print(
    "REMEDY - CHEMBL TARGET → UNIPROT MAPPING"
)
print("=" * 70)

df = pd.read_csv(
    INPUT_FILE,
    low_memory=False
)

target_df = (
    df[
        [
            "target_id",
            "target_name",
            "target_type"
        ]
    ]
    .dropna(
        subset=["target_id"]
    )
    .drop_duplicates(
        subset=["target_id"]
    )
)

target_ids = (
    target_df["target_id"]
    .astype(str)
    .tolist()
)

print(
    f"Unique ChEMBL targets: "
    f"{len(target_ids)}"
)


# ============================================================
# LOAD EXISTING PROGRESS
# ============================================================

completed_targets = set()

if Path(PROGRESS_FILE).exists():

    progress = pd.read_csv(
        PROGRESS_FILE
    )

    completed_targets = set(
        progress[
            progress["status"] == "completed"
        ]["target_id"]
        .astype(str)
    )

print(
    f"Already mapped: "
    f"{len(completed_targets)}"
)

print(
    f"Remaining: "
    f"{len(target_ids) - len(completed_targets)}"
)


# ============================================================
# LOAD EXISTING OUTPUT
# ============================================================

if Path(OUTPUT_FILE).exists():

    mapping_out = pd.read_csv(
        OUTPUT_FILE,
        low_memory=False
    )

else:

    mapping_out = pd.DataFrame()


# ============================================================
# PROCESS IN BATCHES
# ============================================================

remaining_ids = [
    x for x in target_ids
    if x not in completed_targets
]

total_batches = (
    (len(remaining_ids) + BATCH_SIZE - 1)
    // BATCH_SIZE
)

failed_targets = []

for batch_number, start in enumerate(
    range(
        0,
        len(remaining_ids),
        BATCH_SIZE
    ),
    start=1
):

    batch_ids = remaining_ids[
        start:start + BATCH_SIZE
    ]

    print()
    print(
        f"Batch {batch_number}/"
        f"{total_batches}"
    )

    print(
        f"Targets: "
        f"{len(batch_ids)}"
    )

    data = get_targets_batch(
        batch_ids
    )

    # --------------------------------------------------------
    # API FAILURE
    # --------------------------------------------------------

    if data is None:

        print(
            "      Batch failed."
        )

        for target_id in batch_ids:

            failed_targets.append({

                "target_id":
                    target_id,

                "status":
                    "api_failed",

                "timestamp":
                    datetime.now().isoformat()
            })

        continue

    targets = data.get(
        "targets",
        []
    )

    print(
        f"      API returned: "
        f"{len(targets)} targets"
    )

    returned_ids = set()

    batch_records = []
    progress_records = []

    # --------------------------------------------------------
    # PROCESS TARGETS
    # --------------------------------------------------------

    for target in targets:

        target_id = str(
            target.get(
                "target_chembl_id"
            )
        )

        returned_ids.add(
            target_id
        )

        target_name = target.get(
            "pref_name"
        )

        target_type = target.get(
            "target_type"
        )

        target_organism = target.get(
            "organism"
        )

        components = (
            target.get(
                "target_components"
            )
            or []
        )

        # ----------------------------------------------------
        # EXTRACT UNIPROT ACCESSIONS
        # ----------------------------------------------------

        accessions = []

        component_relations = []

        for component in components:

            accession = component.get(
                "accession"
            )

            relation = component.get(
                "relationship"
            )

            if accession:

                accessions.append(
                    str(accession)
                )

                component_relations.append(
                    relation
                )

        accessions = list(
            dict.fromkeys(
                accessions
            )
        )

        # ----------------------------------------------------
        # TARGET WITH NO ACCESSION
        # ----------------------------------------------------

        if not accessions:

            batch_records.append({

                "target_id":
                    target_id,

                "target_name":
                    target_name,

                "target_type":
                    target_type,

                "target_organism":
                    target_organism,

                "uniprot_accession":
                    None,

                "component_count":
                    len(components),

                "mapping_status":
                    "no_uniprot_accession",

                "human_target":
                    (
                        str(
                            target_organism
                            or ""
                        ).strip().lower()
                        == "homo sapiens"
                    ),

                "source":
                    "ChEMBL"
            })

        else:

            for accession in accessions:

                human_target = (
                    str(
                        target_organism
                        or ""
                    ).strip().lower()
                    == "homo sapiens"
                )

                batch_records.append({

                    "target_id":
                        target_id,

                    "target_name":
                        target_name,

                    "target_type":
                        target_type,

                    "target_organism":
                        target_organism,

                    "uniprot_accession":
                        accession,

                    "component_count":
                        len(components),

                    "mapping_status":
                        "mapped",

                    "human_target":
                        human_target,

                    "source":
                        "ChEMBL"
                })

        progress_records.append({

            "target_id":
                target_id,

            "status":
                "completed",

            "timestamp":
                datetime.now().isoformat()
        })

    # --------------------------------------------------------
    # TARGETS NOT RETURNED BY API
    # --------------------------------------------------------

    missing_from_response = (
        set(batch_ids)
        - returned_ids
    )

    for target_id in missing_from_response:

        failed_targets.append({

            "target_id":
                target_id,

            "status":
                "not_returned",

            "timestamp":
                datetime.now().isoformat()
        })

        progress_records.append({

            "target_id":
                target_id,

            "status":
                "api_not_returned",

            "timestamp":
                datetime.now().isoformat()
        })

    # --------------------------------------------------------
    # SAVE BATCH
    # --------------------------------------------------------

    if batch_records:

        batch_df = pd.DataFrame(
            batch_records
        )

        if mapping_out.empty:

            mapping_out = batch_df

        else:

            mapping_out = pd.concat(
                [
                    mapping_out,
                    batch_df
                ],
                ignore_index=True
            )

        mapping_out = (
            mapping_out
            .drop_duplicates(
                subset=[
                    "target_id",
                    "uniprot_accession"
                ]
            )
        )

        mapping_out.to_csv(
            OUTPUT_FILE,
            index=False
        )

    # --------------------------------------------------------
    # SAVE PROGRESS
    # --------------------------------------------------------

    progress_batch = pd.DataFrame(
        progress_records
    )

    if Path(PROGRESS_FILE).exists():

        old_progress = pd.read_csv(
            PROGRESS_FILE
        )

        old_progress = old_progress[
            ~old_progress["target_id"]
            .astype(str)
            .isin(
                progress_batch[
                    "target_id"
                ].astype(str)
            )
        ]

        progress_out = pd.concat(
            [
                old_progress,
                progress_batch
            ],
            ignore_index=True
        )

    else:

        progress_out = progress_batch

    progress_out.to_csv(
        PROGRESS_FILE,
        index=False
    )

    print(
        f"      Mapped records so far: "
        f"{len(mapping_out)}"
    )


# ============================================================
# SAVE FAILED TARGETS
# ============================================================

if failed_targets:

    failed_df = pd.DataFrame(
        failed_targets
    )

    failed_file = Path(
        "data/interim/"
        "chembl_unresolved_uniprot_targets.csv"
    )

    if failed_file.exists():

        old_failed = pd.read_csv(
            failed_file
        )

        failed_df = pd.concat(
            [
                old_failed,
                failed_df
            ],
            ignore_index=True
        )

    failed_df = (
        failed_df
        .drop_duplicates(
            subset=["target_id"]
        )
    )

    failed_df.to_csv(
        failed_file,
        index=False
    )


# ============================================================
# FINAL REPORT
# ============================================================

if Path(OUTPUT_FILE).exists():

    final = pd.read_csv(
        OUTPUT_FILE,
        low_memory=False
    )

    unique_targets = (
        final["target_id"]
        .nunique()
    )

    unique_uniprot = (
        final[
            "uniprot_accession"
        ]
        .dropna()
        .nunique()
    )

    human_targets = (
        final[
            final["human_target"]
            == True
        ]
        ["target_id"]
        .nunique()
    )

    unmapped_targets = (
        final[
            final["uniprot_accession"]
            .isna()
        ]
        ["target_id"]
        .nunique()
    )

else:

    final = pd.DataFrame()

    unique_targets = 0
    unique_uniprot = 0
    human_targets = 0
    unmapped_targets = 0


report = f"""# ChEMBL Target → UniProt Mapping Report

## Input

`{INPUT_FILE}`

## Results

- ChEMBL targets in input: {len(target_ids)}
- ChEMBL targets represented in output: {unique_targets}
- Unique UniProt accessions: {unique_uniprot}
- Human ChEMBL targets: {human_targets}
- Targets without UniProt accession: {unmapped_targets}

## Mapping method

ChEMBL target records were requested in batches using the
`target_chembl_id__in` filter.

UniProt accessions were extracted from ChEMBL
`target_components`.

No protein identifier was inferred from target names.

## Reactome eligibility

Human ChEMBL targets with a UniProt accession can be joined
against the Reactome UniProt mapping using the accession.

Non-human targets are retained in this mapping but should not
be used for the human Reactome pathway layer.

## Files

- `chembl_target_uniprot_mapping.csv`
- `chembl_uniprot_mapping_progress.csv`
- `chembl_unresolved_uniprot_targets.csv` (if needed)
"""

Path(REPORT_FILE).write_text(
    report,
    encoding="utf-8"
)

print()
print("=" * 70)
print("MAPPING COMPLETE")
print("=" * 70)

print(
    f"ChEMBL targets: "
    f"{len(target_ids)}"
)

print(
    f"Unique UniProt: "
    f"{unique_uniprot}"
)

print(
    f"Human targets: "
    f"{human_targets}"
)

print(
    f"Unmapped targets: "
    f"{unmapped_targets}"
)

print()
print(
    f"Saved: {OUTPUT_FILE}"
)

print(
    f"Progress: {PROGRESS_FILE}"
)

print(
    f"Report: {REPORT_FILE}"
)

print("=" * 70)