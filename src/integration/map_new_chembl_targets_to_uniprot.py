from pathlib import Path
import time
import requests
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

RELATIONS_FILE = (
    ROOT / "data/interim/chembl_drug_protein_relations.csv"
)

EXISTING_MAPPING_FILE = (
    ROOT / "data/interim/chembl_target_uniprot_mapping.csv"
)

OUTPUT_FILE = (
    ROOT / "data/interim/chembl_target_uniprot_mapping_extended.csv"
)

REPORT_FILE = (
    ROOT / "reports/chembl_new_target_uniprot_report.md"
)

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

NEW_TARGETS = {
    "CHEMBL2995",
    "CHEMBL3065",
    "CHEMBL3751",
    "CHEMBL4295650",
    "CHEMBL5872",
}


def get_targets(session, target_ids):
    """
    Retrieve target metadata in one ChEMBL API request.
    """
    params = {
        "target_chembl_id__in": ",".join(sorted(target_ids)),
        "limit": 100,
    }

    for attempt in range(1, 6):
        try:
            response = session.get(
                f"{BASE_URL}/target.json",
                params=params,
                timeout=60,
            )

            if response.status_code == 200:
                return response.json()

            if response.status_code in {429, 500, 502, 503, 504}:
                wait = min(30, 2 ** attempt)
                print(
                    f"HTTP {response.status_code}; "
                    f"retrying in {wait}s..."
                )
                time.sleep(wait)
                continue

            raise RuntimeError(
                f"ChEMBL API returned HTTP "
                f"{response.status_code}"
            )

        except requests.RequestException as exc:
            wait = min(30, 2 ** attempt)

            print(
                f"Request error: {exc}; "
                f"retrying in {wait}s..."
            )

            time.sleep(wait)

    raise RuntimeError(
        "Unable to retrieve ChEMBL target metadata."
    )


def extract_uniprot(target):
    """
    Extract UniProt accession from ChEMBL target components.
    Preserve all accessions when multiple components exist.
    """

    accessions = []

    components = target.get(
        "target_components",
        []
    ) or []

    for component in components:

        accession = (
            component.get("accession")
            or ""
        ).strip()

        if accession:
            accessions.append(accession)

        # ChEMBL sometimes exposes component-level
        # UniProt mappings inside database identifiers.
        db_ids = component.get(
            "target_component_xrefs",
            []
        ) or []

        for xref in db_ids:

            db_name = (
                xref.get("xref_src_db")
                or ""
            ).strip().lower()

            xref_id = (
                xref.get("xref_id")
                or ""
            ).strip()

            if (
                "uniprot" in db_name
                and xref_id
            ):
                accessions.append(xref_id)

    return sorted(set(accessions))


def main():

    print("=" * 70)
    print("REMEDY - EXTEND CHEMBL TARGET → UNIPROT MAPPING")
    print("=" * 70)

    # ---------------------------------------------------------
    # 1. Verify these targets really are new
    # ---------------------------------------------------------
    existing = pd.read_csv(
        EXISTING_MAPPING_FILE,
        dtype=str
    ).fillna("")

    existing_ids = set(
        existing["target_id"]
        .astype(str)
    )

    still_new = NEW_TARGETS - existing_ids

    print(
        f"Existing target mappings: "
        f"{len(existing_ids):,}"
    )

    print(
        f"New targets requested: "
        f"{len(NEW_TARGETS)}"
    )

    print(
        f"Still unmapped: "
        f"{len(still_new)}"
    )

    if not still_new:
        print("\nNothing to map.")
        return

    # ---------------------------------------------------------
    # 2. Get names from current relation table
    # ---------------------------------------------------------
    relations = pd.read_csv(
        RELATIONS_FILE,
        dtype=str
    ).fillna("")

    target_info = (
        relations[
            relations["target_id"].isin(still_new)
        ][
            [
                "target_id",
                "target_name"
            ]
        ]
        .drop_duplicates()
        .sort_values("target_id")
    )

    print("\nTargets:")
    print(
        target_info.to_string(
            index=False
        )
    )

    # ---------------------------------------------------------
    # 3. Query ChEMBL
    # ---------------------------------------------------------
    session = requests.Session()

    session.headers.update({
        "User-Agent": "ReMedy-Research-Pipeline/1.0"
    })

    data = get_targets(
        session,
        still_new
    )

    targets = data.get(
        "targets",
        []
    )

    print(
        f"\nChEMBL target records returned: "
        f"{len(targets)}"
    )

    # ---------------------------------------------------------
    # 4. Build mapping rows
    # ---------------------------------------------------------
    rows = []

    returned_ids = set()

    for target in targets:

        target_id = (
            target.get(
                "target_chembl_id",
                ""
            )
            or ""
        ).strip()

        returned_ids.add(target_id)

        target_name = (
            target.get(
                "pref_name",
                ""
            )
            or ""
        ).strip()

        target_type = (
            target.get(
                "target_type",
                ""
            )
            or ""
        ).strip()

        organism = (
            target.get(
                "organism",
                ""
            )
            or ""
        ).strip()

        accessions = extract_uniprot(
            target
        )

        if accessions:

            for accession in accessions:
                rows.append({
                    "target_id": target_id,
                    "target_name": target_name,
                    "target_type": target_type,
                    "target_organism": organism,
                    "uniprot_accession": accession,
                    "component_count": len(
                        target.get(
                            "target_components",
                            []
                        ) or []
                    ),
                    "mapping_status": "mapped",
                    "human_target": (
                        organism.lower()
                        == "homo sapiens"
                    ),
                    "source": "ChEMBL",
                })

        else:

            rows.append({
                "target_id": target_id,
                "target_name": target_name,
                "target_type": target_type,
                "target_organism": organism,
                "uniprot_accession": "",
                "component_count": len(
                    target.get(
                        "target_components",
                        []
                    ) or []
                ),
                "mapping_status": "unmapped",
                "human_target": (
                    organism.lower()
                    == "homo sapiens"
                ),
                "source": "ChEMBL",
            })

    # ---------------------------------------------------------
    # 5. Check missing API records
    # ---------------------------------------------------------
    api_missing = (
        still_new - returned_ids
    )

    if api_missing:
        print(
            "\nTargets not returned by API:"
        )
        print(
            "\n".join(
                sorted(api_missing)
            )
        )

        for target_id in sorted(api_missing):

            name = target_info.loc[
                target_info["target_id"]
                == target_id,
                "target_name"
            ]

            target_name = (
                name.iloc[0]
                if len(name)
                else ""
            )

            rows.append({
                "target_id": target_id,
                "target_name": target_name,
                "target_type": "",
                "target_organism": "",
                "uniprot_accession": "",
                "component_count": 0,
                "mapping_status": "api_not_returned",
                "human_target": False,
                "source": "ChEMBL",
            })

    new_mapping = pd.DataFrame(
        rows
    ).drop_duplicates()

    # ---------------------------------------------------------
    # 6. Combine with existing mapping
    # ---------------------------------------------------------
    extended = pd.concat(
        [
            existing,
            new_mapping
        ],
        ignore_index=True
    ).drop_duplicates()

    # ---------------------------------------------------------
    # 7. Validate
    # ---------------------------------------------------------
    duplicate_pairs = (
        extended[
            ["target_id", "uniprot_accession"]
        ]
        .duplicated()
        .sum()
    )

    if duplicate_pairs:
        print(
            f"\nWarning: duplicate target-UniProt pairs: "
            f"{duplicate_pairs}"
        )

    mapped_new = new_mapping[
        (
            new_mapping["mapping_status"]
            == "mapped"
        )
        &
        (
            new_mapping["uniprot_accession"]
            != ""
        )
    ]["target_id"].nunique()

    unresolved_new = len(still_new) - mapped_new

    # ---------------------------------------------------------
    # 8. Save
    # ---------------------------------------------------------
    extended = extended.sort_values(
        ["target_id", "uniprot_accession"]
    ).reset_index(drop=True)

    extended.to_csv(
        OUTPUT_FILE,
        index=False
    )

    # ---------------------------------------------------------
    # 9. Report
    # ---------------------------------------------------------
    report = f"""# ChEMBL New Target → UniProt Mapping Report

## New targets

- Requested: {len(NEW_TARGETS)}
- Still unmapped before API query: {len(still_new)}
- Successfully mapped: {mapped_new}
- Unresolved: {unresolved_new}

## Extended mapping

- Existing mapping rows: {len(existing):,}
- New mapping rows: {len(new_mapping):,}
- Extended mapping rows: {len(extended):,}
- Unique target IDs in extended mapping:
  {extended["target_id"].nunique():,}
- Unique UniProt accessions:
  {extended["uniprot_accession"].replace("", pd.NA).dropna().nunique():,}

## New target results

{new_mapping.to_markdown(index=False)}

## Output

`data/interim/chembl_target_uniprot_mapping_extended.csv`
"""

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8"
    )

    # ---------------------------------------------------------
    # 10. Terminal output
    # ---------------------------------------------------------
    print()
    print("=" * 70)
    print("TARGET MAPPING EXTENSION COMPLETE")
    print("=" * 70)

    print(
        f"New targets mapped: "
        f"{mapped_new}/{len(still_new)}"
    )

    print(
        f"New targets unresolved: "
        f"{unresolved_new}"
    )

    print(
        f"Extended unique targets: "
        f"{extended['target_id'].nunique():,}"
    )

    print(
        f"Extended unique UniProt: "
        f"{extended['uniprot_accession'].replace('', pd.NA).dropna().nunique():,}"
    )

    print("\nNew target mapping:")
    print(
        new_mapping.to_string(
            index=False
        )
    )

    print()
    print(f"Saved: {OUTPUT_FILE}")
    print(f"Saved: {REPORT_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    main()