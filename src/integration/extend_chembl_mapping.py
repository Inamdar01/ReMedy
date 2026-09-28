from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

MAIN_FILE = ROOT / "data/raw/chembl/drugbank_to_chembl_validated.csv"
ADDITIONAL_FILE = ROOT / "data/raw/chembl/chembl_missing_six_verified.csv"
OUTPUT_FILE = ROOT / "data/raw/chembl/drugbank_to_chembl_validated.csv"


def main():
    print("=" * 70)
    print("REMEDY - EXTEND VALIDATED DRUGBANK → CHEMBL MAPPING")
    print("=" * 70)

    # ---------------------------------------------------------
    # Load existing validated mapping
    # ---------------------------------------------------------
    main = pd.read_csv(
        MAIN_FILE,
        dtype=str
    ).fillna("")

    print(f"Existing validated rows: {len(main):,}")

    # ---------------------------------------------------------
    # Load the six supplemental mappings
    # ---------------------------------------------------------
    extra = pd.read_csv(
        ADDITIONAL_FILE,
        dtype=str
    ).fillna("")

    print(f"Additional mappings: {len(extra):,}")

    # Convert supplemental schema to the main schema
    extra = extra.rename(
        columns={
            "drugbank_name": "drug_name",
            "chembl_name": "chembl_pref_name",
            "mapping_method": "mapping_reason",
        }
    )

    extra["mapping_status"] = "validated"

    extra = extra[
        [
            "drugbank_id",
            "drug_name",
            "chembl_id",
            "chembl_pref_name",
            "mapping_status",
            "mapping_reason",
        ]
    ]

    # ---------------------------------------------------------
    # Check conflicts before merging
    # ---------------------------------------------------------
    existing_ids = set(
        main["drugbank_id"]
        .dropna()
        .astype(str)
    )

    overlapping = extra[
        extra["drugbank_id"].isin(existing_ids)
    ]

    if len(overlapping) > 0:
        print("\nWARNING: DrugBank IDs already present:")
        print(
            overlapping[
                ["drugbank_id", "chembl_id"]
            ].to_string(index=False)
        )
        raise ValueError(
            "Supplemental mapping contains existing DrugBank IDs."
        )

    # ---------------------------------------------------------
    # Merge
    # ---------------------------------------------------------
    merged = pd.concat(
        [main, extra],
        ignore_index=True
    )

    # Remove exact duplicates only
    merged = merged.drop_duplicates()

    # Sort by DrugBank ID
    merged = merged.sort_values(
        "drugbank_id"
    ).reset_index(drop=True)

    # ---------------------------------------------------------
    # Validate uniqueness
    # ---------------------------------------------------------
    duplicate_drugbank = (
        merged["drugbank_id"]
        .duplicated()
        .sum()
    )

    duplicate_chembl = (
        merged["chembl_id"]
        .duplicated()
        .sum()
    )

    print(
        f"\nDuplicate DrugBank IDs: {duplicate_drugbank}"
    )

    print(
        f"Duplicate ChEMBL IDs: {duplicate_chembl}"
    )

    if duplicate_drugbank != 0:
        raise ValueError(
            "Duplicate DrugBank IDs detected after merge."
        )

    # Multiple DrugBank IDs can legitimately map to the same
    # ChEMBL parent, so duplicate ChEMBL IDs are reported but
    # are not automatically treated as an error.

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------
    merged.to_csv(
        OUTPUT_FILE,
        index=False
    )

    # ---------------------------------------------------------
    # Final coverage check
    # ---------------------------------------------------------
    final_ids = set(
        merged["drugbank_id"]
        .dropna()
        .astype(str)
    )

    print()
    print("=" * 70)
    print("MAPPING EXTENSION COMPLETE")
    print("=" * 70)
    print(f"Final validated rows: {len(merged):,}")
    print(f"Unique DrugBank IDs: {len(final_ids):,}")

    expected_new = {
        "DB01036",
        "DB01049",
        "DB01230",
        "DB11273",
        "DB12131",
        "DB13522",
    }

    missing = sorted(expected_new - final_ids)

    print(f"Expected new IDs found: {len(expected_new - set(missing))}")
    print(f"Expected new IDs missing: {len(missing)}")

    if missing:
        print("\nMissing:")
        print("\n".join(missing))
        raise ValueError(
            "Not all six supplemental mappings were merged."
        )

    print()
    print(f"Saved: {OUTPUT_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    main()