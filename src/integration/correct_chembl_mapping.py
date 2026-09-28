from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

FILE = ROOT / "data/raw/chembl/drugbank_to_chembl_validated.csv"

REMOVE_IDS = {
    "DB11077",  # Polyethylene glycol 400
    "DB11161",  # Polyethylene glycol 300
}

VERIFIED = {
    "DB00783": {
        "chembl_id": "CHEMBL135",
        "chembl_pref_name": "ESTRADIOL",
        "mapping_reason": "verified_external_cross_reference",
    },
    "DB13396": {
        "chembl_id": "CHEMBL3707389",
        "chembl_pref_name": "NEOCITRULLAMON",
        "mapping_reason": "verified_external_cross_reference",
    },
    "DB08949": {
        "chembl_id": "CHEMBL1094982",
        "chembl_pref_name": "INOSITOL NICOTINATE",
        "mapping_reason": "verified_external_cross_reference",
    },
    "DB12749": {
        "chembl_id": "CHEMBL248594",
        "chembl_pref_name": "BUTYLPHTHALIDE",
        "mapping_reason": "verified_external_cross_reference",
    },
}


def main():

    print("=" * 70)
    print("REMEDY - CORRECT CHEMBL DRUG MAPPING")
    print("=" * 70)

    df = pd.read_csv(FILE, dtype=str).fillna("")

    print(f"Rows before correction: {len(df):,}")

    # ---------------------------------------------------------
    # Remove unsupported PEG mappings
    # ---------------------------------------------------------
    removed = df[df["drugbank_id"].isin(REMOVE_IDS)].copy()

    print("\nRemoving unsupported mappings:")

    if len(removed):
        print(
            removed[
                ["drugbank_id", "drug_name", "chembl_id"]
            ].to_string(index=False)
        )

    df = df[
        ~df["drugbank_id"].isin(REMOVE_IDS)
    ].copy()

    # ---------------------------------------------------------
    # Correct/confirm four externally verified mappings
    # ---------------------------------------------------------
    for drugbank_id, values in VERIFIED.items():

        mask = df["drugbank_id"] == drugbank_id

        if not mask.any():
            raise ValueError(
                f"{drugbank_id} is missing from the mapping file."
            )

        df.loc[mask, "chembl_id"] = values["chembl_id"]
        df.loc[mask, "chembl_pref_name"] = values["chembl_pref_name"]
        df.loc[mask, "mapping_status"] = "validated"
        df.loc[mask, "mapping_reason"] = values["mapping_reason"]

    # ---------------------------------------------------------
    # Validate
    # ---------------------------------------------------------
    df = df.drop_duplicates()

    duplicate_drugbank = df["drugbank_id"].duplicated().sum()

    if duplicate_drugbank:
        raise ValueError(
            f"Duplicate DrugBank IDs found: {duplicate_drugbank}"
        )

    print()
    print("=" * 70)
    print("CORRECTION COMPLETE")
    print("=" * 70)

    print(f"Rows after correction: {len(df):,}")
    print(
        f"Unique DrugBank IDs: "
        f"{df['drugbank_id'].nunique():,}"
    )

    print(
        f"Unique ChEMBL IDs: "
        f"{df['chembl_id'].nunique():,}"
    )

    print("\nVerified four:")
    print(
        df[
            df["drugbank_id"].isin(VERIFIED.keys())
        ].to_string(index=False)
    )

    print("\nConfirmed unresolved:")
    print("DB11077 | Polyethylene glycol 400")
    print("DB11161 | Polyethylene glycol 300")

    df.to_csv(FILE, index=False)

    print()
    print(f"Saved: {FILE}")
    print("=" * 70)


if __name__ == "__main__":
    main()