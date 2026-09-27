import pandas as pd
from pathlib import Path

INPUT = "data/raw/chembl/drugbank_to_chembl.csv"
OUTPUT = "data/raw/chembl/drugbank_to_chembl_validated.csv"
UNRESOLVED = "reports/unresolved_identifiers.csv"

Path("data/raw/chembl").mkdir(parents=True, exist_ok=True)
Path("reports").mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

# Clearly defensible mappings identified from the ChEMBL candidate review.
validated = {
    "DB00245": ("CHEMBL1201203", "BENZTROPINE", "synonym_match"),
    "DB00340": ("CHEMBL1201342", "METHIXENE", "synonym_match"),
    "DB00353": ("CHEMBL1201356", "METHYLERGONOVINE", "synonym_match"),
    "DB00527": ("CHEMBL1086", "DIBUCAINE", "synonym_match"),
    "DB00647": ("CHEMBL1213351", "PROPOXYPHENE", "parent_compound"),
    "DB00717": ("CHEMBL1162", "NORETHINDRONE", "synonym_match"),
    "DB00783": ("", "", "manual_review"),
    "DB00816": ("CHEMBL776", "METAPROTERENOL", "synonym_match"),
    "DB00849": ("CHEMBL45029", "MEPHOBARBITAL", "synonym_match"),
    "DB00882": ("CHEMBL2355051", "CLOMIPHENE", "synonym_match"),
    "DB00945": ("CHEMBL25", "ASPIRIN", "synonym_match"),
    "DB00977": ("CHEMBL691", "ETHINYL ESTRADIOL", "synonym_match"),
    "DB01001": ("CHEMBL714", "ALBUTEROL", "synonym_match"),
    "DB01045": ("CHEMBL374478", "RIFAMPIN", "synonym_match"),
    "DB01246": ("CHEMBL829", "METHYLPROMAZINE", "synonym_match"),
    "DB01253": ("CHEMBL119443", "ERGONOVINE", "synonym_match"),
    "DB01397": ("CHEMBL3989810", "MAGNESIUM SALICYLATE", "exact_preferred_name"),
    "DB01403": ("CHEMBL1764", "LEVOMEPROMAZINE", "synonym_match"),
    "DB01577": ("CHEMBL1201201", "METHAMPHETAMINE", "synonym_match"),
    "DB06691": ("CHEMBL511", "PYRILAMINE", "synonym_match"),
    "DB08949": ("", "", "manual_review"),
    "DB09167": ("CHEMBL1492500", "DOTHIEPIN", "synonym_match"),
    "DB11077": ("", "", "manual_review"),
    "DB11161": ("", "", "manual_review"),
    "DB11236": ("CHEMBL1201479", "POLYETHYLENE GLYCOL 3350", "manual_review"),
    "DB12478": ("CHEMBL1371770", "PIRIBEDIL", "preferred_name_match"),
    "DB12749": ("", "", "manual_review"),
    "DB13396": ("", "", "not_found"),
    "DB14478": ("CHEMBL2108139", "POVIDONE", "preferred_name_match"),
}

rows = []

for _, row in df.iterrows():

    drugbank_id = row["drugbank_id"]

    if drugbank_id in validated:

        chembl_id, chembl_name, reason = validated[drugbank_id]

        if chembl_id:

            rows.append({
                "drugbank_id": drugbank_id,
                "drug_name": row["drug_name"],
                "chembl_id": chembl_id,
                "chembl_pref_name": chembl_name,
                "mapping_status": "validated",
                "mapping_reason": reason
            })

        else:

            rows.append({
                "drugbank_id": drugbank_id,
                "drug_name": row["drug_name"],
                "chembl_id": "",
                "chembl_pref_name": "",
                "mapping_status": "unresolved",
                "mapping_reason": reason
            })

    elif row["status"] in ["mapped", "mapped_exact_name"]:

        rows.append({
            "drugbank_id": drugbank_id,
            "drug_name": row["drug_name"],
            "chembl_id": row["chembl_id"],
            "chembl_pref_name": row["chembl_pref_name"],
            "mapping_status": "validated",
            "mapping_reason": "previous_mapping"
        })

    else:

        rows.append({
            "drugbank_id": drugbank_id,
            "drug_name": row["drug_name"],
            "chembl_id": "",
            "chembl_pref_name": "",
            "mapping_status": "unresolved",
            "mapping_reason": "ambiguous_or_not_found"
        })


validated_df = pd.DataFrame(rows)

validated_df.to_csv(
    OUTPUT,
    index=False
)

# Create unresolved identifier report
unresolved_df = validated_df[
    validated_df["mapping_status"] == "unresolved"
].copy()

unresolved_df.to_csv(
    UNRESOLVED,
    index=False
)

print()
print("======================================")
print("VALIDATED MAPPING COMPLETE")
print("======================================")
print()

print(
    validated_df["mapping_status"]
    .value_counts()
    .to_string()
)

print()
print("Total drugs:", len(validated_df))
print(
    "Validated:",
    (validated_df["mapping_status"] == "validated").sum()
)
print(
    "Unresolved:",
    (validated_df["mapping_status"] == "unresolved").sum()
)

print()
print("Validated mapping:")
print(OUTPUT)

print()
print("Unresolved identifiers:")
print(UNRESOLVED)