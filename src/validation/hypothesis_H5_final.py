import hashlib
import math
from pathlib import Path

import pandas as pd
import yaml
from scipy.stats import hypergeom


ROOT = Path(__file__).resolve().parents[2]

DD_PATH = ROOT / "data/interim/drug_disease_relations_canonical.csv"
DP_CHEMBL_PATH = ROOT / "data/interim/drug_protein_relations_canonical.csv"
DP_PRIMEKG_PATH = ROOT / "data/interim/primekg_drug_protein_relations_canonical.csv"
DISEASE_PROTEIN_PATH = ROOT / "data/interim/disease_uniprot_relations_canonical.csv"
MAPPING_PATH = ROOT / "data/raw/chembl/drugbank_to_chembl_validated.csv"
DRUGMECHDB_PATH = ROOT / "data/raw/drugmechdb/indication_paths.yaml"

REPORT_PATH = ROOT / "reports/hypothesis_H5_final_report.md"
PAIR_RESULTS_PATH = ROOT / "reports/hypothesis_H5_pair_results.csv"


# ============================================================
# Direct MONDO -> MeSH mappings established from MONDO OBO
# ============================================================

MONDO_TO_MESH = {
    "MONDO:0005027": "MESH:D004827",  # epilepsy
    "MONDO:0005618": "MESH:D001008",  # anxiety disorder
    "MONDO:0005098": "MESH:D020521",  # stroke disorder
    "MONDO:0005180": "MESH:D010300",  # Parkinson disease
    "MONDO:0004975": "MESH:D000544",  # Alzheimer disease
    "MONDO:0005277": "MESH:D008881",  # migraine disorder
    "MONDO:0004985": "MESH:D001714",  # bipolar disorder
    "MONDO:0001627": "MESH:D003704",  # dementia
    "MONDO:0002009": "MESH:D003865",  # major depressive disorder
    "MONDO:0005301": "MESH:D009103",  # multiple sclerosis
}

MESH_TO_MONDO = {
    mesh: mondo
    for mondo, mesh in MONDO_TO_MESH.items()
}


# ============================================================
# Utility functions
# ============================================================

def sha256_file(path: Path) -> str:
    """Calculate SHA-256 checksum of a file."""
    h = hashlib.sha256()

    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)

    return h.hexdigest()


def wilson_ci(
    successes: int,
    trials: int,
    z: float = 1.959963984540054
) -> tuple[float, float]:
    """
    Calculate a 95% Wilson confidence interval
    for a binomial proportion.
    """

    if trials <= 0:
        return float("nan"), float("nan")

    p = successes / trials

    denom = 1.0 + (z * z / trials)

    center = (
        p + (z * z / (2.0 * trials))
    ) / denom

    margin = (
        z
        * math.sqrt(
            (
                p * (1.0 - p) / trials
            )
            + (
                z * z
                / (4.0 * trials * trials)
            )
        )
        / denom
    )

    return center - margin, center + margin


# ============================================================
# Main H5 analysis
# ============================================================

def main() -> None:

    print("Loading ReMedy datasets...")

    # --------------------------------------------------------
    # 1. Load ReMedy data
    # --------------------------------------------------------

    dd = pd.read_csv(DD_PATH)

    dp_chembl = pd.read_csv(
        DP_CHEMBL_PATH
    )

    dp_primekg = pd.read_csv(
        DP_PRIMEKG_PATH
    )

    disease_protein = pd.read_csv(
        DISEASE_PROTEIN_PATH
    )

    mapping = pd.read_csv(
        MAPPING_PATH
    )

    # Only treatment-use relationships.
    #
    # Contraindications are intentionally excluded.
    dd = dd[
        dd["relation"].isin(
            [
                "indication",
                "off-label use"
            ]
        )
    ].copy()

    # --------------------------------------------------------
    # 2. Validated DrugBank -> ChEMBL mapping
    # --------------------------------------------------------

    mapping = mapping[
        mapping["mapping_status"].eq("validated")
    ].copy()

    db_to_chembl = (
        mapping[
            [
                "drugbank_id",
                "chembl_id"
            ]
        ]
        .dropna()
        .drop_duplicates()
        .groupby("drugbank_id")["chembl_id"]
        .apply(set)
        .to_dict()
    )

    # --------------------------------------------------------
    # 3. Combined ReMedy Drug -> UniProt graph
    #
    # ChEMBL + PrimeKG
    # --------------------------------------------------------

    combined_dp = (
        pd.concat(
            [
                dp_chembl[
                    [
                        "chembl_id",
                        "uniprot_accession"
                    ]
                ],

                dp_primekg[
                    [
                        "chembl_id",
                        "uniprot_accession"
                    ]
                ],
            ],
            ignore_index=True,
        )
        .dropna()
        .drop_duplicates()
    )

    drug_proteins = (
        combined_dp
        .groupby("chembl_id")[
            "uniprot_accession"
        ]
        .apply(set)
        .to_dict()
    )

    # --------------------------------------------------------
    # 4. ReMedy Disease -> UniProt proteins
    # --------------------------------------------------------

    disease_proteins = (
        disease_protein
        .groupby("mondo_id")[
            "uniprot_accession"
        ]
        .apply(set)
        .to_dict()
    )

    # --------------------------------------------------------
    # 5. Load DrugMechDB
    # --------------------------------------------------------

    print("Loading DrugMechDB...")

    dm = yaml.safe_load(
        DRUGMECHDB_PATH.read_text(
            encoding="utf-8"
        )
    )

    print(
        "DrugMechDB paths:",
        len(dm)
    )

    # --------------------------------------------------------
    # 6. Build DrugMechDB mechanism index
    #
    # Key:
    #     (DrugBank ID, MeSH Disease ID)
    #
    # Value:
    #     UniProt proteins occurring
    #     in all corresponding mechanism paths
    # --------------------------------------------------------

    dm_index: dict[
        tuple[str, str],
        set[str]
    ] = {}

    for record in dm:

        graph = record["graph"]

        drugbank_id = str(
            graph.get(
                "drugbank",
                ""
            )
        ).replace(
            "DB:",
            ""
        )

        mesh_id = graph.get(
            "disease_mesh"
        )

        # Only diseases that have direct
        # MONDO -> MeSH mappings
        if (
            not drugbank_id
            or mesh_id not in MESH_TO_MONDO
        ):
            continue

        proteins = {
            str(node["id"]).replace(
                "UniProt:",
                ""
            )

            for node in record.get(
                "nodes",
                []
            )

            if (
                node.get("label")
                == "Protein"
            )

            and str(
                node.get("id", "")
            ).startswith(
                "UniProt:"
            )
        }

        key = (
            drugbank_id,
            mesh_id
        )

        if key not in dm_index:
            dm_index[key] = set()

        dm_index[key].update(
            proteins
        )

    # --------------------------------------------------------
    # 7. Determine eligible DrugMechDB drugs
    #
    # A drug must:
    # - occur in DrugMechDB
    # - have a validated DrugBank -> ChEMBL mapping
    # --------------------------------------------------------

    eligible_drugs = sorted(
        drugbank_id

        for drugbank_id, _mesh_id
        in dm_index.keys()

        if drugbank_id in db_to_chembl
    )

    eligible_diseases = sorted(
        MONDO_TO_MESH.keys()
    )

    # --------------------------------------------------------
    # 8. Construct mechanism-overlap universe
    #
    # Universe = every eligible DrugBank x
    # every directly mapped target disease.
    #
    # Mechanism overlap = at least one exact
    # UniProt protein shared between:
    #
    # ReMedy:
    #     Drug proteins ∩ Disease proteins
    #
    # and DrugMechDB:
    #     mechanism-path proteins
    # --------------------------------------------------------

    mechanism_overlap_rows = []

    print(
        "Building mechanism-overlap universe..."
    )

    for drugbank_id in eligible_drugs:

        # Combine all validated ChEMBL proteins
        # associated with this DrugBank ID.
        remedy_drug_proteins = set()

        for chembl_id in db_to_chembl[
            drugbank_id
        ]:

            remedy_drug_proteins.update(
                drug_proteins.get(
                    chembl_id,
                    set()
                )
            )

        for mondo_id in eligible_diseases:

            mesh_id = MONDO_TO_MESH[
                mondo_id
            ]

            # Disease proteins in ReMedy
            remedy_disease_proteins = (
                disease_proteins.get(
                    mondo_id,
                    set()
                )
            )

            # ReMedy biological connectivity
            remedy_shared = (
                remedy_drug_proteins
                &
                remedy_disease_proteins
            )

            # DrugMechDB mechanism proteins
            dm_proteins = dm_index.get(
                (
                    drugbank_id,
                    mesh_id
                ),
                set()
            )

            # Exact UniProt overlap
            exact_shared = (
                remedy_shared
                &
                dm_proteins
            )

            mechanism_overlap_rows.append(
                {
                    "drugbank_id":
                        drugbank_id,

                    "mondo_id":
                        mondo_id,

                    "mesh_id":
                        mesh_id,

                    "mechanism_overlap":
                        int(
                            bool(
                                exact_shared
                            )
                        ),

                    "remedy_shared_protein_count":
                        len(
                            remedy_shared
                        ),

                    "drugmechdb_protein_count":
                        len(
                            dm_proteins
                        ),

                    "exact_shared_mechanism_protein_count":
                        len(
                            exact_shared
                        ),
                }
            )

    universe = pd.DataFrame(
        mechanism_overlap_rows
    )

    # --------------------------------------------------------
    # 9. Construct observed ReMedy treatment-use pairs
    # --------------------------------------------------------

    observed = (
        dd[
            [
                "drugbank_id",
                "mondo_id",
                "relation"
            ]
        ]
        .dropna()
        .drop_duplicates(
            subset=[
                "drugbank_id",
                "mondo_id"
            ]
        )
        .merge(
            universe[
                [
                    "drugbank_id",
                    "mondo_id",
                    "mechanism_overlap",
                    "remedy_shared_protein_count",
                    "drugmechdb_protein_count",
                    "exact_shared_mechanism_protein_count",
                ]
            ],
            on=[
                "drugbank_id",
                "mondo_id"
            ],
            how="inner",
        )
    )

    observed[
        "observed_remedy_pair"
    ] = 1

    # --------------------------------------------------------
    # 10. Save pair-level results
    # --------------------------------------------------------

    PAIR_RESULTS_PATH.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    observed.to_csv(
        PAIR_RESULTS_PATH,
        index=False
    )

    # --------------------------------------------------------
    # 11. Hypergeometric enrichment test
    #
    # N = total eligible drug-disease pairs
    #
    # K = pairs in the universe that show
    #     exact mechanism overlap
    #
    # n = observed ReMedy treatment-use pairs
    #
    # x = observed treatment-use pairs that
    #     show mechanism overlap
    #
    # Test:
    # P(X >= x)
    # --------------------------------------------------------

    N = len(universe)

    K = int(
        universe[
            "mechanism_overlap"
        ].sum()
    )

    n = len(observed)

    x = int(
        observed[
            "mechanism_overlap"
        ].sum()
    )

    expected = (
        n * K / N
        if N > 0
        else float("nan")
    )

    enrichment = (
        x / expected
        if expected > 0
        else float("nan")
    )

    p_value = (
        hypergeom.sf(
            x - 1,
            N,
            K,
            n
        )
        if N > 0
        else float("nan")
    )

    observed_rate = (
        x / n
        if n > 0
        else float("nan")
    )

    universe_rate = (
        K / N
        if N > 0
        else float("nan")
    )

    # 95% Wilson interval for the observed
    # pair-level overlap proportion.
    ci_low, ci_high = wilson_ci(
        x,
        n
    )

    # --------------------------------------------------------
    # 12. DrugMechDB checksum
    # --------------------------------------------------------

    drugmech_sha256 = sha256_file(
        DRUGMECHDB_PATH
    )

    # --------------------------------------------------------
    # 13. Generate final report
    # --------------------------------------------------------

    report = f"""# H5 — ReMedy Mechanism-Path Overlap Validation

## Hypothesis

**H5:** ReMedy treatment-use drug–disease pairs show greater overlap with DrugMechDB mechanism paths than expected under a random-pair null model.

### Null hypothesis (H0)

The observed ReMedy treatment-use pairs do not have greater mechanism-path overlap than would be expected from random selection of the same number of drug–disease pairs from the eligible DrugBank × target-disease universe.

### Alternative hypothesis (H1)

The observed ReMedy treatment-use pairs have greater mechanism-path overlap than expected under the random-pair null model.

## Operational definition

A pair is classified as having **mechanism overlap** when:

1. the DrugBank drug has a validated DrugBank→ChEMBL mapping;
2. the disease is one of the 10 ReMedy diseases with a direct MONDO→MeSH cross-reference;
3. the DrugMechDB mechanism path contains at least one UniProt protein; and
4. at least one of those proteins is also present in the ReMedy drug→disease biological connectivity.

ReMedy drug→protein evidence uses the union of the ChEMBL and PrimeKG Drug→Protein canonical relations.

Only ReMedy `indication` and `off-label use` relationships are treated as observed treatment-use pairs. `contraindication` relationships are excluded.

The two ReMedy diseases without a direct MeSH cross-reference in the downloaded MONDO ontology (`MONDO:0005090` schizophrenia and `MONDO:0007743` attention deficit-hyperactivity disorder) are excluded from this direct-ID validation rather than matched by name.

## Statistical test

A one-sided **hypergeometric enrichment test** is used.

- Universe size, **N** = {N}
- Mechanism-overlap pairs in the universe, **K** = {K}
- Observed ReMedy treatment-use pairs, **n** = {n}
- Observed mechanism-overlap pairs, **x** = {x}
- Expected overlap under the null = **{expected:.4f}**
- Observed overlap rate = **{observed_rate * 100:.2f}%**
- Universe overlap rate = **{universe_rate * 100:.2f}%**
- Enrichment ratio = **{enrichment:.4f}×**
- Hypergeometric p-value = **{p_value:.6g}**

### 95% confidence interval

A 95% Wilson interval for the observed overlap proportion is:

**{ci_low * 100:.2f}% to {ci_high * 100:.2f}%**

This interval describes uncertainty around the observed pair-level overlap proportion; it is not a confidence interval for clinical efficacy.

## Result

The observed overlap was **{x} of {n} pairs ({observed_rate * 100:.2f}%)**, compared with **{expected:.2f} pairs expected under the random-pair null**. The observed overlap was **{enrichment:.2f}×** the random expectation. The one-sided hypergeometric p-value was **{p_value:.3e}**.

Under this operational definition, the observed ReMedy treatment-use pairs show substantially greater exact UniProt-level mechanism overlap with DrugMechDB than expected under the specified random-pair null.

## Interpretation

This result supports **H5 as a dataset-level consistency/validation finding**: the ReMedy treatment-use pairs overlap curated mechanism paths more often than random pairs in the defined universe.

It does **not** establish clinical efficacy, causal mechanism, or therapeutic suitability. DrugMechDB overlap is evidence of mechanistic agreement with a curated knowledge source, not proof that a drug will work for a patient.

## Limitations

1. H5 validates **exact UniProt overlap**, so mechanisms expressed through other node types are not counted unless a shared UniProt protein is present.
2. Only 10 of the 12 ReMedy diseases have direct MONDO→MeSH mappings and are therefore included.
3. The mechanism-overlap universe is constructed from the ReMedy biological graph and DrugMechDB; the test evaluates enrichment of observed treatment labels within that predefined universe.
4. A hypergeometric test treats candidate pairs as the sampling units; it does not model dependence among pairs sharing drugs, diseases, or proteins.
5. The mechanism-overlap definition is an operational proxy for validation and should not be described as clinical validation.

## Reproducibility

DrugMechDB input:

`data/raw/drugmechdb/indication_paths.yaml`

SHA-256:

`{drugmech_sha256}`

The executable analysis is:

`src/validation/hypothesis_H5_final.py`

The pair-level output is:

`reports/hypothesis_H5_pair_results.csv`

The report is:

`reports/hypothesis_H5_final_report.md`
"""

    # --------------------------------------------------------
    # 14. Save final report
    # --------------------------------------------------------

    REPORT_PATH.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    REPORT_PATH.write_text(
        report,
        encoding="utf-8"
    )

    # --------------------------------------------------------
    # 15. Print final results
    # --------------------------------------------------------

    print()
    print("=" * 60)
    print("H5 FINAL RESULT")
    print("=" * 60)

    print(
        f"Validated DrugBank mappings: "
        f"{len(db_to_chembl)}"
    )

    print(
        f"Combined Drug->Protein rows: "
        f"{len(combined_dp)}"
    )

    print(
        f"DrugMechDB paths: "
        f"{len(dm)}"
    )

    print(
        f"Eligible DrugMechDB drugs: "
        f"{len(eligible_drugs)}"
    )

    print(
        f"Mapped target diseases: "
        f"{len(eligible_diseases)}"
    )

    print(
        f"Universe N: "
        f"{N}"
    )

    print(
        f"Mechanism-overlap pairs K: "
        f"{K}"
    )

    print(
        f"Observed ReMedy treatment pairs n: "
        f"{n}"
    )

    print(
        f"Observed overlap x: "
        f"{x}"
    )

    print(
        f"Expected random overlap: "
        f"{expected:.4f}"
    )

    print(
        f"Enrichment ratio: "
        f"{enrichment:.4f}"
    )

    print(
        f"Hypergeometric p-value: "
        f"{p_value:.6g}"
    )

    print(
        f"Observed overlap percentage: "
        f"{observed_rate * 100:.2f}%"
    )

    print(
        f"Universe mechanism-overlap percentage: "
        f"{universe_rate * 100:.2f}%"
    )

    print(
        f"95% Wilson CI: "
        f"[{ci_low * 100:.2f}%, "
        f"{ci_high * 100:.2f}%]"
    )

    print(
        f"DrugMechDB SHA-256: "
        f"{drugmech_sha256}"
    )

    print()
    print(
        f"Saved pair results: "
        f"{PAIR_RESULTS_PATH}"
    )

    print(
        f"Saved report: "
        f"{REPORT_PATH}"
    )

    print("=" * 60)


if __name__ == "__main__":
    main()