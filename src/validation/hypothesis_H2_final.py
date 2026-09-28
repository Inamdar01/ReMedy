from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

EDGES_FILE = ROOT / "data/processed/edges.parquet"

REPORT_FILE = (
    ROOT / "reports/hypothesis_H2_final_report.md"
)

RESULT_FILE = (
    ROOT / "reports/hypothesis_H2_final_results.csv"
)

RANDOM_SEED = 42
PERMUTATIONS = 5000
BOOTSTRAP_ITERATIONS = 5000


# =========================================================
# BUILD BIOLOGICAL RELATION SETS
# =========================================================

def build_relation_sets(edges):

    drug_proteins = (
        edges[
            edges["relation"] == "targets"
        ]
        .groupby("source_id")["target_id"]
        .apply(set)
        .to_dict()
    )

    disease_proteins = (
        edges[
            edges["relation"] == "associated_with"
        ]
        .groupby("source_id")["target_id"]
        .apply(set)
        .to_dict()
    )

    protein_pathways = (
        edges[
            edges["relation"] == "participates_in"
        ]
        .groupby("source_id")["target_id"]
        .apply(set)
        .to_dict()
    )

    return (
        drug_proteins,
        disease_proteins,
        protein_pathways,
    )


def derive_pathways(
    entity_proteins,
    protein_pathways,
):

    entity_pathways = {}

    for entity, proteins in entity_proteins.items():

        pathways = set()

        for protein in proteins:

            pathways.update(
                protein_pathways.get(
                    protein,
                    set(),
                )
            )

        entity_pathways[entity] = pathways

    return entity_pathways


def path_support(
    drug,
    disease,
    drug_proteins,
    disease_proteins,
    drug_pathways,
    disease_pathways,
):

    shared_proteins = (
        drug_proteins.get(
            drug,
            set(),
        )
        &
        disease_proteins.get(
            disease,
            set(),
        )
    )

    shared_pathways = (
        drug_pathways.get(
            drug,
            set(),
        )
        &
        disease_pathways.get(
            disease,
            set(),
        )
    )

    return int(
        len(shared_proteins) > 0
        or
        len(shared_pathways) > 0
    )


# =========================================================
# DEGREE-PRESERVING RANDOMIZATION
# =========================================================

def randomize_known_pairs(
    pairs_by_drug,
    drugs,
    diseases,
    rng,
):

    randomized_pairs = set()

    for drug in drugs:

        degree = len(
            pairs_by_drug.get(
                drug,
                [],
            )
        )

        if degree == 0:
            continue

        sampled_diseases = rng.choice(
            diseases,
            size=degree,
            replace=False,
        )

        for disease in sampled_diseases:

            randomized_pairs.add(
                (
                    drug,
                    disease,
                )
            )

    return randomized_pairs


# =========================================================
# SUPPORT RATE
# =========================================================

def support_rate(
    pairs,
    support_lookup,
):

    if not pairs:
        return 0.0

    supported = sum(
        support_lookup.get(
            pair,
            0,
        )
        for pair in pairs
    )

    return supported / len(pairs)


# =========================================================
# CLUSTER BOOTSTRAP
# =========================================================

def bootstrap_difference_ci(
    pairs_by_drug,
    drugs,
    diseases,
    support_lookup,
    iterations,
    seed,
):

    rng = np.random.default_rng(seed)

    bootstrap_differences = np.empty(
        iterations,
        dtype=float,
    )

    drug_count = len(drugs)

    for i in range(iterations):

        sampled_drugs = rng.choice(
            drugs,
            size=drug_count,
            replace=True,
        )

        observed_numerator = 0.0
        observed_denominator = 0.0

        randomized_numerator = 0.0
        randomized_denominator = 0.0

        for drug in sampled_drugs:

            observed_diseases = pairs_by_drug[
                drug
            ]

            degree = len(
                observed_diseases
            )

            if degree == 0:
                continue

            # Weight contribution by bootstrap
            # multiplicity of the sampled drug.
            observed_support = sum(
                support_lookup[
                    (
                        drug,
                        disease,
                    )
                ]
                for disease in observed_diseases
            )

            observed_numerator += (
                observed_support
            )

            observed_denominator += (
                degree
            )

            # Degree-preserving randomized diseases
            randomized_diseases = rng.choice(
                diseases,
                size=degree,
                replace=False,
            )

            randomized_support = sum(
                support_lookup[
                    (
                        drug,
                        disease,
                    )
                ]
                for disease in randomized_diseases
            )

            randomized_numerator += (
                randomized_support
            )

            randomized_denominator += (
                degree
            )

        if (
            observed_denominator > 0
            and randomized_denominator > 0
        ):

            observed_rate = (
                observed_numerator
                /
                observed_denominator
            )

            randomized_rate = (
                randomized_numerator
                /
                randomized_denominator
            )

            bootstrap_differences[i] = (
                observed_rate
                -
                randomized_rate
            )

        else:

            bootstrap_differences[i] = np.nan

    bootstrap_differences = (
        bootstrap_differences[
            np.isfinite(
                bootstrap_differences
            )
        ]
    )

    lower = np.percentile(
        bootstrap_differences,
        2.5,
    )

    upper = np.percentile(
        bootstrap_differences,
        97.5,
    )

    return (
        float(lower),
        float(upper),
    )


# =========================================================
# MAIN
# =========================================================

def main():

    print("=" * 70)
    print("REMEDY - H2 TEST ON FINAL PARQUET KG")
    print("=" * 70)

    # -----------------------------------------------------
    # 1. Load final KG
    # -----------------------------------------------------

    print("\nLoading final KG...")

    edges = pd.read_parquet(
        EDGES_FILE
    )

    print(
        f"Edges: {len(edges):,}"
    )

    # -----------------------------------------------------
    # 2. Extract Drug-Disease relationships
    # -----------------------------------------------------

    drug_disease = edges[
        edges["relation"].isin(
            [
                "indication",
                "off-label use",
                "contraindication",
            ]
        )
    ].copy()

    drugs = sorted(
        drug_disease[
            "source_id"
        ].unique()
    )

    diseases = sorted(
        drug_disease[
            "target_id"
        ].unique()
    )

    known_pairs = set(
        zip(
            drug_disease[
                "source_id"
            ],
            drug_disease[
                "target_id"
            ],
        )
    )

    print(
        f"Unique drugs: {len(drugs):,}"
    )

    print(
        f"Unique diseases: {len(diseases):,}"
    )

    print(
        f"Known Drug-Disease pairs: "
        f"{len(known_pairs):,}"
    )

    # -----------------------------------------------------
    # 3. Drug degree structure
    # -----------------------------------------------------

    pairs_by_drug = {}

    for drug, disease in known_pairs:

        pairs_by_drug.setdefault(
            drug,
            []
        ).append(disease)

    # -----------------------------------------------------
    # 4. Biological relations
    # -----------------------------------------------------

    print(
        "\nBuilding protein/pathway sets..."
    )

    (
        drug_proteins,
        disease_proteins,
        protein_pathways,
    ) = build_relation_sets(
        edges
    )

    drug_pathways = derive_pathways(
        drug_proteins,
        protein_pathways,
    )

    disease_pathways = derive_pathways(
        disease_proteins,
        protein_pathways,
    )

    # -----------------------------------------------------
    # 5. Precompute support for every Drug-Disease pair
    # -----------------------------------------------------

    print(
        "Calculating pair-level biological support..."
    )

    support_lookup = {}

    pair_rows = []

    for drug in drugs:

        for disease in diseases:

            support = path_support(
                drug,
                disease,
                drug_proteins,
                disease_proteins,
                drug_pathways,
                disease_pathways,
            )

            pair = (
                drug,
                disease,
            )

            support_lookup[pair] = support

            pair_rows.append(
                {
                    "chembl_id": drug,
                    "mondo_id": disease,
                    "path_support": support,
                    "known_pair": int(
                        pair in known_pairs
                    ),
                }
            )

    pair_table = pd.DataFrame(
        pair_rows
    )

    # -----------------------------------------------------
    # 6. Observed rate
    # -----------------------------------------------------

    observed_rate = support_rate(
        known_pairs,
        support_lookup,
    )

    observed_supported = sum(
        support_lookup[pair]
        for pair in known_pairs
    )

    print(
        f"\nObserved path-support rate: "
        f"{observed_rate * 100:.2f}%"
    )

    print(
        f"Observed supported pairs: "
        f"{observed_supported:,}/"
        f"{len(known_pairs):,}"
    )

    # -----------------------------------------------------
    # 7. Permutation test
    # -----------------------------------------------------

    rng = np.random.default_rng(
        RANDOM_SEED
    )

    randomized_rates = np.empty(
        PERMUTATIONS,
        dtype=float,
    )

    print(
        f"\nRunning {PERMUTATIONS:,} permutations..."
    )

    for i in range(
        PERMUTATIONS
    ):

        randomized_pairs = (
            randomize_known_pairs(
                pairs_by_drug,
                drugs,
                diseases,
                rng,
            )
        )

        randomized_rates[i] = (
            support_rate(
                randomized_pairs,
                support_lookup,
            )
        )

        if (
            (i + 1) % 1000 == 0
        ):

            print(
                f"Completed "
                f"{i + 1:,}/"
                f"{PERMUTATIONS:,}"
            )

    # -----------------------------------------------------
    # 8. Permutation p-value
    # -----------------------------------------------------

    extreme_count = int(
        np.sum(
            randomized_rates
            >= observed_rate
        )
    )

    p_value = (
        (extreme_count + 1)
        /
        (PERMUTATIONS + 1)
    )

    # -----------------------------------------------------
    # 9. Null distribution
    # -----------------------------------------------------

    null_mean = float(
        np.mean(
            randomized_rates
        )
    )

    null_std = float(
        np.std(
            randomized_rates,
            ddof=1,
        )
    )

    null_q025 = float(
        np.percentile(
            randomized_rates,
            2.5,
        )
    )

    null_q975 = float(
        np.percentile(
            randomized_rates,
            97.5,
        )
    )

    observed_minus_null = (
        observed_rate
        -
        null_mean
    )

    # -----------------------------------------------------
    # 10. Primary effect size
    # -----------------------------------------------------

    primary_effect_size = (
        observed_minus_null
    )

    if null_std > 0:

        standardized_effect = (
            observed_minus_null
            /
            null_std
        )

    else:

        standardized_effect = 0.0

    # -----------------------------------------------------
    # 11. Correct 95% bootstrap CI
    #
    # Cluster = drug
    # -----------------------------------------------------

    print(
        "\nRunning drug-level bootstrap..."
    )

    ci_low, ci_high = (
        bootstrap_difference_ci(
            pairs_by_drug,
            drugs,
            diseases,
            support_lookup,
            BOOTSTRAP_ITERATIONS,
            RANDOM_SEED,
        )
    )

    print(
        "Bootstrap complete."
    )

    # -----------------------------------------------------
    # 12. Permutation exceedance
    # -----------------------------------------------------

    exceedance_rate = (
        extreme_count
        /
        PERMUTATIONS
        *
        100
    )

    # -----------------------------------------------------
    # 13. Save pair-level table
    # -----------------------------------------------------

    RESULT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    pair_table.to_csv(
        RESULT_FILE,
        index=False,
    )

    # -----------------------------------------------------
    # 14. Interpretation
    # -----------------------------------------------------

    if p_value < 0.05:

        interpretation = (
            "Under the specified degree-preserving "
            "permutation null model, the observed "
            "path-support rate is higher than the "
            "randomized rates."
        )

    else:

        interpretation = (
            "The degree-preserving permutation test "
            "does not provide sufficient evidence that "
            "the observed path-support rate exceeds the "
            "randomized null distribution."
        )

    # -----------------------------------------------------
    # 15. Report
    # -----------------------------------------------------

    report = f"""# ReMedy H2 — Final KG

## Hypothesis

**H2:** Drugs and diseases sharing proteins/pathways are more
connected than unrelated pairs.

## Null hypothesis

The observed Drug-Disease pairs do not have a higher
protein/pathway path-support rate than expected when disease
assignments are randomized while preserving each drug's
number of disease links.

## Alternative hypothesis

The observed Drug-Disease pairs have a higher
protein/pathway path-support rate than the degree-preserving
randomized null.

## Final KG

- Total edges: {len(edges):,}
- Unique drugs: {len(drugs):,}
- Unique diseases: {len(diseases):,}
- Known Drug-Disease pairs: {len(known_pairs):,}

## Path-support definition

A Drug-Disease pair has:

`path_support = 1`

when the drug and disease share at least one UniProt protein
or at least one Reactome pathway reachable from their
associated proteins.

Otherwise:

`path_support = 0`

## Permutation procedure

For every drug, its observed number of disease links was
preserved.

The disease identities were randomly reassigned from the
{len(diseases):,}-disease set.

This was repeated {PERMUTATIONS:,} times.

## Observed result

- Supported known pairs: {observed_supported:,}
- Total known pairs: {len(known_pairs):,}
- Observed support rate: {observed_rate * 100:.2f}%

## Permutation result

- Null mean: {null_mean * 100:.2f}%
- Null standard deviation: {null_std:.6f}
- Observed minus null: {observed_minus_null * 100:.2f} percentage points
- Permutation p-value: {p_value:.6g}
- Null statistics >= observed: {extreme_count:,}/{PERMUTATIONS:,}
- Permutation exceedance rate: {exceedance_rate:.2f}%

## Effect size

Primary effect size:

Observed support rate minus randomized-null mean:

`{primary_effect_size:.6f}`

In percentage-point terms:

`{primary_effect_size * 100:.2f}`

Standardized permutation statistic:

`{standardized_effect:.4f}`

## 95% bootstrap confidence interval

A drug-level clustered bootstrap was used because multiple
Drug-Disease observations share the same drug.

95% bootstrap CI for the observed-minus-randomized
support-rate difference:

`[{ci_low:.6f}, {ci_high:.6f}]`

Percentage-point form:

`[{ci_low * 100:.2f}, {ci_high * 100:.2f}]`

## 95% permutation null interval

For descriptive comparison, the central 95% of the randomized
support-rate distribution was:

`[{null_q025:.6f}, {null_q975:.6f}]`

Percentage-point form:

`[{null_q025 * 100:.2f}, {null_q975 * 100:.2f}]`

This is a permutation null interval, not a confidence interval.

## Interpretation

{interpretation}

The 0.01/0.05 significance threshold is not interpreted as
evidence of clinical efficacy. H2 evaluates biological
connectivity in the integrated knowledge graph only.

## Reproducibility

- Random seed: {RANDOM_SEED}
- Permutations: {PERMUTATIONS}
- Bootstrap iterations: {BOOTSTRAP_ITERATIONS}
- Bootstrap unit: drug

## Output

`reports/hypothesis_H2_final_results.csv`
"""

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    # -----------------------------------------------------
    # 16. Terminal output
    # -----------------------------------------------------

    print()
    print("=" * 70)
    print("H2 FINAL-KG TEST COMPLETE")
    print("=" * 70)

    print(
        f"Known pairs: "
        f"{len(known_pairs):,}"
    )

    print(
        f"Observed support rate: "
        f"{observed_rate * 100:.2f}%"
    )

    print(
        f"Randomized null mean: "
        f"{null_mean * 100:.2f}%"
    )

    print(
        f"Difference: "
        f"{observed_minus_null * 100:.2f} "
        f"percentage points"
    )

    print(
        f"p-value: "
        f"{p_value:.6g}"
    )

    print(
        f"Primary effect size: "
        f"{primary_effect_size:.6f}"
    )

    print(
        f"Standardized permutation statistic: "
        f"{standardized_effect:.4f}"
    )

    print(
        f"95% bootstrap CI: "
        f"[{ci_low:.6f}, {ci_high:.6f}]"
    )

    print(
        f"Null 95% interval: "
        f"[{null_q025:.6f}, {null_q975:.6f}]"
    )

    print(
        f"Null statistics >= observed: "
        f"{extreme_count:,}/"
        f"{PERMUTATIONS:,}"
    )

    print()
    print(
        f"Saved: {REPORT_FILE}"
    )

    print(
        f"Saved: {RESULT_FILE}"
    )

    print("=" * 70)


if __name__ == "__main__":
    main()