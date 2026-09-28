from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]

EVIDENCE_FILE = (
    ROOT / "data/interim/chembl_evidence_clean.csv"
)

REPORT_FILE = (
    ROOT / "reports/hypothesis_H4_final_report.md"
)

RESULT_FILE = (
    ROOT / "reports/hypothesis_H4_final_results.csv"
)

RANDOM_SEED = 42
PERMUTATIONS = 5000
BOOTSTRAP_ITERATIONS = 5000


# =========================================================
# BOOTSTRAP CI
# =========================================================

def bootstrap_mean_difference_ci(
    matched_scores,
    mismatched_scores,
    iterations=5000,
    seed=42,
):

    rng = np.random.default_rng(
        seed
    )

    matched = np.asarray(
        matched_scores,
        dtype=float,
    )

    mismatched = np.asarray(
        mismatched_scores,
        dtype=float,
    )

    differences = np.empty(
        iterations,
        dtype=float,
    )

    for i in range(iterations):

        matched_sample = rng.choice(
            matched,
            size=len(matched),
            replace=True,
        )

        mismatched_sample = rng.choice(
            mismatched,
            size=len(mismatched),
            replace=True,
        )

        differences[i] = (
            np.mean(matched_sample)
            -
            np.mean(mismatched_sample)
        )

    lower = np.percentile(
        differences,
        2.5,
    )

    upper = np.percentile(
        differences,
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
    print("REMEDY - H4 CONTEXT-MATCHED EVIDENCE TEST")
    print("=" * 70)

    # -----------------------------------------------------
    # 1. Load evidence
    # -----------------------------------------------------

    print(
        "\nLoading ChEMBL evidence..."
    )

    df = pd.read_csv(
        EVIDENCE_FILE,
        dtype=str,
    ).fillna("")

    print(
        f"Total evidence rows: "
        f"{len(df):,}"
    )

    # -----------------------------------------------------
    # 2. Define context groups
    # -----------------------------------------------------

    df["context_group"] = (
        df["assay_organism_class"]
        .map(
            {
                "human": "matched",
                "non_human": "mismatched",
            }
        )
    )

    # Exclude unspecified context from the primary H4 test.
    df = df[
        df["context_group"].notna()
    ].copy()

    print(
        f"Context-specified evidence rows: "
        f"{len(df):,}"
    )

    # -----------------------------------------------------
    # 3. Build informativeness indicators
    # -----------------------------------------------------

    df["has_pchembl"] = (
        pd.to_numeric(
            df["pchembl_value"],
            errors="coerce",
        )
        .notna()
        .astype(int)
    )

    df["has_publication"] = (
        df["publication"]
        .astype(str)
        .str.strip()
        .ne("")
        .astype(int)
    )

    df["has_cell_line"] = (
        df["cell_line"]
        .astype(str)
        .str.strip()
        .ne("")
        .astype(int)
    )

    # Equal-weight completeness score.
    df["informativeness_score"] = (
        df["has_pchembl"]
        +
        df["has_publication"]
        +
        df["has_cell_line"]
    )

    # -----------------------------------------------------
    # 4. Split groups
    # -----------------------------------------------------

    matched = df[
        df["context_group"] == "matched"
    ].copy()

    mismatched = df[
        df["context_group"] == "mismatched"
    ].copy()

    matched_scores = (
        matched[
            "informativeness_score"
        ]
        .astype(float)
        .values
    )

    mismatched_scores = (
        mismatched[
            "informativeness_score"
        ]
        .astype(float)
        .values
    )

    print(
        f"\nContext-matched rows: "
        f"{len(matched):,}"
    )

    print(
        f"Context-mismatched rows: "
        f"{len(mismatched):,}"
    )

    # -----------------------------------------------------
    # 5. Descriptive statistics
    # -----------------------------------------------------

    matched_mean = float(
        np.mean(
            matched_scores
        )
    )

    mismatched_mean = float(
        np.mean(
            mismatched_scores
        )
    )

    matched_median = float(
        np.median(
            matched_scores
        )
    )

    mismatched_median = float(
        np.median(
            mismatched_scores
        )
    )

    matched_q25 = float(
        np.percentile(
            matched_scores,
            25,
        )
    )

    matched_q75 = float(
        np.percentile(
            matched_scores,
            75,
        )
    )

    mismatched_q25 = float(
        np.percentile(
            mismatched_scores,
            25,
        )
    )

    mismatched_q75 = float(
        np.percentile(
            mismatched_scores,
            75,
        )
    )

    observed_difference = (
        matched_mean
        -
        mismatched_mean
    )

    # -----------------------------------------------------
    # 6. Permutation test
    #
    # H4 alternative:
    # matched informativeness > mismatched
    # -----------------------------------------------------

    rng = np.random.default_rng(
        RANDOM_SEED
    )

    combined_scores = np.concatenate(
        [
            matched_scores,
            mismatched_scores,
        ]
    )

    n_matched = len(
        matched_scores
    )

    randomized_differences = np.empty(
        PERMUTATIONS,
        dtype=float,
    )

    print(
        f"\nRunning {PERMUTATIONS:,} permutations..."
    )

    for i in range(
        PERMUTATIONS
    ):

        shuffled = rng.permutation(
            combined_scores
        )

        randomized_matched = (
            shuffled[
                :n_matched
            ]
        )

        randomized_mismatched = (
            shuffled[
                n_matched:
            ]
        )

        randomized_differences[i] = (
            np.mean(
                randomized_matched
            )
            -
            np.mean(
                randomized_mismatched
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
    # 7. Permutation p-value
    # -----------------------------------------------------

    extreme_count = int(
        np.sum(
            randomized_differences
            >= observed_difference
        )
    )

    p_value = (
        (extreme_count + 1)
        /
        (PERMUTATIONS + 1)
    )

    # -----------------------------------------------------
    # 8. Bootstrap 95% CI
    # -----------------------------------------------------

    print(
        "\nRunning bootstrap..."
    )

    ci_low, ci_high = (
        bootstrap_mean_difference_ci(
            matched_scores,
            mismatched_scores,
            iterations=BOOTSTRAP_ITERATIONS,
            seed=RANDOM_SEED,
        )
    )

    print(
        "Bootstrap complete."
    )

    # -----------------------------------------------------
    # 9. Standardized permutation effect
    # -----------------------------------------------------

    null_mean = float(
        np.mean(
            randomized_differences
        )
    )

    null_std = float(
        np.std(
            randomized_differences,
            ddof=1,
        )
    )

    if null_std > 0:

        standardized_effect = (
            observed_difference
            -
            null_mean
        ) / null_std

    else:

        standardized_effect = 0.0

    # -----------------------------------------------------
    # 10. Component-level rates
    # -----------------------------------------------------

    component_results = []

    for column, label in [
        (
            "has_pchembl",
            "Quantitative pChEMBL",
        ),
        (
            "has_publication",
            "Publication",
        ),
        (
            "has_cell_line",
            "Cell line",
        ),
    ]:

        matched_rate = (
            matched[column]
            .mean()
            * 100
        )

        mismatched_rate = (
            mismatched[column]
            .mean()
            * 100
        )

        component_results.append(
            {
                "component": label,
                "matched_percentage": matched_rate,
                "mismatched_percentage": mismatched_rate,
                "difference_percentage_points": (
                    matched_rate
                    -
                    mismatched_rate
                ),
            }
        )

    component_df = pd.DataFrame(
        component_results
    )

    # -----------------------------------------------------
    # 11. Save results
    # -----------------------------------------------------

    result_summary = pd.DataFrame(
        [
            {
                "metric": "rows",
                "matched": len(matched),
                "mismatched": len(mismatched),
            },
            {
                "metric": "mean_informativeness",
                "matched": matched_mean,
                "mismatched": mismatched_mean,
            },
            {
                "metric": "median_informativeness",
                "matched": matched_median,
                "mismatched": mismatched_median,
            },
            {
                "metric": "observed_difference",
                "matched": observed_difference,
                "mismatched": np.nan,
            },
            {
                "metric": "bootstrap_ci_low",
                "matched": ci_low,
                "mismatched": np.nan,
            },
            {
                "metric": "bootstrap_ci_high",
                "matched": ci_high,
                "mismatched": np.nan,
            },
            {
                "metric": "permutation_p_value",
                "matched": p_value,
                "mismatched": np.nan,
            },
        ]
    )

    component_export = component_df.copy()

    result_export = pd.concat(
        [
            result_summary,
        ],
        ignore_index=True,
    )

    RESULT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    result_export.to_csv(
        RESULT_FILE,
        index=False,
    )

    # -----------------------------------------------------
    # 12. Interpretation
    # -----------------------------------------------------

    if p_value < 0.05:

        interpretation = (
            "Under the specified evidence-completeness "
            "proxy and permutation test, context-matched "
            "human assay records have a higher mean "
            "informativeness score than context-mismatched "
            "non-human assay records."
        )

    else:

        interpretation = (
            "The permutation test does not provide "
            "sufficient evidence that context-matched "
            "human assay records have a higher mean "
            "informativeness score than context-mismatched "
            "non-human assay records."
        )

    # -----------------------------------------------------
    # 13. Report
    # -----------------------------------------------------

    report = f"""# ReMedy H4 — Final KG

## Hypothesis

**H4:** Context-matched evidence is more informative than
context-mismatched evidence.

## Context definition

For the primary analysis:

- **Context-matched:** human assay organism
- **Context-mismatched:** non-human assay organism
- **Unspecified:** excluded from the primary H4 test

The dataset contains evidence records with unspecified assay
organism; these are not forced into either comparison group.

## Informativeness definition

The current dataset does not contain an independent ground-truth
informativeness label.

Therefore, H4 uses an explicit evidence-completeness proxy:

**Informativeness score =**

1. quantitative pChEMBL value present
2. publication information present
3. cell-line information present

Each available component contributes one point.

Score range:

`0–3`

All three components receive equal weight.

## Null hypothesis

Context-matched and context-mismatched records have the same
mean informativeness score.

## Alternative hypothesis

Context-matched records have a higher mean informativeness
score than context-mismatched records.

## Dataset

- Total ChEMBL evidence records: {len(pd.read_csv(EVIDENCE_FILE, dtype=str)):,}
- Context-matched records: {len(matched):,}
- Context-mismatched records: {len(mismatched):,}
- Unspecified records excluded: {
    len(pd.read_csv(EVIDENCE_FILE, dtype=str))
    - len(df)
:,}

## Informativeness results

| Metric | Context-matched | Context-mismatched |
|---|---:|---:|
| Records | {len(matched):,} | {len(mismatched):,} |
| Mean score | {matched_mean:.4f} | {mismatched_mean:.4f} |
| Median score | {matched_median:.4f} | {mismatched_median:.4f} |
| Q1 | {matched_q25:.4f} | {mismatched_q25:.4f} |
| Q3 | {matched_q75:.4f} | {mismatched_q75:.4f} |

## Component availability

{component_df.to_markdown(index=False)}

## Statistical test

Permutation test comparing:

`mean(context-matched score) - mean(context-mismatched score)`

The context labels were randomly permuted while preserving
the observed group sizes.

Number of permutations:

`{PERMUTATIONS:,}`

## Results

- Observed mean difference: {observed_difference:.6f}
- Observed difference in score units: {observed_difference:.4f}
- Null mean difference: {null_mean:.6f}
- Null standard deviation: {null_std:.6f}
- p-value: {p_value:.6g}
- Permutations >= observed: {extreme_count:,}/{PERMUTATIONS:,}

## Effect size

Primary effect size:

`{observed_difference:.6f}`

Standardized permutation statistic:

`{standardized_effect:.4f}`

## 95% bootstrap confidence interval

`[{ci_low:.6f}, {ci_high:.6f}]`

## Interpretation

{interpretation}

## Important limitation

The informativeness score is an operational evidence-completeness
proxy, not a validated measure of scientific or biological
informativeness.

In particular, the availability of a publication, pChEMBL value,
or cell line does not by itself establish that an experiment is
more scientifically informative.

Therefore, H4 should be reported as a test of contextual evidence
completeness under the defined proxy.

## Reproducibility

- Random seed: {RANDOM_SEED}
- Permutations: {PERMUTATIONS:,}
- Bootstrap iterations: {BOOTSTRAP_ITERATIONS:,}

## Output

`reports/hypothesis_H4_final_results.csv`
"""

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    # -----------------------------------------------------
    # 14. Terminal
    # -----------------------------------------------------

    print()
    print("=" * 70)
    print("H4 FINAL-KG TEST COMPLETE")
    print("=" * 70)

    print(
        f"Context-matched rows: "
        f"{len(matched):,}"
    )

    print(
        f"Context-mismatched rows: "
        f"{len(mismatched):,}"
    )

    print(
        f"Matched mean informativeness: "
        f"{matched_mean:.4f}"
    )

    print(
        f"Mismatched mean informativeness: "
        f"{mismatched_mean:.4f}"
    )

    print(
        f"Mean difference: "
        f"{observed_difference:.6f}"
    )

    print(
        f"p-value: "
        f"{p_value:.6g}"
    )

    print(
        f"Standardized permutation statistic: "
        f"{standardized_effect:.4f}"
    )

    print(
        f"95% CI: "
        f"[{ci_low:.6f}, {ci_high:.6f}]"
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