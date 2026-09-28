from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parents[2]

EDGES_FILE = ROOT / "data/processed/edges.parquet"
EVIDENCE_FILE = ROOT / "data/processed/evidence.parquet"

REPORT_FILE = (
    ROOT / "reports/hypothesis_H3_final_report.md"
)

RESULT_FILE = (
    ROOT / "reports/hypothesis_H3_final_results.csv"
)

RANDOM_SEED = 42
BOOTSTRAP_ITERATIONS = 5000


# =========================================================
# BOOTSTRAP MEDIAN DIFFERENCE
# =========================================================

def bootstrap_median_difference_ci(
    multi_scores,
    single_scores,
    iterations=5000,
    seed=42,
):

    rng = np.random.default_rng(seed)

    multi_scores = np.asarray(
        multi_scores,
        dtype=float,
    )

    single_scores = np.asarray(
        single_scores,
        dtype=float,
    )

    differences = np.empty(
        iterations,
        dtype=float,
    )

    for i in range(iterations):

        multi_sample = rng.choice(
            multi_scores,
            size=len(multi_scores),
            replace=True,
        )

        single_sample = rng.choice(
            single_scores,
            size=len(single_scores),
            replace=True,
        )

        differences[i] = (
            np.median(multi_sample)
            -
            np.median(single_sample)
        )

    return (
        float(
            np.percentile(
                differences,
                2.5,
            )
        ),
        float(
            np.percentile(
                differences,
                97.5,
            )
        ),
    )


# =========================================================
# MAIN
# =========================================================

def main():

    print("=" * 70)
    print("REMEDY - H3 MULTI-SOURCE EVIDENCE TEST")
    print("=" * 70)

    # -----------------------------------------------------
    # 1. Load final KG
    # -----------------------------------------------------

    print("\nLoading final KG...")

    edges = pd.read_parquet(
        EDGES_FILE
    )

    evidence = pd.read_parquet(
        EVIDENCE_FILE
    )

    print(
        f"Edges: {len(edges):,}"
    )

    print(
        f"Evidence rows: {len(evidence):,}"
    )

    # -----------------------------------------------------
    # 2. Count distinct evidence sources per edge
    # -----------------------------------------------------

    source_counts = (
        evidence
        .groupby("edge_id")["source"]
        .nunique()
        .rename("source_count")
        .reset_index()
    )

    analysis = edges.merge(
        source_counts,
        on="edge_id",
        how="left",
    )

    analysis["source_count"] = (
        analysis["source_count"]
        .fillna(0)
        .astype(int)
    )

    # Every current canonical edge has evidence,
    # but retain a defensive filter.
    analysis = analysis[
        analysis["source_count"] > 0
    ].copy()

    # -----------------------------------------------------
    # 3. Classify source support
    # -----------------------------------------------------

    analysis["evidence_class"] = np.where(
        analysis["source_count"] >= 2,
        "multi-source",
        "single-source",
    )

    # -----------------------------------------------------
    # 4. Evidence counts
    # -----------------------------------------------------

    analysis["evidence_count"] = (
        analysis["evidence_count"]
        .astype(int)
    )

    # -----------------------------------------------------
    # 5. Extract groups
    # -----------------------------------------------------

    multi = analysis[
        analysis["evidence_class"]
        == "multi-source"
    ].copy()

    single = analysis[
        analysis["evidence_class"]
        == "single-source"
    ].copy()

    multi_scores = (
        multi["evidence_count"]
        .astype(float)
        .values
    )

    single_scores = (
        single["evidence_count"]
        .astype(float)
        .values
    )

    print()
    print(
        f"Single-source edges: "
        f"{len(single):,}"
    )

    print(
        f"Multi-source edges: "
        f"{len(multi):,}"
    )

    # -----------------------------------------------------
    # 6. Descriptive statistics
    # -----------------------------------------------------

    multi_median = float(
        np.median(
            multi_scores
        )
    )

    single_median = float(
        np.median(
            single_scores
        )
    )

    multi_mean = float(
        np.mean(
            multi_scores
        )
    )

    single_mean = float(
        np.mean(
            single_scores
        )
    )

    multi_q25 = float(
        np.percentile(
            multi_scores,
            25,
        )
    )

    multi_q75 = float(
        np.percentile(
            multi_scores,
            75,
        )
    )

    single_q25 = float(
        np.percentile(
            single_scores,
            25,
        )
    )

    single_q75 = float(
        np.percentile(
            single_scores,
            75,
        )
    )

    # -----------------------------------------------------
    # 7. Mann-Whitney U
    #
    # H3 directional alternative:
    # multi-source evidence_count > single-source
    # -----------------------------------------------------

    u_statistic, p_value = (
        mannwhitneyu(
            multi_scores,
            single_scores,
            alternative="greater",
            method="asymptotic",
        )
    )

    n_multi = len(
        multi_scores
    )

    n_single = len(
        single_scores
    )

    # Rank-biserial effect size.
    rank_biserial = (
        (
            2 * u_statistic
        )
        /
        (
            n_multi
            *
            n_single
        )
    ) - 1

    # -----------------------------------------------------
    # 8. Median difference
    # -----------------------------------------------------

    median_difference = (
        multi_median
        -
        single_median
    )

    # -----------------------------------------------------
    # 9. Bootstrap confidence interval
    # -----------------------------------------------------

    ci_low, ci_high = (
        bootstrap_median_difference_ci(
            multi_scores,
            single_scores,
            iterations=BOOTSTRAP_ITERATIONS,
            seed=RANDOM_SEED,
        )
    )

    # -----------------------------------------------------
    # 10. Additional evidence-support metric
    # -----------------------------------------------------

    multi_high_support = (
        multi["evidence_count"]
        >= 2
    ).mean() * 100

    single_high_support = (
        single["evidence_count"]
        >= 2
    ).mean() * 100

    # -----------------------------------------------------
    # 11. Source combinations
    # -----------------------------------------------------

    source_combination = (
        evidence
        .groupby("edge_id")["source"]
        .apply(
            lambda values:
                " + ".join(
                    sorted(
                        set(
                            str(v).strip()
                            for v in values
                            if str(v).strip()
                        )
                    )
                )
        )
        .rename("source_combination")
        .reset_index()
    )

    combination_counts = (
        source_combination[
            "source_combination"
        ]
        .value_counts()
        .rename_axis(
            "source_combination"
        )
        .reset_index(
            name="edge_count"
        )
    )

    # -----------------------------------------------------
    # 12. Save edge-level results
    # -----------------------------------------------------

    result_columns = [
        "edge_id",
        "source_id",
        "relation",
        "target_id",
        "evidence_count",
        "source_count",
        "evidence_class",
    ]

    analysis[
        result_columns
    ].to_csv(
        RESULT_FILE,
        index=False,
    )

    # -----------------------------------------------------
    # 13. Interpretation
    # -----------------------------------------------------
    # -----------------------------------------------------
    # Numerical reporting for extremely small p-values
    # -----------------------------------------------------

    if p_value == 0:
        p_value_report = "<1e-300"
    else:
        p_value_report = f"{p_value:.6g}"
    if p_value < 0.05:

        interpretation = (
            "Multi-source edges have statistically higher "
            "evidence-count distributions than single-source "
            "edges under the specified one-sided Mann-Whitney "
            "U test."
        )

    else:

        interpretation = (
            "The test does not provide sufficient evidence "
            "that multi-source edges have higher evidence-count "
            "distributions than single-source edges."
        )

    # -----------------------------------------------------
    # 14. Important methodological limitation
    # -----------------------------------------------------

    limitation = (
        "Evidence count is an operational measure of evidence "
        "support, not a direct ground-truth measure of biological "
        "reliability. The current final evidence table does not "
        "contain a source-independent validated confidence label "
        "that can establish true reliability. Therefore, H3 is "
        "interpreted as a test of greater evidence support for "
        "multi-source relations."
    )

    # -----------------------------------------------------
    # 15. Report
    # -----------------------------------------------------

    report = f"""# ReMedy H3 — Final KG

## Hypothesis

**H3:** Multi-source relations are more reliable than
single-source relations.

## Operational definition

Because the current evidence layer does not contain a
source-independent ground-truth reliability label, H3 is
operationalized as:

> Multi-source edges have greater evidence support than
> single-source edges.

Evidence support is measured primarily by the number of
evidence records associated with each canonical edge.

## Null hypothesis

Multi-source edges do not have higher evidence-count
distributions than single-source edges.

## Alternative hypothesis

Multi-source edges have higher evidence-count distributions
than single-source edges.

## Final KG

- Total edges with evidence: {len(analysis):,}
- Single-source edges: {len(single):,}
- Multi-source edges: {len(multi):,}

## Source combinations

{combination_counts.to_markdown(index=False)}

## Evidence-count distributions

| Metric | Multi-source | Single-source |
|---|---:|---:|
| Edge count | {n_multi:,} | {n_single:,} |
| Mean evidence count | {multi_mean:.4f} | {single_mean:.4f} |
| Median evidence count | {multi_median:.4f} | {single_median:.4f} |
| Q1 | {multi_q25:.4f} | {single_q25:.4f} |
| Q3 | {multi_q75:.4f} | {single_q75:.4f} |
| Edges with >=2 evidence records | {multi_high_support:.2f}% | {single_high_support:.2f}% |

## Statistical test

One-sided Mann-Whitney U test:

`multi-source evidence support > single-source evidence support`

## Results

- U statistic: {u_statistic:,.2f}
- p-value: {p_value_report}

### Effect size

Rank-biserial correlation:

`{rank_biserial:.4f}`

### Median difference

Multi-source median minus single-source median:

`{median_difference:.4f}`

### 95% bootstrap confidence interval

`[{ci_low:.4f}, {ci_high:.4f}]`

## Interpretation

{interpretation}

## Limitation

{limitation}

Therefore, a significant result should be described as
evidence of greater evidence support, not as direct proof
that multi-source relations are biologically more reliable.

## Reproducibility

- Random seed: {RANDOM_SEED}
- Bootstrap iterations: {BOOTSTRAP_ITERATIONS}

## Output

`reports/hypothesis_H3_final_results.csv`
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
    print("H3 FINAL-KG TEST COMPLETE")
    print("=" * 70)

    print(
        f"Single-source edges: "
        f"{n_single:,}"
    )

    print(
        f"Multi-source edges: "
        f"{n_multi:,}"
    )

    print(
        f"Single-source median evidence count: "
        f"{single_median:.4f}"
    )

    print(
        f"Multi-source median evidence count: "
        f"{multi_median:.4f}"
    )

    print(
        f"Median difference: "
        f"{median_difference:.4f}"
    )

    print(
        f"U statistic: "
        f"{u_statistic:,.2f}"
    )

    print(
        f"p-value: "
        f"{p_value_report}"
    )

    print(
        f"Rank-biserial effect size: "
        f"{rank_biserial:.4f}"
    )

    print(
        f"95% CI: "
        f"[{ci_low:.4f}, {ci_high:.4f}]"
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