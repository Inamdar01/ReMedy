from pathlib import Path
import itertools
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[2]

NODES_FILE = ROOT / "data/processed/nodes.parquet"
EDGES_FILE = ROOT / "data/processed/edges.parquet"

REPORT_FILE = ROOT / "reports/hypothesis_H1_final_report.md"
RESULT_FILE = ROOT / "reports/hypothesis_H1_final_pair_scores.csv"

RANDOM_SEED = 42
BOOTSTRAP_ITERATIONS = 5000


def bootstrap_median_difference(
    known_scores,
    random_scores,
    iterations=5000,
    seed=42,
):
    rng = np.random.default_rng(seed)

    known_array = np.asarray(
        known_scores,
        dtype=float,
    )

    random_array = np.asarray(
        random_scores,
        dtype=float,
    )

    differences = np.empty(
        iterations,
        dtype=float,
    )

    for i in range(iterations):

        known_sample = rng.choice(
            known_array,
            size=len(known_array),
            replace=True,
        )

        random_sample = rng.choice(
            random_array,
            size=len(random_array),
            replace=True,
        )

        differences[i] = (
            np.median(known_sample)
            - np.median(random_sample)
        )

    return (
        float(np.percentile(differences, 2.5)),
        float(np.percentile(differences, 97.5)),
    )


def main():

    print("=" * 70)
    print("REMEDY - H1 TEST ON FINAL PARQUET KG")
    print("=" * 70)

    # ---------------------------------------------------------
    # 1. Load final KG
    # ---------------------------------------------------------

    print("\nLoading final KG...")

    nodes = pd.read_parquet(
        NODES_FILE
    )

    edges = pd.read_parquet(
        EDGES_FILE
    )

    print(
        f"Nodes: {len(nodes):,}"
    )

    print(
        f"Edges: {len(edges):,}"
    )

    # ---------------------------------------------------------
    # 2. Extract relation layers from final edges
    # ---------------------------------------------------------

    drug_protein_edges = edges[
        edges["relation"] == "targets"
    ].copy()

    disease_protein_edges = edges[
        edges["relation"] == "associated_with"
    ].copy()

    protein_pathway_edges = edges[
        edges["relation"] == "participates_in"
    ].copy()

    drug_disease_edges = edges[
        edges["relation"].isin(
            [
                "indication",
                "off-label use",
                "contraindication",
            ]
        )
    ].copy()

    # ---------------------------------------------------------
    # 3. Build biological sets
    # ---------------------------------------------------------

    drug_proteins = (
        drug_protein_edges
        .groupby("source_id")["target_id"]
        .apply(set)
        .to_dict()
    )

    disease_proteins = (
        disease_protein_edges
        .groupby("source_id")["target_id"]
        .apply(set)
        .to_dict()
    )

    protein_pathways = (
        protein_pathway_edges
        .groupby("source_id")["target_id"]
        .apply(set)
        .to_dict()
    )

    # ---------------------------------------------------------
    # 4. Derive Drug → Pathway
    # ---------------------------------------------------------

    drug_pathways = {}

    for drug, proteins in drug_proteins.items():

        pathways = set()

        for protein in proteins:

            pathways.update(
                protein_pathways.get(
                    protein,
                    set(),
                )
            )

        drug_pathways[drug] = pathways

    # ---------------------------------------------------------
    # 5. Derive Disease → Pathway
    # ---------------------------------------------------------

    disease_pathways = {}

    for disease, proteins in disease_proteins.items():

        pathways = set()

        for protein in proteins:

            pathways.update(
                protein_pathways.get(
                    protein,
                    set(),
                )
            )

        disease_pathways[disease] = pathways

    # ---------------------------------------------------------
    # 6. Biological connectivity
    # ---------------------------------------------------------

    def connectivity(
        drug,
        disease,
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

        return {
            "shared_proteins": len(
                shared_proteins
            ),
            "shared_pathways": len(
                shared_pathways
            ),
            "connectivity_score": (
                len(shared_proteins)
                +
                len(shared_pathways)
            ),
        }

    # ---------------------------------------------------------
    # 7. Known pairs
    # ---------------------------------------------------------

    known_pairs = set(
        zip(
            drug_disease_edges["source_id"],
            drug_disease_edges["target_id"],
        )
    )

    drugs = sorted(
        drug_disease_edges["source_id"].unique()
    )

    diseases = sorted(
        drug_disease_edges["target_id"].unique()
    )

    all_pairs = set(
        itertools.product(
            drugs,
            diseases,
        )
    )

    random_pairs = sorted(
        all_pairs - known_pairs
    )

    print(
        f"\nUnique drugs: {len(drugs):,}"
    )

    print(
        f"Unique diseases: {len(diseases):,}"
    )

    print(
        f"Known Drug-Disease pairs: "
        f"{len(known_pairs):,}"
    )

    print(
        f"Non-known pairs: "
        f"{len(random_pairs):,}"
    )

    # ---------------------------------------------------------
    # 8. Score known pairs
    # ---------------------------------------------------------

    rows = []

    print(
        "\nCalculating known-pair connectivity..."
    )

    for drug, disease in sorted(
        known_pairs
    ):

        score = connectivity(
            drug,
            disease,
        )

        rows.append(
            {
                "pair_type": "known",
                "source_id": drug,
                "target_id": disease,
                **score,
            }
        )

    # ---------------------------------------------------------
    # 9. Score non-known pairs
    # ---------------------------------------------------------

    print(
        "Calculating non-known-pair connectivity..."
    )

    for drug, disease in random_pairs:

        score = connectivity(
            drug,
            disease,
        )

        rows.append(
            {
                "pair_type": "random",
                "source_id": drug,
                "target_id": disease,
                **score,
            }
        )

    scores = pd.DataFrame(rows)

    # ---------------------------------------------------------
    # 10. Statistical test
    # ---------------------------------------------------------

    known_scores = scores.loc[
        scores["pair_type"] == "known",
        "connectivity_score",
    ].astype(float)

    random_scores = scores.loc[
        scores["pair_type"] == "random",
        "connectivity_score",
    ].astype(float)

    u_statistic, p_value = mannwhitneyu(
        known_scores,
        random_scores,
        alternative="greater",
        method="asymptotic",
    )

    n_known = len(known_scores)
    n_random = len(random_scores)

    rank_biserial = (
        (2 * u_statistic)
        /
        (n_known * n_random)
    ) - 1

    known_median = float(
        np.median(known_scores)
    )

    random_median = float(
        np.median(random_scores)
    )

    median_difference = (
        known_median
        - random_median
    )

    ci_low, ci_high = (
        bootstrap_median_difference(
            known_scores.values,
            random_scores.values,
            iterations=BOOTSTRAP_ITERATIONS,
            seed=RANDOM_SEED,
        )
    )

    known_nonzero_pct = (
        known_scores.gt(0).mean()
        * 100
    )

    random_nonzero_pct = (
        random_scores.gt(0).mean()
        * 100
    )

    # ---------------------------------------------------------
    # 11. Save pair scores
    # ---------------------------------------------------------

    RESULT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    scores.to_csv(
        RESULT_FILE,
        index=False,
    )

    # ---------------------------------------------------------
    # 12. Interpretation
    # ---------------------------------------------------------

    if p_value < 0.05:

        interpretation = (
            "The known Drug-Disease pairs have higher "
            "biological connectivity scores than the "
            "non-known pairs under the specified "
            "one-sided Mann-Whitney U test."
        )

    else:

        interpretation = (
            "The test does not provide sufficient evidence "
            "that known Drug-Disease pairs have higher "
            "biological connectivity scores than the "
            "non-known pairs."
        )

    # ---------------------------------------------------------
    # 13. Report
    # ---------------------------------------------------------

    report = f"""# ReMedy H1 — Final KG

## Hypothesis

Known drug-disease pairs have stronger biological
connectivity than random pairs.

## Null hypothesis

Known Drug-Disease pairs do not have higher biological
connectivity than non-known Drug-Disease pairs.

## Alternative hypothesis

Known Drug-Disease pairs have higher biological
connectivity than non-known Drug-Disease pairs.

## Final KG used

- Nodes: {len(nodes):,}
- Edges: {len(edges):,}
- Unique drugs: {len(drugs):,}
- Unique diseases: {len(diseases):,}
- Known Drug-Disease pairs: {n_known:,}
- Non-known pairs: {n_random:,}

## Connectivity definition

For each Drug-Disease pair:

`connectivity_score = shared proteins + shared Reactome pathways`

Shared proteins are proteins connected to both the drug
through `targets` and the disease through `associated_with`.

Shared pathways are Reactome pathways reached from those
drug-target and disease-associated proteins.

## Test

One-sided Mann-Whitney U test:

`known connectivity > non-known connectivity`

## Results

| Metric | Known | Non-known |
|---|---:|---:|
| Number of pairs | {n_known:,} | {n_random:,} |
| Median connectivity | {known_median:.4f} | {random_median:.4f} |
| Non-zero connectivity | {known_nonzero_pct:.2f}% | {random_nonzero_pct:.2f}% |

### Statistical result

- U statistic: {u_statistic:,.2f}
- p-value: {p_value:.6g}
- Rank-biserial effect size: {rank_biserial:.4f}
- Median difference: {median_difference:.4f}
- 95% bootstrap CI: [{ci_low:.4f}, {ci_high:.4f}]

## Interpretation

{interpretation}

This result concerns biological structure in the integrated
knowledge graph. It does not establish clinical efficacy,
safety, or treatment effectiveness.

## Reproducibility

- Random seed: {RANDOM_SEED}
- Bootstrap iterations: {BOOTSTRAP_ITERATIONS}

## Output

`reports/hypothesis_H1_final_pair_scores.csv`
"""

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    # ---------------------------------------------------------
    # 14. Terminal output
    # ---------------------------------------------------------

    print()
    print("=" * 70)
    print("H1 FINAL-KG TEST COMPLETE")
    print("=" * 70)

    print(
        f"Known pairs: {n_known:,}"
    )

    print(
        f"Non-known pairs: {n_random:,}"
    )

    print(
        f"Known median connectivity: "
        f"{known_median:.4f}"
    )

    print(
        f"Non-known median connectivity: "
        f"{random_median:.4f}"
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
        f"{p_value:.6g}"
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