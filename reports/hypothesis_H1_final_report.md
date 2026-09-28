# ReMedy H1 — Final KG

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

- Nodes: 5,113
- Edges: 52,371
- Unique drugs: 568
- Unique diseases: 12
- Known Drug-Disease pairs: 1,086
- Non-known pairs: 5,730

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
| Number of pairs | 1,086 | 5,730 |
| Median connectivity | 126.0000 | 58.0000 |
| Non-zero connectivity | 85.08% | 83.44% |

### Statistical result

- U statistic: 3,660,544.50
- p-value: 1.0546e-20
- Rank-biserial effect size: 0.1765
- Median difference: 68.0000
- 95% bootstrap CI: [59.5000, 74.0000]

## Interpretation

The known Drug-Disease pairs have higher biological connectivity scores than the non-known pairs under the specified one-sided Mann-Whitney U test.

This result concerns biological structure in the integrated
knowledge graph. It does not establish clinical efficacy,
safety, or treatment effectiveness.

## Reproducibility

- Random seed: 42
- Bootstrap iterations: 5000

## Output

`reports/hypothesis_H1_final_pair_scores.csv`
