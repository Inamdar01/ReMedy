# ReMedy H4 — Final KG

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

- Total ChEMBL evidence records: 73,594
- Context-matched records: 3,623
- Context-mismatched records: 1,596
- Unspecified records excluded: 68,375

## Informativeness results

| Metric | Context-matched | Context-mismatched |
|---|---:|---:|
| Records | 3,623 | 1,596 |
| Mean score | 1.8625 | 1.8546 |
| Median score | 2.0000 | 2.0000 |
| Q1 | 1.0000 | 2.0000 |
| Q3 | 2.0000 | 2.0000 |

## Component availability

| component            |   matched_percentage |   mismatched_percentage |   difference_percentage_points |
|:---------------------|---------------------:|------------------------:|-------------------------------:|
| Quantitative pChEMBL |              65.3602 |                76.1278  |                      -10.7676  |
| Publication          |             100      |                99.812   |                        0.18797 |
| Cell line            |              20.8943 |                 9.52381 |                       11.3705  |

## Statistical test

Permutation test comparing:

`mean(context-matched score) - mean(context-mismatched score)`

The context labels were randomly permuted while preserving
the observed group sizes.

Number of permutations:

`5,000`

## Results

- Observed mean difference: 0.007908
- Observed difference in score units: 0.0079
- Null mean difference: 0.000152
- Null standard deviation: 0.018222
- p-value: 0.342132
- Permutations >= observed: 1,710/5,000

## Effect size

Primary effect size:

`0.007908`

Standardized permutation statistic:

`0.4257`

## 95% bootstrap confidence interval

`[-0.023155, 0.039552]`

## Interpretation

The permutation test does not provide sufficient evidence that context-matched human assay records have a higher mean informativeness score than context-mismatched non-human assay records.

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

- Random seed: 42
- Permutations: 5,000
- Bootstrap iterations: 5,000

## Output

`reports/hypothesis_H4_final_results.csv`
