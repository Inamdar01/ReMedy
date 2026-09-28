# ReMedy H2 — Final KG

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

- Total edges: 56,208
- Unique drugs: 568
- Unique diseases: 12
- Known Drug-Disease pairs: 1,086

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
12-disease set.

This was repeated 5,000 times.

## Observed result

- Supported known pairs: 1,050
- Total known pairs: 1,086
- Observed support rate: 96.69%

## Permutation result

- Null mean: 96.42%
- Null standard deviation: 0.001469
- Observed minus null: 0.27 percentage points
- Permutation p-value: 0.0653869
- Null statistics >= observed: 326/5,000
- Permutation exceedance rate: 6.52%

## Effect size

Primary effect size:

Observed support rate minus randomized-null mean:

`0.002676`

In percentage-point terms:

`0.27`

Standardized permutation statistic:

`1.8218`

## 95% bootstrap confidence interval

A drug-level clustered bootstrap was used because multiple
Drug-Disease observations share the same drug.

95% bootstrap CI for the observed-minus-randomized
support-rate difference:

`[-0.001845, 0.008115]`

Percentage-point form:

`[-0.18, 0.81]`

## 95% permutation null interval

For descriptive comparison, the central 95% of the randomized
support-rate distribution was:

`[0.961326, 0.966851]`

Percentage-point form:

`[96.13, 96.69]`

This is a permutation null interval, not a confidence interval.

## Interpretation

The degree-preserving permutation test does not provide sufficient evidence that the observed path-support rate exceeds the randomized null distribution.

The 0.01/0.05 significance threshold is not interpreted as
evidence of clinical efficacy. H2 evaluates biological
connectivity in the integrated knowledge graph only.

## Reproducibility

- Random seed: 42
- Permutations: 5000
- Bootstrap iterations: 5000
- Bootstrap unit: drug

## Output

`reports/hypothesis_H2_final_results.csv`
