# ReMedy H3 — Final KG

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

- Total edges with evidence: 56,208
- Single-source edges: 54,216
- Multi-source edges: 1,992

## Source combinations

| source_combination   |   edge_count |
|:---------------------|-------------:|
| ChEMBL               |        28524 |
| Reactome             |        18088 |
| PrimeKG              |         7604 |
| ChEMBL + PrimeKG     |         1992 |

## Evidence-count distributions

| Metric | Multi-source | Single-source |
|---|---:|---:|
| Edge count | 1,992 | 54,216 |
| Mean evidence count | 9.6295 | 1.3035 |
| Median evidence count | 3.0000 | 1.0000 |
| Q1 | 2.0000 | 1.0000 |
| Q3 | 6.0000 | 1.0000 |
| Edges with >=2 evidence records | 100.00% | 14.73% |

## Statistical test

One-sided Mann-Whitney U test:

`multi-source evidence support > single-source evidence support`

## Results

- U statistic: 102,880,671.00
- p-value: <1e-300

### Effect size

Rank-biserial correlation:

`0.9052`

### Median difference

Multi-source median minus single-source median:

`2.0000`

### 95% bootstrap confidence interval

`[2.0000, 3.0000]`

## Interpretation

Multi-source edges have statistically higher evidence-count distributions than single-source edges under the specified one-sided Mann-Whitney U test.

## Limitation

Evidence count is an operational measure of evidence support, not a direct ground-truth measure of biological reliability. The current final evidence table does not contain a source-independent validated confidence label that can establish true reliability. Therefore, H3 is interpreted as a test of greater evidence support for multi-source relations.

Therefore, a significant result should be described as
evidence of greater evidence support, not as direct proof
that multi-source relations are biologically more reliable.

## Reproducibility

- Random seed: 42
- Bootstrap iterations: 5000

## Output

`reports/hypothesis_H3_final_results.csv`
