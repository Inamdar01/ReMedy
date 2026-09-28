# H3 Multi-Source Evidence Check

## Purpose

H3 asks whether multi-source relations are more reliable than
single-source relations.

## Current final KG

- Total edges: 52,371
- Total evidence records: 84,024
- Single-source edges: 52,371
- Multi-source edges: 0
- Multi-source coverage: 0.0000%

## Source combinations

| source_combination   |   edge_count |
|:---------------------|-------------:|
| ChEMBL               |        30516 |
| Reactome             |        18088 |
| PrimeKG              |         3767 |

## Relation/source distribution

| relation         | evidence_class   |   edge_count |
|:-----------------|:-----------------|-------------:|
| associated_with  | single-source    |         2655 |
| contraindication | single-source    |          816 |
| indication       | single-source    |          249 |
| off-label use    | single-source    |           47 |
| participates_in  | single-source    |        18088 |
| targets          | single-source    |        30516 |

## Decision

No canonical edge currently has evidence from
two or more distinct source databases.

Therefore H3 cannot be statistically tested using the current
final evidence table without adding another independent source
for at least some existing canonical relations.

An H3 result should not be fabricated from single-source edges.
