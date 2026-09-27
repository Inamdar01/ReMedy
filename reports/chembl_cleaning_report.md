# ChEMBL Evidence Cleaning Report

## Input

Source:
`data/raw/chembl/evidence/chembl_target_evidence_final.csv`

## Final cleaned evidence

- Evidence records: 72717
- Unique drugs with evidence: 475
- Unique protein targets: 1842
- Exact duplicate activity rows removed: 0

## Experimental organism classification

Classification uses `assay_organism`.

- Human: 3229
- Non-human: 1543
- Organism unspecified: 67945

`target_organism` is retained separately and is not used to classify the experimental assay.

## Value quality

- Usable: 70544
- Flagged: 2173

## Activity types

- AC50: 28791
- IC50: 13105
- Ki: 10641
- Potency: 6123
- Kd: 5498
- Inhibition: 5313
- Activity: 1721
- EC50: 1506
- Relative IC50: 10
- pKi: 6
- IC20: 3

## Important interpretation rule

ChEMBL activity values are preserved with their original
activity type, value, units, and relation. Different activity
types are not collapsed into one numerical score.

Flagged values are retained and identified through the
`value_quality` field.
