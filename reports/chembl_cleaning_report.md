# ChEMBL Evidence Cleaning Report

## Input

Source:
`data/interim/chembl_target_evidence_combined.csv`

## Final cleaned evidence

- Evidence records: 73594
- Unique drugs with evidence: 480
- Unique protein targets: 1847
- Exact duplicate activity rows removed: 0

## Experimental organism classification

Classification uses `assay_organism`.

- Human: 3623
- Non-human: 1596
- Organism unspecified: 68375

`target_organism` is retained separately and is not used to classify the experimental assay.

## Value quality

- Usable: 71396
- Flagged: 2198

## Activity types

- AC50: 29107
- IC50: 13289
- Ki: 10696
- Potency: 6184
- Kd: 5524
- Inhibition: 5364
- Activity: 1784
- EC50: 1627
- Relative IC50: 10
- pKi: 6
- IC20: 3

## Important interpretation rule

ChEMBL activity values are preserved with their original
activity type, value, units, and relation. Different activity
types are not collapsed into one numerical score.

Flagged values are retained and identified through the
`value_quality` field.
