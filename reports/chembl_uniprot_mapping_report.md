# ChEMBL Target → UniProt Mapping Report

## Input

`data/interim/chembl_evidence_clean.csv`

## Results

- ChEMBL targets in input: 1842
- ChEMBL targets represented in output: 1842
- Unique UniProt accessions: 1841
- Human ChEMBL targets: 1208
- Targets without UniProt accession: 0

## Mapping method

ChEMBL target records were requested in batches using the
`target_chembl_id__in` filter.

UniProt accessions were extracted from ChEMBL
`target_components`.

No protein identifier was inferred from target names.

## Reactome eligibility

Human ChEMBL targets with a UniProt accession can be joined
against the Reactome UniProt mapping using the accession.

Non-human targets are retained in this mapping but should not
be used for the human Reactome pathway layer.

## Files

- `chembl_target_uniprot_mapping.csv`
- `chembl_uniprot_mapping_progress.csv`
- `chembl_unresolved_uniprot_targets.csv` (if needed)
