# Canonical Drug → Protein Report

## Inputs

- ChEMBL drug-target summary:
  `data/interim/chembl_drug_protein_relations.csv`
- Extended target-UniProt mapping:
  `data/interim/chembl_target_uniprot_mapping_extended.csv`

## Core human graph

Only ChEMBL targets explicitly classified as
`Homo sapiens` are included in the canonical core
Drug → UniProt layer.

- Canonical rows: 30,516
- Unique ChEMBL drugs: 478
- Unique ChEMBL targets: 1,210
- Unique UniProt accessions: 1,223
- Human assay evidence represented: 3,685
- Total evidence represented on these target edges: 62,169

## Non-human evidence

Non-human target relationships are preserved separately
rather than used to create human biological links.

- Non-human relation rows: 4,559

## Unresolved target mappings

- Unresolved drug-target rows: 0

## Important interpretation rule

The current source table supports a generic ChEMBL
`Drug → target` relation.

We do NOT infer `inhibits` or `activates` solely from the
activity type because doing so would require an explicit
direction/interaction interpretation.

Activity values remain in the ChEMBL evidence layer.

## Outputs

- `data/interim/drug_protein_relations_canonical.csv`
- `reports/nonhuman_drug_protein_relations.csv`
- `reports/unresolved_drug_protein_edges.csv`
