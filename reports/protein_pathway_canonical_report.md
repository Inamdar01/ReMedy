# Canonical Protein → Pathway Report

## Input

`data/interim/chembl_reactome_protein_pathways.csv`

## Results

- Canonical UniProt → Reactome edges: 18,088
- Unique UniProt proteins: 1,045
- Unique Reactome pathways: 2,069

## Relation

Canonical relation:

`protein -- participates_in --> pathway`

Protein identifier:
`UniProt accession`

Pathway identifier:
`Reactome pathway ID`

Source:
`Reactome`

## Evidence

Reactome evidence codes are preserved and aggregated
per protein-pathway edge.

## Output

`data/interim/protein_pathway_relations_canonical.csv`
