# PrimeKG Drug → Protein Canonical Report

## Source relation

PrimeKG relation:

`drug -- drug_protein --> gene/protein`

Only the DrugBank → NCBI orientation was retained.

## Scope

- Selected ChEMBL drugs: 568
- Relevant DrugBank mappings: 569

## PrimeKG extraction

- DrugBank → NCBI rows: 5,829
- DrugBank IDs represented: 520
- Entrez IDs represented: 803

## Canonical mapping

DrugBank → ChEMBL → Entrez → HGNC → UniProt

- Canonical Drug → Protein edges: 5,829
- Unique ChEMBL drugs: 520
- Unique UniProt proteins: 803
- DrugBank IDs represented in canonical layer: 520
- Entrez IDs represented in canonical layer: 803

## Relation

`drug -- targets --> protein`

## Provenance

PrimeKG source identifiers are retained:

- DrugBank ID
- Entrez ID
- HGNC ID
- PrimeKG source fields

These can later be attached to the canonical Drug → Protein
edge as PrimeKG evidence.

## Important policy

The relation is mapped to the generic `targets` edge.

No `inhibits` or `activates` relation is inferred from the
PrimeKG `drug_protein` relation.

## Outputs

`data/interim/primekg_drug_protein_relations_canonical.csv`

`reports/primekg_drug_protein_canonical_report.md`

Unresolved mappings are not silently invented or reassigned.
