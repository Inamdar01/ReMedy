# PrimeKG Protein → HGNC → UniProt Mapping Report

## Input

- PrimeKG disease-protein relations:
  `data/interim/primekg_disease_protein_relations.csv`
- HGNC complete set:
  `data/raw/hgnc/hgnc_complete_set.txt`

## Mapping strategy

PrimeKG protein identifiers are treated as NCBI/Entrez Gene IDs.

The mapping uses exact Entrez ID matching against HGNC.

Only approved protein-coding HGNC records are retained.

All UniProt accessions supplied by HGNC are preserved.

No identifier is guessed from gene name alone.

## Results

- PrimeKG unique Entrez IDs: 1,599
- Entrez IDs mapped to UniProt: 1,543
- HGNC records represented: 1,543
- Unique UniProt accessions: 1,556
- Unresolved Entrez IDs: 56
- Entrez → UniProt coverage: 96.50%

## Ambiguity

- Entrez IDs with multiple HGNC mappings:
  0

## Harmonized disease-protein layer

- Original disease-protein rows:
  2,697
- Harmonized disease-UniProt rows:
  2,655

## Disease coverage

| disease_name                             |   unique_entrez |   unique_uniprot |   relation_rows |
|:-----------------------------------------|----------------:|-----------------:|----------------:|
| Alzheimer disease                        |              91 |               92 |              92 |
| Parkinson disease                        |             100 |              100 |             100 |
| anxiety disorder                         |             473 |              480 |             480 |
| attention deficit-hyperactivity disorder |              24 |               24 |              24 |
| bipolar disorder                         |             487 |              490 |             490 |
| dementia (disease)                       |              17 |               18 |              18 |
| epilepsy                                 |             193 |              193 |             193 |
| major depressive disorder                |             275 |              280 |             280 |
| migraine disorder                        |               9 |               10 |              10 |
| multiple sclerosis                       |              45 |               45 |              45 |
| schizophrenia                            |             857 |              865 |             865 |
| stroke disorder                          |              58 |               58 |              58 |

## Outputs

- `data/interim/primekg_protein_uniprot_mapping.csv`
- `data/interim/primekg_disease_protein_uniprot.csv`
- `reports/unresolved_primekg_proteins.csv`

## Important limitation

The original PrimeKG relation represents a combined
gene/protein concept. This mapping provides identifier
harmonization to HGNC and UniProt; it does not independently
change the biological meaning or strength of the PrimeKG
association.
