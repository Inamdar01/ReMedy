# Disease → Canonical MONDO Integration Report

## Inputs

- Disease-UniProt layer:
  `data/interim/primekg_disease_protein_uniprot.csv`
- Selected MONDO mappings:
  `data/interim/selected_diseases_mondo_mapping.csv`

## Results

- Canonical diseases represented: 12
- Unique UniProt proteins: 1,556
- Canonical disease-UniProt rows: 2,655
- Unresolved disease mappings: 0

## Disease coverage

| mondo_id      | mondo_name                               |   unique_uniprot |   relation_rows |
|:--------------|:-----------------------------------------|-----------------:|----------------:|
| MONDO:0004975 | Alzheimer disease                        |               92 |              92 |
| MONDO:0005180 | Parkinson disease                        |              100 |             100 |
| MONDO:0005618 | anxiety disorder                         |              480 |             480 |
| MONDO:0007743 | attention deficit-hyperactivity disorder |               24 |              24 |
| MONDO:0004985 | bipolar disorder                         |              490 |             490 |
| MONDO:0001627 | dementia                                 |               18 |              18 |
| MONDO:0005027 | epilepsy                                 |              193 |             193 |
| MONDO:0002009 | major depressive disorder                |              280 |             280 |
| MONDO:0005277 | migraine disorder                        |               10 |              10 |
| MONDO:0005301 | multiple sclerosis                       |               45 |              45 |
| MONDO:0005090 | schizophrenia                            |              865 |             865 |
| MONDO:0005098 | stroke disorder                          |               58 |              58 |

## Identifier policy

PrimeKG disease identifiers are preserved as provenance fields.

The canonical disease identifier used by ReMedy is the MONDO ID.

No disease identifier is inferred from the disease name beyond
the explicit reviewed mapping table.

## Outputs

- `data/interim/disease_uniprot_relations_canonical.csv`
- `reports/unresolved_disease_mondo_mapping.csv`
