# Canonical Drug → Disease Edge Report

## Inputs

- PrimeKG:
  `data/raw/primekg/primekg.csv`
- DrugBank → ChEMBL:
  `data/raw/chembl/drugbank_to_chembl_validated.csv`
- Disease → MONDO:
  `data/interim/selected_diseases_mondo_mapping.csv`

## Results

- PrimeKG selected rows before deduplication: 2,230
- Duplicate orientation rows removed: 1,115
    - Canonical duplicate edges collapsed: 1
- Canonical resolved edges: 1,112
- Unique ChEMBL drugs: 568
- Unique MONDO diseases: 12
- Unresolved edges: 2

## Relation counts

| relation         |   rows |
|:-----------------|-------:|
| contraindication |    816 |
| indication       |    249 |
| off-label use    |     47 |

## Disease coverage

| mondo_id      | mondo_name                               |   unique_drugs |   rows |
|:--------------|:-----------------------------------------|---------------:|-------:|
| MONDO:0004975 | Alzheimer disease                        |             64 |     64 |
| MONDO:0005180 | Parkinson disease                        |             72 |     72 |
| MONDO:0005618 | anxiety disorder                         |            235 |    247 |
| MONDO:0007743 | attention deficit-hyperactivity disorder |             19 |     19 |
| MONDO:0004985 | bipolar disorder                         |             57 |     69 |
| MONDO:0001627 | dementia                                 |             31 |     31 |
| MONDO:0005027 | epilepsy                                 |            332 |    334 |
| MONDO:0002009 | major depressive disorder                |             25 |     25 |
| MONDO:0005277 | migraine disorder                        |             63 |     63 |
| MONDO:0005301 | multiple sclerosis                       |             15 |     15 |
| MONDO:0005090 | schizophrenia                            |             52 |     52 |
| MONDO:0005098 | stroke disorder                          |            121 |    121 |

## Relation semantics

PrimeKG relations are preserved separately:

- `indication`
- `off-label use`
- `contraindication`

They are not collapsed into a single treatment relation.

## Identifier policy

- Drug node identifier: ChEMBL ID
- Disease node identifier: MONDO ID
- PrimeKG DrugBank and disease identifiers are retained as provenance.
- No unmapped identifier is guessed.

## Outputs

- `data/interim/drug_disease_relations_canonical.csv`
- `reports/unresolved_drug_disease_edges.csv`
