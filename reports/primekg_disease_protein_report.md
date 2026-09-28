# PrimeKG Disease-Protein Extraction Report

## Input

- Source: `data/raw/primekg/primekg.csv`
- Relation: `disease_protein`

## Selected diseases

- epilepsy
- anxiety disorder
- stroke disorder
- Parkinson disease
- Alzheimer disease
- migraine disorder
- bipolar disorder
- schizophrenia
- dementia (disease)
- major depressive disorder
- attention deficit-hyperactivity disorder
- multiple sclerosis

## Results

- Disease-protein rows before deduplication: 5,394
- Exact duplicates removed: 2,697
- Final disease-protein rows: 2,697
- Selected diseases represented: 12
- Unique protein identifiers: 1,599

## Coverage

| disease_name                             | disease_source   |   protein_count |   relation_rows |
|:-----------------------------------------|:-----------------|----------------:|----------------:|
| Alzheimer disease                        | MONDO_grouped    |             101 |             101 |
| Parkinson disease                        | MONDO_grouped    |             104 |             104 |
| anxiety disorder                         | MONDO_grouped    |             481 |             481 |
| attention deficit-hyperactivity disorder | MONDO            |              24 |              24 |
| bipolar disorder                         | MONDO_grouped    |             497 |             497 |
| dementia (disease)                       | MONDO            |              17 |              17 |
| epilepsy                                 | MONDO            |             193 |             193 |
| major depressive disorder                | MONDO_grouped    |             278 |             278 |
| migraine disorder                        | MONDO            |               9 |               9 |
| multiple sclerosis                       | MONDO            |              45 |              45 |
| schizophrenia                            | MONDO_grouped    |             890 |             890 |
| stroke disorder                          | MONDO            |              58 |              58 |

## Identifier note

PrimeKG protein identifiers are retained as supplied by PrimeKG.
Disease identifiers are also retained as supplied and are not
treated as canonical MONDO identifiers until an explicit mapping
step is completed.

## Output

`data/interim/primekg_disease_protein_relations.csv`