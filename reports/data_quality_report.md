# ReMedy Data Quality Report
This report was generated from the current canonical relation layers.
## 1. Overall graph statistics
- Total unique nodes: 5,113
- Total relation rows: 52,371
- Drugs: 568
- Diseases: 12
- Genes: 0 (no separate canonical gene relation layer yet)
- Proteins: 2,464
- Pathways: 2,069

## 2. Relation-table quality
| Table | Rows | Missing endpoints | Duplicate edges | Provenance coverage |
|---|---:|---:|---:|---:|
| drug_disease | 1,112 | 0 | 0 | 100.00% |
| drug_protein | 30,516 | 0 | 0 | 100.00% |
| disease_protein | 2,655 | 0 | 0 | 100.00% |
| protein_pathway | 18,088 | 0 | 0 | 100.00% |

## 3. Relation distribution
| Table | Relation | Edge count |
|---|---|---:|
| drug_disease | contraindication | 816 |
| drug_disease | indication | 249 |
| drug_disease | off-label use | 47 |
| protein_pathway | participates_in | 18,088 |

## 4. Missing-value percentage by column
| Table | Column | Missing | Missing % |
|---|---|---:|---:|
| drug_disease | chembl_id | 0 | 0.00% |
| drug_disease | mondo_id | 0 | 0.00% |
| drug_disease | relation | 0 | 0.00% |
| drug_disease | chembl_pref_name | 0 | 0.00% |
| drug_disease | drugbank_id | 0 | 0.00% |
| drug_disease | drug_name | 0 | 0.00% |
| drug_disease | mondo_name | 0 | 0.00% |
| drug_disease | disease_name | 0 | 0.00% |
| drug_disease | primekg_disease_id | 0 | 0.00% |
| drug_disease | disease_source | 0 | 0.00% |
| drug_disease | mapping_status | 0 | 0.00% |
| drug_disease | mapping_reason | 0 | 0.00% |
| drug_protein | chembl_id | 0 | 0.00% |
| drug_protein | drugbank_id | 0 | 0.00% |
| drug_protein | drug_name | 0 | 0.00% |
| drug_protein | target_id | 0 | 0.00% |
| drug_protein | target_name | 0 | 0.00% |
| drug_protein | uniprot_accession | 0 | 0.00% |
| drug_protein | target_organism | 0 | 0.00% |
| drug_protein | evidence_count | 0 | 0.00% |
| drug_protein | human_evidence_count | 0 | 0.00% |
| drug_protein | non_human_evidence_count | 0 | 0.00% |
| drug_protein | unspecified_evidence_count | 0 | 0.00% |
| disease_protein | mondo_id | 0 | 0.00% |
| disease_protein | mondo_name | 0 | 0.00% |
| disease_protein | selected_disease | 0 | 0.00% |
| disease_protein | disease_primekg_id | 0 | 0.00% |
| disease_protein | disease_name | 0 | 0.00% |
| disease_protein | disease_source | 0 | 0.00% |
| disease_protein | protein_entrez_id | 0 | 0.00% |
| disease_protein | protein_gene_symbol | 0 | 0.00% |
| disease_protein | protein_source | 0 | 0.00% |
| disease_protein | hgnc_id | 0 | 0.00% |
| disease_protein | hgnc_symbol | 0 | 0.00% |
| disease_protein | hgnc_name | 0 | 0.00% |
| disease_protein | uniprot_accession | 0 | 0.00% |
| protein_pathway | uniprot_accession | 0 | 0.00% |
| protein_pathway | relation | 0 | 0.00% |
| protein_pathway | reactome_pathway_id | 0 | 0.00% |
| protein_pathway | pathway_name | 0 | 0.00% |
| protein_pathway | evidence_count | 0 | 0.00% |
| protein_pathway | evidence_codes | 0 | 0.00% |
| protein_pathway | source | 0 | 0.00% |

## 5. Drug–disease label distribution
| Relation | Count | Percentage |
|---|---:|---:|
| contraindication | 816 | 73.38% |
| indication | 249 | 22.39% |
| off-label use | 47 | 4.23% |

## 6. High-degree nodes

### drug_disease

| Node | Degree |
|---|---:|
| MONDO:0005027 | 334 |
| MONDO:0005618 | 247 |
| MONDO:0005098 | 121 |
| MONDO:0005180 | 72 |
| MONDO:0004985 | 69 |
| MONDO:0004975 | 64 |
| MONDO:0005277 | 63 |
| MONDO:0005090 | 52 |
| MONDO:0001627 | 31 |
| MONDO:0002009 | 25 |

### drug_protein

| Node | Degree |
|---|---:|
| CHEMBL535 | 578 |
| CHEMBL1336 | 556 |
| CHEMBL553 | 552 |
| CHEMBL413 | 413 |
| P31645 | 405 |
| Q01959 | 403 |
| P35372 | 402 |
| P08908 | 401 |
| P21728 | 399 |
| P23975 | 399 |

### disease_protein

| Node | Degree |
|---|---:|
| MONDO:0005090 | 865 |
| MONDO:0004985 | 490 |
| MONDO:0005618 | 480 |
| MONDO:0002009 | 280 |
| MONDO:0005027 | 193 |
| MONDO:0005180 | 100 |
| MONDO:0004975 | 92 |
| MONDO:0005098 | 58 |
| MONDO:0005301 | 45 |
| MONDO:0007743 | 24 |

### protein_pathway

| Node | Degree |
|---|---:|
| R-HSA-162582 | 490 |
| R-HSA-1430728 | 283 |
| R-HSA-1643685 | 245 |
| R-HSA-168256 | 232 |
| R-HSA-372790 | 225 |
| R-HSA-388396 | 219 |
| P20618 | 182 |
| P28074 | 181 |
| P49721 | 181 |
| P28482 | 180 |

## 7. Identifier mapping
- Validated DrugBank → ChEMBL mappings: 569
- PrimeKG disease-protein → UniProt coverage: 96.50% (from the completed mapping stage)
- ChEMBL target → UniProt unresolved mappings: 0 (from the completed mapping stage)
- Reactome unmatched human ChEMBL targets: 165

## 8. Notes
- Raw source files were not modified.
- Canonical identifiers are used in the current relation layers.
- Drug–disease relations remain separated as indication, off-label use, and contraindication.
- Non-human ChEMBL evidence is preserved separately rather than being treated as human evidence.
- A separate canonical gene relation layer has not yet been created, so genes are reported as 0 rather than inferred from protein identifiers.
