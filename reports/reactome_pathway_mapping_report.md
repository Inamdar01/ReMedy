# Reactome Protein → Pathway Mapping Report

## Input

- ChEMBL → UniProt mapping:
  `data/interim/chembl_target_uniprot_mapping_extended.csv`
- Reactome bulk mapping:
  `data/raw/reactome/UniProt2Reactome_All_Levels.txt`

## Results

- ChEMBL target mappings: 1869
- Human target mappings: 1223
- Unique human UniProt accessions:
  1223

- Human target-pathway mapping rows:
  18088

- Unique ChEMBL targets matched:
  1045

- Unique UniProt proteins matched:
  1045

- Unique Reactome pathways:
  2069

- Human ChEMBL targets without a Reactome pathway:
  165

## Filtering rule

Only Reactome entries with:

`species = Homo sapiens`

were used for the primary ReMedy human pathway layer.

## Output

`data/interim/chembl_reactome_protein_pathways.csv`

Unmatched targets:

`reports/reactome_unmatched_human_targets.csv`
