# Reactome Pathway Hierarchy Integration Report

## Input

- Protein-pathway mapping: `data/interim/chembl_reactome_protein_pathways.csv`
- Reactome pathway definitions: `data/raw/reactome/ReactomePathways.txt`
- Reactome hierarchy: `data/raw/reactome/ReactomePathwaysRelation.txt`

## Results

- Unique pathways connected to ReMedy proteins: 2,069
- Pathways found in Reactome definitions: 2,069
- Total Reactome hierarchy relations: 23,717
- Relevant hierarchy relations: 2,684
- Unique parent pathways: 764
- Unique child pathways: 2,641

## Output

`data/interim/reactome_pathway_hierarchy.csv`

The hierarchy is retained as explicit parent-to-child pathway relationships.
No disease relationships are inferred from the pathway hierarchy.
