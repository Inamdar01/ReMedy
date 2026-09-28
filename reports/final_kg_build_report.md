# ReMedy Final Knowledge Graph Build Report

## Final files

- `data/processed/nodes.parquet`
- `data/processed/edges.parquet`
- `data/processed/evidence.parquet`

## Nodes

- Total nodes: 5,113
- Drugs: 568
- Diseases: 12
- Proteins: 2,464
- Pathways: 2,069
- Genes: 0

## Edges

- Total canonical edges: 52,371

### Relation distribution

| relation         |   edge_count |
|:-----------------|-------------:|
| associated_with  |         2655 |
| contraindication |          816 |
| indication       |          249 |
| off-label use    |           47 |
| participates_in  |        18088 |
| targets          |        30516 |

## Evidence

- Total evidence rows: 84,024
- Edges with >=1 evidence row: 52,371
- Edges with 0 evidence rows: 0

### Evidence source distribution

| source   |   evidence_rows |
|:---------|----------------:|
| ChEMBL   |           62169 |
| PrimeKG  |            3767 |
| Reactome |           18088 |

## Canonical identifier policy

- Drug: ChEMBL ID
- Disease: MONDO ID
- Protein: UniProt accession
- Pathway: Reactome ID

## Relation semantics

- Drug → Disease relations remain separate:
  - indication
  - off-label use
  - contraindication
- Drug → Protein uses the generic `targets` relation.
- Disease → Protein uses the generic `associated_with` relation.
- Protein → Pathway uses `participates_in`.

Mechanistic relations such as `inhibits` or `activates`
were not inferred from activity types.

## Evidence policy

- ChEMBL assay-level evidence is linked to canonical
  human ChEMBL-target → UniProt edges.
- Assay organism is preserved as recorded.
- Missing context is represented as `unknown`.
- PrimeKG provenance is retained in the evidence context.
- Reactome evidence codes are retained in the evidence context.
- Non-human assay evidence is not relabeled as human evidence.

## Important scope note

No separate gene relation layer is constructed because the
current canonical disease-protein source combines gene/protein
information. Gene edges are not inferred from protein records.

## Integrity

The final edge table uses stable deterministic edge IDs
generated from:

`source_id | relation | target_id`

The final evidence table references these edge IDs.
