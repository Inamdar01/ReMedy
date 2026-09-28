# ReMedy Final Knowledge Graph Build Report

## Final files

- `data/processed/nodes.parquet`
- `data/processed/edges.parquet`
- `data/processed/evidence.parquet`

## Nodes

- Total nodes: 5,401
- Drugs: 568
- Diseases: 12
- Proteins: 2,752
- Pathways: 2,069
- Genes: 0

## Edges

- Total canonical edges: 56,208

### Relation distribution

| relation         |   edge_count |
|:-----------------|-------------:|
| associated_with  |         2655 |
| contraindication |          816 |
| indication       |          249 |
| off-label use    |           47 |
| participates_in  |        18088 |
| targets          |        34353 |

## Evidence

- Total evidence rows: 89,853
- Edges with >=1 evidence record: 56,208
- Edges without evidence: 0

### Evidence source distribution

| source   |   evidence_rows |
|:---------|----------------:|
| ChEMBL   |           62169 |
| PrimeKG  |            9596 |
| Reactome |           18088 |

## Multi-source coverage

- Single-source edges: 54,216
- Multi-source edges: 1,992

### Drug → Protein specifically

- Drug → Protein edges: 34,353
- Single-source Drug → Protein edges: 32,361
- Multi-source Drug → Protein edges: 1,992
- Multi-source Drug → Protein coverage: 5.80%

## Canonical identifiers

- Drug: ChEMBL ID
- Disease: MONDO ID
- Protein: UniProt accession
- Pathway: Reactome ID

## Relation semantics

- Drug → Disease:
  - indication
  - off-label use
  - contraindication
- Drug → Protein:
  - targets
- Disease → Protein:
  - associated_with
- Protein → Pathway:
  - participates_in

Mechanistic `inhibits` or `activates` relations are not inferred
from activity values.

## Multi-source Drug → Protein policy

ChEMBL and PrimeKG records that refer to the same canonical
ChEMBL → UniProt edge are represented as one edge with
multiple evidence sources.

Different source records are retained in `evidence.parquet`.

## Evidence policy

- ChEMBL assay-level evidence is preserved.
- PrimeKG provenance is preserved.
- Reactome provenance and evidence codes are preserved.
- Assay organism is preserved as recorded.
- Missing context is represented as `unknown`.
- Non-human assay evidence is not relabeled as human evidence.

## Scope note

No separate gene relation layer is constructed because the
current canonical disease-protein layer combines gene/protein
information. Gene edges are not inferred.

## Integrity

Edge IDs are deterministic hashes of:

`source_id | relation | target_id`

The evidence table references those edge IDs.
