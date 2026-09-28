# ReMedy Final KG Validation Report

Validation performed directly against the final Parquet KG files.

## 1. Files

- `data/processed/nodes.parquet`
- `data/processed/edges.parquet`
- `data/processed/evidence.parquet`

## 2. Counts

- Nodes: 5,401
- Edges: 56,208
- Evidence records: 89,853
- Edges with evidence: 56,208
- Evidence coverage: 100.00%

## 3. Node types

| Node type | Count |
|---|---:|
| disease | 12 |
| drug | 568 |
| pathway | 2,069 |
| protein | 2,752 |

## 4. Relations

| Relation | Count |
|---|---:|
| associated_with | 2,655 |
| contraindication | 816 |
| indication | 249 |
| off-label use | 47 |
| participates_in | 18,088 |
| targets | 34,353 |

## 5. Evidence sources

| Source | Records |
|---|---:|
| ChEMBL | 62,169 |
| PrimeKG | 9,596 |
| Reactome | 18,088 |

## 6. Integrity checks

| Check | Count |
|---|---:|
| Duplicate node IDs | 0 |
| Duplicate edge IDs | 0 |
| Duplicate canonical edges | 0 |
| Orphan edge source IDs | 0 |
| Orphan edge target IDs | 0 |
| Orphan evidence edge IDs | 0 |
| Edges without evidence | 0 |
| Evidence-count mismatches | 0 |
| Invalid relation endpoint types | 0 |
| Duplicate evidence rows | 0 |

## 7. Missing-value counts

### Nodes

| Column | Missing/blank |
|---|---:|
| node_id | 0 |
| node_type | 0 |
| label | 0 |
| source_database | 0 |

### Edges

| Column | Missing/blank |
|---|---:|
| edge_id | 0 |
| source_id | 0 |
| relation | 0 |
| target_id | 0 |
| direction | 0 |
| evidence_count | 0 |

### Evidence

| Column | Missing/blank |
|---|---:|
| edge_id | 0 |
| source | 0 |
| publication | 0 |
| assay | 0 |
| organism | 0 |
| cell_line | 0 |
| activity_value | 0 |
| context | 0 |

## 8. Validation status

**PASSED**

All structural integrity checks passed.
