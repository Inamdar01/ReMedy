# ReMedy-KG v1.0 Freeze Manifest

## Freeze Status

Version: ReMedy-KG v1.0
Status: FROZEN
Freeze date: 2026-09-29

## Final Knowledge Graph

| File | SHA-256 |
|---|---|
| data/processed/nodes.parquet | 6AB0F858AC516F678E2A038518307FD4622AED3AE3D46C15E3FFBAB55C11B549 |
| data/processed/edges.parquet | 94AB5A3A52F03DA14CB5DCF0A28359432D1D223BDD86AF3AFD02636FE6D9D530 |
| data/processed/evidence.parquet | 845E927F7C0E8F0C3918382EDECED8085299774764A9682E42F8E6C6B228C6D9 |

## Final KG Statistics

- Nodes: 5,401
- Edges: 56,208
- Evidence records: 89,853
- Evidence coverage: 100.00%
- Duplicate node IDs: 0
- Duplicate edge IDs: 0
- Duplicate canonical edges: 0
- Orphan edge source IDs: 0
- Orphan edge target IDs: 0
- Orphan evidence edge IDs: 0
- Edges without evidence: 0
- Evidence-count mismatches: 0
- Invalid relation endpoint types: 0
- Duplicate evidence rows: 0

## Validation

Final KG validation: PASSED

H1: finalized
H2: finalized
H3: finalized
H4: finalized
H5: finalized

## H5 Reference

DrugMechDB SHA-256:
B973118FF3A43450AAC9FACC9957EC5CE1310EFA049ED11B6E94B42807050EE7

## Freeze Rule

The files in data/processed/ constitute ReMedy-KG v1.0.

These files must not be modified during model evaluation.

Any future changes to preprocessing, source data, mappings, filters, thresholds, or graph construction must produce a new graph version rather than modifying ReMedy-KG v1.0.

## Purpose

ReMedy-KG v1.0 is the frozen research dataset for downstream model evaluation, including graph-based drug-disease association prediction and explainability experiments.

This dataset is intended for research and hypothesis generation, not clinical recommendation.
