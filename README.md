# ReMedy

Explainable drug-repurposing research on a frozen biomedical knowledge graph. The interface takes a disease name and draws the supported biological path for a few recorded treatments and ranked drugs.

ReMedy-KG v1.0 in `data/processed/` stays frozen. The ranker reads the leakage-safe training graph and does not use treatment edges as features.

## Run

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements-app.txt
.venv/bin/python -m remedy.train
cd frontend && npm install && npm run build && cd ..
.venv/bin/uvicorn remedy.api:app --host 127.0.0.1 --port 8000
```

Open http://127.0.0.1:8000. Type a disease such as Parkinson, epilepsy, or ADHD. Tests for edge direction and path validity:

```bash
.venv/bin/pytest -q tests
```

## How to read a disease graph

Columns run left to right: drug, target protein, pathway, disease-associated protein (only on an indirect bridge), disease. A protein used by several drugs is drawn once.

Arrows keep the recorded direction:

| Relation | Direction in ReMedy-KG | What it means |
| --- | --- | --- |
| `targets` | drug → protein | Documented drug-target interaction |
| `participates_in` | protein → pathway | Documented pathway membership |
| `associated_with` | disease → protein | Documented disease association. The arrow points back toward the protein |
| `indication`, `off-label use` | drug → disease | Recorded use. Not inferred from the pathway |
| `contraindication` | drug → disease | Safety record. Never a treatment |
| predicted association | drug → disease, dashed | Neighborhood rank score. Not a graph edge and not a probability |

Click a drug to dim everything outside its path. The side panel lists the steps, shared proteins, shared pathways, and the evidence stored on those edges (source database, publication or document, assay, organism, activity, and the original node ids).

Turn on **Show indirect pathway bridges** only for drugs that share no protein with the disease. That view draws drug → its target → a pathway ← a different disease-associated protein ← the disease. It does not add a pathway–disease edge.

## What this graph cannot show

ReMedy-KG v1.0 has no pathway–disease edges, no gene layer, and no `inhibits` or `activates` relations. A shared protein or pathway is a graph connection, not a causal mechanism and not evidence of efficacy.

The dashed score is the neighborhood logistic ranker in `models/remedy_ranker.joblib` (held-out accuracy 65.9%, AUPRC 0.694). The earlier R-GCN checkpoint is not the score on the graph: its test scores sit near 0.49 and do not separate candidates. Those scores remain in `reports/model_evaluation/top20_explainable_candidate_report.csv`.

Held-out metrics are written to `reports/model_evaluation/ranker_metrics.json`. The header shows test accuracy, AUPRC, and AUROC at a cutoff of 0.50 chosen before the test split was scored.
