"""Research API for disease graphs and held-out ranking scores."""

from __future__ import annotations

import json
from collections import defaultdict

import joblib
import pandas as pd
from fastapi import FastAPI, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from remedy.diseases import DISEASE_ALIASES, resolve_disease
from remedy.features import BiologicalIndex
from remedy.graph_view import GraphLibrary
from remedy.paths import (
    EDGES_PATH,
    EVIDENCE_PATH,
    FRONTEND_DIST,
    METRICS_PATH,
    MODEL_PATH,
    NODES_PATH,
    TRAINING_EDGES_PATH,
)
from remedy.train import train


TREATMENT_RELATIONS = {"indication", "off-label use"}


def _evidence_records(path) -> dict[str, list[dict]]:
    frame = pd.read_parquet(
        path,
        columns=["edge_id", "source", "publication", "assay", "organism", "activity_value"],
    )
    records: dict[str, list[dict]] = defaultdict(list)
    for row in frame.itertuples(index=False):
        bucket = records[row.edge_id]
        if len(bucket) >= 3:
            continue
        bucket.append(
            {
                "source": str(row.source),
                "publication": str(row.publication),
                "assay": str(row.assay),
                "organism": str(row.organism),
                "activity": str(row.activity_value),
            }
        )
    return dict(records)


class Catalog:
    def __init__(self):
        if not MODEL_PATH.exists() or not METRICS_PATH.exists():
            train()
        nodes = pd.read_parquet(NODES_PATH)
        edges = pd.read_parquet(EDGES_PATH)
        training_edges = pd.read_parquet(TRAINING_EDGES_PATH)
        self.biology = BiologicalIndex(nodes, training_edges)
        self.library = GraphLibrary(
            nodes,
            edges,
            self.biology,
            _evidence_records(EVIDENCE_PATH),
        )
        bundle = joblib.load(MODEL_PATH)
        self.model = bundle["model"]
        self.feature_columns = bundle["feature_columns"]
        self.threshold = float(bundle["threshold"])
        self.metrics = json.loads(METRICS_PATH.read_text(encoding="utf-8"))
        self.diseases = [
            {
                "id": disease_id,
                "label": self.biology.labels[disease_id],
                "aliases": DISEASE_ALIASES.get(disease_id, []),
            }
            for disease_id in self.biology.disease_ids
        ]
        self._ranks: dict[str, list[dict]] = {}
        blocked = edges.loc[
            edges["relation"].isin(TREATMENT_RELATIONS | {"contraindication"}),
            ["source_id", "target_id"],
        ]
        self.blocked = {
            disease_id: set(group["source_id"].astype(str))
            for disease_id, group in blocked.groupby(blocked["target_id"].astype(str))
        }

    def rank(self, disease_id: str) -> list[dict]:
        cached = self._ranks.get(disease_id)
        if cached is not None:
            return cached
        excluded = self.blocked.get(disease_id, set())
        drug_ids = [drug_id for drug_id in self.biology.drug_ids if drug_id not in excluded]
        pairs = pd.DataFrame({"drug_id": drug_ids, "disease_id": disease_id})
        features = self.biology.feature_frame(pairs)
        scores = self.model.predict_proba(features[self.feature_columns])[:, 1]
        ranked = [
            {"drug_id": drug_id, "score": round(float(score), 6)}
            for drug_id, score in zip(drug_ids, scores)
        ]
        ranked.sort(key=lambda item: (-item["score"], item["drug_id"]))
        self._ranks[disease_id] = ranked
        return ranked

    def focus(self, disease_id: str, show_contraindications: bool, show_indirect: bool) -> dict:
        payload = self.library.focus(
            disease_id,
            self.rank(disease_id),
            show_contraindications=show_contraindications,
            show_indirect=show_indirect,
        )
        payload["threshold"] = self.threshold
        return payload


def create_app() -> FastAPI:
    catalog_box: dict[str, Catalog] = {}

    app = FastAPI(title="ReMedy", version="1.0")
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_methods=["*"],
        allow_headers=["*"],
    )

    @app.on_event("startup")
    def load_catalog() -> None:
        catalog_box["catalog"] = Catalog()

    def catalog() -> Catalog:
        return catalog_box["catalog"]

    @app.get("/api/health")
    def health() -> dict:
        return {"status": "ok"}

    @app.get("/api/metrics")
    def metrics() -> dict:
        return catalog().metrics

    @app.get("/api/diseases")
    def diseases() -> dict:
        return {
            "diseases": [
                {"id": disease["id"], "label": disease["label"]}
                for disease in catalog().diseases
            ]
        }

    @app.get("/api/view")
    def view(
        q: str = "",
        disease_id: str = "",
        contraindications: bool = Query(False),
        indirect: bool = Query(False),
    ) -> dict:
        current = catalog()
        if disease_id:
            known = {disease["id"] for disease in current.diseases}
            if disease_id not in known:
                return {
                    "mode": "none",
                    "summary": "That disease identifier is not in ReMedy-KG v1.0.",
                    "matches": [],
                    "graph": {"nodes": [], "edges": []},
                }
            return current.focus(disease_id, contraindications, indirect)

        if q.strip():
            resolved = resolve_disease(q, current.diseases)
            if resolved["status"] == "match":
                payload = current.focus(resolved["disease_id"], contraindications, indirect)
                payload["query"] = q
                return payload
            if resolved["status"] == "ambiguous":
                return {
                    "mode": "ambiguous",
                    "summary": "More than one disease matches that text.",
                    "matches": [
                        {"id": disease["id"], "label": disease["label"]}
                        for disease in resolved["matches"]
                    ],
                    "graph": {"nodes": [], "edges": []},
                }
            return {
                "mode": "none",
                "summary": "No disease in this graph matches that text. The twelve named diseases are listed beside the search.",
                "matches": [],
                "graph": {"nodes": [], "edges": []},
            }

        overview = current.library.overview()
        overview["threshold"] = current.threshold
        return overview

    if FRONTEND_DIST.exists():
        app.mount("/", StaticFiles(directory=FRONTEND_DIST, html=True), name="ui")

    return app


app = create_app()
