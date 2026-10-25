"""Disease graphs whose edges are copied from ReMedy-KG, never invented.

Preferred reading order, left to right:

    drug → protein it targets → pathway that protein belongs to → disease

The disease sits on the right, but `associated_with` keeps its recorded
direction: disease → protein. There is no pathway→disease relation in
ReMedy-KG v1.0, so that hop is named as missing instead of drawn.
"""

from __future__ import annotations

import pandas as pd

from remedy.features import BiologicalIndex
from remedy.path_check import PREDICTED_RELATION, review_edges


TREATMENT_RELATIONS = ("indication", "off-label use")
CONTRAINDICATION = "contraindication"

COLUMN_X = {
    "drug": 0,
    "protein": 340,
    "pathway": 680,
    "disease_protein": 1020,
    "disease": 1360,
}
ROW_GAP = 168

OMISSIONS = [
    "ReMedy-KG v1.0 has no pathway–disease edges. A pathway is only a protein membership.",
    "Disease–protein links are stored as disease → protein. The arrow points back toward the protein.",
    "There is no gene layer and no inhibits or activates relation.",
    "A dashed curve is a neighborhood rank score, not an R-GCN probability and not a recorded treatment.",
]


def short_label(value: str, limit: int = 28) -> str:
    text = " ".join(str(value).split())
    if len(text) <= limit:
        return text
    return text[: limit - 1] + "…"


def _public_evidence(edge_id: str, evidence_by_edge: dict[str, list]) -> list:
    published = []
    for item in evidence_by_edge.get(edge_id, [])[:3]:
        if isinstance(item, str):
            published.append(item)
            continue
        published.append(
            {
                "source": str(item.get("source", "")),
                "publication": str(item.get("publication", "")),
                "assay": str(item.get("assay", "")),
                "organism": str(item.get("organism", "")),
                "activity": str(item.get("activity", "")),
            }
        )
    return published


class GraphLibrary:
    def __init__(
        self,
        nodes: pd.DataFrame,
        edges: pd.DataFrame,
        biology: BiologicalIndex,
        evidence_by_edge: dict[str, list] | None = None,
    ):
        self.biology = biology
        self.evidence_by_edge = evidence_by_edge or {}
        self.labels = biology.labels
        self.node_types = biology.node_types
        self.sources = dict(
            zip(nodes["node_id"].astype(str), nodes["source_database"].astype(str))
        )
        frame = edges.copy()
        frame["source_id"] = frame["source_id"].astype(str)
        frame["target_id"] = frame["target_id"].astype(str)
        self.edges = frame
        self.lookup = {
            (row.source_id, row.relation, row.target_id): row
            for row in frame.itertuples(index=False)
        }

    def label(self, node_id: str) -> str:
        return self.labels.get(node_id, node_id)

    def overview(self) -> dict:
        diseases = self.biology.disease_ids
        nodes = []
        for index, disease_id in enumerate(diseases):
            column = index % 3
            row = index // 3
            nodes.append(
                self._entity(
                    disease_id,
                    "disease",
                    {"x": column * 320, "y": row * 170},
                    detail="Open this disease to read the supported path.",
                )
            )
        return {
            "mode": "overview",
            "title": "Choose a disease",
            "summary": (
                "Each disease opens as a left-to-right path: drug, target protein, "
                "pathway membership, then the disease. Missing hops stay missing."
            ),
            "omissions": OMISSIONS,
            "graph": {"nodes": nodes, "edges": []},
            "candidates": [],
            "known": [],
            "warnings": [],
        }

    def focus(
        self,
        disease_id: str,
        ranked_candidates: list[dict],
        *,
        known_limit: int = 3,
        candidate_limit: int = 3,
        protein_limit: int = 2,
        pathway_limit: int = 1,
        show_contraindications: bool = False,
        show_indirect: bool = False,
        contraindication_limit: int = 2,
    ) -> dict:
        disease_name = self.label(disease_id)
        disease_edges = self.edges.loc[self.edges["target_id"] == disease_id]
        known_by_drug: dict[str, dict] = {}
        for row in disease_edges.loc[
            disease_edges["relation"].isin(TREATMENT_RELATIONS)
        ].itertuples(index=False):
            slot = known_by_drug.setdefault(
                row.source_id,
                {"drug_id": row.source_id, "relations": set(), "evidence_count": 0},
            )
            slot["relations"].add(row.relation)
            slot["evidence_count"] += int(row.evidence_count)
        known_sorted = sorted(
            known_by_drug.values(),
            key=lambda item: (-item["evidence_count"], self.label(item["drug_id"])),
        )
        known_shown = known_sorted[:known_limit]

        contra_rows = disease_edges.loc[disease_edges["relation"] == CONTRAINDICATION]
        contra_evidence = {
            row.source_id: int(row.evidence_count) for row in contra_rows.itertuples(index=False)
        }
        treatment_pairs = {(item["drug_id"], disease_id) for item in known_sorted}
        contraindication_pairs = {(drug_id, disease_id) for drug_id in contra_evidence}

        candidates = []
        warnings: list[dict] = []
        for item in ranked_candidates:
            if item["drug_id"] in contra_evidence:
                warnings.append(
                    {
                        "level": "safety",
                        "edge_id": "",
                        "message": (
                            f"{self.label(item['drug_id'])} has a contraindication record "
                            "for this disease, so it is not drawn as a ranked treatment."
                        ),
                    }
                )
                continue
            if item["drug_id"] in known_by_drug:
                continue
            candidates.append(item)
            if len(candidates) >= candidate_limit:
                break

        nodes: list[dict] = []
        raw_edges: list[dict] = []
        node_ids: set[str] = set()
        drug_order: list[str] = []
        protein_order: list[str] = []
        pathway_order: list[str] = []
        disease_protein_order: list[str] = []
        known_payload = []
        candidate_payload = []
        drug_paths: dict[str, list] = {}
        drug_missing: dict[str, list[str]] = {}
        drug_meta: dict[str, dict] = {}
        drawn_evidence: dict[str, list] = {}

        def remember(node_id: str, bucket: list[str]) -> None:
            if node_id not in bucket:
                bucket.append(node_id)

        def add_catalog_edge(source: str, relation: str, target: str, kind: str, note: str, indirect: bool = False):
            row = self.lookup.get((source, relation, target))
            if row is None:
                return []
            if row.edge_id in drawn_evidence:
                return drawn_evidence[row.edge_id]
            evidence = _public_evidence(row.edge_id, self.evidence_by_edge)
            drawn_evidence[row.edge_id] = evidence
            raw_edges.append(
                {
                    "data": {
                        "id": row.edge_id,
                        "source": source,
                        "target": target,
                        "relation": relation,
                        "label": relation.replace("_", " "),
                        "kind": kind,
                        "direction": str(getattr(row, "direction", "forward")),
                        "evidence_count": int(row.evidence_count),
                        "evidence": evidence,
                        "source_id": source,
                        "target_id": target,
                        "source_type": self.node_types.get(source, ""),
                        "target_type": self.node_types.get(target, ""),
                        "note": note,
                        "indirect": indirect,
                    }
                }
            )
            return evidence

        def add_predicted(drug_id: str, score: float):
            raw_edges.append(
                {
                    "data": {
                        "id": f"predicted:{drug_id}:{disease_id}",
                        "source": drug_id,
                        "target": disease_id,
                        "relation": PREDICTED_RELATION,
                        "label": f"rank {score:.2f}",
                        "kind": "predicted",
                        "direction": "score only",
                        "evidence_count": 0,
                        "evidence": [],
                        "source_id": drug_id,
                        "target_id": disease_id,
                        "source_type": "drug",
                        "target_type": "disease",
                        "note": (
                            "Neighborhood rank score for an unobserved pair. "
                            "Not a recorded edge, not a probability, and not an R-GCN result."
                        ),
                        "indirect": False,
                        "score": score,
                    }
                }
            )

        def attach_direct(drug_id: str, proteins: list[str]) -> None:
            steps_protein = []
            steps_pathway = []
            evidence = []
            for protein_id in proteins:
                remember(protein_id, protein_order)
                target_evidence = add_catalog_edge(
                    drug_id,
                    "targets",
                    protein_id,
                    "target",
                    "Documented drug-target interaction. Not evidence that the drug treats the disease.",
                )
                association_evidence = add_catalog_edge(
                    disease_id,
                    "associated_with",
                    protein_id,
                    "association",
                    "Documented disease-protein association. The arrow points from the disease to the protein.",
                )
                evidence.extend(target_evidence)
                evidence.extend(association_evidence)
                steps_protein.append(
                    f"{self.label(drug_id)} —targets→ {self.label(protein_id)} ({protein_id})"
                )
                steps_protein.append(
                    f"{disease_name} —associated with→ {self.label(protein_id)} ({protein_id})"
                )
                protein_pathways = self.biology.protein_pathways.get(protein_id, set())
                disease_pathways = self.biology.disease_pathways.get(disease_id, set())
                ranked_pathways = sorted(
                    protein_pathways & disease_pathways,
                    key=lambda pathway: (-self.biology.pathway_idf.get(pathway, 0.0), pathway),
                )[:pathway_limit]
                if not ranked_pathways:
                    drug_missing[drug_id].append(
                        f"{self.label(protein_id)} has no pathway membership shared with this disease. "
                        "No pathway–disease edge exists, so that hop is not drawn."
                    )
                for pathway_id in ranked_pathways:
                    remember(pathway_id, pathway_order)
                    membership = add_catalog_edge(
                        protein_id,
                        "participates_in",
                        pathway_id,
                        "membership",
                        "Documented pathway membership. This graph has no pathway–disease edge.",
                    )
                    evidence.extend(membership)
                    steps_pathway.append(
                        f"{self.label(protein_id)} —participates in→ {self.label(pathway_id)}"
                    )
            if steps_protein:
                drug_paths[drug_id].append(
                    {
                        "title": "Shared protein",
                        "steps": steps_protein,
                        "note": (
                            "The same protein is a documented target and a documented disease association. "
                            "That is a graph connection, not proof of a mechanism or of efficacy."
                        ),
                        "evidence": evidence[:4],
                    }
                )
            if steps_pathway:
                drug_paths[drug_id].append(
                    {
                        "title": "Pathway membership",
                        "steps": steps_pathway,
                        "note": "No pathway–disease relationship is stored, so the path stops at the pathway.",
                        "evidence": [],
                    }
                )

        def attach_indirect(drug_id: str) -> bool:
            bridge = self._indirect_bridge(drug_id, disease_id)
            if bridge is None:
                drug_missing[drug_id].append(
                    "No shared protein and no shared pathway connect this drug to the disease."
                )
                return False
            protein_id, pathway_id, disease_protein_id = bridge
            if not show_indirect:
                drug_missing[drug_id].append(
                    "No shared protein. An indirect pathway bridge exists and stays hidden "
                    "until indirect pathway bridges are turned on."
                )
                return False
            remember(protein_id, protein_order)
            remember(pathway_id, pathway_order)
            remember(disease_protein_id, disease_protein_order)
            evidence = []
            evidence.extend(
                add_catalog_edge(
                    drug_id,
                    "targets",
                    protein_id,
                    "target",
                    "Documented drug-target interaction on an indirect path.",
                    indirect=True,
                )
            )
            evidence.extend(
                add_catalog_edge(
                    protein_id,
                    "participates_in",
                    pathway_id,
                    "membership",
                    "Both proteins participate in this pathway. This is not a pathway–disease edge.",
                    indirect=True,
                )
            )
            evidence.extend(
                add_catalog_edge(
                    disease_protein_id,
                    "participates_in",
                    pathway_id,
                    "membership",
                    "Disease-associated protein membership in the same pathway.",
                    indirect=True,
                )
            )
            evidence.extend(
                add_catalog_edge(
                    disease_id,
                    "associated_with",
                    disease_protein_id,
                    "association",
                    "Documented disease-protein association on an indirect path.",
                    indirect=True,
                )
            )
            drug_paths[drug_id].append(
                {
                    "title": "Indirect pathway bridge",
                    "steps": [
                        f"{self.label(drug_id)} —targets→ {self.label(protein_id)}",
                        f"{self.label(protein_id)} —participates in→ {self.label(pathway_id)}",
                        f"{self.label(disease_protein_id)} —participates in→ {self.label(pathway_id)}",
                        f"{disease_name} —associated with→ {self.label(disease_protein_id)}",
                    ],
                    "note": (
                        "The drug and the disease do not share a protein. "
                        "Two proteins meet only because both belong to this pathway. "
                        "That is not a causal mechanism and not evidence of efficacy."
                    ),
                    "evidence": evidence[:4],
                }
            )
            drug_missing[drug_id].append(
                "Direct drug-target to disease-protein link is missing. Only the indirect bridge is shown."
            )
            return True

        def prepare_drug(drug_id: str, role: str, **extra) -> None:
            remember(drug_id, drug_order)
            drug_paths[drug_id] = []
            drug_missing[drug_id] = []
            proteins = self.biology.shared_protein_ids(drug_id, disease_id)[:protein_limit]
            pathways = self.biology.shared_pathway_ids(drug_id, disease_id)
            if proteins:
                attach_direct(drug_id, proteins)
            else:
                attach_indirect(drug_id)
            if drug_id in contra_evidence:
                drug_missing[drug_id].append(
                    "A contraindication is also recorded. That safety record is not evidence of treatment."
                )
                warnings.append(
                    {
                        "level": "conflict",
                        "edge_id": "",
                        "message": (
                            f"{self.label(drug_id)} has both a treatment record and a contraindication "
                            "for this disease. The contraindication is not shown as supporting treatment."
                        ),
                    }
                )
            meta = {
                "role": role,
                "shared_proteins": [self.label(protein) for protein in proteins[:4]],
                "shared_pathways": [self.label(pathway) for pathway in pathways[:3]],
                **extra,
            }
            drug_meta[drug_id] = meta

        for item in known_shown:
            relation_label = " · ".join(sorted(item["relations"]))
            prepare_drug(
                item["drug_id"],
                "recorded",
                relation=relation_label,
                evidence_count=item["evidence_count"],
            )
            for relation in sorted(item["relations"]):
                add_catalog_edge(
                    item["drug_id"],
                    relation,
                    disease_id,
                    "recorded",
                    "Recorded treatment label. It is not inferred from the protein or the pathway.",
                )
            known_payload.append(
                {
                    "drug_id": item["drug_id"],
                    "drug_name": self.label(item["drug_id"]),
                    "relation": relation_label,
                    "evidence_count": item["evidence_count"],
                    "shared_proteins": drug_meta[item["drug_id"]]["shared_proteins"],
                    "shared_pathways": drug_meta[item["drug_id"]]["shared_pathways"],
                    "paths": drug_paths[item["drug_id"]],
                    "missing": drug_missing[item["drug_id"]],
                }
            )

        for rank, item in enumerate(candidates, start=1):
            prepare_drug(
                item["drug_id"],
                "candidate",
                score=item["score"],
                rank=rank,
                predicted=item["score"] >= 0.5,
            )
            add_predicted(item["drug_id"], item["score"])
            candidate_payload.append(
                {
                    "rank": rank,
                    "drug_id": item["drug_id"],
                    "drug_name": self.label(item["drug_id"]),
                    "score": item["score"],
                    "predicted": item["score"] >= 0.5,
                    "shared_proteins": drug_meta[item["drug_id"]]["shared_proteins"],
                    "shared_pathways": drug_meta[item["drug_id"]]["shared_pathways"],
                    "paths": drug_paths[item["drug_id"]],
                    "missing": drug_missing[item["drug_id"]],
                }
            )

        contra_shown = []
        if show_contraindications:
            contra_shown = [
                drug_id
                for drug_id in sorted(
                    contra_evidence,
                    key=lambda drug_id: (-contra_evidence[drug_id], self.label(drug_id)),
                )
                if drug_id not in drug_order
            ][:contraindication_limit]
        for drug_id in contra_shown:
            remember(drug_id, drug_order)
            drug_paths[drug_id] = []
            drug_missing[drug_id] = [
                "Safety relationship only. A contraindication is not evidence that the drug treats the disease."
            ]
            drug_meta[drug_id] = {"role": "safety", "shared_proteins": [], "shared_pathways": []}
            add_catalog_edge(
                drug_id,
                CONTRAINDICATION,
                disease_id,
                "contraindication",
                "Recorded contraindication. Excluded from candidate ranking.",
            )

        positions = {}
        for index, drug_id in enumerate(drug_order):
            positions[drug_id] = {"x": COLUMN_X["drug"], "y": 48 + index * ROW_GAP}
        for index, protein_id in enumerate(protein_order):
            positions[protein_id] = {"x": COLUMN_X["protein"], "y": 48 + index * ROW_GAP}
        for index, pathway_id in enumerate(pathway_order):
            positions[pathway_id] = {"x": COLUMN_X["pathway"], "y": 48 + index * ROW_GAP}
        for index, protein_id in enumerate(disease_protein_order):
            positions[protein_id] = {"x": COLUMN_X["disease_protein"], "y": 48 + index * ROW_GAP}
        anchor_rows = [positions[node_id]["y"] for node_id in protein_order] or [
            positions[node_id]["y"] for node_id in drug_order
        ]
        disease_y = sum(anchor_rows) / len(anchor_rows) if anchor_rows else 48
        positions[disease_id] = {"x": COLUMN_X["disease"], "y": disease_y}

        def place_entity(node_id: str, kind: str, **extra) -> None:
            if node_id in node_ids:
                return
            node_ids.add(node_id)
            nodes.append(self._entity(node_id, kind, positions[node_id], **extra))

        place_entity(
            disease_id,
            "disease",
            detail=(
                f"{len(known_sorted)} recorded treatment links. "
                f"{len(contra_evidence)} contraindicated drugs. "
                "Open a drug to read its supported path."
            ),
        )
        for drug_id in drug_order:
            meta = drug_meta[drug_id]
            role = meta["role"]
            kind = {
                "recorded": "known_drug",
                "candidate": "candidate_drug",
                "safety": "contraindication",
            }[role]
            path_ids = self._path_ids(drug_id, raw_edges)
            detail = self._drug_detail(meta)
            place_entity(
                drug_id,
                kind,
                detail=detail,
                path_ids=path_ids,
                paths=drug_paths.get(drug_id, []),
                missing=drug_missing.get(drug_id, []),
                **{key: value for key, value in meta.items() if key != "role"},
            )
        for protein_id in protein_order:
            place_entity(
                protein_id,
                "protein",
                detail=f"UniProt {protein_id}. Shown once, even when several drugs target it.",
                node_type="protein",
            )
        for pathway_id in pathway_order:
            place_entity(
                pathway_id,
                "pathway",
                detail="Pathway membership of a protein already on this path. Not a disease treatment.",
                node_type="pathway",
            )
        for protein_id in disease_protein_order:
            place_entity(
                protein_id,
                "disease_protein",
                detail=(
                    f"UniProt {protein_id}. Disease-associated protein on an indirect pathway bridge. "
                    "It is not a second copy of a drug target."
                ),
                node_type="protein",
            )

        self._column_headers(
            nodes,
            drugs=bool(drug_order),
            proteins=bool(protein_order),
            pathways=bool(pathway_order),
            disease_proteins=bool(disease_protein_order),
        )
        kept, review_warnings = review_edges(
            raw_edges,
            self.lookup,
            self.node_types,
            treatment_pairs=treatment_pairs,
            contraindication_pairs=contraindication_pairs,
        )
        warnings.extend(review_warnings)
        hidden_known = max(len(known_sorted) - len(known_shown), 0)
        return {
            "mode": "focus",
            "disease_id": disease_id,
            "disease_name": disease_name,
            "title": disease_name,
            "summary": (
                "Read left to right: drug, a protein it targets, a pathway that protein belongs to, "
                "then the disease. Arrows keep the recorded direction, so the disease arrow points "
                f"back at the protein. {len(known_shown)} recorded treatments"
                + (f" ({hidden_known} more not drawn)" if hidden_known else "")
                + f" and {len(candidate_payload)} ranked drugs."
            ),
            "omissions": OMISSIONS,
            "known_total": len(known_sorted),
            "contraindication_total": len(contra_evidence),
            "graph": {"nodes": nodes, "edges": kept},
            "candidates": candidate_payload,
            "known": known_payload,
            "warnings": warnings,
            "checks": {
                "displayed_edges": len(kept),
                "invalid": sum(1 for item in warnings if item["level"] == "invalid"),
            },
        }

    def _indirect_bridge(self, drug_id: str, disease_id: str) -> tuple[str, str, str] | None:
        drug_proteins = self.biology.drug_targets.get(drug_id, set())
        disease_proteins = self.biology.disease_proteins.get(disease_id, set())
        if drug_proteins & disease_proteins:
            return None
        for pathway_id in self.biology.shared_pathway_ids(drug_id, disease_id):
            drug_side = sorted(
                (
                    protein
                    for protein in drug_proteins
                    if pathway_id in self.biology.protein_pathways.get(protein, set())
                ),
                key=lambda protein: (-self.biology.protein_idf.get(protein, 0.0), protein),
            )
            disease_side = sorted(
                (
                    protein
                    for protein in disease_proteins
                    if pathway_id in self.biology.protein_pathways.get(protein, set())
                    and protein not in drug_proteins
                ),
                key=lambda protein: (-self.biology.protein_idf.get(protein, 0.0), protein),
            )
            if drug_side and disease_side:
                return drug_side[0], pathway_id, disease_side[0]
        return None

    def _path_ids(self, drug_id: str, edges: list[dict]) -> list[str]:
        """Nodes on this drug's path, without pulling in other drugs or their proteins."""
        allowed = {drug_id}
        changed = True
        while changed:
            changed = False
            for edge in edges:
                source = edge["data"]["source"]
                target = edge["data"]["target"]
                relation = edge["data"]["relation"]
                if relation == "associated_with":
                    if target in allowed and source not in allowed:
                        allowed.add(source)
                        changed = True
                    continue
                if source not in allowed and target not in allowed:
                    continue
                other = target if source in allowed else source
                if self.node_types.get(other) == "drug" and other != drug_id:
                    continue
                if other not in allowed:
                    allowed.add(other)
                    changed = True
        return sorted(allowed)

    def _drug_detail(self, meta: dict) -> str:
        if meta["role"] == "safety":
            return "Recorded contraindication. This is not a therapeutic path."
        if meta["role"] == "recorded":
            return (
                f"Recorded {meta.get('relation', 'treatment')}. "
                f"Evidence records: {meta.get('evidence_count', 0)}. "
                "The protein path is separate from this label."
            )
        score = meta.get("score")
        return (
            f"Rank score {score:.3f}. Unobserved pair. "
            "The score is not a probability and not a recorded treatment."
        )

    def _entity(self, node_id: str, kind: str, position: dict, **extra) -> dict:
        label = extra.pop("full_label", self.label(node_id))
        data = {
            "id": node_id,
            "label": short_label(label, 32 if kind == "disease" else 22),
            "full_label": label,
            "kind": kind,
            "node_type": self.node_types.get(node_id, extra.get("node_type", "")),
            "source_database": self.sources.get(node_id, ""),
        }
        data.update(extra)
        return {"data": data, "position": position}

    def _column_headers(self, nodes: list[dict], **present: bool) -> None:
        headers = [
            ("drug", "1 · Drug", present["drugs"]),
            ("protein", "2 · Target protein", present["proteins"]),
            ("pathway", "3 · Pathway", present["pathways"]),
            ("disease_protein", "4 · Disease protein", present["disease_proteins"]),
            ("disease", "5 · Disease", True),
        ]
        for key, title, show in headers:
            if not show:
                continue
            nodes.append(
                {
                    "data": {
                        "id": f"column:{key}",
                        "label": title,
                        "full_label": title,
                        "kind": "column",
                        "detail": title,
                    },
                    "position": {"x": COLUMN_X[key], "y": -78},
                }
            )
