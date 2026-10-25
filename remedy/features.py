"""Leakage-safe biological features for drug-disease ranking.

Features use only targets, participates_in, and associated_with.
Treatment and contraindication edges are never turned into features.
"""

from __future__ import annotations

import math
from collections import defaultdict

import pandas as pd


FEATURE_COLUMNS = [
    "protein_jaccard",
    "protein_aa",
    "protein_ra",
    "protein_idf",
    "protein_idf_per_target",
    "pathway_jaccard",
    "pathway_aa",
    "pathway_idf",
    "log_shared_proteins",
    "log_shared_pathways",
]

SUPPORT_RELATIONS = {"targets", "participates_in", "associated_with"}


def _idf(document_frequency: int, collection_size: int) -> float:
    return math.log((collection_size + 1) / (document_frequency + 1)) + 1.0


def _jaccard(left: set[str], right: set[str]) -> float:
    union = len(left | right)
    if union == 0:
        return 0.0
    return len(left & right) / union


class BiologicalIndex:
    """Neighborhoods built from the leakage-safe training graph."""

    def __init__(self, nodes: pd.DataFrame, edges: pd.DataFrame):
        self.labels = dict(zip(nodes["node_id"].astype(str), nodes["label"].astype(str)))
        self.node_types = dict(
            zip(nodes["node_id"].astype(str), nodes["node_type"].astype(str))
        )
        self.drug_ids = sorted(
            nodes.loc[nodes["node_type"] == "drug", "node_id"].astype(str)
        )
        self.disease_ids = sorted(
            nodes.loc[nodes["node_type"] == "disease", "node_id"].astype(str)
        )

        drug_targets: dict[str, set[str]] = defaultdict(set)
        disease_proteins: dict[str, set[str]] = defaultdict(set)
        protein_pathways: dict[str, set[str]] = defaultdict(set)
        protein_drugs: dict[str, set[str]] = defaultdict(set)
        pathway_proteins: dict[str, set[str]] = defaultdict(set)

        unexpected = set()
        for row in edges.itertuples(index=False):
            relation = row.relation
            if relation not in SUPPORT_RELATIONS:
                continue
            source = str(row.source_id)
            target = str(row.target_id)
            if relation == "targets":
                drug_targets[source].add(target)
                protein_drugs[target].add(source)
            elif relation == "associated_with":
                disease_proteins[source].add(target)
            elif relation == "participates_in":
                protein_pathways[source].add(target)
                pathway_proteins[target].add(source)
            else:
                unexpected.add(relation)

        if unexpected:
            raise ValueError(f"Unsupported relations reached the index: {unexpected}")

        self.drug_targets = dict(drug_targets)
        self.disease_proteins = dict(disease_proteins)
        self.protein_pathways = dict(protein_pathways)

        self.drug_pathways = {
            drug: set().union(*(protein_pathways.get(protein, set()) for protein in proteins))
            for drug, proteins in self.drug_targets.items()
        }
        self.disease_pathways = {
            disease: set().union(
                *(protein_pathways.get(protein, set()) for protein in proteins)
            )
            for disease, proteins in self.disease_proteins.items()
        }

        drug_count = max(len(self.drug_ids), 1)
        pathway_count = max(len(pathway_proteins), 1)
        self.protein_idf = {
            protein: _idf(len(drugs), drug_count) for protein, drugs in protein_drugs.items()
        }
        self.pathway_idf = {
            pathway: _idf(len(members), pathway_count)
            for pathway, members in pathway_proteins.items()
        }
        self.protein_drug_degree = {protein: len(drugs) for protein, drugs in protein_drugs.items()}
        self.pathway_degree = {
            pathway: len(members) for pathway, members in pathway_proteins.items()
        }

    def _adamic_adar(self, shared: set[str], degree: dict[str, int]) -> float:
        score = 0.0
        for node_id in shared:
            node_degree = degree.get(node_id, 1)
            if node_degree > 1:
                score += 1.0 / math.log(node_degree)
        return score

    def _resource_allocation(self, shared: set[str], degree: dict[str, int]) -> float:
        return sum(1.0 / max(degree.get(node_id, 1), 1) for node_id in shared)

    def shared_protein_ids(self, drug_id: str, disease_id: str) -> list[str]:
        shared = self.drug_targets.get(drug_id, set()) & self.disease_proteins.get(
            disease_id, set()
        )
        return sorted(shared, key=lambda protein: (-self.protein_idf.get(protein, 0.0), protein))

    def shared_pathway_ids(self, drug_id: str, disease_id: str) -> list[str]:
        shared = self.drug_pathways.get(drug_id, set()) & self.disease_pathways.get(
            disease_id, set()
        )
        return sorted(shared, key=lambda pathway: (-self.pathway_idf.get(pathway, 0.0), pathway))

    def features(self, drug_id: str, disease_id: str) -> dict[str, float]:
        drug_proteins = self.drug_targets.get(drug_id, set())
        disease_proteins = self.disease_proteins.get(disease_id, set())
        shared_proteins = drug_proteins & disease_proteins
        drug_pathways = self.drug_pathways.get(drug_id, set())
        disease_pathways = self.disease_pathways.get(disease_id, set())
        shared_pathways = drug_pathways & disease_pathways
        target_count = max(len(drug_proteins), 1)
        protein_idf = sum(self.protein_idf.get(protein, 0.0) for protein in shared_proteins)
        return {
            "protein_jaccard": _jaccard(drug_proteins, disease_proteins),
            "protein_aa": self._adamic_adar(shared_proteins, self.protein_drug_degree),
            "protein_ra": self._resource_allocation(shared_proteins, self.protein_drug_degree),
            "protein_idf": protein_idf,
            "protein_idf_per_target": protein_idf / target_count,
            "pathway_jaccard": _jaccard(drug_pathways, disease_pathways),
            "pathway_aa": self._adamic_adar(shared_pathways, self.pathway_degree),
            "pathway_idf": sum(self.pathway_idf.get(pathway, 0.0) for pathway in shared_pathways),
            "log_shared_proteins": math.log1p(len(shared_proteins)),
            "log_shared_pathways": math.log1p(len(shared_pathways)),
        }

    def feature_frame(self, pairs: pd.DataFrame) -> pd.DataFrame:
        rows = [
            self.features(str(row.drug_id), str(row.disease_id))
            for row in pairs.itertuples(index=False)
        ]
        return pd.DataFrame(rows, columns=FEATURE_COLUMNS)
