import pandas as pd

from remedy.diseases import resolve_disease
from remedy.features import BiologicalIndex
from remedy.graph_view import GraphLibrary
from remedy.path_check import review_edges


DISEASES = [
    {"id": "MONDO:0005180", "label": "Parkinson disease", "aliases": ["parkinson", "parkinson's disease"]},
    {"id": "MONDO:0004975", "label": "Alzheimer disease", "aliases": ["alzheimer"]},
    {"id": "MONDO:0005098", "label": "stroke disorder", "aliases": ["stroke"]},
]


def _library(evidence=None) -> GraphLibrary:
    nodes = pd.DataFrame(
        [
            {"node_id": "MONDO:1", "node_type": "disease", "label": "Example disease", "source_database": "MONDO"},
            {"node_id": "CHEMBL1", "node_type": "drug", "label": "Known drug", "source_database": "ChEMBL"},
            {"node_id": "CHEMBL2", "node_type": "drug", "label": "Candidate drug", "source_database": "ChEMBL"},
            {"node_id": "CHEMBL3", "node_type": "drug", "label": "Unsafe drug", "source_database": "ChEMBL"},
            {"node_id": "CHEMBL4", "node_type": "drug", "label": "Indirect drug", "source_database": "ChEMBL"},
            {"node_id": "P11111", "node_type": "protein", "label": "Protein one", "source_database": "UniProt"},
            {"node_id": "P22222", "node_type": "protein", "label": "Target protein", "source_database": "UniProt"},
            {"node_id": "P33333", "node_type": "protein", "label": "Disease protein", "source_database": "UniProt"},
            {"node_id": "R-HSA-1", "node_type": "pathway", "label": "Pathway one", "source_database": "Reactome"},
            {"node_id": "R-HSA-2", "node_type": "pathway", "label": "Pathway two", "source_database": "Reactome"},
        ]
    )
    edges = pd.DataFrame(
        [
            {"edge_id": "e1", "source_id": "CHEMBL1", "relation": "indication", "target_id": "MONDO:1", "direction": "forward", "evidence_count": 2},
            {"edge_id": "e2", "source_id": "CHEMBL1", "relation": "targets", "target_id": "P11111", "direction": "forward", "evidence_count": 4},
            {"edge_id": "e3", "source_id": "CHEMBL2", "relation": "targets", "target_id": "P11111", "direction": "forward", "evidence_count": 3},
            {"edge_id": "e4", "source_id": "MONDO:1", "relation": "associated_with", "target_id": "P11111", "direction": "forward", "evidence_count": 5},
            {"edge_id": "e5", "source_id": "P11111", "relation": "participates_in", "target_id": "R-HSA-1", "direction": "forward", "evidence_count": 1},
            {"edge_id": "e6", "source_id": "CHEMBL3", "relation": "contraindication", "target_id": "MONDO:1", "direction": "forward", "evidence_count": 1},
            {"edge_id": "e7", "source_id": "CHEMBL4", "relation": "targets", "target_id": "P22222", "direction": "forward", "evidence_count": 2},
            {"edge_id": "e8", "source_id": "P22222", "relation": "participates_in", "target_id": "R-HSA-2", "direction": "forward", "evidence_count": 1},
            {"edge_id": "e9", "source_id": "P33333", "relation": "participates_in", "target_id": "R-HSA-2", "direction": "forward", "evidence_count": 1},
            {"edge_id": "e10", "source_id": "MONDO:1", "relation": "associated_with", "target_id": "P33333", "direction": "forward", "evidence_count": 2},
        ]
    )
    biology = BiologicalIndex(nodes, edges)
    return GraphLibrary(nodes, edges, biology, evidence or {"e2": ["ChEMBL · Ki = 10 nM"]})


def test_parkinson_alias_resolves_to_one_disease():
    resolved = resolve_disease("Parkinson's", DISEASES)
    assert resolved["status"] == "match"
    assert resolved["disease_id"] == "MONDO:0005180"


def test_short_query_does_not_match_every_disease():
    resolved = resolve_disease("di", DISEASES)
    assert resolved["status"] == "none"


def test_focus_keeps_real_direction_and_one_protein():
    view = _library().focus(
        "MONDO:1",
        [
            {"drug_id": "CHEMBL2", "score": 0.82},
            {"drug_id": "CHEMBL3", "score": 0.91},
        ],
    )
    by_id = {node["data"]["id"]: node for node in view["graph"]["nodes"]}
    assert [node["data"]["id"] for node in view["graph"]["nodes"]].count("P11111") == 1
    assert by_id["CHEMBL1"]["position"]["x"] < by_id["P11111"]["position"]["x"]
    assert by_id["P11111"]["position"]["x"] < by_id["R-HSA-1"]["position"]["x"]
    assert by_id["R-HSA-1"]["position"]["x"] < by_id["MONDO:1"]["position"]["x"]

    relations = {edge["data"]["relation"] for edge in view["graph"]["edges"]}
    assert "targets" in relations
    assert "participates_in" in relations
    assert "associated_with" in relations
    assert "indication" in relations
    assert "predicted association" in relations
    assert not any(
        edge["data"]["target"] == "MONDO:1" and edge["data"]["source"].startswith("R-HSA")
        for edge in view["graph"]["edges"]
    )

    association = next(edge for edge in view["graph"]["edges"] if edge["data"]["id"] == "e4")
    assert association["data"]["source"] == "MONDO:1"
    assert association["data"]["target"] == "P11111"
    assert association["data"]["kind"] == "association"

    target = next(edge for edge in view["graph"]["edges"] if edge["data"]["id"] == "e2")
    assert target["data"]["source"] == "CHEMBL1"
    assert target["data"]["target"] == "P11111"
    assert target["data"]["evidence"] == ["ChEMBL · Ki = 10 nM"]

    predicted = next(edge for edge in view["graph"]["edges"] if edge["data"]["kind"] == "predicted")
    assert predicted["data"]["source"] == "CHEMBL2"
    assert predicted["data"]["target"] == "MONDO:1"
    assert "not a recorded edge" in predicted["data"]["note"].lower()
    assert view["candidates"][0]["predicted"] is True
    assert all(candidate["drug_id"] != "CHEMBL3" for candidate in view["candidates"])
    assert any("contraindication" in warning["message"] for warning in view["warnings"])
    assert view["checks"]["invalid"] == 0
    assert any("provenance" == warning["level"] for warning in view["warnings"])
    assert "P11111" in by_id["CHEMBL1"]["data"]["path_ids"]
    assert "CHEMBL2" not in by_id["CHEMBL1"]["data"]["path_ids"]


def test_indirect_bridge_is_hidden_until_requested_and_never_invents_a_disease_edge():
    hidden = _library().focus("MONDO:1", [{"drug_id": "CHEMBL4", "score": 0.4}])
    hidden_ids = {node["data"]["id"] for node in hidden["graph"]["nodes"]}
    assert "P22222" not in hidden_ids
    assert any("indirect pathway bridge" in item for item in hidden["candidates"][0]["missing"])

    shown = _library().focus(
        "MONDO:1",
        [{"drug_id": "CHEMBL4", "score": 0.4}],
        show_indirect=True,
    )
    shown_ids = {node["data"]["id"] for node in shown["graph"]["nodes"]}
    assert {"CHEMBL4", "P22222", "R-HSA-2", "P33333", "MONDO:1"} <= shown_ids
    disease_protein = next(node for node in shown["graph"]["nodes"] if node["data"]["id"] == "P33333")
    pathway = next(node for node in shown["graph"]["nodes"] if node["data"]["id"] == "R-HSA-2")
    disease = next(node for node in shown["graph"]["nodes"] if node["data"]["id"] == "MONDO:1")
    assert pathway["position"]["x"] < disease_protein["position"]["x"] < disease["position"]["x"]
    memberships = [
        edge
        for edge in shown["graph"]["edges"]
        if edge["data"]["relation"] == "participates_in" and edge["data"]["target"] == "R-HSA-2"
    ]
    assert {edge["data"]["source"] for edge in memberships} == {"P22222", "P33333"}
    assert not any(edge["data"]["target"] == "MONDO:1" and edge["data"]["kind"] == "membership" for edge in shown["graph"]["edges"])


def test_contraindication_stays_a_safety_edge():
    view = _library().focus(
        "MONDO:1",
        [{"drug_id": "CHEMBL2", "score": 0.2}],
        show_contraindications=True,
    )
    contra = next(edge for edge in view["graph"]["edges"] if edge["data"]["relation"] == "contraindication")
    assert contra["data"]["kind"] == "contraindication"
    assert contra["data"]["source"] == "CHEMBL3"
    assert contra["data"]["target"] == "MONDO:1"
    assert not any(
        edge["data"]["relation"] == "contraindication" and edge["data"]["kind"] == "recorded"
        for edge in view["graph"]["edges"]
    )


def test_review_drops_a_fabricated_path_and_a_prediction_on_a_contraindication():
    nodes = {"CHEMBL9": "drug", "MONDO:1": "disease", "P11111": "protein"}
    fabricated = {
        "data": {
            "id": "fake",
            "source": "R-HSA-9",
            "target": "MONDO:1",
            "relation": "participates_in",
            "kind": "membership",
            "evidence_count": 0,
            "evidence": [],
        }
    }
    unsafe_score = {
        "data": {
            "id": "predicted:CHEMBL9:MONDO:1",
            "source": "CHEMBL9",
            "target": "MONDO:1",
            "relation": "predicted association",
            "kind": "predicted",
            "evidence_count": 0,
            "evidence": [],
        }
    }
    kept, warnings = review_edges(
        [fabricated, unsafe_score],
        {},
        nodes,
        treatment_pairs=set(),
        contraindication_pairs={("CHEMBL9", "MONDO:1")},
    )
    assert kept == []
    assert any(warning["level"] == "invalid" for warning in warnings)
    assert any("contraindication" in warning["message"] for warning in warnings)
