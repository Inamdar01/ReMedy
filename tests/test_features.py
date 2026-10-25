import pandas as pd

from remedy.features import BiologicalIndex


def _toy():
    nodes = pd.DataFrame(
        [
            {"node_id": "DRUG", "node_type": "drug", "label": "Example drug", "source_database": "test"},
            {"node_id": "DISEASE", "node_type": "disease", "label": "Example disease", "source_database": "test"},
            {"node_id": "PROT", "node_type": "protein", "label": "Shared protein", "source_database": "test"},
            {"node_id": "OTHER", "node_type": "protein", "label": "Other protein", "source_database": "test"},
            {"node_id": "PATH", "node_type": "pathway", "label": "Shared pathway", "source_database": "test"},
        ]
    )
    edges = pd.DataFrame(
        [
            {"source_id": "DRUG", "relation": "targets", "target_id": "PROT"},
            {"source_id": "DRUG", "relation": "indication", "target_id": "DISEASE"},
            {"source_id": "DISEASE", "relation": "associated_with", "target_id": "PROT"},
            {"source_id": "DISEASE", "relation": "associated_with", "target_id": "OTHER"},
            {"source_id": "PROT", "relation": "participates_in", "target_id": "PATH"},
            {"source_id": "OTHER", "relation": "participates_in", "target_id": "PATH"},
        ]
    )
    return BiologicalIndex(nodes, edges)


def test_treatment_edges_do_not_become_features():
    biology = _toy()
    features = biology.features("DRUG", "DISEASE")
    assert features["log_shared_proteins"] > 0
    assert biology.shared_protein_ids("DRUG", "DISEASE") == ["PROT"]
    assert "indication" not in biology.features("DRUG", "DISEASE")


def test_unconnected_drug_has_zero_shared_proteins():
    biology = _toy()
    features = biology.features("MISSING", "DISEASE")
    assert features["protein_jaccard"] == 0
    assert features["protein_idf"] == 0
    assert features["log_shared_proteins"] == 0
