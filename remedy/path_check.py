"""Check that a drawn path uses real edges and the recorded direction.

Displayed edges must already exist in the selected graph, except a predicted
drug-disease score, which is not a graph edge. A contraindication is never
accepted as a treatment edge.
"""

from __future__ import annotations


ENDPOINT_TYPES = {
    "targets": ("drug", "protein"),
    "participates_in": ("protein", "pathway"),
    "associated_with": ("disease", "protein"),
    "indication": ("drug", "disease"),
    "off-label use": ("drug", "disease"),
    "contraindication": ("drug", "disease"),
}

PREDICTED_RELATION = "predicted association"

# Kinds that are allowed to leave the knowledge graph.
SYNTHETIC_KINDS = {"predicted"}


def endpoint_issue(node_types: dict[str, str], source: str, relation: str, target: str) -> str | None:
    expected = ENDPOINT_TYPES.get(relation)
    if expected is None:
        return f"Unknown relation {relation}."
    if source not in node_types:
        return f"Missing node {source}."
    if target not in node_types:
        return f"Missing node {target}."
    actual = (node_types.get(source), node_types.get(target))
    if actual != expected:
        return (
            f"{relation} expects {expected[0]} → {expected[1]}, "
            f"got {actual[0]} → {actual[1]} ({source} → {target})."
        )
    return None


def review_edges(
    edges: list[dict],
    lookup: dict,
    node_types: dict[str, str],
    *,
    treatment_pairs: set[tuple[str, str]],
    contraindication_pairs: set[tuple[str, str]],
) -> tuple[list[dict], list[dict]]:
    """Drop edges that are missing, mistyped, duplicated, or unsafe.

    Provenance gaps are reported and the edge is kept, because the relationship
    itself is still in the graph.
    """

    kept: list[dict] = []
    warnings: list[dict] = []
    seen: set[tuple[str, str, str]] = set()
    for edge in edges:
        data = edge["data"]
        source = str(data["source"])
        target = str(data["target"])
        relation = str(data["relation"])
        kind = str(data["kind"])
        key = (source, relation, target)
        if key in seen:
            warnings.append(
                {
                    "level": "invalid",
                    "edge_id": data.get("id", ""),
                    "message": f"Duplicate edge {source} | {relation} | {target}.",
                }
            )
            continue
        seen.add(key)

        if relation == "contraindication" and kind != "contraindication":
            warnings.append(
                {
                    "level": "invalid",
                    "edge_id": data.get("id", ""),
                    "message": (
                        f"{source} → {target} is a contraindication and cannot be shown as a treatment."
                    ),
                }
            )
            continue

        if kind in SYNTHETIC_KINDS or relation == PREDICTED_RELATION:
            if (source, target) in contraindication_pairs:
                warnings.append(
                    {
                        "level": "invalid",
                        "edge_id": data.get("id", ""),
                        "message": (
                            f"Predicted score for {source} → {target} contradicts a contraindication record."
                        ),
                    }
                )
                continue
            if (source, target) in treatment_pairs:
                warnings.append(
                    {
                        "level": "invalid",
                        "edge_id": data.get("id", ""),
                        "message": (
                            f"Predicted score for {source} → {target} repeats a recorded treatment."
                        ),
                    }
                )
                continue
            if key in lookup or endpoint_issue(node_types, source, "indication", target):
                issue = endpoint_issue(node_types, source, "indication", target)
                if issue:
                    warnings.append(
                        {
                            "level": "invalid",
                            "edge_id": data.get("id", ""),
                            "message": issue,
                        }
                    )
                    continue
            kept.append(edge)
            continue

        row = lookup.get(key)
        if row is None:
            warnings.append(
                {
                    "level": "invalid",
                    "edge_id": data.get("id", ""),
                    "message": f"No graph edge for {source} | {relation} | {target}.",
                }
            )
            continue
        issue = endpoint_issue(node_types, source, relation, target)
        if issue:
            warnings.append({"level": "invalid", "edge_id": data.get("id", ""), "message": issue})
            continue
        evidence_count = int(data.get("evidence_count") or 0)
        if evidence_count > 0 and not data.get("evidence"):
            warnings.append(
                {
                    "level": "provenance",
                    "edge_id": str(data.get("id", "")),
                    "message": (
                        f"Edge {data.get('id', '')} records {evidence_count} evidence "
                        "rows, but none were attached for display."
                    ),
                }
            )
        kept.append(edge)
    return kept, warnings
