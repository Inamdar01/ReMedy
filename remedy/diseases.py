"""Resolve a typed disease name onto the 12 ReMedy diseases."""

from __future__ import annotations

import re


DISEASE_ALIASES = {
    "MONDO:0001627": ["dementia"],
    "MONDO:0002009": [
        "major depressive disorder",
        "major depression",
        "depression",
        "mdd",
    ],
    "MONDO:0004975": [
        "alzheimer disease",
        "alzheimer's disease",
        "alzheimers disease",
        "alzheimers",
        "alzheimer",
    ],
    "MONDO:0004985": ["bipolar disorder", "bipolar"],
    "MONDO:0005027": ["epilepsy"],
    "MONDO:0005090": ["schizophrenia"],
    "MONDO:0005098": ["stroke disorder", "stroke"],
    "MONDO:0005180": [
        "parkinson disease",
        "parkinson's disease",
        "parkinsons disease",
        "parkinsons",
        "parkinson",
    ],
    "MONDO:0005277": ["migraine disorder", "migraine"],
    "MONDO:0005301": ["multiple sclerosis"],
    "MONDO:0005618": ["anxiety disorder", "anxiety"],
    "MONDO:0007743": [
        "attention deficit-hyperactivity disorder",
        "attention deficit hyperactivity disorder",
        "adhd",
    ],
}


def normalize_query(value: str) -> str:
    text = value.casefold().strip()
    text = text.replace("’", "'")
    text = re.sub(r"[^a-z0-9:+'\- ]+", " ", text)
    text = text.replace("-", " ")
    return re.sub(r"\s+", " ", text).strip()


def resolve_disease(query: str, diseases: list[dict]) -> dict:
    """Match free text to one disease, or return the competing matches."""
    cleaned = normalize_query(query)
    if not cleaned:
        return {"status": "empty", "matches": []}

    by_id = {disease["id"]: disease for disease in diseases}
    if cleaned.upper().startswith("MONDO:") or cleaned.startswith("mondo:"):
        mondo_id = cleaned.upper()
        if mondo_id in by_id:
            return {"status": "match", "disease_id": mondo_id, "matches": [by_id[mondo_id]]}
        return {"status": "none", "matches": []}

    exact = []
    partial = []
    for disease in diseases:
        names = [disease["label"], *disease.get("aliases", [])]
        normalized = [normalize_query(name) for name in names]
        if cleaned in normalized:
            exact.append(disease)
        elif len(cleaned) >= 4 and any(cleaned in name for name in normalized):
            partial.append(disease)

    if len(exact) == 1:
        return {"status": "match", "disease_id": exact[0]["id"], "matches": exact}
    if len(exact) > 1:
        return {"status": "ambiguous", "matches": exact}
    if len(partial) == 1:
        return {"status": "match", "disease_id": partial[0]["id"], "matches": partial}
    if len(partial) > 1:
        return {"status": "ambiguous", "matches": partial}
    return {"status": "none", "matches": []}
