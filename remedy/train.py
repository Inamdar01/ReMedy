"""Train the drug-disease ranker and write held-out metrics.

Hyperparameters are chosen on the validation split by AUPRC.
The decision threshold is 0.5, fixed before the test split is scored.
The classifier is then refit on train and validation and applied once to test.
"""

from __future__ import annotations

import json

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from remedy.features import FEATURE_COLUMNS, BiologicalIndex
from remedy.paths import (
    METRICS_PATH,
    MODEL_PATH,
    NODES_PATH,
    PREDICTIONS_PATH,
    SPLITS_DIR,
    TRAINING_EDGES_PATH,
)


SEED = 42
C_GRID = [0.03, 0.1, 0.3, 1.0, 3.0]
THRESHOLD = 0.5


def _load_split(name: str) -> pd.DataFrame:
    frame = pd.read_csv(
        SPLITS_DIR / f"{name}_pairs.csv",
        dtype={"drug_id": str, "disease_id": str},
    )
    frame["label"] = frame["label"].astype(int)
    return frame


def _ranking_metrics(labels: np.ndarray, scores: np.ndarray) -> dict:
    labels = np.asarray(labels).astype(int)
    scores = np.asarray(scores).astype(float)
    order = np.argsort(-scores)
    k = min(10, len(labels))
    top = labels[order[:k]]
    positives = int(labels.sum())
    oracle = max(float(((scores >= cutoff) == labels).mean()) for cutoff in np.unique(scores))
    predicted = scores >= THRESHOLD
    correct = int((predicted == labels).sum())
    return {
        "auroc": float(roc_auc_score(labels, scores)),
        "auprc": float(average_precision_score(labels, scores)),
        "accuracy": float(correct / len(labels)),
        "correct": correct,
        "total": int(len(labels)),
        "precision_at_10": float(top.mean()),
        "recall_at_10": float(top.sum() / positives) if positives else 0.0,
        "oracle_accuracy": oracle,
    }


def _pipeline(c_value: float) -> Pipeline:
    return Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "classifier",
                LogisticRegression(
                    C=c_value,
                    max_iter=4000,
                    random_state=SEED,
                ),
            ),
        ]
    )


def train() -> dict:
    nodes = pd.read_parquet(NODES_PATH)
    edges = pd.read_parquet(TRAINING_EDGES_PATH)
    biology = BiologicalIndex(nodes, edges)

    splits = {name: _load_split(name) for name in ("train", "validation", "test")}
    features = {name: biology.feature_frame(frame) for name, frame in splits.items()}
    labels = {name: frame["label"].to_numpy() for name, frame in splits.items()}

    selection = []
    best_c = None
    best_ap = -1.0
    for c_value in C_GRID:
        candidate = _pipeline(c_value)
        candidate.fit(features["train"][FEATURE_COLUMNS], labels["train"])
        validation_scores = candidate.predict_proba(features["validation"][FEATURE_COLUMNS])[:, 1]
        validation_ap = float(average_precision_score(labels["validation"], validation_scores))
        selection.append({"C": c_value, "validation_auprc": validation_ap})
        if validation_ap > best_ap:
            best_ap = validation_ap
            best_c = c_value

    train_validation = pd.concat(
        [features["train"], features["validation"]],
        ignore_index=True,
    )
    train_validation_labels = np.concatenate([labels["train"], labels["validation"]])
    model = _pipeline(best_c)
    model.fit(train_validation[FEATURE_COLUMNS], train_validation_labels)

    test_scores = model.predict_proba(features["test"][FEATURE_COLUMNS])[:, 1]
    validation_model = _pipeline(best_c)
    validation_model.fit(features["train"][FEATURE_COLUMNS], labels["train"])
    validation_scores = validation_model.predict_proba(features["validation"][FEATURE_COLUMNS])[:, 1]

    similarity = (
        features["test"]["protein_jaccard"] + features["test"]["pathway_jaccard"]
    ).to_numpy()
    validation_similarity = (
        features["validation"]["protein_jaccard"] + features["validation"]["pathway_jaccard"]
    ).to_numpy()

    metrics = {
        "model_name": "Specific protein and pathway logistic regression",
        "dataset_version": "ReMedy-KG v1.0",
        "seed": SEED,
        "selected_C": best_c,
        "threshold": THRESHOLD,
        "threshold_rule": (
            "Fixed at 0.5 before test evaluation. "
            "C was selected by validation AUPRC on a model fit to training pairs only. "
            "The shipped model was then refit on training plus validation."
        ),
        "features": FEATURE_COLUMNS,
        "treatment_labels_used_as_features": False,
        "feature_relations": ["targets", "participates_in", "associated_with"],
        "validation_selection": selection,
        "validation": _ranking_metrics(labels["validation"], validation_scores),
        "test": _ranking_metrics(labels["test"], test_scores),
        "similarity_baseline_test": _ranking_metrics(labels["test"], similarity),
        "notes": (
            "Oracle accuracy is the best cutoff on that split's scores. "
            "It is reported so the ceiling is visible. The shipped cutoff remains 0.5."
        ),
    }
    # Similarity accuracy at 0.5 is not meaningful because Jaccard is not a probability.
    # Replace it with a cutoff chosen on validation similarity only.
    best_similarity_cutoff = 0.0
    best_similarity_accuracy = -1.0
    for cutoff in np.unique(validation_similarity):
        accuracy = float(((validation_similarity >= cutoff) == labels["validation"]).mean())
        if accuracy > best_similarity_accuracy:
            best_similarity_accuracy = accuracy
            best_similarity_cutoff = float(cutoff)
    similarity_predicted = similarity >= best_similarity_cutoff
    similarity_correct = int((similarity_predicted == labels["test"]).sum())
    metrics["similarity_baseline_test"]["accuracy"] = similarity_correct / len(labels["test"])
    metrics["similarity_baseline_test"]["correct"] = similarity_correct
    metrics["similarity_baseline_test"]["threshold"] = best_similarity_cutoff
    metrics["similarity_baseline_test"]["threshold_rule"] = (
        "Cutoff chosen on validation similarity, then applied once to test."
    )

    MODEL_PATH.parent.mkdir(parents=True, exist_ok=True)
    METRICS_PATH.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": model,
            "feature_columns": FEATURE_COLUMNS,
            "selected_C": best_c,
            "threshold": THRESHOLD,
            "dataset_version": "ReMedy-KG v1.0",
        },
        MODEL_PATH,
    )
    METRICS_PATH.write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    predictions = splits["test"][["drug_id", "disease_id", "label", "label_source"]].copy()
    predictions["score"] = test_scores
    predictions["predicted_label"] = (test_scores >= THRESHOLD).astype(int)
    predictions.sort_values(["score", "drug_id"], ascending=[False, True]).to_csv(
        PREDICTIONS_PATH,
        index=False,
    )
    return metrics


def main() -> None:
    metrics = train()
    test = metrics["test"]
    print(
        f"Selected C={metrics['selected_C']}  "
        f"validation AUPRC={metrics['validation']['auprc']:.3f}"
    )
    print(
        f"Test accuracy={test['accuracy']:.3f} "
        f"({test['correct']}/{test['total']})  "
        f"AUPRC={test['auprc']:.3f}  AUROC={test['auroc']:.3f}  "
        f"oracle accuracy={test['oracle_accuracy']:.3f}"
    )
    print(f"Saved {MODEL_PATH}")
    print(f"Saved {METRICS_PATH}")


if __name__ == "__main__":
    main()
