from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

NODES_PATH = ROOT / "data" / "processed" / "nodes.parquet"
EDGES_PATH = ROOT / "data" / "processed" / "edges.parquet"
EVIDENCE_PATH = ROOT / "data" / "processed" / "evidence.parquet"
TRAINING_EDGES_PATH = ROOT / "data" / "model" / "training_graph_edges.parquet"
SPLITS_DIR = ROOT / "data" / "splits"
MODEL_PATH = ROOT / "models" / "remedy_ranker.joblib"
METRICS_PATH = ROOT / "reports" / "model_evaluation" / "ranker_metrics.json"
PREDICTIONS_PATH = ROOT / "reports" / "model_evaluation" / "ranker_test_predictions.csv"
FRONTEND_DIST = ROOT / "frontend" / "dist"
