"""Phase 5 baselines evaluation harness — scGen-style CELL-LEVEL prediction.

REWRITTEN 2026-05-13 after user-flagged degeneracy in response-vector-grain
version. Field-level finding (Ahlmann-Eltze 2025 Nat Methods, Diversity by
Design arXiv 2506.22641): mean baseline competitive with deep learning models
is a well-known pattern. Standard fix = cell-level grain per scGen (Lotfollahi
2019 Nat Methods), CPA (Lotfollahi 2023 MSB), CellOT (Bunne 2023). None pair
donors randomly.

CELL-LEVEL EVAL DESIGN:

  Each Lee cell has its own log-normalized HVG expression vector. Mock cells =
  baseline pool. Infected cells (sars_cov_2 or iav) = target pool.

  Baselines predict counterfactual "if this mock cell were infected with v":

  predict_mean (constant per (bucket, virus)):
    predicted = mean(train infected cells | virus=v, bucket=b)
    Same for all test mock cells in (bucket, virus).

  linear_delta (scGen-style δ in gene space):
    δ_v[b] = mean(train_infected | v, b) - mean(train_mock | b)
    predicted = test_mock + δ_v[bucket]
    Varies per cell (each mock has own baseline).

  knn (k=5, cosine):
    For each test mock cell: K nearest INFECTED training cells where virus=v.
    predicted = mean(K nearest infected cells' expression).

  Cross-virus per Issue 15: hold out virus v; train on other viruses. knn +
  linear_delta use OTHER virus' infected pool as their reference.

METRICS:
  Per (bucket, virus, split_type): Pearson r + R² + L1 error between predicted
  cell-mean (over mock query cells) and observed cell-mean (over held-out
  infected cells).

CAVEAT (Ahlmann-Eltze 2025 / Diversity by Design 2025): even with cell-level
grain, mean baseline may stay competitive. Headline differences depend on
within-cell-type signal heterogeneity. Phase 7 GATE 2 calibration may need
interpolated-duplicate positive control + dynamic-range-fraction (open
methodology debate per Oct 2025 rebuttal bioRxiv 2025.10.20.683304).

OUTPUT:
  results/tables/phase5_baselines_eval.csv (wide format per row)
  results/tables/phase5_baselines_eval_long.csv (long format per baseline)
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"

LEE_FILE = PROC / "lee_2020_reannotated.h5ad"
HARMONY_GLOBAL = PROC / "harmony_global_embedding.h5ad"

BUCKETS = ("monocyte", "B", "NK", "CD4T", "CD8T")
VIRUSES = ("sars_cov_2", "iav")
SEED = 42
KNN_K = 5  # Issue 24 default; sensitivity at k=10/25 in Phase 5 supplementary
WITHIN_VIRUS_TRAIN_FRAC = 0.8


def load_lee_cell_level() -> tuple[np.ndarray, pd.DataFrame, list[str]]:
    """Load Lee 2020 → cell × HVG matrix log-normalized on harmony_global 4000 HVG."""
    logger.info("loading harmony_global HVG list...")
    hg = ad.read_h5ad(HARMONY_GLOBAL, backed="r")
    hvg_genes = list(hg.var.index)

    logger.info("loading harmony_global for Lee bucket assignment...")
    hg_lee = hg.obs[hg.obs["study_id"] == "lee_2020"][["coarse"]].copy()
    hg_lee["positional_idx"] = hg_lee.index.str.replace("-lee_2020", "", regex=False).astype(int)
    pos_to_coarse = dict(zip(hg_lee["positional_idx"], hg_lee["coarse"], strict=False))

    logger.info("loading Lee reannotated h5ad...")
    lee = ad.read_h5ad(LEE_FILE)
    lee.obs["coarse"] = pd.Series(
        [pos_to_coarse.get(i) for i in range(lee.shape[0])],
        index=lee.obs.index,
        dtype="object",
    )
    keep_mask = (lee.obs["coarse"].isin(BUCKETS) & lee.obs["virus"].isin(["mock", *VIRUSES])).values
    lee = lee[keep_mask].copy()
    logger.info("Lee in-bucket valid-virus: %d cells", lee.shape[0])

    var_sym = lee.var.index.astype(str)
    common = pd.Index(hvg_genes).intersection(var_sym)
    lee = lee[:, common].copy()

    sc.pp.normalize_total(lee, target_sum=1e4)
    sc.pp.log1p(lee)

    X = lee.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)

    # Reindex to harmony_global HVG order (zero-fill missing)
    n_genes = len(hvg_genes)
    sym_to_idx = {g: i for i, g in enumerate(common)}
    X_aligned = np.zeros((X.shape[0], n_genes), dtype=np.float32)
    for j, g in enumerate(hvg_genes):
        if g in sym_to_idx:
            X_aligned[:, j] = X[:, sym_to_idx[g]]

    obs_df = lee.obs[["donor_id", "coarse", "virus"]].copy().reset_index(drop=True)
    logger.info("cell-level X: %s; obs: %s", X_aligned.shape, obs_df.shape)
    return X_aligned, obs_df, hvg_genes


def pearson_r(pred: np.ndarray, true: np.ndarray) -> float:
    if pred.std() < 1e-12 or true.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(pred, true)[0, 1])


def r_squared(pred: np.ndarray, true: np.ndarray) -> float:
    ss_res = float(np.sum((true - pred) ** 2))
    ss_tot = float(np.sum((true - true.mean()) ** 2))
    if ss_tot < 1e-12:
        return float("nan")
    return 1.0 - ss_res / ss_tot


def l1_error(pred: np.ndarray, true: np.ndarray) -> float:
    return float(np.mean(np.abs(pred - true)))


def predict_mean_pred(
    X_train_inf: np.ndarray,
    virus_train: np.ndarray,
    bucket_train: np.ndarray,
    bucket_test: str,
    virus_test: str,
    n_test: int,
    n_genes: int,
) -> np.ndarray:
    mask = (virus_train == virus_test) & (bucket_train == bucket_test)
    if mask.sum() == 0:
        mask = virus_train == virus_test
    if mask.sum() == 0:
        return np.zeros((n_test, n_genes), dtype=np.float32)
    m = X_train_inf[mask].mean(axis=0)
    return np.broadcast_to(m, (n_test, n_genes)).copy()


def linear_delta_pred(
    X_test_mock: np.ndarray,
    X_train_inf: np.ndarray,
    X_train_mock: np.ndarray,
    virus_train_inf: np.ndarray,
    bucket_train_inf: np.ndarray,
    bucket_train_mock: np.ndarray,
    bucket_test: str,
    virus_test: str,
) -> np.ndarray:
    inf_mask = (virus_train_inf == virus_test) & (bucket_train_inf == bucket_test)
    if inf_mask.sum() == 0:
        inf_mask = virus_train_inf == virus_test
    mock_mask = bucket_train_mock == bucket_test
    if mock_mask.sum() == 0:
        mock_mask = np.ones(len(bucket_train_mock), dtype=bool)
    if inf_mask.sum() == 0 or mock_mask.sum() == 0:
        return X_test_mock.copy()
    delta = X_train_inf[inf_mask].mean(axis=0) - X_train_mock[mock_mask].mean(axis=0)
    return X_test_mock + delta[None, :]


def knn_pred(
    X_test_mock: np.ndarray,
    X_train_inf: np.ndarray,
    virus_train_inf: np.ndarray,
    bucket_train_inf: np.ndarray,
    bucket_test: str,
    virus_test: str,
    k: int = KNN_K,
) -> np.ndarray:
    pool_mask = (virus_train_inf == virus_test) & (bucket_train_inf == bucket_test)
    if pool_mask.sum() < k:
        pool_mask = virus_train_inf == virus_test
    if pool_mask.sum() < k:
        pool_mask = np.ones(len(virus_train_inf), dtype=bool)
    pool = X_train_inf[pool_mask]
    knn = NearestNeighbors(n_neighbors=min(k, len(pool)), metric="cosine").fit(pool)
    _, neighbor_idx = knn.kneighbors(X_test_mock)
    return pool[neighbor_idx].mean(axis=1)


def evaluate(pred_cells: np.ndarray, observed_cells: np.ndarray, name: str) -> dict:
    pred_mean = pred_cells.mean(axis=0)
    obs_mean = observed_cells.mean(axis=0)
    return {
        f"{name}_r_pearson": round(pearson_r(pred_mean, obs_mean), 4),
        f"{name}_r_squared": round(r_squared(pred_mean, obs_mean), 4),
        f"{name}_l1_error": round(l1_error(pred_mean, obs_mean), 4),
    }


def eval_within_virus(X: np.ndarray, obs: pd.DataFrame) -> list[dict]:
    """80/20 random infected-cell split per (bucket, virus)."""
    rng = np.random.default_rng(SEED)
    rows = []
    for bucket in BUCKETS:
        mock_idx = ((obs["virus"] == "mock") & (obs["coarse"] == bucket)).values.nonzero()[0]
        if len(mock_idx) < 5:
            continue
        for virus in VIRUSES:
            inf_idx = ((obs["virus"] == virus) & (obs["coarse"] == bucket)).values.nonzero()[0]
            if len(inf_idx) < 50:
                continue
            shuffled = inf_idx.copy()
            rng.shuffle(shuffled)
            n_train = int(WITHIN_VIRUS_TRAIN_FRAC * len(shuffled))
            train_inf = shuffled[:n_train]
            test_inf = shuffled[n_train:]
            if len(train_inf) < KNN_K or len(test_inf) < 5:
                continue

            X_train_inf = X[train_inf]
            X_train_mock = X[mock_idx]
            X_query_mock = X[mock_idx]  # query = all mock cells in bucket
            X_observed = X[test_inf]
            bucket_train_inf = np.full(len(train_inf), bucket)
            virus_train_inf = np.full(len(train_inf), virus)
            bucket_train_mock = np.full(len(mock_idx), bucket)
            n_test, n_genes = X_query_mock.shape

            pred_pm = predict_mean_pred(
                X_train_inf, virus_train_inf, bucket_train_inf, bucket, virus, n_test, n_genes
            )
            pred_ld = linear_delta_pred(
                X_query_mock,
                X_train_inf,
                X_train_mock,
                virus_train_inf,
                bucket_train_inf,
                bucket_train_mock,
                bucket,
                virus,
            )
            pred_knn = knn_pred(
                X_query_mock, X_train_inf, virus_train_inf, bucket_train_inf, bucket, virus
            )

            row = {
                "split_type": "within_virus",
                "bucket": bucket,
                "virus": virus,
                "n_train_infected": len(train_inf),
                "n_test_infected_observed": len(test_inf),
                "n_mock_query": len(mock_idx),
            }
            row.update(evaluate(pred_pm, X_observed, "predict_mean"))
            row.update(evaluate(pred_ld, X_observed, "linear_delta"))
            row.update(evaluate(pred_knn, X_observed, "knn"))
            rows.append(row)
            logger.info(
                "  within %s %s: pm_r=%.4f ld_r=%.4f knn_r=%.4f (n_inf=%d→%d)",
                bucket,
                virus,
                row["predict_mean_r_pearson"],
                row["linear_delta_r_pearson"],
                row["knn_r_pearson"],
                len(train_inf),
                len(test_inf),
            )
    return rows


def eval_cross_virus(X: np.ndarray, obs: pd.DataFrame) -> list[dict]:
    """Leave-one-virus-out: train on N-1 viruses, predict held virus per Issue 15."""
    rows = []
    for bucket in BUCKETS:
        mock_idx = ((obs["virus"] == "mock") & (obs["coarse"] == bucket)).values.nonzero()[0]
        if len(mock_idx) < 5:
            continue
        X_mock = X[mock_idx]
        bucket_mock = np.full(len(mock_idx), bucket)
        for held_virus in VIRUSES:
            train_viruses = [v for v in VIRUSES if v != held_virus]
            train_inf_idx = (
                obs["virus"].isin(train_viruses) & (obs["coarse"] == bucket)
            ).values.nonzero()[0]
            held_inf_idx = (
                (obs["virus"] == held_virus) & (obs["coarse"] == bucket)
            ).values.nonzero()[0]
            if len(train_inf_idx) < KNN_K or len(held_inf_idx) < 5:
                continue
            X_train_inf = X[train_inf_idx]
            X_observed = X[held_inf_idx]
            bucket_train_inf = np.full(len(train_inf_idx), bucket)
            virus_train_inf = obs.iloc[train_inf_idx]["virus"].values
            n_test, n_genes = X_mock.shape

            # NB: for cross-virus, held_virus absent from train so we use TRAIN virus pool
            # to compute mean / delta / knn. Predicts "shift seen in train" applied to mock.
            train_virus_for_pred = train_viruses[0]  # e.g., predict iav from sars-trained δ
            pred_pm = predict_mean_pred(
                X_train_inf,
                virus_train_inf,
                bucket_train_inf,
                bucket,
                train_virus_for_pred,
                n_test,
                n_genes,
            )
            pred_ld = linear_delta_pred(
                X_mock,
                X_train_inf,
                X_mock,
                virus_train_inf,
                bucket_train_inf,
                bucket_mock,
                bucket,
                train_virus_for_pred,
            )
            pred_knn = knn_pred(
                X_mock,
                X_train_inf,
                virus_train_inf,
                bucket_train_inf,
                bucket,
                train_virus_for_pred,
            )

            row = {
                "split_type": "cross_virus",
                "bucket": bucket,
                "virus": held_virus,
                "n_train_infected": len(train_inf_idx),
                "n_test_infected_observed": len(held_inf_idx),
                "n_mock_query": len(mock_idx),
            }
            row.update(evaluate(pred_pm, X_observed, "predict_mean"))
            row.update(evaluate(pred_ld, X_observed, "linear_delta"))
            row.update(evaluate(pred_knn, X_observed, "knn"))
            rows.append(row)
            logger.info(
                "  cross %s held=%s: pm_r=%.4f ld_r=%.4f knn_r=%.4f (n_train=%d n_held=%d)",
                bucket,
                held_virus,
                row["predict_mean_r_pearson"],
                row["linear_delta_r_pearson"],
                row["knn_r_pearson"],
                len(train_inf_idx),
                len(held_inf_idx),
            )
    return rows


def main() -> int:
    X, obs, _hvg = load_lee_cell_level()
    logger.info("=== within-virus eval ===")
    within_rows = eval_within_virus(X, obs)
    logger.info("=== cross-virus eval ===")
    cross_rows = eval_cross_virus(X, obs)

    df = pd.DataFrame(within_rows + cross_rows)
    out = TABLES / "phase5_baselines_eval.csv"
    df.to_csv(out, index=False)
    logger.info("wrote %s (%d rows)", out.name, len(df))

    long_rows = []
    for _, r in df.iterrows():
        for name in ("predict_mean", "linear_delta", "knn"):
            long_rows.append(
                {
                    "baseline": name,
                    "split_type": r["split_type"],
                    "bucket": r["bucket"],
                    "virus": r["virus"],
                    "n_train_infected": r["n_train_infected"],
                    "n_test_infected_observed": r["n_test_infected_observed"],
                    "r_pearson": r[f"{name}_r_pearson"],
                    "r_squared": r[f"{name}_r_squared"],
                    "l1_error": r[f"{name}_l1_error"],
                }
            )
    long_df = pd.DataFrame(long_rows)
    long_out = TABLES / "phase5_baselines_eval_long.csv"
    long_df.to_csv(long_out, index=False)
    logger.info("wrote %s (%d rows)", long_out.name, len(long_df))

    logger.info("\n=== Phase 5 cell-level headline ===")
    summary = (
        long_df.groupby(["baseline", "split_type"], observed=True)[
            ["r_pearson", "r_squared", "l1_error"]
        ]
        .agg(["mean", "min", "max"])
        .round(4)
    )
    logger.info("\n%s", summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
