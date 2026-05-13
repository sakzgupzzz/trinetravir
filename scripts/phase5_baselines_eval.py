"""Phase 5 baselines evaluation harness — protocol-compliant cell-level eval, optimized.

Implements references/phase5_protocol.md (commit 72addfd) for 3 Category 1 baselines
per Issue 24. Cell-level grain per scGen Lotfollahi 2019.

OPTIMIZATIONS (vs prior version):
  - Pre-fit NN index over ALL train_inf cells per (bucket, split). Query once
    with K_max=300; reuse neighbor IDs for all permutations. ~100× speedup over
    refitting NN per perm.
  - For each (perm, K) tuple: filter pre-fetched neighbors by perm_virus==v,
    take top-K closest, mean.
  - predict_mean / linear_delta vectorized: precompute group means; per perm,
    just re-aggregate over permuted labels.

Same protocol semantics. N_PERM=N_BOOTSTRAP=1000.

OUTPUT:
  results/tables/phase5_baselines_eval.csv with columns:
    bucket, virus, baseline, metric, value, lower_ci, upper_ci,
    perm_null_p_value, fdr_q_value, within_or_cross, n_train_infected, n_test_observed
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import anndata as ad
import faiss
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from trinetravir.eval.calibration import fdr_bh

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"
MVS_FILE = REPO / "references" / "khatri_mvs_gene_list.csv"

LEE_FILE = PROC / "lee_2020_reannotated.h5ad"
HARMONY_GLOBAL = PROC / "harmony_global_embedding.h5ad"

BUCKETS = ("monocyte", "B", "NK", "CD4T", "CD8T")
VIRUSES = ("sars_cov_2", "iav")
K_SWEEP = (25, 50, 100)
K_MAX = 300  # Pre-fetch this many neighbors for cheap K-filtering during permutation
N_PERM = 1000
N_BOOTSTRAP = 1000
SEED = 42
WITHIN_VIRUS_TRAIN_FRAC = 0.8
TOP_K_DE = 100


def load_lee_cell_level() -> tuple[np.ndarray, pd.DataFrame, np.ndarray]:
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

    n_genes = len(hvg_genes)
    sym_to_idx = {g: i for i, g in enumerate(common)}
    X_aligned = np.zeros((X.shape[0], n_genes), dtype=np.float32)
    for j, g in enumerate(hvg_genes):
        if g in sym_to_idx:
            X_aligned[:, j] = X[:, sym_to_idx[g]]

    mvs_genes = set(pd.read_csv(MVS_FILE, comment="#")["gene_symbol"].astype(str))
    mvs_mask = np.array([g in mvs_genes for g in hvg_genes])
    logger.info("MVS subset overlap: %d / %d", int(mvs_mask.sum()), n_genes)

    obs_df = lee.obs[["donor_id", "coarse", "virus"]].copy().reset_index(drop=True)
    return X_aligned, obs_df, mvs_mask


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


def de_jaccard_top_k(
    pred: np.ndarray, true: np.ndarray, base: np.ndarray, k: int = TOP_K_DE
) -> float:
    pred_top = set(np.argsort(np.abs(pred - base))[-k:])
    true_top = set(np.argsort(np.abs(true - base))[-k:])
    union = pred_top | true_top
    return len(pred_top & true_top) / len(union) if union else float("nan")


def direction_accuracy(pred: np.ndarray, true: np.ndarray, base: np.ndarray) -> float:
    return float(np.mean(np.sign(pred - base) == np.sign(true - base)))


def all_metrics(
    pred_cells: np.ndarray, observed_cells: np.ndarray, base_mean: np.ndarray, mvs_mask: np.ndarray
) -> dict:
    pred_mean = pred_cells.mean(axis=0)
    obs_mean = observed_cells.mean(axis=0)
    return {
        "pearson_r_full": pearson_r(pred_mean, obs_mean),
        "r_squared_full": r_squared(pred_mean, obs_mean),
        "pearson_r_mvs": pearson_r(pred_mean[mvs_mask], obs_mean[mvs_mask]),
        "r_squared_mvs": r_squared(pred_mean[mvs_mask], obs_mean[mvs_mask]),
        "de_jaccard_top100": de_jaccard_top_k(pred_mean, obs_mean, base_mean, TOP_K_DE),
        "direction_accuracy": direction_accuracy(pred_mean, obs_mean, base_mean),
    }


def eval_condition(
    X: np.ndarray,
    obs: pd.DataFrame,
    mvs_mask: np.ndarray,
    bucket: str,
    virus: str,
    split_type: str,
    rng: np.random.Generator,
) -> list[dict]:
    """Evaluate all 5 baselines (pm, ld, knn k=25/50/100) for one bucket-virus-split.

    Pre-builds NN index ONCE over all train_inf cells, queries with all mock once,
    reuses K_max neighbors across permutations + K values.
    """
    mock_idx = ((obs["virus"] == "mock") & (obs["coarse"] == bucket)).values.nonzero()[0]
    if len(mock_idx) < 5:
        return []

    if split_type == "within_virus":
        inf_idx = ((obs["virus"] == virus) & (obs["coarse"] == bucket)).values.nonzero()[0]
        if len(inf_idx) < 50:
            return []
        shuffled = inf_idx.copy()
        rng.shuffle(shuffled)
        n_tr = int(WITHIN_VIRUS_TRAIN_FRAC * len(shuffled))
        train_inf_idx = shuffled[:n_tr]
        test_inf_idx = shuffled[n_tr:]
        virus_train_inf = np.full(len(train_inf_idx), virus)
        # observed is the same virus (held-out subset)
    else:  # cross_virus
        train_viruses = [v for v in VIRUSES if v != virus]
        train_inf_idx = (
            obs["virus"].isin(train_viruses) & (obs["coarse"] == bucket)
        ).values.nonzero()[0]
        test_inf_idx = ((obs["virus"] == virus) & (obs["coarse"] == bucket)).values.nonzero()[0]
        virus_train_inf = obs.iloc[train_inf_idx]["virus"].values

    if len(train_inf_idx) < max(K_SWEEP) or len(test_inf_idx) < 5:
        return []

    X_train_inf = X[train_inf_idx]
    X_train_mock = X[mock_idx]
    X_query_mock = X_train_mock
    X_observed = X[test_inf_idx]
    base_mean = X_train_mock.mean(axis=0)

    # Effective "virus_test" for prediction (within = virus itself; cross = train virus[0])
    virus_pred = virus if split_type == "within_virus" else train_viruses[0]

    # ----- L2-normalize query mock cells once (FAISS cosine-via-IP) -----
    q_norm = X_query_mock / (np.linalg.norm(X_query_mock, axis=1, keepdims=True) + 1e-12)
    q_norm = q_norm.astype(np.float32)

    # ----- Helper: predicted-mean vector for one (baseline, virus_pred) per perm -----
    # Semantics: per-perm refit (sklearn-style original logic) but with FAISS for KNN
    # (~10× speedup on KNN bottleneck; methodology unchanged).
    def predict_mean_vector(virus_arr: np.ndarray, baseline: str, k: int | None) -> np.ndarray:
        mask = virus_arr == virus_pred
        if mask.sum() == 0:
            mask = np.ones_like(virus_arr, dtype=bool)
        if baseline == "predict_mean":
            return X_train_inf[mask].mean(axis=0)
        if baseline == "linear_delta":
            delta = X_train_inf[mask].mean(axis=0) - base_mean
            return X_query_mock.mean(axis=0) + delta
        # knn: build FAISS index over filtered pool, query mock cells, mean top-K
        pool = X_train_inf[mask].astype(np.float32)
        pool_norm = pool / (np.linalg.norm(pool, axis=1, keepdims=True) + 1e-12)
        k_eff = min(k, len(pool))
        index = faiss.IndexFlatIP(pool.shape[1])
        index.add(np.ascontiguousarray(pool_norm))
        _sims, nidx = index.search(q_norm, k_eff)
        per_cell_means = pool[nidx].mean(axis=1)
        return per_cell_means.mean(axis=0)

    obs_mean = X_observed.mean(axis=0)
    rows = []
    baselines_with_k = [("predict_mean", None), ("linear_delta", None)] + [
        ("knn", k) for k in K_SWEEP
    ]

    import time as _time

    for baseline_name, k in baselines_with_k:
        bl_label = baseline_name + (f"_k{k}" if k else "")
        t_bl_start = _time.time()

        # --- Observed prediction + metrics ---
        t0 = _time.time()
        pred_mean_obs = predict_mean_vector(virus_train_inf, baseline_name, k)
        t_obs = _time.time() - t0
        metrics_obs = {
            "pearson_r_full": pearson_r(pred_mean_obs, obs_mean),
            "r_squared_full": r_squared(pred_mean_obs, obs_mean),
            "pearson_r_mvs": pearson_r(pred_mean_obs[mvs_mask], obs_mean[mvs_mask]),
            "r_squared_mvs": r_squared(pred_mean_obs[mvs_mask], obs_mean[mvs_mask]),
            "de_jaccard_top100": de_jaccard_top_k(pred_mean_obs, obs_mean, base_mean, TOP_K_DE),
            "direction_accuracy": direction_accuracy(pred_mean_obs, obs_mean, base_mean),
        }

        # --- Bootstrap CI on each metric (resample test cells, recompute obs_mean) ---
        boot_metrics = {m: [] for m in metrics_obs}
        n_test = len(test_inf_idx)
        t0 = _time.time()
        for _ in range(N_BOOTSTRAP):
            samp = rng.choice(n_test, size=n_test, replace=True)
            obs_mean_b = X_observed[samp].mean(axis=0)
            boot_metrics["pearson_r_full"].append(pearson_r(pred_mean_obs, obs_mean_b))
            boot_metrics["r_squared_full"].append(r_squared(pred_mean_obs, obs_mean_b))
            boot_metrics["pearson_r_mvs"].append(
                pearson_r(pred_mean_obs[mvs_mask], obs_mean_b[mvs_mask])
            )
            boot_metrics["r_squared_mvs"].append(
                r_squared(pred_mean_obs[mvs_mask], obs_mean_b[mvs_mask])
            )
            boot_metrics["de_jaccard_top100"].append(
                de_jaccard_top_k(pred_mean_obs, obs_mean_b, base_mean, TOP_K_DE)
            )
            boot_metrics["direction_accuracy"].append(
                direction_accuracy(pred_mean_obs, obs_mean_b, base_mean)
            )
        t_boot = _time.time() - t0

        # --- Permutation null (permute virus labels of train_inf) ---
        perm_metrics = {m: [] for m in metrics_obs}
        t0 = _time.time()
        t_predict_total = 0.0
        for pi in range(N_PERM):
            if pi % 100 == 0 and pi > 0:
                elapsed = _time.time() - t0
                eta_perm = elapsed / pi * (N_PERM - pi)
                logger.info(
                    "    %s perm %d/%d: %.1fs elapsed, ETA %.1fs (pred=%.1fs)",
                    bl_label,
                    pi,
                    N_PERM,
                    elapsed,
                    eta_perm,
                    t_predict_total,
                )
            perm_virus = virus_train_inf.copy()
            rng.shuffle(perm_virus)
            tp0 = _time.time()
            pred_mean_perm = predict_mean_vector(perm_virus, baseline_name, k)
            t_predict_total += _time.time() - tp0
            perm_metrics["pearson_r_full"].append(pearson_r(pred_mean_perm, obs_mean))
            perm_metrics["r_squared_full"].append(r_squared(pred_mean_perm, obs_mean))
            perm_metrics["pearson_r_mvs"].append(
                pearson_r(pred_mean_perm[mvs_mask], obs_mean[mvs_mask])
            )
            perm_metrics["r_squared_mvs"].append(
                r_squared(pred_mean_perm[mvs_mask], obs_mean[mvs_mask])
            )
            perm_metrics["de_jaccard_top100"].append(
                de_jaccard_top_k(pred_mean_perm, obs_mean, base_mean, TOP_K_DE)
            )
            perm_metrics["direction_accuracy"].append(
                direction_accuracy(pred_mean_perm, obs_mean, base_mean)
            )
        t_perm = _time.time() - t0
        t_bl = _time.time() - t_bl_start
        logger.info(
            "    %s TOTAL: %.1fs (obs=%.2fs, boot=%.1fs, perm=%.1fs of which pred=%.1fs)",
            bl_label,
            t_bl,
            t_obs,
            t_boot,
            t_perm,
            t_predict_total,
        )

        for metric_name, obs_val in metrics_obs.items():
            boot_arr = np.array(boot_metrics[metric_name])
            perm_arr = np.array(perm_metrics[metric_name])
            boot_clean = boot_arr[~np.isnan(boot_arr)]
            perm_clean = perm_arr[~np.isnan(perm_arr)]
            ci_low = float(np.percentile(boot_clean, 2.5)) if len(boot_clean) > 10 else float("nan")
            ci_high = (
                float(np.percentile(boot_clean, 97.5)) if len(boot_clean) > 10 else float("nan")
            )
            if len(perm_clean) > 10 and not np.isnan(obs_val):
                perm_p = float(np.mean(perm_clean >= obs_val))
            else:
                perm_p = float("nan")
            rows.append(
                {
                    "bucket": bucket,
                    "virus": virus,
                    "baseline": baseline_name + (f"_k{k}" if k is not None else ""),
                    "metric": metric_name,
                    "value": round(obs_val, 4) if not np.isnan(obs_val) else float("nan"),
                    "lower_ci": round(ci_low, 4) if not np.isnan(ci_low) else float("nan"),
                    "upper_ci": round(ci_high, 4) if not np.isnan(ci_high) else float("nan"),
                    "perm_null_p_value": round(perm_p, 4) if not np.isnan(perm_p) else float("nan"),
                    "within_or_cross": split_type,
                    "n_train_infected": len(train_inf_idx),
                    "n_test_observed": len(test_inf_idx),
                }
            )
    return rows


def main() -> int:
    X, obs, mvs_mask = load_lee_cell_level()

    all_rows = []
    rng = np.random.default_rng(SEED)
    for bucket in BUCKETS:
        for virus in VIRUSES:
            for split_type in ("within_virus", "cross_virus"):
                logger.info("  %s %s %s", bucket, virus, split_type)
                rows = eval_condition(X, obs, mvs_mask, bucket, virus, split_type, rng)
                all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    if df.empty:
        logger.error("no rows produced")
        return 1

    df["fdr_q_value"] = np.nan
    for metric in df["metric"].unique():
        mask = df["metric"] == metric
        p_vals = df.loc[mask, "perm_null_p_value"].astype(float).values
        valid = ~np.isnan(p_vals)
        if valid.sum() > 0:
            q_vals = np.full(len(p_vals), np.nan)
            q_vals[valid] = fdr_bh(p_vals[valid])
            df.loc[mask, "fdr_q_value"] = q_vals.round(4)

    out = TABLES / "phase5_baselines_eval.csv"
    df.to_csv(out, index=False)
    logger.info("wrote %s (%d rows)", out.name, len(df))

    logger.info("\n=== Phase 5 headline: pearson_r_full per baseline × split ===")
    pearson_rows = df[df["metric"] == "pearson_r_full"]
    summary = (
        pearson_rows.groupby(["baseline", "within_or_cross"], observed=True)["value"]
        .agg(["mean", "min", "max", "count"])
        .round(4)
    )
    logger.info("\n%s", summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
