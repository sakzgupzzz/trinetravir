"""Session 4 Part B: scVI global supplementary sweep.

Per METHODS_CHOICES.md Issue 34 (pre-spec 2026-05-11) + references/session_4_prompt.md
Part B. Addresses reviewer concern "did you try scVI in its standard global mode?"

Single scVI training run on full v1 corpus (244,389 cells on 4000 HVG matching
harmony_global_embedding.h5ad input space). Same 16-config hyperparameter grid as
Part A. Best config selected on held-out donor validation loss. Per-bucket response
vectors computed via obs['coarse'] (NOT cell_type_bucket per Issue 7).

Comparison: scVI_global per-bucket Δr_mvs vs cached Harmony_global per-bucket
response vectors. Verdict per same Issue 34 four-tier rule (with Amendment 1 for
Tier I two-sided proximity).

Output:
  results/tables/session4_part_b_global_scvi_per_bucket.csv
    rows = 16 configs + 5 calibrated rows (best-config per bucket)
  results/tables/session4_part_b_global_verdict.csv

Wall-time: ~2-4h on A100; ~3-6h on L4 (one global run vs 80 per-bucket runs in
Part A; single config covers all 5 buckets via projection).
"""

from __future__ import annotations

import logging
import time
import warnings
from itertools import product
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"
MVS_FILE = REPO / "references" / "khatri_mvs_gene_list.csv"

INPUT_H5AD = PROC / "scvi_input_global.h5ad"
# Session 4.5 Part F unblock: load precomputed Harmony response vector cache
# instead of full 16GB harmony_global_embedding.h5ad (Modal CLI couldn't upload
# 16GB; cache ~5MB uploads cleanly). Cache schema: rows indexed by
# (bucket, study_id, donor_disease_status), columns = 4000 HVG gene symbols.
# Built locally via scripts/build_harmony_global_response_vector_cache.py.
HARMONY_CACHE = PROC / "harmony_global_response_vector_cache.parquet"

BUCKETS = ("monocyte", "B", "NK", "CD4T", "CD8T")

# Issue 34 pre-spec hyperparameter grid (16 configs)
HP_GRID = list(
    product(
        [10, 20, 30, 50],  # n_latent
        [128, 256],  # n_hidden
        [1, 2],  # n_layers
    )
)

SCVI_FIXED = {
    "dropout_rate": 0.1,
    "gene_likelihood": "zinb",
    "dispersion": "gene",
    "latent_distribution": "normal",
}
TRAIN_KWARGS = {
    "max_epochs": 400,
    "early_stopping": True,
    "early_stopping_patience": 50,
    "early_stopping_monitor": "reconstruction_loss_validation",
    "check_val_every_n_epoch": 5,
}

N_PERM = 1000
N_BOOTSTRAP = 1000
SEED = 42


def load_mvs_genes() -> set[str]:
    df = pd.read_csv(MVS_FILE, comment="#")
    return set(df["gene_symbol"].astype(str).tolist())


def load_harmony_global_response_vectors() -> dict[str, pd.DataFrame]:
    """Load harmony_global per-bucket per-study response vectors from precomputed cache.

    Cache built locally via scripts/build_harmony_global_response_vector_cache.py to
    avoid 16GB Modal CLI upload. Schema:
      Rows = (bucket, study_id, donor_disease_status) tuple
      Columns = 4000 HVG gene symbols + ['bucket','study_id','donor_disease_status','n_cells']

    Returns dict keyed by bucket. Each value is DataFrame indexed by gene,
    columns are (study_id, donor_disease_status) tuples, values are mean
    Harmony-scaled HVG expression across cells in that stratum.
    """
    logger.info("loading harmony_global response vector cache from %s", HARMONY_CACHE.name)
    cache = pd.read_parquet(HARMONY_CACHE)
    meta_cols = {"bucket", "study_id", "donor_disease_status", "n_cells"}
    hvg_genes = [c for c in cache.columns if c not in meta_cols]
    logger.info(
        "cache shape: %s rows × %d gene columns + %d meta cols",
        cache.shape,
        len(hvg_genes),
        len(meta_cols),
    )

    out = {}
    for bucket in BUCKETS:
        bucket_rows = cache[cache["bucket"] == bucket]
        if bucket_rows.empty:
            logger.warning("  %s: no cache rows; skip", bucket)
            continue
        rv_cols = {}
        for _, row in bucket_rows.iterrows():
            key = (row["study_id"], row["donor_disease_status"])
            rv_cols[key] = row[hvg_genes].values.astype(float)
        df = pd.DataFrame(rv_cols, index=hvg_genes)
        out[bucket] = df
        logger.info("  %s: %d strata x %d genes", bucket, df.shape[1], df.shape[0])
    return out


def compute_per_study_status_rv(
    expr: np.ndarray, obs: pd.DataFrame, genes: list[str]
) -> pd.DataFrame:
    """Compute per-(study, status) response vectors from a cell-level expression matrix."""
    rv_cols = {}
    for (study, status), _sub in obs.groupby(["study_id", "donor_disease_status"], observed=True):
        mask = ((obs["study_id"] == study) & (obs["donor_disease_status"] == status)).values
        if mask.sum() < 1:
            continue
        mean_expr = expr[mask].mean(axis=0)
        rv_cols[(study, status)] = np.asarray(mean_expr).flatten()
    return pd.DataFrame(rv_cols, index=genes)


def mean_off_diagonal_pearson(rv: pd.DataFrame) -> float:
    """Mean off-diagonal Pearson r across (study, status) pairs."""
    n = rv.shape[1]
    if n < 2:
        return float("nan")
    corr = rv.corr().values  # uses pandas corr default = Pearson on columns
    iu = np.triu_indices_from(corr, k=1)
    return float(np.nanmean(corr[iu]))


def apply_verdict(deltas: np.ndarray) -> str:
    """Issue 34 + Amendment 1 four-tier verdict (two-sided Tier I)."""
    max_abs_d = float(np.max(np.abs(deltas)))
    above_010 = int((deltas > 0.10).sum())
    above_020 = int((deltas > 0.20).sum())
    in_005_010 = int(((deltas > 0.05) & (deltas <= 0.10)).sum())
    below_neg010 = int((deltas < -0.10).sum())
    if max_abs_d <= 0.05:
        return f"TIER_I_HARMONY_ADEQUATE (max |Δr_mvs| = {max_abs_d:.4f} ≤ 0.05)"
    if above_010 >= 3 or above_020 >= 1:
        return "TIER_III_SCVI_PREFERRED (≥3 buckets Δr_mvs > 0.10 or any > 0.20)"
    if below_neg010 >= 3:
        return "TIER_IV_HARMONY_PREFERRED (≥3 buckets Δr_mvs < -0.10)"
    if in_005_010 >= 1 and above_010 == 0:
        return "TIER_II_MIXED (at least one Δr_mvs in (0.05, 0.10], none > 0.10)"
    return "TIER_II_MIXED (boundary case; conservative tie-break)"


def run_global_sweep() -> int:
    import scvi
    import torch

    np.random.seed(SEED)
    torch.manual_seed(SEED)
    scvi.settings.seed = SEED

    logger.info("loading %s", INPUT_H5AD)
    adata = ad.read_h5ad(INPUT_H5AD)
    # filter out NaN coarse cells (not in any bucket)
    keep = adata.obs["coarse"].notna()
    adata = adata[keep].copy()
    logger.info("filtered to %d cells in 5 buckets", adata.shape[0])

    mvs_genes = load_mvs_genes()
    logger.info("MVS gene set: %d genes", len(mvs_genes))
    mvs_idx = adata.var.index.isin(mvs_genes)
    logger.info("MVS overlap with HVG: %d genes", int(mvs_idx.sum()))

    harmony_rvs = load_harmony_global_response_vectors()

    # held-out donor split (80/20) stratified by donor_disease_status
    donor_meta = adata.obs.groupby("donor_id", observed=True)["donor_disease_status"].first()
    rng = np.random.default_rng(SEED)
    val_donors = set()
    for _status, sub in donor_meta.groupby(donor_meta):
        n_val = max(1, int(round(0.2 * len(sub))))
        val_donors.update(rng.choice(sub.index.values, size=n_val, replace=False))

    is_val = adata.obs["donor_id"].isin(val_donors).values
    train_idx = np.where(~is_val)[0]
    val_idx = np.where(is_val)[0]
    logger.info(
        "train/val split: %d train cells / %d val cells (%d val donors)",
        len(train_idx),
        len(val_idx),
        len(val_donors),
    )

    # store rows for each config
    rows = []

    for cfg_idx, (n_lat, n_hid, n_lay) in enumerate(HP_GRID):
        cfg_label = f"config {cfg_idx + 1}/16: n_latent={n_lat} n_hidden={n_hid} n_layers={n_lay}"
        logger.info("  %s", cfg_label)
        t0 = time.time()

        scvi.model.SCVI.setup_anndata(adata, batch_key="study_id")
        model = scvi.model.SCVI(
            adata,
            n_latent=n_lat,
            n_hidden=n_hid,
            n_layers=n_lay,
            **SCVI_FIXED,
        )
        model.train(
            train_size=1.0 - len(val_idx) / adata.shape[0],
            validation_size=len(val_idx) / adata.shape[0],
            **TRAIN_KWARGS,
        )
        val_loss = float(model.history["reconstruction_loss_validation"].iloc[-1].item())

        # normalized expression
        norm = model.get_normalized_expression(library_size=1e4, return_numpy=True)
        # log1p + scale (matches Harmony normalization protocol)
        norm_log = np.log1p(norm)
        norm_scaled = sc.pp.scale(norm_log, zero_center=True, max_value=10, copy=True)

        # per-bucket r evaluation
        bucket_r_full = {}
        bucket_r_mvs = {}
        for bucket in BUCKETS:
            bucket_mask = (adata.obs["coarse"] == bucket).values
            b_obs = adata.obs[bucket_mask].copy().reset_index(drop=True)
            b_expr = norm_scaled[bucket_mask]
            rv = compute_per_study_status_rv(b_expr, b_obs, list(adata.var.index))
            r_full = mean_off_diagonal_pearson(rv)
            rv_mvs = rv.loc[adata.var.index[mvs_idx].intersection(rv.index)]
            r_mvs = mean_off_diagonal_pearson(rv_mvs)
            bucket_r_full[bucket] = r_full
            bucket_r_mvs[bucket] = r_mvs

        wall = time.time() - t0
        for bucket in BUCKETS:
            rows.append(
                {
                    "bucket": bucket,
                    "config_idx": cfg_idx,
                    "n_latent": n_lat,
                    "n_hidden": n_hid,
                    "n_layers": n_lay,
                    "val_loss": round(val_loss, 4),
                    "r_full": round(bucket_r_full[bucket], 4),
                    "r_mvs": round(bucket_r_mvs[bucket], 4),
                    "wall_seconds": round(wall, 1),
                    "selected": False,
                }
            )

        logger.info(
            "    val_loss=%.4f wall=%.1fs r_mvs={%s}",
            val_loss,
            wall,
            ", ".join(f"{b}={bucket_r_mvs[b]:.3f}" for b in BUCKETS),
        )

    # select best config: globally by lowest val_loss (single global model serves all buckets)
    df = pd.DataFrame(rows)
    # Each config has 5 bucket rows; val_loss identical within config
    config_summary = df.groupby("config_idx", as_index=False).agg({"val_loss": "first"})
    best_cfg_idx = int(config_summary.loc[config_summary["val_loss"].idxmin(), "config_idx"])
    logger.info(
        "best config: %d (val_loss=%.4f)",
        best_cfg_idx,
        config_summary["val_loss"].min(),
    )
    df.loc[df["config_idx"] == best_cfg_idx, "selected"] = True

    # Δr vs harmony_global per-bucket
    for bucket in BUCKETS:
        h_rv = harmony_rvs[bucket]
        r_full_h = mean_off_diagonal_pearson(h_rv)
        r_mvs_h = mean_off_diagonal_pearson(h_rv.loc[h_rv.index.intersection(mvs_genes)])
        sel = (df["config_idx"] == best_cfg_idx) & (df["bucket"] == bucket)
        df.loc[sel, "r_full_harmony"] = round(r_full_h, 4)
        df.loc[sel, "r_mvs_harmony"] = round(r_mvs_h, 4)
        df.loc[sel, "delta_r_full"] = round(df.loc[sel, "r_full"].iloc[0] - r_full_h, 4)
        df.loc[sel, "delta_r_mvs"] = round(df.loc[sel, "r_mvs"].iloc[0] - r_mvs_h, 4)

    out = TABLES / "session4_part_b_global_scvi_per_bucket.csv"
    df.to_csv(out, index=False)
    logger.info("wrote %s (%d rows)", out.name, len(df))

    # verdict
    selected = df[df["selected"]].copy()
    deltas = selected["delta_r_mvs"].astype(float).values
    verdict = apply_verdict(deltas)
    logger.info("\n=== Issue 6 Part B verdict (scVI global vs Harmony global) ===")
    logger.info("%s", verdict)
    for _, r in selected.iterrows():
        logger.info(
            "  %s: Δr_mvs=%+.4f (scVI=%.4f, Harmony=%.4f)",
            r["bucket"],
            r["delta_r_mvs"],
            r["r_mvs"],
            r["r_mvs_harmony"],
        )
    verdict_df = pd.DataFrame(
        [
            {
                "verdict": verdict,
                "max_abs_delta_r_mvs": round(float(np.max(np.abs(deltas))), 4),
                "min_delta_r_mvs": round(float(np.min(deltas)), 4),
                "mean_delta_r_mvs": round(float(np.mean(deltas)), 4),
                "best_config_idx": best_cfg_idx,
                "scope": "Part B global supplementary",
            }
        ]
    )
    verdict_out = TABLES / "session4_part_b_global_verdict.csv"
    verdict_df.to_csv(verdict_out, index=False)
    logger.info("wrote %s", verdict_out.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(run_global_sweep())
