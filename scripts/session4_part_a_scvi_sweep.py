"""Session 4 Part A: scVI per-bucket sensitivity sweep.

Per METHODS_CHOICES.md Issue 34 (pre-spec 2026-05-11).

For each v1 bucket (monocyte, B, NK, CD4T, CD8T):
  1. Load data/processed/scvi_input_<bucket>.h5ad (raw counts on 4000 HVG).
  2. Sweep 16 hyperparameter configurations (n_latent × n_hidden × n_layers).
  3. Per config: train scVI, get_normalized_expression(library_size=1e4),
     log1p + scale, compute per-(study, status) response vectors,
     mean off-diagonal cross-study Pearson r (full HVG + MVS subset).
  4. Select best config per bucket by max MVS Pearson r.
  5. Calibration framework v2 on best config:
     - permutation null N=1000 (donor-label shuffle stratified by study)
     - bootstrap CI N=1000 on observed r (donor resample with replacement)
  6. Δr (scVI minus per-bucket Harmony) on full + MVS metrics.

Output per bucket:
  results/tables/session4_scvi_per_bucket_<bucket>.csv
    rows = 16 configs + 1 best-config calibrated row
    cols = config params, r_full, r_mvs, delta_r_full, delta_r_mvs,
           perm_p_mvs, ci_low_mvs, ci_high_mvs, selected

Combined output:
  results/tables/session4_scvi_per_bucket_combined.csv (5 best-config rows)
  results/tables/session4_issue6_verdict.csv (Tier I/II/III/IV per Issue 34)

Wall-time: ~6.5-7h on A100; ~10-12h on A10G g5.xlarge.

Reproducibility:
  random_state=42 (numpy + torch + scvi)
  max_epochs=400, patience=50, monitor='reconstruction_loss_validation'
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
    """Khatri MVS canonical 86-gene list per Issue 18."""
    df = pd.read_csv(MVS_FILE, comment="#")
    return set(df["gene_symbol"].astype(str).tolist())


def load_harmony_response_vectors() -> dict[str, pd.Series]:
    """Cached per-bucket Harmony response vectors (gene-indexed, training_consensus = mean across studies)."""
    out = {}
    for bucket in BUCKETS:
        p = PROC / f"phase3_response_vectors_{bucket}.parquet"
        if not p.exists():
            logger.warning("missing %s", p)
            continue
        df = pd.read_parquet(p)
        out[bucket] = df.mean(axis=1)
    return out


def train_scvi(adata: ad.AnnData, n_latent: int, n_hidden: int, n_layers: int):
    """Set up + train scVI on the bucket's raw counts with batch_key=study_id."""
    import scvi
    import torch

    scvi.settings.seed = SEED
    torch.manual_seed(SEED)
    np.random.seed(SEED)

    scvi.model.SCVI.setup_anndata(adata, batch_key="study_id")
    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_hidden=n_hidden,
        n_layers=n_layers,
        **SCVI_FIXED,
    )
    model.train(**TRAIN_KWARGS)
    return model


def scvi_normalized_log_scaled(model, adata: ad.AnnData) -> np.ndarray:
    """A.3 recipe: get_normalized_expression(library_size=1e4) + log1p + sc.pp.scale."""
    scvi_norm = model.get_normalized_expression(adata, library_size=1e4, return_numpy=True)
    scvi_log = np.log1p(scvi_norm)
    a_scaled = ad.AnnData(X=scvi_log.astype(np.float32), obs=adata.obs.copy(), var=adata.var.copy())
    sc.pp.scale(a_scaled, zero_center=True, max_value=10)
    return np.asarray(a_scaled.X)


def per_study_status_response_vectors(
    X_scaled: np.ndarray, obs: pd.DataFrame
) -> dict[str, dict[str, np.ndarray]]:
    """Compute mean expression per (study, donor_disease_status) cell group.

    Uses boolean masks (positional) rather than groupby().groups indices (label-based)
    to avoid IndexError when obs.index is non-integer (e.g., cell_id strings).
    """
    out: dict[str, dict[str, np.ndarray]] = {}
    study_arr = obs["study_id"].astype(str).values
    status_arr = obs["donor_disease_status"].astype(str).values
    for study in pd.unique(study_arr):
        for status in pd.unique(status_arr):
            mask = (study_arr == study) & (status_arr == status)
            if not mask.any():
                continue
            out.setdefault(study, {})[status] = X_scaled[mask].mean(axis=0)
    return out


def per_study_response_diff(
    rvs: dict[str, dict[str, np.ndarray]],
) -> dict[str, np.ndarray]:
    """Per study: diseased - healthy_control response vector (gene-space)."""
    out = {}
    for study, by_status in rvs.items():
        d = by_status.get("diseased")
        h = by_status.get("healthy_control")
        if d is None or h is None:
            continue
        out[study] = d - h
    return out


def mean_off_diag_pearson(
    rv_dict: dict[str, np.ndarray], gene_names: list[str], restrict: set | None = None
) -> tuple[float, int]:
    """Mean off-diagonal Pearson r across studies (genes × studies matrix)."""
    if restrict is not None:
        mask = np.array([g in restrict for g in gene_names])
        n_genes = int(mask.sum())
        if n_genes < 10:
            return float("nan"), n_genes
        mat = np.column_stack([rv[mask] for rv in rv_dict.values()])
    else:
        mat = np.column_stack(list(rv_dict.values()))
        n_genes = mat.shape[0]
    if mat.shape[1] < 2:
        return float("nan"), n_genes
    corr = np.corrcoef(mat.T)
    n_studies = corr.shape[0]
    off_diag = corr[~np.eye(n_studies, dtype=bool)]
    return float(np.nanmean(off_diag)), n_genes


def permutation_null(
    X_scaled: np.ndarray,
    obs: pd.DataFrame,
    gene_names: list[str],
    restrict: set | None,
    n_perm: int = N_PERM,
    seed: int = SEED,
) -> np.ndarray:
    """Donor-level shuffle of disease label within each study; compute null r distribution."""
    rng = np.random.default_rng(seed)
    null_rs = []
    obs_arr = obs[["study_id", "donor_id", "donor_disease_status"]].copy()
    studies = obs_arr["study_id"].unique()
    # Per-study donor-status mapping
    donor_status_per_study = {}
    for s in studies:
        ds = obs_arr[obs_arr["study_id"] == s][
            ["donor_id", "donor_disease_status"]
        ].drop_duplicates()
        donor_status_per_study[s] = ds.reset_index(drop=True)

    obs_idx = obs.reset_index(drop=True)
    for _ in range(n_perm):
        # Permute donor → status within each study
        perm_status_map = {}
        for s, ds in donor_status_per_study.items():
            perm_labels = rng.permutation(ds["donor_disease_status"].values)
            for d, ls in zip(ds["donor_id"].values, perm_labels, strict=False):
                perm_status_map[(s, d)] = ls
        # Apply to cells
        perm_labels = np.array(
            [
                perm_status_map[(s, d)]
                for s, d in zip(obs_idx["study_id"], obs_idx["donor_id"], strict=False)
            ]
        )
        # Recompute response vectors per (study, status_perm)
        rvs: dict[str, dict[str, np.ndarray]] = {}
        for s in studies:
            for ds_label in ("diseased", "healthy_control"):
                mask = (obs_idx["study_id"].values == s) & (perm_labels == ds_label)
                if mask.sum() == 0:
                    continue
                rvs.setdefault(s, {})[ds_label] = X_scaled[mask].mean(axis=0)
        diffs = per_study_response_diff(rvs)
        if len(diffs) < 2:
            continue
        r, _ = mean_off_diag_pearson(diffs, gene_names, restrict=restrict)
        if not np.isnan(r):
            null_rs.append(r)
    return np.asarray(null_rs)


def bootstrap_ci(
    X_scaled: np.ndarray,
    obs: pd.DataFrame,
    gene_names: list[str],
    restrict: set | None,
    n_boot: int = N_BOOTSTRAP,
    seed: int = SEED + 1,
) -> tuple[float, float]:
    """Bootstrap CI on observed r: resample donors within each study with replacement."""
    rng = np.random.default_rng(seed)
    obs_idx = obs.reset_index(drop=True)
    studies = obs_idx["study_id"].unique()
    boot_rs = []
    for _ in range(n_boot):
        sampled_idx = []
        for s in studies:
            study_mask = obs_idx["study_id"].values == s
            study_donors = obs_idx[study_mask]["donor_id"].unique()
            sampled_donors = rng.choice(study_donors, size=len(study_donors), replace=True)
            for d in sampled_donors:
                donor_idx = np.where(study_mask & (obs_idx["donor_id"].values == d))[0]
                sampled_idx.extend(donor_idx.tolist())
        sub_obs = obs_idx.iloc[sampled_idx]
        sub_X = X_scaled[sampled_idx]
        rvs = per_study_status_response_vectors(sub_X, sub_obs)
        diffs = per_study_response_diff(rvs)
        if len(diffs) < 2:
            continue
        r, _ = mean_off_diag_pearson(diffs, gene_names, restrict=restrict)
        if not np.isnan(r):
            boot_rs.append(r)
    if not boot_rs:
        return float("nan"), float("nan")
    return float(np.percentile(boot_rs, 2.5)), float(np.percentile(boot_rs, 97.5))


def run_bucket(bucket: str, mvs_genes: set[str], harmony_rvs: dict[str, pd.Series]) -> pd.DataFrame:
    """Sweep 16 configs on a bucket; select best by MVS r; calibrate best."""
    h5 = PROC / f"scvi_input_{bucket}.h5ad"
    if not h5.exists():
        raise FileNotFoundError(h5)
    logger.info("=== %s ===  loading %s", bucket, h5.name)
    adata = ad.read_h5ad(h5)
    gene_names = list(adata.var_names)
    logger.info(
        "  %d cells × %d HVG; %d studies; donors=%d",
        adata.n_obs,
        adata.n_vars,
        adata.obs["study_id"].nunique(),
        adata.obs["donor_id"].nunique(),
    )

    sweep_rows = []
    for cfg_idx, (n_lat, n_hid, n_lay) in enumerate(HP_GRID):
        t0 = time.time()
        logger.info(
            "  config %d/16: n_latent=%d n_hidden=%d n_layers=%d", cfg_idx + 1, n_lat, n_hid, n_lay
        )
        model = train_scvi(adata, n_lat, n_hid, n_lay)
        X_scaled = scvi_normalized_log_scaled(model, adata)
        rvs = per_study_status_response_vectors(X_scaled, adata.obs)
        diffs = per_study_response_diff(rvs)
        r_full, n_full = mean_off_diag_pearson(diffs, gene_names)
        r_mvs, n_mvs = mean_off_diag_pearson(diffs, gene_names, restrict=mvs_genes)
        wall = time.time() - t0
        sweep_rows.append(
            {
                "bucket": bucket,
                "config_idx": cfg_idx,
                "n_latent": n_lat,
                "n_hidden": n_hid,
                "n_layers": n_lay,
                "r_full": round(r_full, 4),
                "r_mvs": round(r_mvs, 4),
                "n_genes_full": n_full,
                "n_genes_mvs": n_mvs,
                "wall_seconds": round(wall, 1),
                "selected": False,
            }
        )
        logger.info("    r_full=%.4f r_mvs=%.4f (%.1fs)", r_full, r_mvs, wall)
        del model

    df = pd.DataFrame(sweep_rows)

    # Selection: max r_mvs (primary metric per Issue 34)
    best_idx = df["r_mvs"].astype(float).idxmax()
    df.loc[best_idx, "selected"] = True
    best = df.loc[best_idx]
    logger.info(
        "  selected best config: idx=%d n_latent=%d n_hidden=%d n_layers=%d r_mvs=%.4f",
        best["config_idx"],
        best["n_latent"],
        best["n_hidden"],
        best["n_layers"],
        best["r_mvs"],
    )

    # Re-train best config; full calibration
    model = train_scvi(adata, int(best["n_latent"]), int(best["n_hidden"]), int(best["n_layers"]))
    X_scaled = scvi_normalized_log_scaled(model, adata)
    rvs = per_study_status_response_vectors(X_scaled, adata.obs)
    diffs = per_study_response_diff(rvs)
    r_full_obs, _ = mean_off_diag_pearson(diffs, gene_names)
    r_mvs_obs, _ = mean_off_diag_pearson(diffs, gene_names, restrict=mvs_genes)

    # Δr vs Harmony reference: project Harmony rv onto scVI gene space, compute mean off-diag
    # Harmony rv is gene-indexed (Series); the cached file is per-study × bucket, mean'd.
    # Compute Harmony per-study off-diag r the same way for fair comparison.
    harmony_df = pd.read_parquet(PROC / f"phase3_response_vectors_{bucket}.parquet")
    harmony_per_study = {s: harmony_df[s].reindex(gene_names).values for s in harmony_df.columns}
    # Drop NaN rows for fair comparison
    harmony_mat = np.column_stack(list(harmony_per_study.values()))
    valid_mask = ~np.isnan(harmony_mat).any(axis=1)
    h_gene_names = [g for g, v in zip(gene_names, valid_mask, strict=False) if v]
    harmony_filtered = {s: harmony_per_study[s][valid_mask] for s in harmony_per_study}
    r_full_harmony, _ = mean_off_diag_pearson(harmony_filtered, h_gene_names)
    r_mvs_harmony, _ = mean_off_diag_pearson(harmony_filtered, h_gene_names, restrict=mvs_genes)
    delta_r_full = r_full_obs - r_full_harmony
    delta_r_mvs = r_mvs_obs - r_mvs_harmony

    # Calibration: perm null + bootstrap CI on MVS
    logger.info("  running %d permutations + %d bootstrap (best config)...", N_PERM, N_BOOTSTRAP)
    null_r_mvs = permutation_null(X_scaled, adata.obs, gene_names, restrict=mvs_genes)
    p_mvs = (
        float(((null_r_mvs >= r_mvs_obs).sum() + 1) / (len(null_r_mvs) + 1))
        if len(null_r_mvs)
        else float("nan")
    )
    ci_low_mvs, ci_high_mvs = bootstrap_ci(X_scaled, adata.obs, gene_names, restrict=mvs_genes)

    # Append calibrated best-config row
    calib_row = {
        "bucket": bucket,
        "config_idx": int(best["config_idx"]),
        "n_latent": int(best["n_latent"]),
        "n_hidden": int(best["n_hidden"]),
        "n_layers": int(best["n_layers"]),
        "r_full": round(r_full_obs, 4),
        "r_mvs": round(r_mvs_obs, 4),
        "r_full_harmony": round(r_full_harmony, 4),
        "r_mvs_harmony": round(r_mvs_harmony, 4),
        "delta_r_full": round(delta_r_full, 4),
        "delta_r_mvs": round(delta_r_mvs, 4),
        "perm_p_mvs": round(p_mvs, 4),
        "ci_low_mvs": round(ci_low_mvs, 4) if not np.isnan(ci_low_mvs) else float("nan"),
        "ci_high_mvs": round(ci_high_mvs, 4) if not np.isnan(ci_high_mvs) else float("nan"),
        "n_perm": len(null_r_mvs),
        "selected": True,
        "calibrated": True,
    }
    df["calibrated"] = False
    df = pd.concat([df, pd.DataFrame([calib_row])], ignore_index=True)
    logger.info(
        "  CALIBRATED: r_mvs=%.4f Δr_mvs=%.4f p=%.4f CI=[%.4f, %.4f]",
        r_mvs_obs,
        delta_r_mvs,
        p_mvs,
        ci_low_mvs,
        ci_high_mvs,
    )
    return df


def apply_verdict(combined: pd.DataFrame) -> str:
    """Issue 34 four-tier verdict from per-bucket Δr_mvs."""
    deltas = combined["delta_r_mvs"].astype(float).values
    max_d = float(np.max(deltas))
    if max_d <= 0.05:
        return f"TIER_I_HARMONY_ADEQUATE (max Δr_mvs = {max_d:.4f} ≤ 0.05)"
    above_010 = int((deltas > 0.10).sum())
    above_020 = int((deltas > 0.20).sum())
    in_005_010 = int(((deltas > 0.05) & (deltas <= 0.10)).sum())
    below_neg010 = int((deltas < -0.10).sum())
    if above_010 >= 3 or above_020 >= 1:
        return "TIER_III_SCVI_PREFERRED (≥3 buckets Δr_mvs > 0.10 or any > 0.20)"
    if below_neg010 >= 3:
        return "TIER_IV_HARMONY_PREFERRED (≥3 buckets Δr_mvs < -0.10)"
    if in_005_010 >= 1 and above_010 == 0:
        return "TIER_II_MIXED (at least one Δr_mvs in (0.05, 0.10], none > 0.10)"
    # Boundary case → conservative
    return "TIER_II_MIXED (boundary case; conservative tie-break)"


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    mvs_genes = load_mvs_genes()
    logger.info("MVS gene set: %d genes", len(mvs_genes))
    harmony_rvs = load_harmony_response_vectors()
    logger.info("Harmony reference rvs loaded for: %s", list(harmony_rvs.keys()))

    all_dfs = []
    for bucket in BUCKETS:
        try:
            df = run_bucket(bucket, mvs_genes, harmony_rvs)
        except Exception as e:
            logger.error("%s failed: %s", bucket, e)
            raise
        out = TABLES / f"session4_scvi_per_bucket_{bucket}.csv"
        df.to_csv(out, index=False)
        logger.info("wrote %s: %d rows", out.name, len(df))
        all_dfs.append(df)

    # Combined: just the calibrated best-config rows
    combined = pd.concat([df[df["calibrated"]] for df in all_dfs], ignore_index=True)
    combined_out = TABLES / "session4_scvi_per_bucket_combined.csv"
    combined.to_csv(combined_out, index=False)
    logger.info("wrote %s: %d rows", combined_out.name, len(combined))

    verdict = apply_verdict(combined)
    logger.info("\n=== Issue 6 verdict (Part A primary; per-bucket scVI vs per-bucket Harmony) ===")
    logger.info("%s", verdict)
    for _, r in combined.iterrows():
        logger.info(
            "  %s: Δr_mvs=%+.4f (scVI=%.4f, Harmony=%.4f, p=%.4f, CI=[%.4f, %.4f])",
            r["bucket"],
            r["delta_r_mvs"],
            r["r_mvs"],
            r["r_mvs_harmony"],
            r["perm_p_mvs"],
            r["ci_low_mvs"],
            r["ci_high_mvs"],
        )
    verdict_df = pd.DataFrame(
        [
            {
                "verdict": verdict,
                "max_delta_r_mvs": round(float(combined["delta_r_mvs"].max()), 4),
                "min_delta_r_mvs": round(float(combined["delta_r_mvs"].min()), 4),
                "mean_delta_r_mvs": round(float(combined["delta_r_mvs"].mean()), 4),
            }
        ]
    )
    verdict_out = TABLES / "session4_issue6_verdict.csv"
    verdict_df.to_csv(verdict_out, index=False)
    logger.info("wrote %s", verdict_out.name)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
