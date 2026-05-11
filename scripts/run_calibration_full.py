"""Session 3 full calibration sweep.

Produces per-dataset calibration tables, one row per (bucket, metric):

  results/tables/calibration_phase3.csv               (Phase 3 original labels)
  results/tables/calibration_phase35_low.csv          (Phase 3.5 Low unified)
  results/tables/calibration_phase35_high.csv         (Phase 3.5 High unified, 3 buckets)
  results/tables/calibration_phase35_subbucket.csv    (Phase 3.5 Low sub-bucket)
  results/tables/calibration_phase3_global_harmony.csv (global Harmony, obsm only)

Each row: bucket, metric, observed, perm_p95, perm_p99, perm_p_value,
          split_half_ceiling, sh_ci_low, sh_ci_high, in_split_half_ci_alpha05,
          calibrated_pass_p95_alpha05, calibrated_pass_p99_alpha05,
          summary_mean, summary_median, summary_min, n_perm, n_splits.

MMD-RBF: observed-value-only sensitivity column (no permutation null;
documented limitation for v1).

Memory: loads only the per-bucket gene-space matrix or obsm slice required;
never loads the full 4000-HVG layer of harmony_global_embedding.h5ad.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from trinetravir.eval.calibration import (
    calibrated_gate_verdict,
    permutation_null_with_metric,
    split_half_with_metric,
    summary_stats_off_diag,
)
from trinetravir.eval.metrics import (
    de_jaccard_off_diag,
    mmd_rbf_off_diag,
    pearson_off_diag,
    spearman_off_diag,
)

REPO = Path(__file__).resolve().parents[1]
PROCESSED = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"
CACHE = PROCESSED / "calibration_cache"

# Response-vector metrics (operate on dict[study -> response_vector])
RV_METRICS = {
    "pearson": pearson_off_diag,
    "spearman": spearman_off_diag,
    "de_jaccard_top100": de_jaccard_off_diag,
}

DEFAULT_N_PERM = 1000
DEFAULT_N_SPLITS = 50
DEFAULT_SEED = 42


def load_phase3_per_celltype(bucket: str) -> tuple[np.ndarray, pd.DataFrame]:
    """Load Phase 3 per-cell-type Harmony h5ad for one bucket."""
    p = PROCESSED / f"harmony_per_celltype_{bucket}.h5ad"
    a = ad.read_h5ad(p)
    x = np.asarray(a.obsm["X_harmony_scaled_hvg"])
    obs_df = pd.DataFrame(
        {
            "study_id": a.obs["study_id"].astype(str).values,
            "donor_id": a.obs["donor_id"].astype(str).values,
            "donor_disease_status": a.obs["donor_disease_status"].astype(str).values,
        }
    )
    return x, obs_df


def load_phase35_consolidated(suffix: str) -> tuple[ad.AnnData, str]:
    """Load a Phase 3.5 consolidated h5ad. Returns (adata, bucket_column)."""
    p = PROCESSED / f"harmony_per_celltype_phase35_{suffix}.h5ad"
    if suffix == "subbucket_low":
        p = PROCESSED / "harmony_subbucket_phase35_low.h5ad"
    a = ad.read_h5ad(p)
    return a, a.uns["bucket_column"]


def split_phase35_bucket(
    a: ad.AnnData, bucket: str, bucket_col: str
) -> tuple[np.ndarray, pd.DataFrame]:
    """Slice a consolidated Phase 3.5 h5ad to one bucket's (X, obs).

    The consolidated h5ad stores the bucket label in ``obs['bucket']`` (set
    by the phase3_5_harmony.py writer). The bucket_col argument from uns
    refers to the source column on the per-study reannotated h5ads
    (cell_type_bucket_unified or cell_type_subbucket_unified) — we don't
    use it here.
    """
    key = f"X_corrected_{bucket}"
    if key not in a.uns:
        raise KeyError(f"missing uns key {key}")
    x = np.asarray(a.uns[key])
    mask = a.obs["bucket"].astype(str) == bucket
    obs_df = pd.DataFrame(
        {
            "study_id": a.obs.loc[mask, "study_id"].astype(str).values,
            "donor_id": a.obs.loc[mask, "donor_id"].astype(str).values,
            "donor_disease_status": a.obs.loc[mask, "donor_disease_status"].astype(str).values,
        }
    )
    if len(obs_df) != len(x):
        # Fall back: reconstruct by contiguous bucket block ordering.
        bucket_order = list(a.uns.get("per_bucket_r", {}).keys())
        if bucket in bucket_order:
            idx_so_far = 0
            for b in bucket_order:
                if b == bucket:
                    break
                idx_so_far += int(np.sum(a.obs["bucket"].astype(str) == b))
            n = len(x)
            sub = a.obs.iloc[idx_so_far : idx_so_far + n]
            obs_df = pd.DataFrame(
                {
                    "study_id": sub["study_id"].astype(str).values,
                    "donor_id": sub["donor_id"].astype(str).values,
                    "donor_disease_status": sub["donor_disease_status"].astype(str).values,
                }
            )
    return x, obs_df


def load_global_harmony() -> tuple[np.ndarray, pd.DataFrame]:
    """Load global Harmony embedding (obsm only, NOT full layer)."""
    p = PROCESSED / "harmony_global_embedding.h5ad"
    a = ad.read_h5ad(p, backed="r")
    obs_df = pd.DataFrame(
        {
            "study_id": a.obs["study_id"].astype(str).values,
            "donor_id": a.obs["donor_id"].astype(str).values,
            "donor_disease_status": a.obs["donor_disease_status"].astype(str).values,
            "coarse": a.obs["coarse"].astype(str).values,
        }
    )
    # X_harmony is small (~50 dim); copy and close file.
    x = np.asarray(a.obsm["X_harmony"])
    a.file.close()
    return x, obs_df


def filter_sars_mock(obs_df: pd.DataFrame, x: np.ndarray) -> tuple[np.ndarray, pd.DataFrame]:
    """Match the harmony_per_bucket SARS+mock filter."""
    mask = obs_df["donor_disease_status"].isin(["diseased", "healthy_control"]).values
    return x[mask], obs_df.loc[mask].reset_index(drop=True)


def _observed_value(
    x: np.ndarray, obs: pd.DataFrame, metric_name: str, metric_fn
) -> tuple[float, dict]:
    """Compute observed metric on per-study response vectors. Returns (value, summary_stats)."""
    from trinetravir.eval.calibration import _per_study_donor_status, _response_vector

    per_study = _per_study_donor_status(obs)
    cell_study = obs["study_id"].astype(str).to_numpy()
    cell_donor = obs["donor_id"].astype(str).to_numpy()
    rvs = {}
    for sid in sorted(per_study.keys()):
        ds = per_study[sid]
        donor_to_label = ds.to_dict()
        m = cell_study == sid
        rv = _response_vector(x[m], cell_donor[m], donor_to_label)
        if rv is not None:
            rvs[sid] = rv
    if len(rvs) < 2:
        return float("nan"), {"mean": float("nan"), "median": float("nan"), "min": float("nan")}
    summary = summary_stats_off_diag(rvs)
    return float(metric_fn(rvs)), summary


def _observed_mmd(x: np.ndarray, obs: pd.DataFrame) -> float:
    """MMD-RBF: per-study (diseased-cell) vs (healthy-cell) centered embedding."""
    from trinetravir.eval.calibration import _per_study_donor_status

    per_study = _per_study_donor_status(obs)
    cell_study = obs["study_id"].astype(str).to_numpy()
    cell_donor = obs["donor_id"].astype(str).to_numpy()
    embeddings_per_study: dict[str, np.ndarray] = {}
    for sid in sorted(per_study.keys()):
        ds = per_study[sid].to_dict()
        m = cell_study == sid
        sub_x = x[m]
        sub_donor = cell_donor[m]
        labels = np.array([ds.get(d, "unknown") for d in sub_donor])
        d_mask = labels == "diseased"
        h_mask = labels == "healthy_control"
        if not d_mask.any() or not h_mask.any():
            continue
        # Center each cell on its study's healthy mean -> diseased-distribution residual.
        h_mean = sub_x[h_mask].mean(axis=0)
        embeddings_per_study[sid] = sub_x[d_mask] - h_mean
    if len(embeddings_per_study) < 2:
        return float("nan")
    return float(mmd_rbf_off_diag(embeddings_per_study))


def calibrate_bucket(
    x: np.ndarray,
    obs: pd.DataFrame,
    bucket: str,
    dataset_label: str,
    *,
    n_perm: int = DEFAULT_N_PERM,
    n_splits: int = DEFAULT_N_SPLITS,
    seed: int = DEFAULT_SEED,
) -> list[dict]:
    """Run calibration for one (X, obs, bucket). Returns list of per-metric rows."""
    out_rows = []
    for metric_name, metric_fn in RV_METRICS.items():
        observed, summary = _observed_value(x, obs, metric_name, metric_fn)
        if np.isnan(observed):
            out_rows.append(
                {
                    "bucket": bucket,
                    "metric": metric_name,
                    "observed": float("nan"),
                    "error": "insufficient_data",
                }
            )
            continue
        cache_label = f"{dataset_label}_{metric_name}"
        perm = permutation_null_with_metric(
            x,
            obs,
            bucket,
            metric_fn=metric_fn,
            n_permutations=n_perm,
            seed=seed,
            cache_dir=CACHE,
            cache_label=cache_label,
        )
        sh = split_half_with_metric(
            x,
            obs,
            bucket,
            metric_fn=metric_fn,
            n_splits=n_splits,
            seed=seed,
        )
        v95 = calibrated_gate_verdict(
            observed, perm["null"], sh["split_half_distribution"], percentile=95, alpha=0.05
        )
        v99 = calibrated_gate_verdict(
            observed, perm["null"], sh["split_half_distribution"], percentile=99, alpha=0.05
        )
        out_rows.append(
            {
                "bucket": bucket,
                "metric": metric_name,
                "observed": round(observed, 4),
                "summary_mean": round(summary["mean"], 4),
                "summary_median": round(summary["median"], 4),
                "summary_min": round(summary["min"], 4),
                "perm_p95": round(v95["null_threshold"], 4),
                "perm_p99": round(v99["null_threshold"], 4),
                "perm_p_value": round(v99["p_value"], 4),
                "perm_n_actual": perm["n_actual"],
                "split_half_ceiling": round(sh["ceiling"], 4)
                if not np.isnan(sh["ceiling"])
                else float("nan"),
                "sh_ci_low_alpha05": round(v99["ci_low"], 4)
                if not np.isnan(v99["ci_low"])
                else float("nan"),
                "sh_ci_high_alpha05": round(v99["ci_high"], 4)
                if not np.isnan(v99["ci_high"])
                else float("nan"),
                "in_split_half_ci_alpha05": bool(v99["in_split_half_ci"]),
                "calibrated_pass_p95_alpha05": bool(v95["pass"]),
                "calibrated_pass_p99_alpha05": bool(v99["pass"]),
                "n_perm": n_perm,
                "n_splits": n_splits,
            }
        )
    # MMD observed-only (no permutation; documented limitation)
    mmd_obs = _observed_mmd(x, obs)
    out_rows.append(
        {
            "bucket": bucket,
            "metric": "mmd_rbf",
            "observed": round(mmd_obs, 4) if not np.isnan(mmd_obs) else float("nan"),
            "note": "observed_only_no_perm_null_v1_limitation",
        }
    )
    return out_rows


def run_phase3() -> pd.DataFrame:
    rows = []
    for bucket in ("monocyte", "B", "NK", "CD4T", "CD8T"):
        print(f"[phase3] bucket {bucket}")
        x, obs = load_phase3_per_celltype(bucket)
        rows.extend(calibrate_bucket(x, obs, bucket, "phase3"))
    return pd.DataFrame(rows)


def run_phase35_consolidated(suffix: str, dataset_label: str) -> pd.DataFrame:
    print(f"[{dataset_label}] loading consolidated")
    a, bucket_col = load_phase35_consolidated(suffix)
    buckets = sorted(
        {k.removeprefix("X_corrected_") for k in a.uns if k.startswith("X_corrected_")}
    )
    rows = []
    for bucket in buckets:
        print(f"[{dataset_label}] bucket {bucket}")
        try:
            x, obs = split_phase35_bucket(a, bucket, bucket_col)
        except KeyError as e:
            print(f"  skip: {e}")
            continue
        rows.extend(calibrate_bucket(x, obs, bucket, dataset_label))
    return pd.DataFrame(rows)


def run_global_harmony() -> pd.DataFrame:
    print("[global_harmony] loading obsm only")
    x_all, obs_all = load_global_harmony()
    x_all, obs_all = filter_sars_mock(obs_all, x_all)
    rows = []
    for bucket in ("monocyte", "B", "NK", "CD4T", "CD8T"):
        print(f"[global_harmony] bucket {bucket}")
        mask = (obs_all["coarse"] == bucket).values
        x = x_all[mask]
        obs = obs_all.loc[mask].reset_index(drop=True)
        if len(obs) < 100:
            print(f"  skip bucket {bucket}: too few cells")
            continue
        rows.extend(calibrate_bucket(x, obs, bucket, "global_harmony"))
    return pd.DataFrame(rows)


def main() -> None:
    targets = (
        sys.argv[1:]
        if len(sys.argv) > 1
        else [
            "phase3",
            "phase35_low",
            "phase35_high",
            "phase35_subbucket",
            "global_harmony",
        ]
    )
    TABLES.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    for tag in targets:
        out = TABLES / f"calibration_{tag}.csv"
        if tag == "phase3":
            df = run_phase3()
        elif tag == "phase35_low":
            df = run_phase35_consolidated("low", "phase35_low")
        elif tag == "phase35_high":
            df = run_phase35_consolidated("high", "phase35_high")
        elif tag == "phase35_subbucket":
            df = run_phase35_consolidated("subbucket_low", "phase35_subbucket")
        elif tag == "global_harmony":
            df = run_global_harmony()
            out = TABLES / "calibration_phase3_global_harmony.csv"
        else:
            print(f"unknown target {tag}")
            continue
        df.to_csv(out, index=False)
        print(f"wrote {out}: {len(df)} rows")


if __name__ == "__main__":
    main()
