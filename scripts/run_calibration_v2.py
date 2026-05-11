"""Session 5 v2 calibration: re-run verdicts with corrected bootstrap CI direction
(Part A1) + add observed-r bootstrap CI (Part A2) + add FDR-BH correction (Part A3).

Reads existing permutation null caches from data/processed/calibration_cache/ so
no perm re-compute. Writes _v2 tables alongside the v1 tables (v1 preserved for
audit trail per Issue 17 atomic schema rule).

Output:
  results/tables/calibration_phase3_v2.csv
  results/tables/calibration_phase3_global_harmony_v2.csv
  results/tables/calibration_phase35_high_v2.csv
  results/tables/calibration_phase35_low_v2.csv
  results/tables/calibration_phase35_subbucket_v2.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from trinetravir.eval.calibration import (
    _per_study_donor_status,
    _response_vector,
    bootstrap_observed_r,
    calibrated_gate_verdict,
    fdr_bh,
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

RV_METRICS = {
    "pearson": pearson_off_diag,
    "spearman": spearman_off_diag,
    "de_jaccard_top100": de_jaccard_off_diag,
}

N_PERM = 1000
N_SPLITS = 50
N_BOOTSTRAP = 200  # reduced from 1000 for wall-time; 200 gives stable 95% percentile to ~0.02
SEED = 42


def load_phase3_per_celltype(bucket: str) -> tuple[np.ndarray, pd.DataFrame]:
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


def load_phase35_consolidated(suffix: str) -> ad.AnnData:
    p = PROCESSED / f"harmony_per_celltype_phase35_{suffix}.h5ad"
    if suffix == "subbucket_low":
        p = PROCESSED / "harmony_subbucket_phase35_low.h5ad"
    return ad.read_h5ad(p)


def split_phase35_bucket(a: ad.AnnData, bucket: str) -> tuple[np.ndarray, pd.DataFrame]:
    key = f"X_corrected_{bucket}"
    if key not in a.uns:
        raise KeyError(key)
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
        bucket_order = list(a.uns.get("per_bucket_r", {}).keys())
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
    x = np.asarray(a.obsm["X_harmony"])
    a.file.close()
    return x, obs_df


def load_cached_null(label: str, bucket: str, metric: str) -> np.ndarray:
    p = CACHE / f"perm_null_{label}_{metric}_{bucket}_{N_PERM}_{SEED}.npz"
    if not p.exists():
        return np.array([])
    return np.load(p)["null"]


def observed_value(x, obs, metric_fn):
    per_study = _per_study_donor_status(obs)
    cs = obs["study_id"].astype(str).to_numpy()
    cd = obs["donor_id"].astype(str).to_numpy()
    rvs = {}
    for sid in sorted(per_study.keys()):
        ds = per_study[sid]
        m = cs == sid
        rv = _response_vector(x[m], cd[m], ds.to_dict())
        if rv is not None:
            rvs[sid] = rv
    if len(rvs) < 2:
        return (
            float("nan"),
            {"mean": float("nan"), "median": float("nan"), "min": float("nan")},
            rvs,
        )
    summary = summary_stats_off_diag(rvs)
    return float(metric_fn(rvs)), summary, rvs


def observed_mmd(x, obs):
    per_study = _per_study_donor_status(obs)
    cs = obs["study_id"].astype(str).to_numpy()
    cd = obs["donor_id"].astype(str).to_numpy()
    embeddings = {}
    for sid in sorted(per_study.keys()):
        ds = per_study[sid].to_dict()
        m = cs == sid
        sub_x = x[m]
        sub_donor = cd[m]
        labels = np.array([ds.get(d, "unknown") for d in sub_donor])
        d_mask = labels == "diseased"
        h_mask = labels == "healthy_control"
        if not d_mask.any() or not h_mask.any():
            continue
        h_mean = sub_x[h_mask].mean(axis=0)
        embeddings[sid] = sub_x[d_mask] - h_mean
    if len(embeddings) < 2:
        return float("nan")
    return float(mmd_rbf_off_diag(embeddings, seed=SEED))


def calibrate_bucket_v2(x, obs, bucket: str, dataset_label: str) -> list[dict]:
    rows = []
    for metric_name, metric_fn in RV_METRICS.items():
        observed, summary, _ = observed_value(x, obs, metric_fn)
        if np.isnan(observed):
            rows.append(
                {
                    "bucket": bucket,
                    "metric": metric_name,
                    "observed": float("nan"),
                    "error": "insufficient_data",
                }
            )
            continue
        null = load_cached_null(dataset_label, bucket, metric_name)
        if len(null) == 0:
            rows.append(
                {
                    "bucket": bucket,
                    "metric": metric_name,
                    "observed": round(observed, 4),
                    "error": "no_cached_null",
                }
            )
            continue
        sh = split_half_with_metric(
            x, obs, bucket, metric_fn=metric_fn, n_splits=N_SPLITS, seed=SEED
        )
        # Observed-r bootstrap CI (Part A2)
        boot = bootstrap_observed_r(x, obs, metric_fn=metric_fn, n_bootstrap=N_BOOTSTRAP, seed=SEED)
        v95 = calibrated_gate_verdict(
            observed,
            null,
            sh["split_half_distribution"],
            percentile=95,
            alpha=0.05,
            use_corrected_ci=True,
        )
        v99 = calibrated_gate_verdict(
            observed,
            null,
            sh["split_half_distribution"],
            percentile=99,
            alpha=0.05,
            use_corrected_ci=True,
        )
        # Legacy verdict for audit comparison
        v99_legacy = calibrated_gate_verdict(
            observed,
            null,
            sh["split_half_distribution"],
            percentile=99,
            alpha=0.05,
            use_corrected_ci=False,
        )
        rows.append(
            {
                "bucket": bucket,
                "metric": metric_name,
                "observed": round(observed, 4),
                "observed_ci_low": round(boot["observed_ci_low"], 4)
                if not np.isnan(boot["observed_ci_low"])
                else float("nan"),
                "observed_ci_high": round(boot["observed_ci_high"], 4)
                if not np.isnan(boot["observed_ci_high"])
                else float("nan"),
                "summary_mean": round(summary["mean"], 4),
                "summary_median": round(summary["median"], 4),
                "summary_min": round(summary["min"], 4),
                "perm_p95": round(v95["null_threshold"], 4),
                "perm_p99": round(v99["null_threshold"], 4),
                "perm_p_value": round(v99["p_value"], 4),
                "perm_n_actual": len(null),
                "split_half_ceiling": round(sh["ceiling"], 4)
                if not np.isnan(sh["ceiling"])
                else float("nan"),
                "sh_ci_low_alpha05": round(v99["ci_low"], 4)
                if not np.isnan(v99["ci_low"])
                else float("nan"),
                "sh_ci_high_alpha05": round(v99["ci_high"], 4)
                if not np.isnan(v99["ci_high"])
                else float("nan"),
                "at_or_above_sh_ci_low": bool(v99["at_or_above_ci_low"]),
                "in_sh_ci_legacy": bool(v99["in_split_half_ci"]),
                "calibrated_pass_p95_corrected": bool(v95["pass"]),
                "calibrated_pass_p99_corrected": bool(v99["pass"]),
                "calibrated_pass_p99_legacy_in_ci": bool(v99_legacy["pass"]),
                "verdict_changed_vs_v1": bool(v99["pass"]) != bool(v99_legacy["pass"]),
                "n_perm": N_PERM,
                "n_splits": N_SPLITS,
                "n_bootstrap": N_BOOTSTRAP,
            }
        )
    # MMD observed-only
    mmd_obs = observed_mmd(x, obs)
    rows.append(
        {
            "bucket": bucket,
            "metric": "mmd_rbf",
            "observed": round(mmd_obs, 4) if not np.isnan(mmd_obs) else float("nan"),
            "note": "observed_only_no_perm_null_v1_limitation",
        }
    )
    return rows


def apply_fdr(rows: list[dict]) -> list[dict]:
    """Add fdr_corrected_p column to all rows that have perm_p_value."""
    ps = np.array([r.get("perm_p_value", np.nan) for r in rows], dtype=np.float64)
    adj = fdr_bh(ps)
    for r, q in zip(rows, adj, strict=False):
        r["fdr_corrected_p"] = round(float(q), 4) if not np.isnan(q) else float("nan")
        # Update verdicts to use FDR-corrected p < 0.01 instead of raw p99.
        if "perm_p_value" in r and not np.isnan(q):
            r["calibrated_pass_p99_corrected_fdr"] = bool(
                q < 0.01 and r.get("at_or_above_sh_ci_low", False)
            )
        else:
            r["calibrated_pass_p99_corrected_fdr"] = False
    return rows


def write_csv(rows: list[dict], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"wrote {out}: {len(df)} rows")


def run_phase3():
    rows = []
    for bucket in ("monocyte", "B", "NK", "CD4T", "CD8T"):
        print(f"[phase3 v2] {bucket}")
        x, obs = load_phase3_per_celltype(bucket)
        rows.extend(calibrate_bucket_v2(x, obs, bucket, "phase3"))
    rows = apply_fdr(rows)
    write_csv(rows, TABLES / "calibration_phase3_v2.csv")


def run_phase35_consolidated(suffix: str, label: str, out_name: str):
    a = load_phase35_consolidated(suffix)
    buckets = sorted(
        {k.removeprefix("X_corrected_") for k in a.uns if k.startswith("X_corrected_")}
    )
    rows = []
    for bucket in buckets:
        print(f"[{label} v2] {bucket}")
        try:
            x, obs = split_phase35_bucket(a, bucket)
        except KeyError:
            continue
        rows.extend(calibrate_bucket_v2(x, obs, bucket, label))
    rows = apply_fdr(rows)
    write_csv(rows, TABLES / out_name)


def run_global_harmony():
    x_all, obs_all = load_global_harmony()
    mask = obs_all["donor_disease_status"].isin(["diseased", "healthy_control"]).values
    x_all, obs_all = x_all[mask], obs_all.loc[mask].reset_index(drop=True)
    rows = []
    for bucket in ("monocyte", "B", "NK", "CD4T", "CD8T"):
        print(f"[global_harmony v2] {bucket}")
        m = (obs_all["coarse"] == bucket).values
        x = x_all[m]
        obs = obs_all.loc[m].reset_index(drop=True)
        if len(obs) < 100:
            continue
        rows.extend(calibrate_bucket_v2(x, obs, bucket, "global_harmony"))
    rows = apply_fdr(rows)
    write_csv(rows, TABLES / "calibration_phase3_global_harmony_v2.csv")


def main():
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
    for tag in targets:
        if tag == "phase3":
            run_phase3()
        elif tag == "phase35_low":
            run_phase35_consolidated("low", "phase35_low", "calibration_phase35_low_v2.csv")
        elif tag == "phase35_high":
            run_phase35_consolidated("high", "phase35_high", "calibration_phase35_high_v2.csv")
        elif tag == "phase35_subbucket":
            run_phase35_consolidated(
                "subbucket_low", "phase35_subbucket", "calibration_phase35_subbucket_v2.csv"
            )
        elif tag == "global_harmony":
            run_global_harmony()
        else:
            print(f"unknown target {tag}")


if __name__ == "__main__":
    main()
