"""Build harmonization_protocol_calibrated_comparison.csv (Issue 7).

For each (bucket, metric), join the per-cell-type and global Harmony
calibration tables to produce per-bucket calibrated verdicts and a
difference-significance test.

Difference significance: bootstrap CI overlap of the two per-bucket
permutation null distributions (per-cell-type vs global). If observed
delta exceeds 95% of bootstrap deltas under the null, flag significant.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
TABLES = REPO / "results" / "tables"
CACHE = REPO / "data" / "processed" / "calibration_cache"


def load_null(
    label: str, metric: str, bucket: str, n_perm: int = 1000, seed: int = 42
) -> np.ndarray:
    p = CACHE / f"perm_null_{label}_{metric}_{bucket}_{n_perm}_{seed}.npz"
    if not p.exists():
        return np.array([])
    return np.load(p)["null"]


def diff_significance(
    obs_a: float,
    null_a: np.ndarray,
    obs_b: float,
    null_b: np.ndarray,
    *,
    alpha: float = 0.05,
    n_bootstrap: int = 1000,
    seed: int = 42,
) -> dict:
    """Bootstrap test: is observed delta (obs_a - obs_b) significantly different from 0?

    Bootstraps a null delta distribution by resampling values from null_a and null_b
    independently, computing their delta. Significant iff observed delta falls outside
    the (1-alpha) CI of the bootstrap null deltas.
    """
    if len(null_a) == 0 or len(null_b) == 0:
        return {
            "observed_delta": obs_a - obs_b,
            "ci_low": float("nan"),
            "ci_high": float("nan"),
            "significant": False,
        }
    rng = np.random.default_rng(seed)
    deltas = rng.choice(null_a, size=n_bootstrap, replace=True) - rng.choice(
        null_b, size=n_bootstrap, replace=True
    )
    lo = float(np.percentile(deltas, 100 * alpha / 2))
    hi = float(np.percentile(deltas, 100 * (1 - alpha / 2)))
    obs_delta = float(obs_a - obs_b)
    return {
        "observed_delta": obs_delta,
        "ci_low": lo,
        "ci_high": hi,
        "significant": bool(obs_delta < lo or obs_delta > hi),
    }


def main() -> None:
    p_pc = TABLES / "calibration_phase3.csv"
    p_gl = TABLES / "calibration_phase3_global_harmony.csv"
    if not p_pc.exists() or not p_gl.exists():
        raise SystemExit("missing one of the calibration tables")

    pc = pd.read_csv(p_pc)
    gl = pd.read_csv(p_gl)
    pc = pc[pc["metric"] != "mmd_rbf"]
    gl = gl[gl["metric"] != "mmd_rbf"]

    out_rows = []
    for (bucket, metric), pc_row in pc.set_index(["bucket", "metric"]).iterrows():
        try:
            gl_row = gl.set_index(["bucket", "metric"]).loc[(bucket, metric)]
        except KeyError:
            continue
        null_pc = load_null("phase3", metric, bucket)
        null_gl = load_null("global_harmony", metric, bucket)
        sig = diff_significance(
            float(pc_row["observed"]),
            null_pc,
            float(gl_row["observed"]),
            null_gl,
        )
        out_rows.append(
            {
                "bucket": bucket,
                "metric": metric,
                "per_cell_type_observed_r": pc_row["observed"],
                "per_cell_type_null_p99": pc_row["perm_p99"],
                "per_cell_type_split_half_ceiling": pc_row["split_half_ceiling"],
                "per_cell_type_calibrated_verdict_p99": pc_row["calibrated_pass_p99_alpha05"],
                "global_observed_r": gl_row["observed"],
                "global_null_p99": gl_row["perm_p99"],
                "global_split_half_ceiling": gl_row["split_half_ceiling"],
                "global_calibrated_verdict_p99": gl_row["calibrated_pass_p99_alpha05"],
                "per_bucket_difference": round(sig["observed_delta"], 4),
                "diff_ci_low_alpha05": round(sig["ci_low"], 4)
                if not np.isnan(sig["ci_low"])
                else float("nan"),
                "diff_ci_high_alpha05": round(sig["ci_high"], 4)
                if not np.isnan(sig["ci_high"])
                else float("nan"),
                "difference_significant_alpha05": sig["significant"],
            }
        )

    out = TABLES / "harmonization_protocol_calibrated_comparison.csv"
    pd.DataFrame(out_rows).to_csv(out, index=False)
    print(f"wrote {out}: {len(out_rows)} rows")


if __name__ == "__main__":
    main()
