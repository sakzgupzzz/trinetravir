"""Session 4.5 Part C: Wilk mito-fraction z-score measurement.

Per METHODS_CHOICES.md Issue 39 Amendment (Session 4.5 Part C pre-spec).
Decision rule pre-committed before this measurement runs. This script
computes the z-score and applies the rule mechanically.

Protocol:
  - For each of 4 v1 training studies, compute per-cell mito% (MT- gene
    counts / total counts × 100) on the reannotated h5ad (raw counts).
  - Aggregate to study-level mean.
  - Wilk z-score = (Wilk_mean - mean_others) / SD_others (other-3 reference).
  - Apply pre-committed tier rule (Tier I/II/III).

Output:
  results/tables/session4_5_wilk_mito_sensitivity.csv

Companion script for Tier II / Tier III escalation: write per-bucket
response-vector recompute with MT- genes excluded if needed (defer to
secondary script if measurement triggers escalation).
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"

STUDY_FILES = {
    "wilk_2020": PROC / "wilk_2020_reannotated.h5ad",
    "lee_2020": PROC / "lee_2020_reannotated.h5ad",
    "arunachalam_2020": PROC / "arunachalam_2020_reannotated.h5ad",
    "schulte_schrepping_2020": PROC / "schulte_schrepping_2020_reannotated.h5ad",
}

# Sample size per study for compute envelope; 10K cells per study sufficient for
# stable mean estimate of mito%.
SAMPLE_SIZE = 10_000
SEED = 42


def compute_study_mito_mean(study_path: Path, study_id: str) -> dict:
    """Compute per-cell mito% on sample, aggregate to study mean."""
    a = ad.read_h5ad(study_path, backed="r")
    n_cells = a.shape[0]

    # Resolve symbol column
    sym_col = next(
        (
            c
            for c in a.var.columns
            if c.lower() in ("feature_name", "gene_symbol", "name", "symbol")
        ),
        None,
    )
    symbols = a.var[sym_col] if sym_col else a.var.index.astype(str)
    mt_mask = np.array([str(s).startswith("MT-") for s in symbols])
    n_mt = int(mt_mask.sum())

    rng = np.random.default_rng(SEED)
    sample_idx = rng.choice(n_cells, size=min(SAMPLE_SIZE, n_cells), replace=False)
    sample_idx.sort()

    Xs = a.X[sample_idx]
    if sp.issparse(Xs):
        Xs = Xs.toarray()
    Xs = np.asarray(Xs, dtype=np.float64)

    total = Xs.sum(axis=1)
    mt_sum = Xs[:, mt_mask].sum(axis=1)
    pct = 100.0 * mt_sum / np.maximum(total, 1e-9)

    return {
        "study_id": study_id,
        "n_mt_genes": n_mt,
        "n_cells_total": n_cells,
        "n_cells_sampled": len(sample_idx),
        "mito_pct_mean": float(np.mean(pct)),
        "mito_pct_std": float(np.std(pct, ddof=1)),
        "mito_pct_median": float(np.median(pct)),
        "mito_pct_max": float(np.max(pct)),
    }


def main() -> int:
    rows = []
    for study_id, study_path in STUDY_FILES.items():
        logger.info("=== %s ===", study_id)
        r = compute_study_mito_mean(study_path, study_id)
        logger.info(
            "  %s: mito_pct mean=%.3f sd=%.3f median=%.3f max=%.3f (n_mt=%d, n_sampled=%d)",
            study_id,
            r["mito_pct_mean"],
            r["mito_pct_std"],
            r["mito_pct_median"],
            r["mito_pct_max"],
            r["n_mt_genes"],
            r["n_cells_sampled"],
        )
        rows.append(r)

    df = pd.DataFrame(rows)
    # Wilk z-score vs other-3 reference
    wilk = df[df["study_id"] == "wilk_2020"].iloc[0]
    others = df[df["study_id"] != "wilk_2020"]
    others_mean = float(others["mito_pct_mean"].mean())
    others_sd = float(others["mito_pct_mean"].std(ddof=1))  # SD across 3 study means
    z_wilk = (wilk["mito_pct_mean"] - others_mean) / others_sd if others_sd > 0 else float("nan")

    logger.info("\n=== Wilk z-score analysis ===")
    logger.info("  other 3 studies mean(mito_pct_mean) = %.3f", others_mean)
    logger.info("  other 3 studies SD(mito_pct_mean)   = %.3f", others_sd)
    logger.info("  Wilk mito_pct_mean                  = %.3f", wilk["mito_pct_mean"])
    logger.info("  Wilk z-score                        = %.3f", z_wilk)

    # Apply pre-committed Tier rule (Issue 39 Amendment, commit 1ea9dc9)
    abs_z = abs(z_wilk)
    if abs_z <= 1.0:
        tier = "TIER_I_REDUNDANT_DEFER"
        action = "Close Issue 39 as redundant given Part A.5b watchpoint clean. Document as v1 limitation."
    elif abs_z <= 2.0:
        tier = "TIER_II_RUN_SENSITIVITY"
        action = "Run mito-excluded response-vector recompute; verify |Δr_no_mito| ≤ 0.05 across 5 buckets."
    else:
        tier = "TIER_III_ESCALATE_NOW"
        action = (
            "Run Tier II sensitivity AND consider Issue 6 Tier IV verdict Wilk-specific caveat."
        )

    # Boundary tie-break note
    boundary_note = ""
    if 0.9 <= abs_z <= 1.1 or 1.9 <= abs_z <= 2.1:
        boundary_note = (
            "BOUNDARY case (|z| within ±0.1 of threshold); conservative tier applied per pre-spec."
        )
        logger.info("  %s", boundary_note)

    logger.info("\n=== Pre-committed Tier rule outcome ===")
    logger.info("  |z| = %.3f", abs_z)
    logger.info("  Tier: %s", tier)
    logger.info("  Action: %s", action)

    # Build output table
    df["wilk_z_score"] = np.where(df["study_id"] == "wilk_2020", z_wilk, np.nan)
    df["others_mean"] = others_mean
    df["others_sd"] = others_sd
    df["pre_committed_tier"] = tier
    df["pre_committed_action"] = action
    df["boundary_note"] = boundary_note

    out = TABLES / "session4_5_wilk_mito_sensitivity.csv"
    df.to_csv(out, index=False)
    logger.info("\nwrote %s", out.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
