"""Session 4.5 Part C Tier II: mito-excluded response-vector recompute.

Triggered by Issue 39 Amendment pre-spec (commit 1ea9dc9) Tier II (1.0 < |z| ≤ 2.0).
Wilk measured z=1.859 → Tier II RUN_SENSITIVITY.

Protocol per pre-spec:
  - Load per-bucket Harmony response vectors (4-study columns)
  - Drop MT- prefixed genes from gene set
  - Recompute cross-study mean off-diagonal Pearson r per bucket
  - Compare to original r (with mito genes) → compute Δr_no_mito
  - PASS if |Δr_no_mito| ≤ 0.05 across all 5 buckets
  - FAIL if any bucket exceeds 0.05

Output:
  results/tables/session4_5_wilk_mito_tier2_sensitivity.csv

Verdict:
  PASS  → Issue 39 resolved as "no mito-driven artifact"
  FAIL  → Escalate to Phase 5 supplementary with explicit manuscript disclosure
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"

BUCKETS = ("monocyte", "B", "NK", "CD4T", "CD8T")
PASS_THRESHOLD = 0.05


def mean_off_diagonal_pearson(df: pd.DataFrame) -> float:
    """Mean off-diagonal Pearson r across columns (studies)."""
    corr = df.corr().values
    iu = np.triu_indices_from(corr, k=1)
    return float(np.nanmean(corr[iu]))


def main() -> int:
    rows = []
    for bucket in BUCKETS:
        p = PROC / f"phase3_response_vectors_{bucket}.parquet"
        if not p.exists():
            logger.warning("missing %s", p.name)
            continue
        rv = pd.read_parquet(p)

        # Identify MT- genes
        mt_mask = rv.index.astype(str).str.startswith("MT-")
        n_mt = int(mt_mask.sum())
        n_genes_total = rv.shape[0]

        # Cross-study Pearson r WITH MT
        r_with_mito = mean_off_diagonal_pearson(rv)

        # Cross-study Pearson r WITHOUT MT
        rv_no_mito = rv[~mt_mask]
        r_no_mito = mean_off_diagonal_pearson(rv_no_mito)

        delta_r = r_no_mito - r_with_mito

        rows.append(
            {
                "bucket": bucket,
                "n_genes_total": n_genes_total,
                "n_mt_genes": n_mt,
                "n_genes_no_mito": int((~mt_mask).sum()),
                "r_with_mito": round(r_with_mito, 4),
                "r_no_mito": round(r_no_mito, 4),
                "delta_r_no_mito": round(delta_r, 4),
                "abs_delta_r": round(abs(delta_r), 4),
                "pass_005_threshold": bool(abs(delta_r) <= PASS_THRESHOLD),
            }
        )
        logger.info(
            "  %s: r_with_mito=%.4f r_no_mito=%.4f Δr=%+.4f (|Δr|=%.4f) %s (n_mt=%d/%d)",
            bucket,
            r_with_mito,
            r_no_mito,
            delta_r,
            abs(delta_r),
            "PASS" if abs(delta_r) <= PASS_THRESHOLD else "FAIL",
            n_mt,
            n_genes_total,
        )

    df = pd.DataFrame(rows)
    out = TABLES / "session4_5_wilk_mito_tier2_sensitivity.csv"
    df.to_csv(out, index=False)
    logger.info("\nwrote %s", out.name)

    # Aggregate verdict
    all_pass = bool(df["pass_005_threshold"].all())
    max_abs = float(df["abs_delta_r"].max())
    logger.info("\n=== Tier II verdict ===")
    logger.info("  max |Δr_no_mito| across 5 buckets = %.4f (threshold 0.05)", max_abs)
    if all_pass:
        logger.info("  ALL 5 BUCKETS PASS → Issue 39 resolved: NO_MITO_DRIVEN_ARTIFACT")
        logger.info("  Mito-gene inclusion does not materially affect cross-study coherence.")
    else:
        failing = df[~df["pass_005_threshold"]]["bucket"].tolist()
        logger.info("  FAIL: %s exceed |Δr|=0.05 threshold", failing)
        logger.info("  → Escalate to Phase 5 supplementary; explicit manuscript disclosure.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
