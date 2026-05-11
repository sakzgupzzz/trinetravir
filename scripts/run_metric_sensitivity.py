"""Phase 3 metric sensitivity analysis (METHODS_CHOICES Issue 3).

Computes mean off-diagonal cross-study coherence per bucket under three
metrics — Pearson r, Spearman r, top-100 DE Jaccard — using the response
vectors cached by notebook 04 (data/processed/phase3_response_vectors_*.parquet).

Output: results/tables/metric_sensitivity_phase3.csv with one row per
(bucket, metric) plus the verdict vs heuristic threshold. Used to assess
whether the Phase 3 gate verdict is robust to choice of cross-study
coherence metric.

MMD (RBF kernel, median heuristic) is deferred to v1.5 because it requires
per-cell distributions in the Harmony-corrected embedding space; v1's
calibration cache stores summary statistics only. Documented in
METHODS_CHOICES Issue 3 resolution.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd

from trinetravir.data.harmonize import COARSE_BUCKETS
from trinetravir.eval.metrics import (
    load_phase3_response_vectors,
    metric_sensitivity_table,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

REPO = Path(__file__).resolve().parents[1]
PROCESSED = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)


def main() -> int:
    rvs = load_phase3_response_vectors(PROCESSED, list(COARSE_BUCKETS))
    if not rvs:
        print(
            f"No response vector parquets under {PROCESSED}. Run notebook 04 first.",
            file=sys.stderr,
        )
        return 1
    print(f"Loaded response vectors for {len(rvs)} buckets: {sorted(rvs.keys())}")

    table = metric_sensitivity_table(rvs)
    out_csv = TABLES / "metric_sensitivity_phase3.csv"
    table.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")
    print("\n--- Phase 3 metric sensitivity table ---")
    pd.set_option("display.width", 220)
    print(table.to_string(index=False))

    # Verdict consistency check: per bucket, count how many metrics flag exceeds_threshold=True
    print("\n--- per-bucket verdict consistency ---")
    pivot = table.pivot(index="bucket", columns="metric", values="exceeds_threshold")
    print(pivot.to_string())
    consistent = pivot.apply(lambda r: r.sum() >= 2, axis=1)
    print("\nbuckets where >= 2/3 metrics agree on PASS:")
    print(consistent.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
