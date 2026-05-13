"""Regenerate Session 4 Part A combined CSV + Issue 6 verdict CSV locally.

Triggered when the modal sweep does NOT complete the full 5-bucket loop
in a single run (e.g., 16h timeout cut off CD8T → CD8T solo resume in
a separate run → 5 per-bucket CSVs exist independently but no canonical
combined CSV).

Reads:
  results/tables/session4_scvi_per_bucket_<bucket>.csv × 5

Writes (overwriting any bucket-incomplete versions from partial sweeps):
  results/tables/session4_scvi_per_bucket_combined.csv
  results/tables/session4_issue6_verdict.csv

Verdict logic per METHODS_CHOICES.md Issue 34 (pre-spec):
  Tier I HARMONY_ADEQUATE: max(Δr_mvs) ≤ 0.05
  Tier II MIXED: ≥1 bucket Δr_mvs ∈ (0.05, 0.10], none > 0.10
  Tier III SCVI_PREFERRED: ≥3 buckets Δr_mvs > 0.10 OR any > 0.20
  Tier IV HARMONY_PREFERRED: ≥3 buckets Δr_mvs < -0.10
  Boundary ties → conservative tier per Issue 34.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
TABLES = REPO / "results" / "tables"

BUCKETS = ("monocyte", "B", "NK", "CD4T", "CD8T")


def apply_verdict(combined: pd.DataFrame) -> str:
    """Issue 34 four-tier verdict from per-bucket Δr_mvs.

    Tier I CORRECTED per Issue 34 Amendment 1 (2026-05-12, METHODS_CHOICES.md):
    two-sided proximity (max(|Δr|) ≤ 0.05) rather than one-sided (max(Δr) ≤ 0.05).
    The one-sided rule fired when Harmony arbitrarily-strongly beat scVI,
    contradicting the label "HARMONY_ADEQUATE: scVI ≈ Harmony, no difference."
    Amendment is data-direction-independent: same correction applies regardless of
    which method dominated.
    """
    deltas = combined["delta_r_mvs"].astype(float).values
    max_abs_d = float(np.max(np.abs(deltas)))
    above_010 = int((deltas > 0.10).sum())
    above_020 = int((deltas > 0.20).sum())
    below_neg010 = int((deltas < -0.10).sum())
    in_005_010 = int(((deltas > 0.05) & (deltas <= 0.10)).sum())

    # Tier I (Amendment 1): two-sided proximity
    if max_abs_d <= 0.05:
        return f"TIER_I_HARMONY_ADEQUATE (max |Δr_mvs| = {max_abs_d:.4f} ≤ 0.05)"
    # Tier III SCVI_PREFERRED (one-sided)
    if above_010 >= 3 or above_020 >= 1:
        return (
            f"TIER_III_SCVI_PREFERRED (≥3 buckets Δr_mvs > 0.10 [n={above_010}] "
            f"OR any > 0.20 [n={above_020}])"
        )
    # Tier IV HARMONY_PREFERRED (one-sided)
    if below_neg010 >= 3:
        return f"TIER_IV_HARMONY_PREFERRED (≥3 buckets Δr_mvs < -0.10 [n={below_neg010}])"
    # Tier II MIXED (scVI marginal)
    if in_005_010 >= 1 and above_010 == 0:
        return "TIER_II_MIXED (at least one Δr_mvs in (0.05, 0.10], none > 0.10)"
    # Boundary case → conservative
    return "TIER_II_MIXED (boundary case; conservative tie-break)"


def main() -> int:
    missing = []
    per_bucket_dfs = []
    for bucket in BUCKETS:
        p = TABLES / f"session4_scvi_per_bucket_{bucket}.csv"
        if not p.exists():
            missing.append(bucket)
            continue
        df = pd.read_csv(p)
        per_bucket_dfs.append(df)

    if missing:
        logger.error("missing per-bucket CSVs: %s. Cannot regenerate.", missing)
        logger.error("Expected files: results/tables/session4_scvi_per_bucket_<bucket>.csv × 5")
        return 1

    logger.info("loaded all 5 per-bucket CSVs: %s", list(BUCKETS))

    # Combined: only the calibrated best-config rows
    combined = pd.concat(
        [df[df["calibrated"] == True] for df in per_bucket_dfs],  # noqa: E712
        ignore_index=True,
    )

    combined_out = TABLES / "session4_scvi_per_bucket_combined.csv"
    combined.to_csv(combined_out, index=False)
    logger.info("wrote %s: %d rows", combined_out.name, len(combined))

    # Verdict
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
                "n_buckets": int(len(combined)),
                "regenerated_locally": True,
                "source": "scripts/session4_part_a_regen_combined.py",
            }
        ]
    )
    verdict_out = TABLES / "session4_issue6_verdict.csv"
    verdict_df.to_csv(verdict_out, index=False)
    logger.info("wrote %s", verdict_out.name)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
