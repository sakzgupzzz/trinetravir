"""Session 6B Step 4: combine all 4 cohort v2 calibration tables + mechanical decision-rule application.

Reads:
  results/tables/heldout_v2_calibration_<cohort>.csv (4 cohorts)

Outputs:
  results/tables/heldout_v2_calibration_combined.csv  (all 14 rows + FDR + Issue 27-30 verdict columns)
  results/tables/heldout_issue_verdicts.csv  (per-Issue side-by-side rule + observed + verdict)

Pre-committed decision rules (mechanically applied, no interpretation):
  Issue 27 (Randolph monocyte): MVS r ≥ 0.40 supports H1; < 0.20 challenges; in [0.20, 0.40] inconclusive.
  Issue 28 (Yoshida monocyte cross-age): MVS r ≥ 0.30 supports H1; < 0.10 challenges; in [0.10, 0.30] partial.
  Issue 29 (Allen Atlas monocyte chronic-vs-acute): MVS r in [0.10, 0.40] appropriate; > 0.50 over-prediction concerning; < 0.05 no shared biology concerning.
  Issue 30 (GSE157829 CD4T retrovirus): MVS r in [0.00, 0.20] expected; > 0.40 surprising; < -0.10 anti-correlation interpretable.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from trinetravir.eval.calibration import fdr_bh

REPO = Path(__file__).resolve().parents[1]
TABLES = REPO / "results" / "tables"

COHORTS = ["randolph_2021", "yoshida_2022", "allen_atlas_monocyte", "gse157829"]


def issue_27_verdict(r_mvs: float) -> str:
    """Randolph monocyte cross-context conserved-component."""
    if np.isnan(r_mvs):
        return "n/a"
    if r_mvs >= 0.40:
        return "SUPPORTS_H1 (r_mvs ≥ 0.40)"
    if r_mvs < 0.20:
        return "CHALLENGES_H1 (r_mvs < 0.20)"
    return "INCONCLUSIVE (r_mvs in [0.20, 0.40])"


def issue_28_verdict(r_mvs: float) -> str:
    """Yoshida monocyte cross-age."""
    if np.isnan(r_mvs):
        return "n/a"
    if r_mvs >= 0.30:
        return "SUPPORTS_H1 (r_mvs ≥ 0.30; transfer across age groups)"
    if r_mvs < 0.10:
        return "CHALLENGES_H1 (r_mvs < 0.10; conserved component does NOT transfer to pediatric)"
    return "PARTIAL_TRANSFER (r_mvs in [0.10, 0.30])"


def issue_29_verdict(r_mvs: float) -> str:
    """Allen Atlas monocyte chronic-latent CMV vs naive."""
    if np.isnan(r_mvs):
        return "n/a"
    if 0.10 <= r_mvs <= 0.40:
        return "APPROPRIATE_DISCRIMINATION (r_mvs in [0.10, 0.40])"
    if r_mvs > 0.50:
        return "CONCERNING_OVER_PREDICTION (r_mvs > 0.50; conserved component is acute-disease-non-specific)"
    if r_mvs < 0.05:
        return "CONCERNING_NO_SHARED_BIOLOGY (r_mvs < 0.05; conserved component is chronic-naive non-discriminative or anti-correlated)"
    return f"BORDERLINE (r_mvs = {r_mvs:.3f}; between ranges)"


def issue_30_verdict(r_mvs: float) -> str:
    """GSE157829 CD4T retrovirus context."""
    if np.isnan(r_mvs):
        return "n/a"
    if 0.00 <= r_mvs <= 0.20:
        return "EXPECTED_RETROVIRUS_DISTINCTNESS (r_mvs in [0.00, 0.20])"
    if r_mvs > 0.40:
        return "SURPRISING_HIGH (r_mvs > 0.40; framework not capturing retrovirus-specific biology)"
    if r_mvs < -0.10:
        return "ANTI_CORRELATION (r_mvs < -0.10; HIV CD4T response opposite to acute RNA virus)"
    return f"BORDERLINE (r_mvs = {r_mvs:.3f})"


def main() -> int:
    parts = []
    for cohort in COHORTS:
        p = TABLES / f"heldout_v2_calibration_{cohort}.csv"
        if not p.exists():
            print(f"missing {p.name}; skip")
            continue
        df = pd.read_csv(p)
        parts.append(df)
    combined = pd.concat(parts, ignore_index=True)

    # FDR-BH across all rows on MVS p-value (where present)
    p_mvs = combined["perm_p_value_mvs"].astype(float).values
    combined["fdr_corrected_p_mvs"] = np.round(fdr_bh(p_mvs), 4)
    combined["calibrated_pass_p99_mvs_fdr"] = combined["fdr_corrected_p_mvs"] < 0.01

    combined.to_csv(TABLES / "heldout_v2_calibration_combined.csv", index=False)
    print(f"wrote heldout_v2_calibration_combined.csv: {len(combined)} rows")

    # Per-Issue side-by-side
    issue_rows = []

    # Issue 27 (Randolph monocyte)
    r27 = combined[(combined["cohort"] == "randolph_2021") & (combined["bucket"] == "monocyte")]
    if not r27.empty:
        row = r27.iloc[0]
        issue_rows.append(
            {
                "issue": "27 Randolph ex vivo IAV",
                "primary_bucket": "monocyte",
                "rule": "r_mvs >= 0.40 SUPPORTS_H1; < 0.20 CHALLENGES; in [0.20, 0.40] INCONCLUSIVE",
                "observed_r_mvs": row["observed_r_mvs"],
                "observed_r_full": row["observed_r_full"],
                "perm_p_value_mvs": row["perm_p_value_mvs"],
                "fdr_corrected_p_mvs": combined.loc[r27.index[0], "fdr_corrected_p_mvs"],
                "verdict": issue_27_verdict(row["observed_r_mvs"]),
            }
        )

    # Issue 28 (Yoshida monocyte cross-age)
    r28 = combined[(combined["cohort"] == "yoshida_2022") & (combined["bucket"] == "monocyte")]
    if not r28.empty:
        row = r28.iloc[0]
        issue_rows.append(
            {
                "issue": "28 Yoshida pediatric vs adult cross-age",
                "primary_bucket": "monocyte",
                "rule": "r_mvs >= 0.30 SUPPORTS_H1; < 0.10 CHALLENGES; in [0.10, 0.30] PARTIAL",
                "observed_r_mvs": row["observed_r_mvs"],
                "observed_r_full": row["observed_r_full"],
                "perm_p_value_mvs": row["perm_p_value_mvs"],
                "fdr_corrected_p_mvs": combined.loc[r28.index[0], "fdr_corrected_p_mvs"],
                "verdict": issue_28_verdict(row["observed_r_mvs"]),
            }
        )

    # Issue 29 (Allen Atlas monocyte chronic-vs-naive)
    r29 = combined[
        (combined["cohort"] == "allen_atlas_monocyte") & (combined["bucket"] == "monocyte")
    ]
    if not r29.empty:
        row = r29.iloc[0]
        issue_rows.append(
            {
                "issue": "29 Allen Atlas chronic-latent CMV vs naive",
                "primary_bucket": "monocyte",
                "rule": "r_mvs in [0.10, 0.40] APPROPRIATE; > 0.50 OVER_PREDICTION; < 0.05 NO_SHARED_BIOLOGY",
                "observed_r_mvs": row["observed_r_mvs"],
                "observed_r_full": row["observed_r_full"],
                "perm_p_value_mvs": row["perm_p_value_mvs"],
                "fdr_corrected_p_mvs": combined.loc[r29.index[0], "fdr_corrected_p_mvs"],
                "verdict": issue_29_verdict(row["observed_r_mvs"]),
            }
        )

    # Issue 30 (GSE157829 CD4T retrovirus)
    r30 = combined[(combined["cohort"] == "gse157829") & (combined["bucket"] == "CD4T")]
    if not r30.empty:
        row = r30.iloc[0]
        issue_rows.append(
            {
                "issue": "30 GSE157829 HIV retrovirus context",
                "primary_bucket": "CD4T",
                "rule": "r_mvs in [0.00, 0.20] EXPECTED; > 0.40 SURPRISING; < -0.10 ANTI_CORRELATION",
                "observed_r_mvs": row["observed_r_mvs"],
                "observed_r_full": row["observed_r_full"],
                "perm_p_value_mvs": row["perm_p_value_mvs"],
                "fdr_corrected_p_mvs": combined.loc[r30.index[0], "fdr_corrected_p_mvs"],
                "verdict": issue_30_verdict(row["observed_r_mvs"]),
            }
        )

    df_issues = pd.DataFrame(issue_rows)
    df_issues.to_csv(TABLES / "heldout_issue_verdicts.csv", index=False)
    print("\n=== Issue verdicts (mechanical application) ===\n")
    for _, r in df_issues.iterrows():
        print(f"{r['issue']} ({r['primary_bucket']}):")
        print(f"  rule: {r['rule']}")
        print(
            f"  observed r_mvs = {r['observed_r_mvs']:.4f}  (r_full = {r['observed_r_full']:.4f}; perm p={r['perm_p_value_mvs']:.3f}; FDR-corrected p={r['fdr_corrected_p_mvs']:.3f})"
        )
        print(f"  verdict: {r['verdict']}")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
