"""Session 6B post-bg consolidator: combine N=1000 results + apply Issues 27-30 + 31 verdicts.

Reads:
  results/tables/heldout_v2_calibration_<cohort>_n1000.csv (yoshida, allen, gse157829)
  results/tables/heldout_v2_calibration_randolph_2021_n1000_issue31.csv (Randolph with
      Issue 31 cross-bucket healthy reference; supersedes pre-Issue-31 Randolph N=1000).

Outputs:
  results/tables/heldout_v2_calibration_combined_n1000.csv (overwrites; 15 rows total)
  results/tables/heldout_issue_verdicts_n1000.csv (mechanical Issue 27-30 verdicts)

Pre-committed decision rules (mechanically applied, no interpretation):
  Issue 27 (Randolph monocyte_infected, primary post-Issue-31):
    MVS r ≥ 0.40 SUPPORTS_H1; < 0.20 CHALLENGES; in [0.20, 0.40] INCONCLUSIVE.
  Issue 27-sensitivity (Randolph monocyte bystander):
    Same rule; reported as supplementary sensitivity.
  Issue 28 (Yoshida monocyte cross-age):
    MVS r ≥ 0.30 SUPPORTS_H1; < 0.10 CHALLENGES; in [0.10, 0.30] PARTIAL.
  Issue 29 (Allen Atlas monocyte CMV chronic-vs-naive):
    MVS r in [0.10, 0.40] APPROPRIATE; > 0.50 OVER_PREDICTION; < 0.05 NO_SHARED_BIOLOGY.
  Issue 30 (GSE157829 CD4T retrovirus):
    MVS r in [0.00, 0.20] EXPECTED; > 0.40 SURPRISING; < -0.10 ANTI_CORRELATION.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from trinetravir.eval.calibration import fdr_bh

REPO = Path(__file__).resolve().parents[1]
TABLES = REPO / "results" / "tables"

# Cohort -> source CSV filename for N=1000 combine.
COHORT_SOURCES = {
    "yoshida_2022": "heldout_v2_calibration_yoshida_2022_n1000.csv",
    "allen_atlas_monocyte": "heldout_v2_calibration_allen_atlas_monocyte_n1000.csv",
    "gse157829": "heldout_v2_calibration_gse157829_n1000.csv",
    "randolph_2021": "heldout_v2_calibration_randolph_2021_n1000_issue31.csv",
}


def issue_27_verdict(r_mvs: float) -> str:
    if np.isnan(r_mvs):
        return "n/a"
    if r_mvs >= 0.40:
        return "SUPPORTS_H1 (r_mvs ≥ 0.40)"
    if r_mvs < 0.20:
        return "CHALLENGES_H1 (r_mvs < 0.20)"
    return "INCONCLUSIVE (r_mvs in [0.20, 0.40])"


def issue_28_verdict(r_mvs: float) -> str:
    if np.isnan(r_mvs):
        return "n/a"
    if r_mvs >= 0.30:
        return "SUPPORTS_H1 (r_mvs ≥ 0.30; transfer across age groups)"
    if r_mvs < 0.10:
        return "CHALLENGES_H1 (r_mvs < 0.10; conserved component does NOT transfer to pediatric)"
    return "PARTIAL_TRANSFER (r_mvs in [0.10, 0.30])"


def issue_29_verdict(r_mvs: float) -> str:
    if np.isnan(r_mvs):
        return "n/a"
    if 0.10 <= r_mvs <= 0.40:
        return "APPROPRIATE_DISCRIMINATION (r_mvs in [0.10, 0.40])"
    if r_mvs > 0.50:
        return "CONCERNING_OVER_PREDICTION (r_mvs > 0.50)"
    if r_mvs < 0.05:
        return "CONCERNING_NO_SHARED_BIOLOGY (r_mvs < 0.05)"
    return f"BORDERLINE (r_mvs = {r_mvs:.3f})"


def issue_30_verdict(r_mvs: float) -> str:
    if np.isnan(r_mvs):
        return "n/a"
    if 0.00 <= r_mvs <= 0.20:
        return "EXPECTED_RETROVIRUS_DISTINCTNESS (r_mvs in [0.00, 0.20])"
    if r_mvs > 0.40:
        return "SURPRISING_HIGH (r_mvs > 0.40)"
    if r_mvs < -0.10:
        return "ANTI_CORRELATION (r_mvs < -0.10)"
    return f"BORDERLINE (r_mvs = {r_mvs:.3f})"


def main() -> int:
    parts = []
    for cohort, fname in COHORT_SOURCES.items():
        p = TABLES / fname
        if not p.exists():
            print(f"missing {p.name}; skip")
            continue
        df = pd.read_csv(p)
        if "cohort" not in df.columns:
            df["cohort"] = cohort
        parts.append(df)
    combined = pd.concat(parts, ignore_index=True)

    # FDR-BH on perm_p_value_mvs across all 15 tests
    p_mvs = combined["perm_p_value_mvs"].astype(float).values
    combined["fdr_corrected_p_mvs"] = np.round(fdr_bh(p_mvs), 4)
    combined["fdr_corrected_p_full"] = np.round(
        fdr_bh(combined["perm_p_value_full"].astype(float).values), 4
    )
    combined["calibrated_pass_mvs_fdr_05"] = combined["fdr_corrected_p_mvs"] < 0.05
    combined["calibrated_pass_mvs_fdr_01"] = combined["fdr_corrected_p_mvs"] < 0.01

    out_combined = TABLES / "heldout_v2_calibration_combined_n1000.csv"
    combined.to_csv(out_combined, index=False)
    print(f"wrote {out_combined.name}: {len(combined)} rows")

    # Per-Issue side-by-side
    issue_rows = []

    def add_row(issue_label: str, primary_bucket: str, cohort: str, rule: str, verdict_fn):
        sel = combined[(combined["cohort"] == cohort) & (combined["bucket"] == primary_bucket)]
        if sel.empty:
            return
        row = sel.iloc[0]
        idx = sel.index[0]
        issue_rows.append(
            {
                "issue": issue_label,
                "primary_bucket": primary_bucket,
                "rule": rule,
                "observed_r_full": row["observed_r_full"],
                "observed_r_mvs": row["observed_r_mvs"],
                "perm_p_value_full": row["perm_p_value_full"],
                "perm_p_value_mvs": row["perm_p_value_mvs"],
                "fdr_corrected_p_mvs": combined.loc[idx, "fdr_corrected_p_mvs"],
                "fdr_corrected_p_full": combined.loc[idx, "fdr_corrected_p_full"],
                "ci_low_mvs": row.get("ci_low_mvs", float("nan")),
                "ci_high_mvs": row.get("ci_high_mvs", float("nan")),
                "verdict": verdict_fn(row["observed_r_mvs"]),
                "n_d_cells": row.get("n_d_cells", float("nan")),
                "n_h_cells": row.get("n_h_cells", float("nan")),
                "healthy_reference_bucket": row.get("healthy_reference_bucket", "self"),
            }
        )

    # Issue 27 PRIMARY: Randolph monocyte_infected (corrected per Issue 31)
    add_row(
        "27 Randolph cross-context IAV (PRIMARY, post-Issue-31)",
        "monocyte_infected",
        "randolph_2021",
        "r_mvs >= 0.40 SUPPORTS_H1; < 0.20 CHALLENGES; in [0.20, 0.40] INCONCLUSIVE",
        issue_27_verdict,
    )
    # Issue 27 SENSITIVITY: Randolph bystander monocyte
    add_row(
        "27 Randolph cross-context IAV (SENSITIVITY, bystander)",
        "monocyte",
        "randolph_2021",
        "r_mvs >= 0.40 SUPPORTS_H1; < 0.20 CHALLENGES; in [0.20, 0.40] INCONCLUSIVE",
        issue_27_verdict,
    )
    add_row(
        "28 Yoshida pediatric vs adult cross-age",
        "monocyte",
        "yoshida_2022",
        "r_mvs >= 0.30 SUPPORTS_H1; < 0.10 CHALLENGES; in [0.10, 0.30] PARTIAL",
        issue_28_verdict,
    )
    add_row(
        "29 Allen Atlas chronic-latent CMV vs naive",
        "monocyte",
        "allen_atlas_monocyte",
        "r_mvs in [0.10, 0.40] APPROPRIATE; > 0.50 OVER_PREDICTION; < 0.05 NO_SHARED_BIOLOGY",
        issue_29_verdict,
    )
    add_row(
        "30 GSE157829 HIV retrovirus context",
        "CD4T",
        "gse157829",
        "r_mvs in [0.00, 0.20] EXPECTED; > 0.40 SURPRISING; < -0.10 ANTI_CORRELATION",
        issue_30_verdict,
    )

    df_issues = pd.DataFrame(issue_rows)
    out_verdicts = TABLES / "heldout_issue_verdicts_n1000.csv"
    df_issues.to_csv(out_verdicts, index=False)
    print(f"\nwrote {out_verdicts.name}\n")
    print("=== Issue verdicts (mechanical application, N=1000) ===\n")
    for _, r in df_issues.iterrows():
        print(f"{r['issue']} ({r['primary_bucket']}, healthy_ref={r['healthy_reference_bucket']}):")
        print(f"  rule: {r['rule']}")
        print(
            f"  observed r_mvs = {r['observed_r_mvs']:.4f}  "
            f"(r_full = {r['observed_r_full']:.4f}; perm p_mvs = {r['perm_p_value_mvs']:.3f}; "
            f"FDR p_mvs = {r['fdr_corrected_p_mvs']:.3f})"
        )
        if not (pd.isna(r["ci_low_mvs"]) or pd.isna(r["ci_high_mvs"])):
            print(f"  bootstrap CI r_mvs = [{r['ci_low_mvs']:.3f}, {r['ci_high_mvs']:.3f}]")
        print(f"  verdict: {r['verdict']}\n")

    n_fdr_05 = int(combined["calibrated_pass_mvs_fdr_05"].sum())
    n_fdr_01 = int(combined["calibrated_pass_mvs_fdr_01"].sum())
    print(f"FDR<0.05 hits: {n_fdr_05}/{len(combined)}; FDR<0.01 hits: {n_fdr_01}/{len(combined)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
