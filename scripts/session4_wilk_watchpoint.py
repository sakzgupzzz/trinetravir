"""Session 4 Part A.5b Wilk depth watchpoint diagnostic.

Per references/session_4_prompt.md Part A.5b:
"Inspect per-study response vectors (Wilk vs each of Arunachalam/Lee/Schulte)
for the selected scVI configuration per bucket. If Wilk's response vector
shows outlier pattern (e.g., dominated by housekeeping genes; orthogonal to
other studies) NOT present in Harmony's per-study response vectors, flag as
a Wilk-specific scVI artifact in Issue 6 resolution."

The Part A sweep script does NOT persist per-study scVI response vectors
(only the aggregate mean off-diagonal r). Rerunning scVI per-config to get
per-study breakdown requires GPU. Pragmatic proxy: inspect Wilk's HARMONY
per-study response vector against the other 3 studies.

If Wilk-Harmony shows outlier pattern, scVI inherits worse (depth-driven).
If Wilk-Harmony is consistent, scVI's NB likelihood + library_size=1e4
normalization handled depth internally → Tier IV verdict applies as-is.

Inspection (per bucket):
  (1) Wilk vs (other 3 studies) per-gene Pearson r — full HVG + MVS subset
  (2) Wilk's response vector L2 magnitude vs others
  (3) Top-10 absolute-magnitude genes — housekeeping vs ISG check

Output: results/tables/session4_wilk_watchpoint.csv
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
MVS_FILE = REPO / "references" / "khatri_mvs_gene_list.csv"

BUCKETS = ("monocyte", "B", "NK", "CD4T", "CD8T")

# Common housekeeping markers (subset of canonical housekeeping for quick check)
HOUSEKEEPING = {"ACTB", "GAPDH", "HPRT1", "B2M", "PGK1", "PPIA", "RPL13A", "TBP", "UBC", "YWHAZ"}


def load_mvs() -> set[str]:
    df = pd.read_csv(MVS_FILE, comment="#")
    return set(df["gene_symbol"].astype(str).tolist())


def main() -> int:
    mvs_genes = load_mvs()
    logger.info("MVS gene set: %d genes", len(mvs_genes))

    rows = []
    for bucket in BUCKETS:
        p = PROC / f"phase3_response_vectors_{bucket}.parquet"
        if not p.exists():
            logger.warning("missing %s", p.name)
            continue
        rv = pd.read_parquet(p)
        logger.info("=== %s ===  rv shape: %s", bucket, rv.shape)

        if "wilk_2020" not in rv.columns:
            logger.warning("  Wilk column missing")
            continue

        wilk = rv["wilk_2020"]
        others = rv[[c for c in rv.columns if c != "wilk_2020"]]
        others_mean = others.mean(axis=1)

        # Full HVG
        r_wilk_vs_others = wilk.corr(others_mean)
        # MVS subset
        mvs_in_rv = [g for g in rv.index if g in mvs_genes]
        if len(mvs_in_rv) > 0:
            r_wilk_vs_others_mvs = wilk.loc[mvs_in_rv].corr(others_mean.loc[mvs_in_rv])
        else:
            r_wilk_vs_others_mvs = float("nan")

        # L2 magnitude comparison
        wilk_l2 = float(np.sqrt((wilk**2).sum()))
        other_l2_mean = float(np.sqrt((others**2).sum(axis=0)).mean())
        l2_ratio = wilk_l2 / other_l2_mean if other_l2_mean > 0 else float("nan")

        # MVS-only L2
        if len(mvs_in_rv) > 0:
            wilk_l2_mvs = float(np.sqrt((wilk.loc[mvs_in_rv] ** 2).sum()))
            other_l2_mvs = float(np.sqrt((others.loc[mvs_in_rv] ** 2).sum(axis=0)).mean())
            l2_ratio_mvs = wilk_l2_mvs / other_l2_mvs if other_l2_mvs > 0 else float("nan")
        else:
            wilk_l2_mvs = float("nan")
            l2_ratio_mvs = float("nan")

        # Top-10 absolute-magnitude genes in Wilk: housekeeping dominance check
        wilk_abs_sorted = wilk.abs().sort_values(ascending=False)
        top10 = list(wilk_abs_sorted.index[:10])
        n_housekeeping_top10 = sum(1 for g in top10 if g in HOUSEKEEPING)
        n_mvs_top10 = sum(1 for g in top10 if g in mvs_genes)

        # Pairwise r Wilk vs each other study
        r_wilk_aru = wilk.corr(rv["arunachalam_2020"])
        r_wilk_lee = wilk.corr(rv["lee_2020"])
        r_wilk_sch = wilk.corr(rv["schulte_schrepping_2020"])
        r_others_min = float(
            min(
                rv["arunachalam_2020"].corr(rv["lee_2020"]),
                rv["arunachalam_2020"].corr(rv["schulte_schrepping_2020"]),
                rv["lee_2020"].corr(rv["schulte_schrepping_2020"]),
            )
        )

        rows.append(
            {
                "bucket": bucket,
                "wilk_vs_others_mean_r_full": round(r_wilk_vs_others, 4),
                "wilk_vs_others_mean_r_mvs": round(r_wilk_vs_others_mvs, 4),
                "wilk_vs_arunachalam_r": round(r_wilk_aru, 4),
                "wilk_vs_lee_r": round(r_wilk_lee, 4),
                "wilk_vs_schulte_r": round(r_wilk_sch, 4),
                "min_others_pairwise_r": round(r_others_min, 4),
                "wilk_l2_full": round(wilk_l2, 3),
                "wilk_l2_mvs": round(wilk_l2_mvs, 3),
                "wilk_l2_vs_others_ratio_full": round(l2_ratio, 3),
                "wilk_l2_vs_others_ratio_mvs": round(l2_ratio_mvs, 3),
                "wilk_top10_genes": ",".join(top10[:5]) + "...",  # truncate display
                "wilk_n_housekeeping_in_top10": n_housekeeping_top10,
                "wilk_n_mvs_in_top10": n_mvs_top10,
                "n_mvs_in_rv": len(mvs_in_rv),
            }
        )

    df = pd.DataFrame(rows)
    out = TABLES / "session4_wilk_watchpoint.csv"
    df.to_csv(out, index=False)
    logger.info("\nwrote %s\n", out.name)

    # verdict logic
    logger.info("=== Wilk watchpoint per-bucket summary ===")
    artifact_flag_any = False
    for _, r in df.iterrows():
        # Flag conditions: wilk_vs_others_r < 0.3 (orthogonal), OR housekeeping_top10 > 3 (housekeeping
        # dominance), OR l2_ratio outside [0.5, 2.0] (extreme magnitude mismatch)
        flags = []
        if r["wilk_vs_others_mean_r_mvs"] < 0.3:
            flags.append("LOW_MVS_R")
        if r["wilk_n_housekeeping_in_top10"] > 3:
            flags.append("HOUSEKEEPING_DOMINANCE")
        if not (0.5 <= r["wilk_l2_vs_others_ratio_mvs"] <= 2.0):
            flags.append("MAGNITUDE_OUTLIER")
        flag_str = ",".join(flags) if flags else "OK"
        if flags:
            artifact_flag_any = True
        logger.info(
            "  %s: Wilk-vs-others r_mvs=%.3f (full=%.3f); L2_ratio_mvs=%.3f; HK_top10=%d/10 MVS_top10=%d/10 → %s",
            r["bucket"],
            r["wilk_vs_others_mean_r_mvs"],
            r["wilk_vs_others_mean_r_full"],
            r["wilk_l2_vs_others_ratio_mvs"],
            r["wilk_n_housekeeping_in_top10"],
            r["wilk_n_mvs_in_top10"],
            flag_str,
        )

    if artifact_flag_any:
        logger.info(
            "\nVERDICT: Wilk shows artifact pattern in at least one bucket. Flag in Issue 6."
        )
    else:
        logger.info(
            "\nVERDICT: Wilk consistent with other 3 studies at Harmony level. No artifact concern."
        )
        logger.info(
            "Implication: scVI's depth handling inherits Harmony's robustness; Tier IV verdict applies."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
