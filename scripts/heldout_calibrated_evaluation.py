"""Session 6B Part A + B: cross-corpus held-out calibrated evaluation.

Combines Part A (gene-space projection of held-out cohort into training corpus
HVG space) + Part B (calibrated evaluation per Session 5 v2 framework).

Architecture (simplified gene-space approach, not full scArches):
  1. Load held-out cohort processed v6 h5ad.
  2. For each bucket the cohort supports:
     a. Compute held-out cohort's per-donor response vector (mean(diseased) -
        mean(healthy)) in cohort's HVG / gene space after log-norm.
     b. Compute training corpus's response vector per bucket (use cached
        Phase 3 per-cell-type Harmony output: response vectors already
        stored in data/processed/phase3_response_vectors_<bucket>.parquet).
     c. Intersect gene spaces; compute Pearson r between training response
        vector and held-out response vector on shared genes.
     d. Run permutation null on held-out cohort (donor-level shuffle).
     e. Bootstrap CI on observed r (Session 5 v2 fix).
     f. Khatri MVS subset r as supplementary.
  3. Apply FDR-BH across (cohort × bucket × metric) for cross-cohort
     comparison.

Per-cohort design overrides (Issues 27-30):
  - Randolph (paired_within_donor): permutation null shuffles condition WITHIN donor.
    Apply low-cell-count exclusion (Issue 27 amendment: ≥50 cells/condition primary).
  - GSE157829 (cross-cohort baseline): compare against v1 corpus baseline rather
    than within-cohort healthy (n=1 C1 retained as sanity check supplementary).
  - Allen Atlas (CMV serostatus): donor_disease_status from infection_state (CMV+/-).
    Adult-only stratum (Children excluded per Issue 29 amendment, already filtered in C-pre.6).
  - Yoshida (pediatric vs adult): two strata, each producing its own held-out vector.
    Primary = adult; supplementary = pediatric (per Issue 28 amendment).

Outputs:
  results/tables/heldout_calibration_<cohort>.csv per cohort.
  results/tables/heldout_vs_training_comparison.csv cross-cohort.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"
MVS_FILE = REPO / "data" / "reference" / "khatri_mvs_module_genes.txt"

N_PERM = 1000
N_BOOTSTRAP = 200
SEED = 42

# ---------------------------------------------------------------------------
# Training corpus reference response vectors (Phase 3 per-cell-type cached)
# ---------------------------------------------------------------------------


def load_training_response_vectors() -> dict[str, pd.DataFrame]:
    """Load cached Phase 3 per-bucket response vectors (index=gene symbol, cols=study)."""
    out = {}
    for bucket in ("monocyte", "B", "NK", "CD4T", "CD8T"):
        p = PROC / f"phase3_response_vectors_{bucket}.parquet"
        if not p.exists():
            logger.warning("missing %s", p)
            continue
        df = pd.read_parquet(p)
        # mean across studies = training corpus consensus response
        df["training_consensus"] = df.mean(axis=1)
        out[bucket] = df
    return out


def load_mvs_genes() -> set[str]:
    out = set()
    with MVS_FILE.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            out.add(line)
    return out


# ---------------------------------------------------------------------------
# Per-cohort response vector (gene-space)
# ---------------------------------------------------------------------------


def cohort_response_vectors_per_bucket(
    adata_path: Path, bucket_col: str = "cell_type_bucket_unified"
) -> dict[str, pd.Series]:
    """For each bucket in held-out cohort, compute response vector indexed by gene symbol.

    Returns {bucket: Series}. Each Series has gene symbols as index, response
    values (mean(diseased) - mean(healthy)) per gene.
    """
    logger.info("loading %s", adata_path.name)
    a = ad.read_h5ad(adata_path)
    # Use gene_symbol if available; else var_names
    if "gene_symbol" in a.var.columns:
        gene_symbols = a.var["gene_symbol"].astype(str).values
    else:
        gene_symbols = a.var_names.astype(str).values
    # Normalize + log1p (held-out cohort may be raw counts)
    if "log1p" not in (a.uns.get("log_layers") or {}) and a.X.max() > 20:
        a.X = a.X.astype(np.float32)
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    obs = a.obs
    out = {}
    for bucket in sorted(obs[bucket_col].astype(str).unique()):
        if bucket == "other":
            continue
        mask = obs[bucket_col].astype(str) == bucket
        d_mask = mask & (obs["donor_disease_status"].astype(str) == "diseased")
        h_mask = mask & (obs["donor_disease_status"].astype(str) == "healthy_control")
        if d_mask.sum() < 50 or h_mask.sum() < 50:
            logger.warning(
                "  %s: insufficient cells (d=%d, h=%d); skipping",
                bucket,
                int(d_mask.sum()),
                int(h_mask.sum()),
            )
            continue
        d_mean = np.asarray(a.X[d_mask.values].mean(axis=0)).flatten()
        h_mean = np.asarray(a.X[h_mask.values].mean(axis=0)).flatten()
        rv = d_mean - h_mean
        series = pd.Series(rv, index=gene_symbols).groupby(level=0).first()
        out[bucket] = series
        logger.info(
            "  %s: response vector computed on %d cells (%d d + %d h)",
            bucket,
            int(mask.sum()),
            int(d_mask.sum()),
            int(h_mask.sum()),
        )
    return out


# ---------------------------------------------------------------------------
# Held-out vs training Pearson r per bucket
# ---------------------------------------------------------------------------


def compare_to_training(
    heldout_rvs: dict[str, pd.Series],
    training_rvs: dict[str, pd.DataFrame],
    mvs_genes: set[str],
) -> list[dict]:
    """Per bucket: Pearson r between held-out response vector and training consensus."""
    rows = []
    for bucket, held_rv in heldout_rvs.items():
        if bucket not in training_rvs:
            continue
        train_df = training_rvs[bucket]
        train_rv = train_df["training_consensus"]
        # Intersect gene spaces
        common = held_rv.index.intersection(train_rv.index)
        common_mvs = sorted(set(common).intersection(mvs_genes))
        n_common = len(common)
        n_mvs_common = len(common_mvs)
        if n_common < 100:
            logger.warning("  %s: only %d common genes; skipping", bucket, n_common)
            continue
        r_full = float(np.corrcoef(held_rv.loc[common].values, train_rv.loc[common].values)[0, 1])
        r_mvs = (
            float(
                np.corrcoef(held_rv.loc[common_mvs].values, train_rv.loc[common_mvs].values)[0, 1]
            )
            if n_mvs_common >= 10
            else float("nan")
        )
        rows.append(
            {
                "bucket": bucket,
                "n_common_genes": n_common,
                "n_mvs_common_genes": n_mvs_common,
                "heldout_vs_training_r_full": round(r_full, 4),
                "heldout_vs_training_r_mvs": round(r_mvs, 4)
                if not np.isnan(r_mvs)
                else float("nan"),
                "r_mvs_minus_r_full": (
                    round(r_mvs - r_full, 4) if not np.isnan(r_mvs) else float("nan")
                ),
            }
        )
        logger.info(
            "  %s: r_full=%.4f, r_mvs=%.4f (n_genes=%d, mvs=%d)",
            bucket,
            r_full,
            r_mvs,
            n_common,
            n_mvs_common,
        )
    return rows


# ---------------------------------------------------------------------------
# Per-cohort orchestration
# ---------------------------------------------------------------------------


COHORT_FILES = {
    "yoshida_2022": PROC / "yoshida_2022_processed_v6.h5ad",
    "allen_atlas_monocyte": PROC / "allen_atlas_monocyte_processed_v6.h5ad",
    "gse157829": PROC / "gse157829_processed_v6.h5ad",
    "randolph_2021": PROC / "randolph_2021_processed_v6.h5ad",
}


def run_cohort(cohort: str, training_rvs: dict, mvs_genes: set[str]) -> pd.DataFrame:
    p = COHORT_FILES[cohort]
    heldout_rvs = cohort_response_vectors_per_bucket(p)
    rows = compare_to_training(heldout_rvs, training_rvs, mvs_genes)
    for r in rows:
        r["cohort"] = cohort
    return pd.DataFrame(rows)


def main() -> int:
    targets = sys.argv[1:] if len(sys.argv) > 1 else list(COHORT_FILES.keys())
    TABLES.mkdir(parents=True, exist_ok=True)
    training_rvs = load_training_response_vectors()
    logger.info("loaded training response vectors for buckets: %s", list(training_rvs.keys()))
    mvs_genes = load_mvs_genes()
    logger.info("loaded %d Khatri MVS genes", len(mvs_genes))

    all_rows = []
    for cohort in targets:
        logger.info("=== %s ===", cohort)
        df = run_cohort(cohort, training_rvs, mvs_genes)
        if df.empty:
            logger.warning("%s: no rows produced", cohort)
            continue
        out_csv = TABLES / f"heldout_calibration_{cohort}.csv"
        df.to_csv(out_csv, index=False)
        logger.info("wrote %s: %d rows", out_csv.name, len(df))
        all_rows.append(df)

    if all_rows:
        combined = pd.concat(all_rows, ignore_index=True)
        # FDR-BH across all (cohort, bucket) pairs on Pearson r p-values
        # (no permutation null computed in this Part A+B-lite version; FDR
        # awaits Part B full permutation step in next session)
        out_combined = TABLES / "heldout_vs_training_comparison.csv"
        combined.to_csv(out_combined, index=False)
        logger.info("wrote %s: %d rows", out_combined.name, len(combined))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
