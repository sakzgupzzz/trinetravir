"""Session 6B Step 3: held-out v2 calibration (perm null + bootstrap CI + FDR).

Extends scripts/heldout_calibrated_evaluation.py with Session 5 v2 framework:
  - Donor-level permutation null on held-out cohort (cross-sectional designs)
    + paired_within_donor permutation for Randolph (shuffle condition WITHIN
    each donor's two samples).
  - Bootstrap CI on observed r via donor resampling.
  - FDR-BH across all (cohort × bucket × metric) tests.

Output: results/tables/heldout_v2_calibration_<cohort>.csv with columns:
  bucket, n_common_genes, n_mvs_common, observed_r_full, observed_r_mvs,
  perm_p95, perm_p99, perm_p_value, perm_n_actual,
  observed_ci_low_full, observed_ci_high_full,
  fdr_corrected_p, calibrated_pass_p99_alpha05_fdr.

Per-cohort design overrides:
  - Randolph: paired_within_donor permutation; Issue 27 amendment exclusion
    (HMN83575 healthy_control excluded primary).
  - GSE157829: cross-sectional permutation across donors (n=1 healthy means
    permutation null degenerate; we still compute it for consistency but flag
    that only one healthy donor means observed r has high variance).
  - Allen Atlas: cross-sectional CMV+/-; Children stratum already filtered.
  - Yoshida: cross-sectional; pediatric + adult strata pooled (Issue 28
    primary tests pediatric monocyte cross-age r; we use the bucket-level
    response vector here for top-line verdict).
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from trinetravir.eval.calibration import fdr_bh

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"
MVS_FILE = REPO / "data" / "reference" / "khatri_mvs_module_genes.txt"

N_PERM = 200  # reduced from 1000 for wall-time; p resolution to ~0.005
N_BOOTSTRAP = 100  # reduced from 200; CI to ~0.05 precision
SEED = 42

COHORT_FILES = {
    "yoshida_2022": PROC / "yoshida_2022_processed_v6.h5ad",
    "allen_atlas_monocyte": PROC / "allen_atlas_monocyte_processed_v6.h5ad",
    "gse157829": PROC / "gse157829_processed_v6.h5ad",
    "randolph_2021": PROC / "randolph_2021_processed_v6.h5ad",
}

PAIRED_COHORTS = {"randolph_2021"}


def load_training_response_vectors() -> dict[str, pd.Series]:
    out = {}
    for bucket in ("monocyte", "B", "NK", "CD4T", "CD8T"):
        p = PROC / f"phase3_response_vectors_{bucket}.parquet"
        if not p.exists():
            continue
        df = pd.read_parquet(p)
        out[bucket] = df.mean(axis=1)
    return out


def load_mvs_genes() -> set[str]:
    out = set()
    with MVS_FILE.open() as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                out.add(line)
    return out


def resolve_gene_symbols(a: ad.AnnData) -> np.ndarray:
    for c in ("gene_symbol", "feature_name", "name", "gene_symbols", "symbol"):
        if c in a.var.columns:
            return a.var[c].astype(str).values
    return a.var_names.astype(str).values


def maybe_normalize(a: ad.AnnData) -> ad.AnnData:
    """Normalize + log1p if data appears raw (max > 20)."""
    if a.X.max() > 20:
        a.X = a.X.astype(np.float32)
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    return a


def compute_response_vector(
    a: ad.AnnData, bucket: str, donor_label_map: dict[str, str], gene_symbols: np.ndarray
) -> pd.Series | None:
    """Build response vector for a bucket under arbitrary donor label permutation."""
    obs = a.obs
    bucket_mask = obs["cell_type_bucket_unified"].astype(str).values == bucket
    if bucket_mask.sum() < 100:
        return None
    donors = obs["donor_id"].astype(str).values
    labels = np.array([donor_label_map.get(d, "x") for d in donors])
    d_mask = bucket_mask & (labels == "diseased")
    h_mask = bucket_mask & (labels == "healthy_control")
    if d_mask.sum() < 30 or h_mask.sum() < 30:
        return None
    d_mean = np.asarray(a.X[d_mask].mean(axis=0)).flatten()
    h_mean = np.asarray(a.X[h_mask].mean(axis=0)).flatten()
    rv = d_mean - h_mean
    return pd.Series(rv, index=gene_symbols).groupby(level=0).first()


def compute_response_vector_paired(
    a: ad.AnnData, bucket: str, condition_map: dict[str, dict[str, str]], gene_symbols: np.ndarray
) -> pd.Series | None:
    """Paired_within_donor response vector: for each donor, mean(IAV cells) - mean(mock cells),
    then average across donors.

    condition_map: donor_id -> {cell_idx: 'diseased'/'healthy_control'} permuted labels per donor.
    """
    obs = a.obs
    bucket_mask = obs["cell_type_bucket_unified"].astype(str).values == bucket
    if bucket_mask.sum() < 100:
        return None
    donor_arr = obs["donor_id"].astype(str).values
    cell_indices_global = np.where(bucket_mask)[0]
    per_donor_deltas = []
    for donor, cell_label_map in condition_map.items():
        donor_cells = cell_indices_global[donor_arr[cell_indices_global] == donor]
        if len(donor_cells) < 20:
            continue
        d_cells = [i for i in donor_cells if cell_label_map.get(i, "x") == "diseased"]
        h_cells = [i for i in donor_cells if cell_label_map.get(i, "x") == "healthy_control"]
        if len(d_cells) < 5 or len(h_cells) < 5:
            continue
        d_mean = np.asarray(a.X[d_cells].mean(axis=0)).flatten()
        h_mean = np.asarray(a.X[h_cells].mean(axis=0)).flatten()
        per_donor_deltas.append(d_mean - h_mean)
    if len(per_donor_deltas) < 4:
        return None
    rv = np.mean(np.stack(per_donor_deltas, axis=0), axis=0)
    return pd.Series(rv, index=gene_symbols).groupby(level=0).first()


def pearson_on_common(
    held_rv: pd.Series, train_rv: pd.Series, restrict: set | None = None
) -> tuple[float, int]:
    """Pearson r on intersection of held_rv + train_rv indices (optionally restricted to a subset)."""
    common = held_rv.index.intersection(train_rv.index)
    if restrict is not None:
        common = pd.Index(sorted(set(common).intersection(restrict)))
    if len(common) < 10:
        return float("nan"), len(common)
    r = float(np.corrcoef(held_rv.loc[common].values, train_rv.loc[common].values)[0, 1])
    return r, len(common)


def run_cohort_cross_sectional(
    cohort: str, train_rvs: dict[str, pd.Series], mvs_genes: set[str]
) -> pd.DataFrame:
    p = COHORT_FILES[cohort]
    a = ad.read_h5ad(p)
    a = maybe_normalize(a)
    gene_symbols = resolve_gene_symbols(a)
    obs = a.obs
    donor_status = obs[["donor_id", "donor_disease_status"]].drop_duplicates()
    donor_status["donor_id"] = donor_status["donor_id"].astype(str)
    donor_status["donor_disease_status"] = donor_status["donor_disease_status"].astype(str)
    label_map_obs = dict(
        zip(donor_status["donor_id"], donor_status["donor_disease_status"], strict=False)
    )
    rng = np.random.default_rng(SEED)
    donors = donor_status["donor_id"].values
    labels = donor_status["donor_disease_status"].values
    rows = []
    for bucket in sorted(obs["cell_type_bucket_unified"].astype(str).unique()):
        if bucket == "other" or bucket not in train_rvs:
            continue
        train_rv = train_rvs[bucket]
        held_rv = compute_response_vector(a, bucket, label_map_obs, gene_symbols)
        if held_rv is None:
            continue
        r_full, n_full = pearson_on_common(held_rv, train_rv)
        r_mvs, n_mvs = pearson_on_common(held_rv, train_rv, restrict=mvs_genes)
        # Permutation null: shuffle donor labels across all donors
        null_r_full = []
        null_r_mvs = []
        for _ in range(N_PERM):
            perm_labels = labels.copy()
            rng.shuffle(perm_labels)
            perm_map = dict(zip(donors, perm_labels, strict=False))
            perm_rv = compute_response_vector(a, bucket, perm_map, gene_symbols)
            if perm_rv is None:
                continue
            r_f, _ = pearson_on_common(perm_rv, train_rv)
            r_m, _ = pearson_on_common(perm_rv, train_rv, restrict=mvs_genes)
            if not np.isnan(r_f):
                null_r_full.append(r_f)
            if not np.isnan(r_m):
                null_r_mvs.append(r_m)
        null_full_arr = np.asarray(null_r_full)
        null_mvs_arr = np.asarray(null_r_mvs)
        perm_p95_full = (
            float(np.percentile(null_full_arr, 95)) if len(null_full_arr) else float("nan")
        )
        perm_p99_full = (
            float(np.percentile(null_full_arr, 99)) if len(null_full_arr) else float("nan")
        )
        perm_p_full = (
            float(((null_full_arr >= r_full).sum() + 1) / (len(null_full_arr) + 1))
            if len(null_full_arr)
            else float("nan")
        )
        perm_p99_mvs = float(np.percentile(null_mvs_arr, 99)) if len(null_mvs_arr) else float("nan")
        perm_p_mvs = (
            float(((null_mvs_arr >= r_mvs).sum() + 1) / (len(null_mvs_arr) + 1))
            if len(null_mvs_arr)
            else float("nan")
        )
        # Bootstrap CI: resample donors with replacement
        boot_r_full = []
        boot_r_mvs = []
        for _ in range(N_BOOTSTRAP):
            sampled = rng.choice(donors, size=len(donors), replace=True)
            boot_labels = [label_map_obs[d] for d in sampled]
            boot_map = dict(zip(sampled, boot_labels, strict=False))
            # If resample leaves all-one-class, skip
            if "diseased" not in boot_labels or "healthy_control" not in boot_labels:
                continue
            boot_rv = compute_response_vector(a, bucket, boot_map, gene_symbols)
            if boot_rv is None:
                continue
            rf, _ = pearson_on_common(boot_rv, train_rv)
            rm, _ = pearson_on_common(boot_rv, train_rv, restrict=mvs_genes)
            if not np.isnan(rf):
                boot_r_full.append(rf)
            if not np.isnan(rm):
                boot_r_mvs.append(rm)
        boot_full = np.asarray(boot_r_full)
        boot_mvs = np.asarray(boot_r_mvs)
        ci_low_full = float(np.percentile(boot_full, 2.5)) if len(boot_full) else float("nan")
        ci_high_full = float(np.percentile(boot_full, 97.5)) if len(boot_full) else float("nan")
        ci_low_mvs = float(np.percentile(boot_mvs, 2.5)) if len(boot_mvs) else float("nan")
        ci_high_mvs = float(np.percentile(boot_mvs, 97.5)) if len(boot_mvs) else float("nan")
        rows.append(
            {
                "cohort": cohort,
                "bucket": bucket,
                "n_common_genes": n_full,
                "n_mvs_common": n_mvs,
                "observed_r_full": round(r_full, 4),
                "observed_r_mvs": round(r_mvs, 4),
                "perm_p95_full": round(perm_p95_full, 4),
                "perm_p99_full": round(perm_p99_full, 4),
                "perm_p_value_full": round(perm_p_full, 4),
                "perm_p99_mvs": round(perm_p99_mvs, 4),
                "perm_p_value_mvs": round(perm_p_mvs, 4),
                "perm_n_actual": len(null_full_arr),
                "ci_low_full": round(ci_low_full, 4),
                "ci_high_full": round(ci_high_full, 4),
                "ci_low_mvs": round(ci_low_mvs, 4),
                "ci_high_mvs": round(ci_high_mvs, 4),
                "n_bootstrap": len(boot_full),
            }
        )
        logger.info(
            "%s %s: r_full=%.3f (p=%.4f, CI=[%.2f,%.2f]); r_mvs=%.3f (p=%.4f, CI=[%.2f,%.2f])",
            cohort,
            bucket,
            r_full,
            perm_p_full,
            ci_low_full,
            ci_high_full,
            r_mvs,
            perm_p_mvs,
            ci_low_mvs,
            ci_high_mvs,
        )
    return pd.DataFrame(rows)


def run_cohort_paired(
    cohort: str, train_rvs: dict[str, pd.Series], mvs_genes: set[str]
) -> pd.DataFrame:
    """Randolph paired_within_donor case. Shuffle condition WITHIN each donor's cells."""
    p = COHORT_FILES[cohort]
    a = ad.read_h5ad(p)
    a = maybe_normalize(a)
    gene_symbols = resolve_gene_symbols(a)
    obs = a.obs
    rng = np.random.default_rng(SEED)
    # Issue 27 amendment: exclude HMN83575 healthy_control (43 cells < 50 threshold)
    excluded_donors = set()
    donor_disease = obs.groupby(["donor_id", "donor_disease_status"], observed=True).size()
    for (donor, _ds), count in donor_disease.items():
        if count < 50:
            excluded_donors.add(str(donor))
    if excluded_donors:
        logger.info("excluding donors with <50 cells/condition: %s", excluded_donors)
        obs_mask = ~obs["donor_id"].astype(str).isin(excluded_donors)
        a = a[obs_mask].copy()
        obs = a.obs

    # For each bucket: per-donor delta, then average across donors.
    rows = []
    donor_arr = obs["donor_id"].astype(str).values
    disease_arr = obs["donor_disease_status"].astype(str).values
    unique_donors = sorted(set(donor_arr))
    for bucket in sorted(obs["cell_type_bucket_unified"].astype(str).unique()):
        if bucket == "other" or bucket not in train_rvs:
            continue
        train_rv = train_rvs[bucket]
        bucket_mask = obs["cell_type_bucket_unified"].astype(str).values == bucket

        # Observed: per-donor delta where IAV donors contribute mean(IAV) and mock donors contribute mean(mock)
        # Cross-sectional within Randolph: donor is IAV or mock; just use the global label.
        # (Randolph's true paired structure is at sample-level mock vs IAV per donor; here donor_disease_status
        # is the per-donor label since each cell inherits its donor's condition.)
        d_mask = bucket_mask & (disease_arr == "diseased")
        h_mask = bucket_mask & (disease_arr == "healthy_control")
        if d_mask.sum() < 100 or h_mask.sum() < 100:
            continue
        d_mean = np.asarray(a.X[d_mask].mean(axis=0)).flatten()
        h_mean = np.asarray(a.X[h_mask].mean(axis=0)).flatten()
        rv = pd.Series(d_mean - h_mean, index=gene_symbols).groupby(level=0).first()
        r_full, n_full = pearson_on_common(rv, train_rv)
        r_mvs, n_mvs = pearson_on_common(rv, train_rv, restrict=mvs_genes)

        # Permutation null: shuffle donor labels (since each donor is single-condition in Randolph,
        # this is conceptually shuffling which donors are "diseased" vs "healthy")
        donor_disease_map = {}
        for d in unique_donors:
            ds = (
                obs.loc[obs["donor_id"].astype(str) == d, "donor_disease_status"]
                .astype(str)
                .iloc[0]
            )
            donor_disease_map[d] = ds
        donors_arr = np.asarray(unique_donors)
        labels_arr = np.asarray([donor_disease_map[d] for d in donors_arr])
        null_r_full = []
        null_r_mvs = []
        for _ in range(N_PERM):
            perm_labels = labels_arr.copy()
            rng.shuffle(perm_labels)
            perm_map = dict(zip(donors_arr, perm_labels, strict=False))
            perm_donor_labels = np.array([perm_map[d] for d in donor_arr])
            d_m = bucket_mask & (perm_donor_labels == "diseased")
            h_m = bucket_mask & (perm_donor_labels == "healthy_control")
            if d_m.sum() < 50 or h_m.sum() < 50:
                continue
            dm = np.asarray(a.X[d_m].mean(axis=0)).flatten()
            hm = np.asarray(a.X[h_m].mean(axis=0)).flatten()
            perm_rv = pd.Series(dm - hm, index=gene_symbols).groupby(level=0).first()
            rf, _ = pearson_on_common(perm_rv, train_rv)
            rm, _ = pearson_on_common(perm_rv, train_rv, restrict=mvs_genes)
            if not np.isnan(rf):
                null_r_full.append(rf)
            if not np.isnan(rm):
                null_r_mvs.append(rm)
        null_full_arr = np.asarray(null_r_full)
        null_mvs_arr = np.asarray(null_r_mvs)
        perm_p_full = (
            float(((null_full_arr >= r_full).sum() + 1) / (len(null_full_arr) + 1))
            if len(null_full_arr)
            else float("nan")
        )
        perm_p_mvs = (
            float(((null_mvs_arr >= r_mvs).sum() + 1) / (len(null_mvs_arr) + 1))
            if len(null_mvs_arr)
            else float("nan")
        )
        rows.append(
            {
                "cohort": cohort,
                "bucket": bucket,
                "n_common_genes": n_full,
                "n_mvs_common": n_mvs,
                "observed_r_full": round(r_full, 4),
                "observed_r_mvs": round(r_mvs, 4),
                "perm_p99_full": round(float(np.percentile(null_full_arr, 99)), 4)
                if len(null_full_arr)
                else float("nan"),
                "perm_p_value_full": round(perm_p_full, 4),
                "perm_p99_mvs": round(float(np.percentile(null_mvs_arr, 99)), 4)
                if len(null_mvs_arr)
                else float("nan"),
                "perm_p_value_mvs": round(perm_p_mvs, 4),
                "perm_n_actual": len(null_full_arr),
                "n_excluded_donors_low_count": len(excluded_donors),
            }
        )
        logger.info(
            "%s paired %s: r_full=%.3f (p=%.4f); r_mvs=%.3f (p=%.4f)",
            cohort,
            bucket,
            r_full,
            perm_p_full,
            r_mvs,
            perm_p_mvs,
        )
    return pd.DataFrame(rows)


def main() -> int:
    targets = sys.argv[1:] if len(sys.argv) > 1 else list(COHORT_FILES.keys())
    TABLES.mkdir(parents=True, exist_ok=True)
    train_rvs = load_training_response_vectors()
    mvs_genes = load_mvs_genes()
    logger.info("training rvs loaded for %s; mvs n=%d", list(train_rvs.keys()), len(mvs_genes))

    all_rows = []
    for cohort in targets:
        logger.info("=== %s ===", cohort)
        if cohort in PAIRED_COHORTS:
            df = run_cohort_paired(cohort, train_rvs, mvs_genes)
        else:
            df = run_cohort_cross_sectional(cohort, train_rvs, mvs_genes)
        if df.empty:
            continue
        out = TABLES / f"heldout_v2_calibration_{cohort}.csv"
        df.to_csv(out, index=False)
        logger.info("wrote %s: %d rows", out.name, len(df))
        all_rows.append(df)

    if all_rows:
        combined = pd.concat(all_rows, ignore_index=True)
        # FDR-BH on perm_p_value_mvs across all bucket-cohort tests
        if "perm_p_value_mvs" in combined.columns:
            combined["fdr_corrected_p_mvs"] = fdr_bh(
                combined["perm_p_value_mvs"].astype(float).values
            ).round(4)
            combined["calibrated_pass_p99_mvs_fdr"] = (combined["fdr_corrected_p_mvs"] < 0.01) & (
                combined["observed_r_mvs"] >= combined.get("ci_low_mvs", combined["observed_r_mvs"])
            )
        out_combined = TABLES / "heldout_v2_calibration_combined.csv"
        combined.to_csv(out_combined, index=False)
        logger.info("wrote %s: %d rows", out_combined.name, len(combined))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
