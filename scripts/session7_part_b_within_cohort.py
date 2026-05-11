"""Session 7 Part B: within-cohort-only sensitivity analysis.

Per Issue 33 (METHODS_CHOICES.md, pre-spec 2026-05-11).

For each v1 cohort (Wilk, Lee, Arunachalam, Schulte-Schrepping):
  1. Load raw normalized counts (no Harmony) from per-study reannotated h5ad.
  2. For each bucket: per-donor pseudobulk (diseased - healthy) → per-bucket
     response vector in gene space.
  3. Pairwise Pearson r BETWEEN BUCKETS within-cohort (bucket_i × bucket_j).
  4. v2 framework: permutation null (N=500, shuffle donor labels within cohort),
     bootstrap CI (N=200, resample donors with replacement).
  5. Write `results/tables/sensitivity_within_cohort.csv`.

Then aggregate vs cross-study harmonized:
  - Cross-study harmonized bucket-pair r: from cached
    phase3_response_vectors_<bucket>.parquet, take mean across studies per
    bucket, then pairwise r between bucket means.
  - Within-cohort mean bucket-pair r: average across 4 cohorts.
  - Sign concordance + magnitude alignment per bucket pair.
  - Write `results/tables/sensitivity_within_vs_cross.csv`.

Apply Issue 33 mechanical decision rule:
  Sign concordance ≥80% with magnitude divergence ≤0.20 → biology consistent
  Sign concordance 50-80%, magnitude divergence 0.20-0.50 → partial
  Sign concordance <50% OR systematic magnitude reversal → artifact
"""

from __future__ import annotations

import logging
import warnings
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"
MVS_FILE = REPO / "data" / "reference" / "khatri_mvs_module_genes.txt"

BUCKETS = ("monocyte", "B", "NK", "CD4T", "CD8T")
STUDY_FILE = {
    "wilk_2020": PROC / "wilk_2020_reannotated.h5ad",
    "lee_2020": PROC / "lee_2020_reannotated.h5ad",
    "arunachalam_2020": PROC / "arunachalam_2020_reannotated.h5ad",
    "schulte_schrepping_2020": PROC / "schulte_schrepping_2020_reannotated.h5ad",
}
N_PERM = 500
N_BOOTSTRAP = 200
SEED = 42


def load_mvs_genes() -> set[str]:
    out = set()
    with MVS_FILE.open() as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                out.add(s)
    return out


def resolve_symbols(a: ad.AnnData) -> np.ndarray:
    for c in ("gene_symbol", "feature_name", "name", "gene_symbols", "symbol"):
        if c in a.var.columns:
            return a.var[c].astype(str).values
    return a.var_names.astype(str).values


def harmony_donor_lookup(bucket: str) -> pd.DataFrame:
    """Per-study, per-cell donor_id + donor_disease_status from harmony h5ad obs.

    Returns DataFrame indexed by (study_id, row_idx) (row_idx = cell idx in study h5ad).
    """
    h = ad.read_h5ad(PROC / f"harmony_per_celltype_{bucket}.h5ad", backed="r")
    obs = h.obs[["study_id", "donor_id", "donor_disease_status"]].copy()
    obs.index = obs.index.astype(str)
    studies = obs["study_id"].astype(str).values
    row_idx = np.array(
        [int(idx.replace(f"-{s}", "")) for idx, s in zip(obs.index, studies, strict=False)]
    )
    obs["row_idx"] = row_idx
    return obs


def per_bucket_response_vectors_within_cohort(
    cohort: str, hvg_dict: dict[str, list[str]]
) -> tuple[dict[str, pd.Series], dict[str, list[str]], dict[str, list[str]]]:
    """For one cohort: per-bucket response vector + donor lists.

    Returns:
      bucket_rv: {bucket: gene-indexed response vector}
      donor_means_d: {bucket: list of (donor_id, gene_vector) for diseased}
      donor_means_h: {bucket: list of (donor_id, gene_vector) for healthy}
    """
    a = ad.read_h5ad(STUDY_FILE[cohort])
    sym = resolve_symbols(a)
    a.var_names = sym
    a.var_names_make_unique()
    # normalize once on full study, then slice per bucket
    a.X = a.X.astype(np.float32)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    bucket_rv = {}
    donor_means_d_all = {}
    donor_means_h_all = {}
    for bucket in BUCKETS:
        hvg = hvg_dict[bucket]
        donor_obs = harmony_donor_lookup(bucket)
        donor_cohort = donor_obs[donor_obs["study_id"] == cohort]
        if donor_cohort.empty:
            continue
        cell_idx = donor_cohort["row_idx"].values
        sub = a[cell_idx, :].copy()
        # restrict to HVG present
        hvg_in = [g for g in hvg if g in set(sub.var_names)]
        sub = sub[:, hvg_in]
        sub.obs["donor_id_h"] = donor_cohort["donor_id"].astype(str).values
        sub.obs["donor_disease_status_h"] = donor_cohort["donor_disease_status"].astype(str).values
        X = np.asarray(sub.X.todense() if hasattr(sub.X, "todense") else sub.X)
        donor_means_d = []
        donor_means_h = []
        donor_df = pd.DataFrame(
            {
                "donor_id": sub.obs["donor_id_h"].values,
                "donor_disease_status": sub.obs["donor_disease_status_h"].values,
            }
        )
        for (donor, ds), grp_idx in donor_df.groupby(
            ["donor_id", "donor_disease_status"], observed=True
        ).groups.items():
            mean_vec = X[list(grp_idx)].mean(axis=0)
            if ds == "diseased":
                donor_means_d.append((donor, mean_vec, hvg_in))
            elif ds == "healthy_control":
                donor_means_h.append((donor, mean_vec, hvg_in))
        if not donor_means_d or not donor_means_h:
            continue
        d_mean = np.mean(np.stack([v for _, v, _ in donor_means_d]), axis=0)
        h_mean = np.mean(np.stack([v for _, v, _ in donor_means_h]), axis=0)
        rv = pd.Series(d_mean - h_mean, index=hvg_in)
        rv = rv.reindex(hvg)
        bucket_rv[bucket] = rv
        donor_means_d_all[bucket] = donor_means_d
        donor_means_h_all[bucket] = donor_means_h
        logger.info(
            "  %s %s: rv (%d d donors, %d h donors, %d HVG)",
            cohort,
            bucket,
            len(donor_means_d),
            len(donor_means_h),
            len(hvg_in),
        )
    return bucket_rv, donor_means_d_all, donor_means_h_all


def pearson_on_common(a: pd.Series, b: pd.Series, restrict: set | None = None) -> float:
    common = a.index.intersection(b.index)
    if restrict is not None:
        common = pd.Index(sorted(set(common).intersection(restrict)))
    av = a.loc[common].dropna()
    bv = b.loc[common].dropna()
    both = av.index.intersection(bv.index)
    if len(both) < 10:
        return float("nan")
    return float(np.corrcoef(av.loc[both].values, bv.loc[both].values)[0, 1])


def perm_null_bucket_pair(
    donor_means_d_a: list,
    donor_means_h_a: list,
    donor_means_d_b: list,
    donor_means_h_b: list,
    hvg_a: list[str],
    hvg_b: list[str],
    restrict: set | None = None,
    n_perm: int = N_PERM,
    seed: int = SEED,
) -> tuple[float, np.ndarray]:
    """Permutation null: shuffle donor labels within each bucket; recompute r between buckets."""
    rng = np.random.default_rng(seed)

    def make_rv(d_list, h_list, hvg):
        d_means = np.stack([v for _, v, _ in d_list])
        h_means = np.stack([v for _, v, _ in h_list])
        d_mean = d_means.mean(axis=0)
        h_mean = h_means.mean(axis=0)
        return pd.Series(d_mean - h_mean, index=hvg)

    observed_r = pearson_on_common(
        make_rv(donor_means_d_a, donor_means_h_a, hvg_a),
        make_rv(donor_means_d_b, donor_means_h_b, hvg_b),
        restrict=restrict,
    )
    # For permutation: combine all donors from both conditions within each bucket,
    # then shuffle disease labels.
    null_r = []
    for _ in range(n_perm):
        # Bucket A: shuffle disease label across A's donor pool
        all_a = donor_means_d_a + donor_means_h_a
        n_d_a = len(donor_means_d_a)
        rng.shuffle(all_a)
        perm_d_a = all_a[:n_d_a]
        perm_h_a = all_a[n_d_a:]
        rv_a = make_rv(perm_d_a, perm_h_a, hvg_a)
        # Bucket B: same
        all_b = donor_means_d_b + donor_means_h_b
        n_d_b = len(donor_means_d_b)
        rng.shuffle(all_b)
        perm_d_b = all_b[:n_d_b]
        perm_h_b = all_b[n_d_b:]
        rv_b = make_rv(perm_d_b, perm_h_b, hvg_b)
        r = pearson_on_common(rv_a, rv_b, restrict=restrict)
        if not np.isnan(r):
            null_r.append(r)
    return observed_r, np.asarray(null_r)


def bootstrap_bucket_pair(
    donor_means_d_a: list,
    donor_means_h_a: list,
    donor_means_d_b: list,
    donor_means_h_b: list,
    hvg_a: list[str],
    hvg_b: list[str],
    restrict: set | None = None,
    n_boot: int = N_BOOTSTRAP,
    seed: int = SEED + 1,
) -> np.ndarray:
    """Bootstrap CI: resample donors per condition per bucket with replacement; recompute r."""
    rng = np.random.default_rng(seed)

    def make_rv(d_list, h_list, hvg):
        d_means = np.stack([v for _, v, _ in d_list])
        h_means = np.stack([v for _, v, _ in h_list])
        return pd.Series(d_means.mean(axis=0) - h_means.mean(axis=0), index=hvg)

    boot_r = []
    for _ in range(n_boot):
        b_d_a = [
            donor_means_d_a[i] for i in rng.integers(0, len(donor_means_d_a), len(donor_means_d_a))
        ]
        b_h_a = [
            donor_means_h_a[i] for i in rng.integers(0, len(donor_means_h_a), len(donor_means_h_a))
        ]
        b_d_b = [
            donor_means_d_b[i] for i in rng.integers(0, len(donor_means_d_b), len(donor_means_d_b))
        ]
        b_h_b = [
            donor_means_h_b[i] for i in rng.integers(0, len(donor_means_h_b), len(donor_means_h_b))
        ]
        rv_a = make_rv(b_d_a, b_h_a, hvg_a)
        rv_b = make_rv(b_d_b, b_h_b, hvg_b)
        r = pearson_on_common(rv_a, rv_b, restrict=restrict)
        if not np.isnan(r):
            boot_r.append(r)
    return np.asarray(boot_r)


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    mvs_genes = load_mvs_genes()
    # Load HVG per bucket from harmony h5ads
    hvg_dict = {}
    for bucket in BUCKETS:
        h = ad.read_h5ad(PROC / f"harmony_per_celltype_{bucket}.h5ad", backed="r")
        hvg_dict[bucket] = [str(g) for g in h.uns["hvg_genes"]]

    rows = []
    for cohort in STUDY_FILE:
        logger.info("=== %s ===", cohort)
        bucket_rv, dm_d, dm_h = per_bucket_response_vectors_within_cohort(cohort, hvg_dict)
        # Pairwise bucket combinations
        for b_a, b_b in combinations(sorted(bucket_rv.keys()), 2):
            for gene_set, restrict in (("full", None), ("MVS", mvs_genes)):
                obs_r, null = perm_null_bucket_pair(
                    dm_d[b_a],
                    dm_h[b_a],
                    dm_d[b_b],
                    dm_h[b_b],
                    hvg_dict[b_a],
                    hvg_dict[b_b],
                    restrict=restrict,
                )
                if np.isnan(obs_r):
                    continue
                p = (
                    float(((null >= obs_r).sum() + 1) / (len(null) + 1))
                    if len(null)
                    else float("nan")
                )
                boot = bootstrap_bucket_pair(
                    dm_d[b_a],
                    dm_h[b_a],
                    dm_d[b_b],
                    dm_h[b_b],
                    hvg_dict[b_a],
                    hvg_dict[b_b],
                    restrict=restrict,
                )
                ci_low = float(np.percentile(boot, 2.5)) if len(boot) else float("nan")
                ci_high = float(np.percentile(boot, 97.5)) if len(boot) else float("nan")
                rows.append(
                    {
                        "cohort": cohort,
                        "bucket_pair": f"{b_a}_vs_{b_b}",
                        "gene_set": gene_set,
                        "observed_r": round(obs_r, 4),
                        "perm_p_raw": round(p, 4),
                        "ci_low": round(ci_low, 4),
                        "ci_high": round(ci_high, 4),
                        "n_perm": len(null),
                        "n_bootstrap": len(boot),
                    }
                )
                logger.info(
                    "  %s %s %s: r=%.4f p=%.4f CI=[%.3f, %.3f]",
                    cohort,
                    f"{b_a}_vs_{b_b}",
                    gene_set,
                    obs_r,
                    p,
                    ci_low,
                    ci_high,
                )

    df = pd.DataFrame(rows)
    out = TABLES / "sensitivity_within_cohort.csv"
    df.to_csv(out, index=False)
    logger.info("wrote %s: %d rows", out.name, len(df))

    # Aggregate vs cross-study harmonized.
    # Cross-study harmonized bucket-pair r: load cached parquets, average per
    # bucket across studies (training_consensus), pairwise r between buckets.
    train_rv = {}
    for bucket in BUCKETS:
        p = PROC / f"phase3_response_vectors_{bucket}.parquet"
        if not p.exists():
            continue
        d = pd.read_parquet(p)
        train_rv[bucket] = d.mean(axis=1)

    agg_rows = []
    for gene_set, restrict in (("full", None), ("MVS", mvs_genes)):
        for b_a, b_b in combinations(sorted(BUCKETS), 2):
            if b_a not in train_rv or b_b not in train_rv:
                continue
            cross_r = pearson_on_common(train_rv[b_a], train_rv[b_b], restrict=restrict)
            pair_str = f"{b_a}_vs_{b_b}"
            within = df[(df["bucket_pair"] == pair_str) & (df["gene_set"] == gene_set)]
            if within.empty:
                continue
            within_rs = within["observed_r"].astype(float).values
            mean_within = float(np.mean(within_rs))
            # Sign concordance: fraction of cohorts where sign matches cross_r sign
            if not np.isnan(cross_r):
                sign_cross = np.sign(cross_r)
                sign_within = np.sign(within_rs)
                concordance = float(np.mean(sign_within == sign_cross))
            else:
                concordance = float("nan")
            magnitude = float(np.mean(np.abs(within_rs - cross_r)))
            agg_rows.append(
                {
                    "bucket_pair": pair_str,
                    "gene_set": gene_set,
                    "n_cohorts_observed": len(within_rs),
                    "mean_within_cohort_r": round(mean_within, 4),
                    "cross_study_harmonized_r": round(cross_r, 4)
                    if not np.isnan(cross_r)
                    else float("nan"),
                    "sign_concordance": round(concordance, 4),
                    "magnitude_alignment": round(magnitude, 4),
                }
            )

    agg_df = pd.DataFrame(agg_rows)
    out2 = TABLES / "sensitivity_within_vs_cross.csv"
    agg_df.to_csv(out2, index=False)
    logger.info("wrote %s: %d rows", out2.name, len(agg_df))
    print("\n=== Aggregate (Issue 33) ===\n")
    print(agg_df.to_string(index=False))

    # Issue 33 mechanical verdict
    print("\n=== Issue 33 verdict (mechanical, per gene_set) ===")
    for gene_set in ("full", "MVS"):
        sub = agg_df[agg_df["gene_set"] == gene_set]
        if sub.empty:
            continue
        mean_concordance = float(sub["sign_concordance"].mean())
        mean_magnitude = float(sub["magnitude_alignment"].mean())
        if mean_concordance >= 0.80 and mean_magnitude <= 0.20:
            verdict = "BIOLOGY_CONSISTENT (sign concordance ≥80%, mag div ≤0.20)"
        elif mean_concordance >= 0.50 and mean_magnitude <= 0.50:
            verdict = "PARTIAL_ALIGNMENT (sign 50-80%, mag 0.20-0.50)"
        else:
            verdict = "ARTIFACT (sign <50% OR systematic magnitude reversal)"
        print(
            f"  {gene_set}: mean_concordance={mean_concordance:.3f}, "
            f"mean_magnitude_alignment={mean_magnitude:.3f} → {verdict}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
