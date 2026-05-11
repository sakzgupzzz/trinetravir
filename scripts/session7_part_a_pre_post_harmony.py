"""Session 7 Part A: pre/post-Harmony response vector Δr quantification.

Per Issue 32 (METHODS_CHOICES.md, pre-spec 2026-05-11) and
references/session_7_prompt.md.

For each v1 bucket (monocyte, B, NK, CD4T, CD8T):
  1. Load `harmony_per_celltype_<bucket>.h5ad` to extract:
     - HVG gene list (4000 genes, Harmony's input space)
     - Cell membership: (cell_id, study_id, donor_id, donor_disease_status)
  2. For pre-Harmony: re-extract those exact cells from per-study reannotated
     h5ads, normalize_total + log1p, restrict to bucket's HVG, compute
     per-donor pseudobulk, then per-study response vector (mean diseased -
     mean healthy across donors).
  3. For post-Harmony: load cached `phase3_response_vectors_<bucket>.parquet`
     (gene-space x_corrected mean per study × condition difference).
  4. For each gene_set in (full HVG, MVS subset):
     - Compute mean off-diagonal cross-study Pearson r for pre-Harmony
       response-vector matrix (genes × studies).
     - Same for post-Harmony.
     - Δr = r_post − r_pre.
  5. Write `results/tables/sensitivity_pre_post_harmony.csv` with rows
     (bucket, gene_set, r_pre, r_post, delta_r, n_studies, n_donors_total).

Apply Issue 32 mechanical decision rule to Δr per bucket per gene_set:
  Δr ≤ 0.10 → biology dominant
  Δr ∈ (0.10, 0.30] → mixed
  Δr > 0.30 → Harmony dominant
"""

from __future__ import annotations

import logging
import warnings
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


def compute_pre_harmony_rv_per_study(
    bucket: str, harmony_obs: pd.DataFrame, hvg_genes: list[str]
) -> pd.DataFrame:
    """Per-study pre-Harmony response vector on bucket's HVG.

    Returns DataFrame indexed by gene symbol; one column per study; values are
    mean(diseased_donors) - mean(healthy_donors) per gene after normalize_total
    + log1p on raw counts.
    """
    studies = sorted(harmony_obs["study_id"].unique())
    series_per_study: dict[str, pd.Series] = {}
    for study in studies:
        # Harmony obs index format: '<int>-<study_id>'. Strip suffix → positional row idx.
        h_study_obs = harmony_obs.loc[harmony_obs["study_id"] == study].copy()
        suffix = f"-{study}"
        h_study_obs["row_idx"] = (
            h_study_obs.index.astype(str).str.replace(suffix, "", regex=False).astype(int)
        )
        cell_idx_keep = h_study_obs["row_idx"].values
        a = ad.read_h5ad(STUDY_FILE[study])
        # Resolve gene symbols from this study's var.
        sym = resolve_symbols(a)
        a.var_names = sym
        a.var_names_make_unique()
        # Filter to bucket cells that survived to harmony (positional index).
        if len(cell_idx_keep) == 0:
            logger.warning("  %s %s: no cells overlap; skipping", study, bucket)
            continue
        a = a[cell_idx_keep].copy()
        # Bring donor_disease_status / donor_id from harmony_obs.
        # h_study_obs aligned with cell_idx_keep order via .loc on row_idx-sorted.
        h_aligned = h_study_obs.set_index("row_idx").loc[cell_idx_keep]
        a.obs["donor_id_h"] = h_aligned["donor_id"].astype(str).values
        a.obs["donor_disease_status_h"] = h_aligned["donor_disease_status"].astype(str).values
        # normalize + log1p (raw counts)
        a.X = a.X.astype(np.float32)
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        # restrict to HVG (intersect)
        hvg_set = [g for g in hvg_genes if g in set(a.var_names)]
        if not hvg_set:
            logger.warning("  %s %s: no HVG overlap", study, bucket)
            continue
        a = a[:, hvg_set]
        # per-donor mean per condition
        donor_df = pd.DataFrame(
            {
                "donor_id": a.obs["donor_id_h"].astype(str).values,
                "donor_disease_status": a.obs["donor_disease_status_h"].astype(str).values,
            }
        )
        X = np.asarray(a.X.todense() if hasattr(a.X, "todense") else a.X)
        donor_means: dict[tuple[str, str], np.ndarray] = {}
        for (donor, ds), grp_idx in donor_df.groupby(
            ["donor_id", "donor_disease_status"], observed=True
        ).groups.items():
            donor_means[(donor, ds)] = X[list(grp_idx)].mean(axis=0)
        d_donors = [k for k in donor_means if k[1] == "diseased"]
        h_donors = [k for k in donor_means if k[1] == "healthy_control"]
        if not d_donors or not h_donors:
            logger.warning("  %s %s: missing one condition; skipping", study, bucket)
            continue
        d_mean = np.mean(np.stack([donor_means[k] for k in d_donors]), axis=0)
        h_mean = np.mean(np.stack([donor_means[k] for k in h_donors]), axis=0)
        rv = pd.Series(d_mean - h_mean, index=hvg_set)
        # Ensure series indexed on canonical HVG order (NaN for missing)
        rv = rv.reindex(hvg_genes)
        series_per_study[study] = rv
        logger.info(
            "  %s %s: pre-Harmony rv done (%d donors d, %d donors h, %d HVG)",
            study,
            bucket,
            len(d_donors),
            len(h_donors),
            int((~rv.isna()).sum()),
        )
    return pd.DataFrame(series_per_study)


def mean_off_diag_pearson(df: pd.DataFrame, restrict: set | None = None) -> tuple[float, int]:
    """Mean off-diagonal cross-study Pearson r over a gene index optionally restricted to a subset."""
    if restrict is not None:
        df = df.loc[df.index.intersection(pd.Index(sorted(restrict)))]
    df = df.dropna(how="any")
    if len(df) < 10 or df.shape[1] < 2:
        return float("nan"), len(df)
    corr = df.corr().values  # studies × studies
    n = corr.shape[0]
    off = corr[~np.eye(n, dtype=bool)]
    return float(np.nanmean(off)), len(df)


def run() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    mvs_genes = load_mvs_genes()
    logger.info("MVS gene set: %d", len(mvs_genes))
    rows = []
    for bucket in BUCKETS:
        h5 = PROC / f"harmony_per_celltype_{bucket}.h5ad"
        if not h5.exists():
            logger.warning("missing %s; skip", h5.name)
            continue
        a = ad.read_h5ad(h5)
        hvg_genes = [str(g) for g in a.uns["hvg_genes"]]
        harmony_obs = a.obs[["study_id", "donor_id", "donor_disease_status"]].copy()
        harmony_obs.index = harmony_obs.index.astype(str)
        logger.info("=== %s ===", bucket)
        # Pre-Harmony per-study response vectors
        pre_df = compute_pre_harmony_rv_per_study(bucket, harmony_obs, hvg_genes)
        # Post-Harmony from cached parquet
        post_path = PROC / f"phase3_response_vectors_{bucket}.parquet"
        if not post_path.exists():
            logger.warning("missing %s; skip", post_path.name)
            continue
        post_df = pd.read_parquet(post_path)
        # Align to same HVG and same study columns
        common_studies = sorted(set(pre_df.columns).intersection(set(post_df.columns)))
        pre_df = pre_df[common_studies]
        post_df = post_df[common_studies]
        # Index alignment: use HVG list from harmony (canonical order)
        pre_df = pre_df.reindex(hvg_genes)
        post_df = post_df.reindex(hvg_genes)
        n_donors_total = harmony_obs[["study_id", "donor_id"]].drop_duplicates().shape[0]
        # Full HVG
        r_pre_full, n_full = mean_off_diag_pearson(pre_df)
        r_post_full, _ = mean_off_diag_pearson(post_df)
        # MVS subset
        r_pre_mvs, n_mvs = mean_off_diag_pearson(pre_df, restrict=mvs_genes)
        r_post_mvs, _ = mean_off_diag_pearson(post_df, restrict=mvs_genes)
        for gene_set, r_pre, r_post, n_genes in (
            ("full", r_pre_full, r_post_full, n_full),
            ("MVS", r_pre_mvs, r_post_mvs, n_mvs),
        ):
            delta_r = (
                round(r_post - r_pre, 4)
                if not (np.isnan(r_pre) or np.isnan(r_post))
                else float("nan")
            )
            rows.append(
                {
                    "bucket": bucket,
                    "gene_set": gene_set,
                    "n_genes": n_genes,
                    "n_studies": len(common_studies),
                    "n_donors_total": n_donors_total,
                    "r_pre": round(r_pre, 4) if not np.isnan(r_pre) else float("nan"),
                    "r_post": round(r_post, 4) if not np.isnan(r_post) else float("nan"),
                    "delta_r": delta_r,
                }
            )
            logger.info(
                "  %s %s: r_pre=%.4f r_post=%.4f Δr=%.4f (n_genes=%d, n_studies=%d)",
                bucket,
                gene_set,
                r_pre,
                r_post,
                delta_r
                if not isinstance(delta_r, float) or not np.isnan(delta_r)
                else float("nan"),
                n_genes,
                len(common_studies),
            )
    df = pd.DataFrame(rows)
    out = TABLES / "sensitivity_pre_post_harmony.csv"
    df.to_csv(out, index=False)
    logger.info("wrote %s: %d rows", out.name, len(df))
    print()
    print(df.to_string(index=False))
    # Issue 32 mechanical verdict
    print("\n=== Issue 32 verdict (mechanical, per bucket × gene_set) ===")
    for _, r in df.iterrows():
        dr = r["delta_r"]
        if np.isnan(dr):
            verdict = "n/a"
        elif dr <= 0.10:
            verdict = "BIOLOGY_DOMINANT (Δr ≤ 0.10)"
        elif dr <= 0.30:
            verdict = "MIXED (Δr in (0.10, 0.30])"
        else:
            verdict = "HARMONY_DOMINANT (Δr > 0.30)"
        print(f"  {r['bucket']:10s} {r['gene_set']:5s}: Δr={dr:+.4f} → {verdict}")
    # Aggregate: most-buckets verdict per gene_set
    for gs in ("full", "MVS"):
        sub = df[df["gene_set"] == gs]
        bio = (sub["delta_r"] <= 0.10).sum()
        mix = ((sub["delta_r"] > 0.10) & (sub["delta_r"] <= 0.30)).sum()
        har = (sub["delta_r"] > 0.30).sum()
        total = len(sub)
        if bio >= 3:
            agg = "BIOLOGY_DOMINANT"
        elif har >= 3:
            agg = "HARMONY_DOMINANT"
        elif mix + har >= 3:
            agg = "MIXED"
        else:
            agg = "PER_BUCKET_DISCLOSURE"
        print(f"\n  Aggregate {gs}: bio={bio}, mix={mix}, har={har} of {total} → {agg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
