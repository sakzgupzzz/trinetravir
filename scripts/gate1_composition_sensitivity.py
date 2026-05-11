"""Gate 1 composition-correction sensitivity (Part F, Issues 5 + 16).

Lee 2020 SARS-vs-IAV cross-virus Pearson r under three approaches:
  1. Bulk PBMC (composition-confounded baseline).
  2. Per-stratum (chosen primary; per-bucket then averaged).
  3. Bulk PBMC with composition correction (reweight IAV cells so cell-type
     proportions match SARS donors; then bulk response vector).

Output: results/tables/gate1_composition_sensitivity.csv
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
PROCESSED = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"

BUCKETS = ("monocyte", "B", "NK", "CD4T", "CD8T")


def normalize_log(a: ad.AnnData, n_top_genes: int = 4000) -> ad.AnnData:
    """log1p + HVG (batch_key=virus) to get a consistent gene set."""
    a = a.copy()
    a.X = a.X.astype(np.float32)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=n_top_genes, flavor="seurat", subset=True)
    return a


def response_vectors_for_lee(a: ad.AnnData, bucket: str | None) -> dict[str, np.ndarray]:
    """Return {'sars_cov_2': rv, 'iav': rv} where rv = mean(virus) - mean(healthy)."""
    obs = a.obs
    if bucket is not None:
        mask_bucket = obs["cell_type_bucket_unified"].astype(str) == bucket
    else:
        # 'other' included as bulk; user spec says "bulk PBMC".
        mask_bucket = np.ones(len(obs), dtype=bool)
    healthy = mask_bucket & (obs["donor_disease_status"].astype(str) == "healthy_control")
    sars = mask_bucket & (obs["virus"].astype(str) == "sars_cov_2")
    iav = mask_bucket & (obs["virus"].astype(str) == "iav")
    if healthy.sum() < 50 or sars.sum() < 50 or iav.sum() < 50:
        return {}
    X = a.X
    h_mean = np.asarray(X[healthy.values].mean(axis=0)).flatten()
    s_mean = np.asarray(X[sars.values].mean(axis=0)).flatten()
    i_mean = np.asarray(X[iav.values].mean(axis=0)).flatten()
    return {"sars_cov_2": s_mean - h_mean, "iav": i_mean - h_mean}


def pearson_pair(rvs: dict[str, np.ndarray]) -> float:
    if "sars_cov_2" not in rvs or "iav" not in rvs:
        return float("nan")
    return float(np.corrcoef(rvs["sars_cov_2"], rvs["iav"])[0, 1])


def composition_reweighted_rv(a: ad.AnnData) -> dict[str, np.ndarray]:
    """Reweight IAV cells so cell-type proportions match SARS, then bulk response."""
    obs = a.obs
    bucket_col = obs["cell_type_bucket_unified"].astype(str).values
    virus = obs["virus"].astype(str).values
    disease = obs["donor_disease_status"].astype(str).values

    target_buckets = list(BUCKETS)  # match SARS composition over the 5 buckets
    sars_mask = virus == "sars_cov_2"
    iav_mask = virus == "iav"
    healthy_mask = disease == "healthy_control"

    sars_counts = pd.Series({b: int(((bucket_col == b) & sars_mask).sum()) for b in target_buckets})
    iav_counts = pd.Series({b: int(((bucket_col == b) & iav_mask).sum()) for b in target_buckets})
    sars_prop = sars_counts / sars_counts.sum()
    iav_prop = iav_counts / iav_counts.sum()

    # Per-IAV-cell weight = sars_prop[bucket] / iav_prop[bucket]; cells outside
    # target buckets get weight 0 (excluded).
    weights = np.zeros(len(obs), dtype=np.float64)
    for b in target_buckets:
        mask_b = (bucket_col == b) & iav_mask
        if iav_prop[b] > 0:
            weights[mask_b] = sars_prop[b] / iav_prop[b]

    # Weighted IAV mean (per-bucket reweighted to sars composition).
    X = a.X
    iav_in_targets = iav_mask & np.isin(bucket_col, target_buckets)
    w = weights[iav_in_targets]
    Xi = np.asarray(
        X[iav_in_targets].toarray() if hasattr(X[iav_in_targets], "toarray") else X[iav_in_targets]
    )
    iav_mean = (Xi * w[:, None]).sum(axis=0) / w.sum()

    # SARS unweighted mean (bulk over the same buckets).
    sars_in_targets = sars_mask & np.isin(bucket_col, target_buckets)
    Xs = np.asarray(
        X[sars_in_targets].toarray()
        if hasattr(X[sars_in_targets], "toarray")
        else X[sars_in_targets]
    )
    sars_mean = Xs.mean(axis=0)

    # Healthy unweighted mean (bulk over the same buckets).
    healthy_in_targets = healthy_mask & np.isin(bucket_col, target_buckets)
    Xh = np.asarray(
        X[healthy_in_targets].toarray()
        if hasattr(X[healthy_in_targets], "toarray")
        else X[healthy_in_targets]
    )
    healthy_mean = Xh.mean(axis=0)

    return {"sars_cov_2": sars_mean - healthy_mean, "iav": iav_mean - healthy_mean}


def main() -> None:
    print("loading Lee reannotated_low.h5ad")
    a_raw = ad.read_h5ad(PROCESSED / "lee_2020_reannotated_low.h5ad")
    print(f"  {a_raw.shape}")

    # Normalize + HVG on Lee alone so all three approaches share gene set.
    print("normalize + HVG")
    a = normalize_log(a_raw, n_top_genes=4000)
    print(f"  after HVG: {a.n_vars} genes")

    rows = []

    # 1. Bulk
    print("approach 1: bulk PBMC")
    rvs_bulk = response_vectors_for_lee(a, bucket=None)
    r_bulk = pearson_pair(rvs_bulk)
    rows.append(
        {
            "approach": "bulk_pbmc_confounded_baseline",
            "bucket": "bulk",
            "pearson_r": round(r_bulk, 4),
            "interpretation": "Confounded by Lee IAV vs SARS composition differences.",
        }
    )
    print(f"  bulk r = {r_bulk:.4f}")

    # 2. Per-stratum: per-bucket Pearson + mean
    print("approach 2: per-stratum")
    per_bucket_rs = {}
    for b in BUCKETS:
        rvs = response_vectors_for_lee(a, bucket=b)
        r = pearson_pair(rvs)
        per_bucket_rs[b] = r
        rows.append(
            {
                "approach": "per_stratum_primary",
                "bucket": b,
                "pearson_r": round(r, 4),
                "interpretation": "Within-bucket SARS vs IAV; biology, not composition.",
            }
        )
        print(f"  bucket {b}: r = {r:.4f}")
    mean_per_stratum = float(np.mean([v for v in per_bucket_rs.values() if not np.isnan(v)]))
    rows.append(
        {
            "approach": "per_stratum_primary_mean",
            "bucket": "mean_across_buckets",
            "pearson_r": round(mean_per_stratum, 4),
            "interpretation": "Headline cross-virus value under per-stratum primary.",
        }
    )
    print(f"  per-stratum mean = {mean_per_stratum:.4f}")

    # 3. Bulk with composition reweighting
    print("approach 3: bulk with composition correction")
    rvs_corrected = composition_reweighted_rv(a)
    r_corrected = pearson_pair(rvs_corrected)
    rows.append(
        {
            "approach": "bulk_with_composition_correction",
            "bucket": "bulk_reweighted",
            "pearson_r": round(r_corrected, 4),
            "interpretation": "Bulk Pearson with IAV cell-type proportions reweighted to match SARS.",
        }
    )
    print(f"  composition-corrected bulk r = {r_corrected:.4f}")

    out = TABLES / "gate1_composition_sensitivity.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
