"""Phase 3 harmonization protocol sensitivity (METHODS_CHOICES Issue 7).

Compares per-cell-type Harmony (v1 chosen protocol) against global Harmony
on the combined AnnData. For each protocol, computes per-bucket cross-study
mean off-diagonal Pearson r and tabulates the verdict.

Per-cell-type protocol: run Harmony separately on each of 5 buckets with
  ``study_id`` as the sole batch key (the v1 protocol; results already in
  ``data/processed/phase3_response_vectors_<bucket>.parquet`` from notebook
  04).

Global protocol: run Harmony once on all cells (combined AnnData) with
  ``study_id`` as the sole batch key, then compute per-bucket response
  vectors in the shared embedding by sub-setting cells per coarse bucket
  from the SAME corrected embedding.

The comparison answers Issue 7's question: does the per-cell-type protocol
materially change the cross-study coherence verdict, or is global Harmony
adequate?

Wall time: per-cell-type results are cached. Global Harmony on ~244k cells
takes ~20-30 min on laptop CPU. Sequential.

Output: results/tables/harmonization_protocol_sensitivity.csv
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import harmonypy
import numpy as np
import pandas as pd
import scanpy as sc

from trinetravir.data.harmonize import (
    COARSE_BUCKETS,
    GATE_THRESHOLDS,
    concat_clean_studies,
    load_clean_studies,
)
from trinetravir.eval.metrics import pearson_off_diag

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw"
PROCESSED = REPO / "data" / "processed"
CFG = REPO / "configs" / "datasets.yaml"
TABLES = REPO / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)


def global_harmony_response_vectors(
    combined: ad.AnnData,
    *,
    n_top_genes: int = 4000,
    n_pcs: int = 50,
    random_state: int = 42,
    persist_embedding_to: Path | None = None,
) -> dict[str, dict[str, np.ndarray]]:
    """Run one global Harmony pass then derive per-bucket response vectors.

    Mirrors the per-cell-type pipeline (normalize -> log -> HVG -> scale ->
    PCA -> Harmony -> inverse-PCA) but on ALL cells together, then per-bucket
    response vectors are computed by sub-setting cells in the SHARED
    Harmony-corrected scaled-HVG space.
    """
    work = combined.copy()
    work.X = work.X.astype(np.float32)
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(
        work, n_top_genes=n_top_genes, batch_key="study_id", flavor="seurat", subset=True
    )
    hvg_genes = list(work.var_names)
    logging.info("global Harmony: %d HVG", len(hvg_genes))
    sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=n_pcs, random_state=random_state)
    pca_embedding = work.obsm["X_pca"].astype(np.float64)
    harmony_out = harmonypy.run_harmony(
        pca_embedding, work.obs, ["study_id"], random_state=random_state
    )
    z_corr = np.asarray(harmony_out.Z_corr)
    if z_corr.shape == (n_pcs, work.n_obs):
        z_corr = z_corr.T
    work.obsm["X_pca_harmony"] = z_corr
    loadings = work.varm["PCs"]
    x_corrected = work.obsm["X_pca_harmony"] @ loadings.T  # (n_cells, n_hvg)

    if persist_embedding_to is not None:
        # Persist the full integrated AnnData with the global-Harmony embedding
        # in obsm['X_harmony'] (the PCA-corrected coordinates) and x_corrected
        # mirrored as a layer for direct cell-level computation in downstream
        # calibration. See METHODS_CHOICES Issue 7.
        persist_embedding_to.parent.mkdir(parents=True, exist_ok=True)
        out_adata = work.copy()
        out_adata.obsm["X_harmony"] = z_corr
        out_adata.layers["X_harmony_scaled_hvg"] = x_corrected
        out_adata.uns["harmonization_protocol"] = "global_harmony_study_id_only"
        out_adata.uns["harmony_n_top_genes"] = n_top_genes
        out_adata.uns["harmony_n_pcs"] = n_pcs
        out_adata.uns["harmony_random_state"] = random_state
        out_adata.write_h5ad(persist_embedding_to)
        logging.info("persisted global Harmony embedding -> %s", persist_embedding_to)

    out: dict[str, dict[str, np.ndarray]] = {}
    for bucket in COARSE_BUCKETS:
        bucket_mask = (work.obs["coarse"] == bucket).to_numpy()
        sars_mask = work.obs["virus"].isin(["sars_cov_2", "mock"]).to_numpy()
        keep_cells = bucket_mask & sars_mask
        rvs: dict[str, np.ndarray] = {}
        for sid in sorted(work.obs["study_id"].astype(str).unique()):
            study_mask = (work.obs["study_id"].astype(str) == sid).to_numpy()
            sel = keep_cells & study_mask
            obs = work.obs.loc[sel]
            cells_d = (obs["donor_disease_status"] == "diseased").to_numpy()
            cells_h = (obs["donor_disease_status"] == "healthy_control").to_numpy()
            if cells_d.sum() < 50 or cells_h.sum() < 50:
                continue
            x_sel = x_corrected[sel]
            rv = x_sel[cells_d].mean(axis=0) - x_sel[cells_h].mean(axis=0)
            rvs[sid] = rv
        if len(rvs) >= 2:
            out[bucket] = rvs
        else:
            logging.warning(
                "global Harmony bucket %s: only %d studies with >=50 cells per class; skipping",
                bucket,
                len(rvs),
            )
    return out


def load_per_celltype_rvs(parquet_dir: Path) -> dict[str, dict[str, np.ndarray]]:
    """Load cached per-cell-type Phase 3 response vectors."""
    out: dict[str, dict[str, np.ndarray]] = {}
    for bucket in COARSE_BUCKETS:
        path = parquet_dir / f"phase3_response_vectors_{bucket}.parquet"
        if not path.exists():
            continue
        df = pd.read_parquet(path)
        out[bucket] = {sid: df[sid].to_numpy() for sid in df.columns}
    return out


def main() -> int:
    # Per-cell-type (cached)
    per_ct = load_per_celltype_rvs(PROCESSED)
    if not per_ct:
        print("Per-cell-type cache missing. Run notebook 04 first.", file=sys.stderr)
        return 1
    print(f"Loaded per-cell-type response vectors for {len(per_ct)} buckets")

    # Global Harmony (compute)
    print("Loading studies + concatenating for global Harmony pass...")
    studies = load_clean_studies(CFG, RAW)
    combined = concat_clean_studies(studies)
    print(f"Combined: {combined.n_obs} cells x {combined.n_vars} genes")
    print("Running global Harmony (one pass, study_id only)...")
    embedding_path = PROCESSED / "harmony_global_embedding.h5ad"
    global_rvs = global_harmony_response_vectors(combined, persist_embedding_to=embedding_path)
    print(f"Global Harmony produced response vectors for {len(global_rvs)} buckets")

    rows = []
    for bucket in COARSE_BUCKETS:
        th = GATE_THRESHOLDS[bucket]
        per_ct_r = pearson_off_diag(per_ct.get(bucket, {})) if per_ct.get(bucket) else float("nan")
        global_r = (
            pearson_off_diag(global_rvs.get(bucket, {})) if global_rvs.get(bucket) else float("nan")
        )
        rows.append(
            {
                "bucket": bucket,
                "per_celltype_r": round(float(per_ct_r), 3)
                if not np.isnan(per_ct_r)
                else float("nan"),
                "global_r": round(float(global_r), 3) if not np.isnan(global_r) else float("nan"),
                "delta_global_minus_perct": (
                    round(float(global_r - per_ct_r), 3)
                    if not (np.isnan(global_r) or np.isnan(per_ct_r))
                    else float("nan")
                ),
                "threshold": th,
                "per_celltype_pass": (per_ct_r >= th) if not np.isnan(per_ct_r) else False,
                "global_pass": (global_r >= th) if not np.isnan(global_r) else False,
                "verdict_matches": (
                    (per_ct_r >= th) == (global_r >= th)
                    if not (np.isnan(global_r) or np.isnan(per_ct_r))
                    else False
                ),
            }
        )
    table = pd.DataFrame(rows)
    out_csv = TABLES / "harmonization_protocol_sensitivity.csv"
    table.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")
    print("\n--- harmonization protocol sensitivity ---")
    pd.set_option("display.width", 220)
    print(table.to_string(index=False))

    # Persist global RVs so future work can re-use them
    for bucket, rvs in global_rvs.items():
        if not rvs:
            continue
        df = pd.DataFrame(dict(rvs))
        df.to_parquet(PROCESSED / f"phase3_global_response_vectors_{bucket}.parquet")
    return 0


if __name__ == "__main__":
    sys.exit(main())
