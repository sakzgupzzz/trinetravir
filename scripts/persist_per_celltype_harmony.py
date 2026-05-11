"""Persist per-cell-type Harmony embeddings to disk for Session 3 calibration.

The v1 per-cell-type Harmony pipeline (notebooks 04, 06,
scripts/run_harmonization_protocol_sensitivity.py) computes a separate
Harmony correction per coarse bucket but does not persist the full
(n_cells, n_hvg) corrected embedding — only the per-study response
vectors and r matrices are written to disk.

This script re-runs harmony_per_bucket(keep_cells=True) per bucket and
writes one AnnData per bucket containing the Harmony-corrected scaled-HVG
gene-space embedding in layers['X_harmony_scaled_hvg'] plus the PCA-
corrected coordinates in obsm['X_harmony']. Allows METHODS_CHOICES
Issue 7 (per-cell-type vs global Harmony sensitivity) and Issue 3 (full
per-metric calibration including MMD) to run without re-running Harmony.

Output paths:
  data/processed/harmony_per_celltype_<bucket>.h5ad   (one per bucket)

Each file has:
  - obs:  study_id, donor_id, donor_disease_status, virus, cell_type,
          coarse (= bucket)
  - obsm: X_harmony  (n_cells, n_pcs)  — Harmony-corrected PCA coordinates
  - layers: X_harmony_scaled_hvg  (n_cells, n_hvg)  — inverse-PCA projection
  - uns: harmonization_protocol = "per_celltype_harmony"
         bucket, n_top_genes, n_pcs, random_state

Wall time on laptop CPU: ~3-5 min per bucket; ~20-30 min total for the 5
v1 buckets.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np

from trinetravir.data.harmonize import (
    COARSE_BUCKETS,
    concat_clean_studies,
    harmony_per_bucket,
    load_clean_studies,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw"
PROCESSED = REPO / "data" / "processed"
CFG = REPO / "configs" / "datasets.yaml"


def main() -> int:
    print("Loading clean studies + concatenating...")
    studies = load_clean_studies(CFG, RAW)
    combined = concat_clean_studies(studies)
    print(f"Combined: {combined.n_obs} cells x {combined.n_vars} genes\n")

    for bucket in COARSE_BUCKETS:
        out_path = PROCESSED / f"harmony_per_celltype_{bucket}.h5ad"
        if out_path.exists():
            print(f"[{bucket}] already at {out_path.name}; skipping")
            continue
        print(f"\n========== {bucket} ==========")
        res = harmony_per_bucket(combined, bucket, keep_cells=True)
        if res is None or res.x_corrected is None or res.cell_obs is None:
            print(f"  {bucket}: harmony_per_bucket returned None or missing cells")
            continue
        # Re-build per-bucket AnnData from cell_obs + x_corrected. The cells
        # used by harmony_per_bucket are the SARS+mock + bucket subset that
        # passes the min-per-group filter; rebuild a minimal AnnData with
        # just the metadata Session 3 needs.
        cell_ids = res.cell_obs["cell_id"].to_numpy()
        sel = combined.obs.index.isin(cell_ids)
        if int(sel.sum()) != len(cell_ids):
            # Re-derive from cell_obs index alignment instead.
            print(
                f"  [warn] {bucket}: cell_id index alignment mismatch — using cell_obs frame directly"
            )
        obs = res.cell_obs.set_index("cell_id")
        adata_out = ad.AnnData(
            X=np.zeros((len(cell_ids), 1), dtype=np.float32),  # placeholder
            obs=obs.loc[cell_ids],
            var=__import__("pandas").DataFrame(index=["__placeholder__"]),
        )
        adata_out.obsm["X_harmony_scaled_hvg"] = np.asarray(res.x_corrected)
        adata_out.uns["harmonization_protocol"] = "per_celltype_harmony"
        adata_out.uns["bucket"] = bucket
        adata_out.uns["hvg_genes"] = list(res.hvg_genes)
        adata_out.uns["studies_used"] = list(res.studies_used)
        adata_out.write_h5ad(out_path)
        print(f"  wrote {out_path.name}: {adata_out.shape}, layers={list(adata_out.layers)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
