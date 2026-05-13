"""Session 4.5 Part F.1: build precomputed Harmony global response vector cache.

Replaces 16GB harmony_global_embedding.h5ad upload for Part B with a ~5MB
precomputed cache containing only the per-(bucket, study, status) Harmony
response vector means needed for Δr_global comparison.

Cache shape: ~40 rows (5 buckets × 4 studies × 2 statuses, minus empty strata)
× 4000 HVG columns. Storable as parquet with gene symbols as columns.

Per Session 4.5 spec Part F.1: bucket assignment from obs['coarse'] (NOT
cell_type_bucket per Issue 7 status note).

Input:
  data/processed/harmony_global_embedding.h5ad (16 GB, stays local)

Output:
  data/processed/harmony_global_response_vector_cache.parquet (~5 MB)

Schema:
  Rows indexed by (bucket, study_id, donor_disease_status) tuple
  Columns = 4000 HVG gene symbols
  Values = mean Harmony-scaled HVG expression across cells in stratum
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
INPUT = PROC / "harmony_global_embedding.h5ad"
OUTPUT = PROC / "harmony_global_response_vector_cache.parquet"

BUCKETS = ("monocyte", "B", "NK", "CD4T", "CD8T")


def main() -> int:
    logger.info("loading %s...", INPUT.name)
    a = ad.read_h5ad(INPUT)
    logger.info("loaded: shape=%s", a.shape)

    # Use X_harmony_scaled_hvg layer if present, else X
    X = a.layers.get("X_harmony_scaled_hvg", a.X)
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)
    logger.info("X dtype=%s shape=%s", X.dtype, X.shape)

    hvg_genes = list(a.var.index)
    obs = a.obs[["study_id", "donor_disease_status", "coarse"]].copy()

    # Build per-(bucket, study, status) mean response vectors
    rows = []
    for bucket in BUCKETS:
        bucket_mask = (obs["coarse"] == bucket).values
        if bucket_mask.sum() == 0:
            logger.warning("  %s: no cells, skip", bucket)
            continue
        b_obs = obs[bucket_mask].copy().reset_index(drop=True)
        b_X = X[bucket_mask]

        for (study, status), sub in b_obs.groupby(
            ["study_id", "donor_disease_status"], observed=True
        ):
            if len(sub) < 1:
                continue
            sub_mask = (
                (b_obs["study_id"] == study) & (b_obs["donor_disease_status"] == status)
            ).values
            mean_expr = b_X[sub_mask].mean(axis=0)
            row = {
                "bucket": bucket,
                "study_id": str(study),
                "donor_disease_status": str(status),
                "n_cells": int(sub_mask.sum()),
            }
            for i, gene in enumerate(hvg_genes):
                row[gene] = float(mean_expr[i])
            rows.append(row)
            logger.info(
                "  %s | %s | %s: n_cells=%d, mean_expr_norm=%.3f",
                bucket,
                study,
                status,
                int(sub_mask.sum()),
                float(np.linalg.norm(mean_expr)),
            )

    df = pd.DataFrame(rows)
    logger.info("cache shape: %s", df.shape)

    df.to_parquet(OUTPUT, compression="snappy", index=False)
    logger.info("wrote %s (%.1f MB)", OUTPUT.name, OUTPUT.stat().st_size / 1e6)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
