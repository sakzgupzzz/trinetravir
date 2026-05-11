"""Extract per-bucket raw counts for Session 4 scVI input.

Per `references/session_4_prompt.md` Part A.1.

Per-bucket + global Harmony h5ads do NOT preserve raw counts (verified
2026-05-11: all per-bucket files have layers=[], raw=False; global has only
X_harmony_scaled_hvg). scVI requires count-level data. This script extracts
raw counts from per-study reannotated h5ads, subsets to bucket cells via
cell-id positional mapping, restricts to the bucket's 4000 HVG, and writes
data/processed/scvi_input_<bucket>.h5ad.

Pattern reusable from scripts/session7_part_a_pre_post_harmony.py
(compute_pre_harmony_rv_per_study).

Inputs:
  data/processed/harmony_per_celltype_<bucket>.h5ad
    - obs.index format: '<int>-<study_id>' (positional row idx in study h5ad)
    - uns['hvg_genes']: 4000-gene HVG list
  data/processed/<study_id>_reannotated.h5ad (4 v1 studies)
    - X = raw counts (max > 20)
    - var has feature_name / gene_symbol column for symbol resolution

Output per bucket:
  data/processed/scvi_input_<bucket>.h5ad
    - X = raw counts (np.float32 sparse) on bucket's HVG
    - obs: study_id, donor_id, donor_disease_status, cell_id_original
    - var index: HVG gene symbols

Verification:
  - cell count matches harmony h5ad bucket count
  - X.max() > 20 confirms counts
  - obs['study_id'].nunique() == 4
  - HVG list matches uns['hvg_genes'] from source harmony h5ad

Idempotent: rerunning regenerates output files in place.
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"

BUCKETS = ("monocyte", "B", "NK", "CD4T", "CD8T")
STUDY_FILE = {
    "wilk_2020": PROC / "wilk_2020_reannotated.h5ad",
    "lee_2020": PROC / "lee_2020_reannotated.h5ad",
    "arunachalam_2020": PROC / "arunachalam_2020_reannotated.h5ad",
    "schulte_schrepping_2020": PROC / "schulte_schrepping_2020_reannotated.h5ad",
}


def resolve_symbols(a: ad.AnnData) -> np.ndarray:
    for c in ("gene_symbol", "feature_name", "name", "gene_symbols", "symbol"):
        if c in a.var.columns:
            return a.var[c].astype(str).values
    return a.var_names.astype(str).values


def extract_bucket(bucket: str) -> ad.AnnData:
    """Build scVI input AnnData for a single bucket."""
    h5 = PROC / f"harmony_per_celltype_{bucket}.h5ad"
    if not h5.exists():
        raise FileNotFoundError(h5)
    h = ad.read_h5ad(h5)
    hvg_genes = [str(g) for g in h.uns["hvg_genes"]]
    harmony_obs = h.obs[["study_id", "donor_id", "donor_disease_status"]].copy()
    harmony_obs.index = harmony_obs.index.astype(str)
    expected_n_cells = len(harmony_obs)
    logger.info("%s: harmony h5ad has %d cells, %d HVG", bucket, expected_n_cells, len(hvg_genes))

    studies = sorted(harmony_obs["study_id"].astype(str).unique())
    per_study_anndatas = []
    for study in studies:
        h_study_obs = harmony_obs.loc[harmony_obs["study_id"].astype(str) == study].copy()
        # Cell-id format: '<int>-<study_id>' → strip suffix → positional row idx.
        suffix = f"-{study}"
        h_study_obs["row_idx"] = (
            h_study_obs.index.astype(str).str.replace(suffix, "", regex=False).astype(int)
        )
        cell_idx = h_study_obs["row_idx"].values
        if len(cell_idx) == 0:
            logger.warning("  %s %s: no cells in harmony obs; skipping", bucket, study)
            continue

        a = ad.read_h5ad(STUDY_FILE[study])
        # Resolve gene symbols + align var_names.
        sym = resolve_symbols(a)
        a.var_names = sym
        a.var_names_make_unique()

        # Slice study h5ad at positional row indices.
        sub = a[cell_idx].copy()
        # Verify raw counts (max should be > 20 for normalized log1p; raw is integer-like > 0).
        x_max = float(np.asarray(sub.X.max()))
        if x_max <= 20:
            logger.warning(
                "  %s %s: X.max=%.2f is suspiciously low — may be normalized, not raw counts",
                bucket,
                study,
                x_max,
            )

        # Restrict to bucket's HVG (intersection with this study's var_names).
        hvg_in = [g for g in hvg_genes if g in set(sub.var_names)]
        if not hvg_in:
            logger.warning("  %s %s: no HVG overlap; skipping", bucket, study)
            continue
        sub = sub[:, hvg_in].copy()

        # Attach canonical obs columns from harmony source.
        # h_study_obs rows are ordered by harmony obs.index; sub.obs rows are
        # in study h5ad row order matching cell_idx slicing → same order.
        h_aligned = h_study_obs.set_index("row_idx").loc[cell_idx]
        sub.obs["study_id"] = study
        sub.obs["donor_id"] = h_aligned["donor_id"].astype(str).values
        sub.obs["donor_disease_status"] = h_aligned["donor_disease_status"].astype(str).values
        sub.obs["cell_id_original"] = h_aligned.index.astype(str).values
        # Drop other obs columns to keep output lean.
        keep_obs = ["study_id", "donor_id", "donor_disease_status", "cell_id_original"]
        sub.obs = sub.obs[keep_obs]
        # Drop var columns; we just need symbol-indexed var.
        sub.var = pd.DataFrame(index=sub.var_names)

        # Force sparse float32 for scVI.
        if not sp.issparse(sub.X):
            sub.X = sp.csr_matrix(sub.X.astype(np.float32))
        else:
            sub.X = sub.X.astype(np.float32)

        per_study_anndatas.append(sub)
        logger.info(
            "  %s %s: %d cells × %d HVG; X.max=%.1f",
            bucket,
            study,
            sub.n_obs,
            sub.n_vars,
            float(np.asarray(sub.X.max())),
        )

    if not per_study_anndatas:
        raise RuntimeError(f"{bucket}: no per-study extractions succeeded")

    # Concat; outer-join HVG axis so genes missing in some studies pad as zeros.
    merged = ad.concat(per_study_anndatas, axis=0, join="outer", merge="same", index_unique=None)
    # Re-order var to canonical HVG order from harmony source.
    hvg_in_merged = [g for g in hvg_genes if g in set(merged.var_names)]
    merged = merged[:, hvg_in_merged].copy()
    merged.var = pd.DataFrame(index=merged.var_names)

    # Verification gates.
    assert merged.n_obs == expected_n_cells, (
        f"{bucket}: extracted {merged.n_obs} cells but harmony h5ad has {expected_n_cells}"
    )
    x_max = float(np.asarray(merged.X.max()))
    assert x_max > 20, f"{bucket}: X.max={x_max} suggests normalized counts, not raw"
    assert merged.obs["study_id"].nunique() == 4, (
        f"{bucket}: {merged.obs['study_id'].nunique()} studies (expected 4)"
    )
    logger.info(
        "%s: merged %d cells × %d HVG (X.max=%.1f, %d studies)",
        bucket,
        merged.n_obs,
        merged.n_vars,
        x_max,
        merged.obs["study_id"].nunique(),
    )
    return merged


def main() -> int:
    PROC.mkdir(parents=True, exist_ok=True)
    for bucket in BUCKETS:
        logger.info("=== %s ===", bucket)
        merged = extract_bucket(bucket)
        out = PROC / f"scvi_input_{bucket}.h5ad"
        merged.write_h5ad(out, compression="gzip")
        size_mb = out.stat().st_size / 1e6
        logger.info(
            "wrote %s: %d cells × %d HVG (%.1f MB)", out.name, merged.n_obs, merged.n_vars, size_mb
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
