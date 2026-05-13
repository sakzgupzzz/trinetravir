"""Extract raw counts on harmony_global's 4000 HVG for Part B scVI global sweep.

Per references/session_4_prompt.md Part B + user 2026-05-12 spec verification:
input to global scVI must use the SAME HVG list as harmony_global_embedding.h5ad
(4000 genes via var.index) so that Δr_global comparison measures method difference,
not HVG-space difference.

Recipe:
  INPUTS:
    data/processed/harmony_global_embedding.h5ad
      → obs.index format: '<positional_idx>-<study_id>'
      → var.index = 4000 gene symbols (HVG)
      → obs has: study_id, donor_id, donor_disease_status, coarse

    data/processed/<study_id>_reannotated.h5ad × 4
      → wilk_2020, lee_2020, arunachalam_2020, schulte_schrepping_2020
      → X = raw counts on full ~61K genes
      → var.index = gene symbols

  ALGORITHM:
    1. Load harmony_global obs + var.index (4000 HVG gene symbols) + obs columns.
    2. For each study in obs['study_id']:
         a. Load study reannotated h5ad.
         b. Strip suffix '-<study_id>' from harmony_global obs.index → positional idx.
         c. Slice study h5ad at those positional indices.
         d. Subset to intersection of study's var.index with harmony_global's 4000 HVG.
         e. Reindex columns to harmony_global var order (NaN columns filled with 0
            if any genes missing in study).
         f. Append obs columns: study_id, donor_id, donor_disease_status, coarse.
    3. Concat across studies → single AnnData with:
         X = raw counts on 4000 HVG (harmony_global order)
         obs columns: study_id, donor_id, donor_disease_status, coarse
         var.index = harmony_global HVG list (4000 symbols)
    4. Write data/processed/scvi_input_global.h5ad with compression='gzip'.

  VERIFICATION GATE:
    - cell count == 244,389 (harmony_global)
    - X.max() > 20 (raw counts, not log-normalized)
    - 4 unique study_ids
    - 5 coarse buckets (+ NaN dropped) — keep NaN for downstream filter
    - 4000 genes in var
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
HARMONY_GLOBAL = PROC / "harmony_global_embedding.h5ad"
OUTPUT = PROC / "scvi_input_global.h5ad"

STUDY_FILES = {
    "wilk_2020": PROC / "wilk_2020_reannotated.h5ad",
    "lee_2020": PROC / "lee_2020_reannotated.h5ad",
    "arunachalam_2020": PROC / "arunachalam_2020_reannotated.h5ad",
    "schulte_schrepping_2020": PROC / "schulte_schrepping_2020_reannotated.h5ad",
}


def main() -> int:
    logger.info("loading harmony_global reference...")
    hg = ad.read_h5ad(HARMONY_GLOBAL)
    hvg_symbols = list(hg.var.index)
    logger.info("harmony_global: shape=%s, HVG n=%d", hg.shape, len(hvg_symbols))

    keep_obs_cols = ["study_id", "donor_id", "donor_disease_status", "coarse"]
    hg_obs = hg.obs[keep_obs_cols].copy()
    hg_obs["original_id"] = hg.obs.index

    accumulated_X = []
    accumulated_obs = []
    n_total_expected = hg.shape[0]

    for study_id, study_path in STUDY_FILES.items():
        logger.info("=== %s ===", study_id)
        s = ad.read_h5ad(study_path)
        logger.info("  loaded %s: shape=%s", study_path.name, s.shape)

        # cells in harmony_global from this study
        study_mask = hg_obs["study_id"] == study_id
        study_obs = hg_obs[study_mask].copy()
        if study_obs.empty:
            logger.warning("  no cells in harmony_global from %s; skip", study_id)
            continue
        logger.info("  %s cells in harmony_global from %s", len(study_obs), study_id)

        # positional indices: strip '-<study_id>' suffix from harmony_global obs.index
        suffix = f"-{study_id}"
        positional_idx = (
            study_obs["original_id"].str.replace(suffix, "", regex=False).astype(int).values
        )

        # slice study h5ad at these positional indices
        s_sub = s[positional_idx, :]
        logger.info("  sliced study h5ad: shape=%s", s_sub.shape)

        # subset to harmony_global HVG (intersection on var.index = gene symbol)
        study_symbols = s_sub.var.index.astype(str)
        hvg_in_study = pd.Index(hvg_symbols).intersection(study_symbols)
        missing = pd.Index(hvg_symbols).difference(study_symbols)
        logger.info(
            "  HVG: %d in study, %d missing (will be zero-filled)",
            len(hvg_in_study),
            len(missing),
        )

        # subset + reindex to harmony_global HVG order
        s_hvg = s_sub[:, hvg_in_study]
        # build full HVG-ordered matrix (zero-fill missing genes)
        if len(missing) == 0 and list(s_hvg.var.index) == hvg_symbols:
            X_aligned = s_hvg.X
        else:
            n_cells = s_hvg.shape[0]
            X_full = sp.lil_matrix((n_cells, len(hvg_symbols)), dtype=np.float32)
            study_col_lookup = {g: i for i, g in enumerate(s_hvg.var.index)}
            for j, sym in enumerate(hvg_symbols):
                if sym in study_col_lookup:
                    src_col = study_col_lookup[sym]
                    src = s_hvg.X[:, src_col]
                    if sp.issparse(src):
                        src = src.toarray().flatten()
                    X_full[:, j] = src.reshape(-1, 1)
            X_aligned = X_full.tocsr()

        accumulated_X.append(X_aligned)
        accumulated_obs.append(study_obs.reset_index(drop=True))

    # concat
    logger.info("concatenating across %d studies...", len(accumulated_X))
    X_all = (
        sp.vstack(accumulated_X)
        if all(sp.issparse(x) for x in accumulated_X)
        else np.vstack(accumulated_X)
    )
    obs_all = pd.concat(accumulated_obs, ignore_index=True)

    # build output AnnData
    out = ad.AnnData(
        X=X_all,
        obs=obs_all,
        var=pd.DataFrame(index=hvg_symbols),
    )
    out.uns["hvg_genes"] = hvg_symbols
    out.uns["source"] = "scripts/extract_global_counts.py"
    out.uns["harmony_global_reference"] = str(HARMONY_GLOBAL.name)

    # verification
    logger.info("output shape: %s (expected %d cells)", out.shape, n_total_expected)
    assert out.shape[0] == n_total_expected, (
        f"cell count mismatch: {out.shape[0]} vs {n_total_expected}"
    )
    assert out.shape[1] == 4000, f"HVG count mismatch: {out.shape[1]} vs 4000"
    x_max = float(out.X.max())
    logger.info("X.max = %.1f (expect >20 for raw counts)", x_max)
    assert x_max > 20, f"X.max={x_max} ≤ 20 suggests not raw counts"
    n_studies = out.obs["study_id"].nunique()
    assert n_studies == 4, f"study count: {n_studies} vs 4"

    logger.info("writing %s...", OUTPUT)
    out.write_h5ad(OUTPUT, compression="gzip")
    logger.info("DONE: wrote %s (%.1f MB)", OUTPUT.name, OUTPUT.stat().st_size / 1e6)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
