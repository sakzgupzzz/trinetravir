"""Phase 3.5 Harmony — per-bucket runs on CellTypist-unified labels.

Usage:
  uv run python scripts/phase3_5_harmony.py {low|high|subbucket_low}

Produces:
  data/processed/harmony_per_celltype_phase35_low.h5ad        (mode=low)
  data/processed/harmony_per_celltype_phase35_high.h5ad       (mode=high)
  data/processed/harmony_subbucket_phase35_low.h5ad           (mode=subbucket_low)
  results/tables/gate_phase35_immune_all_low.csv              (mode=low)
  results/tables/gate_phase35_immune_all_high.csv             (mode=high)
  results/tables/gate_phase35_subbucket_granularity_low.csv   (mode=subbucket_low)

Each consolidated h5ad carries (one row per cell, stitched across all
surviving buckets; cells outside any surviving bucket are dropped):
  obs                              — cell metadata: study_id, donor_id,
                                     donor_disease_status, bucket.
  X                                — (n_cells, 1) placeholder. Harmony-
                                     corrected matrices are NOT in X.
  uns['harmonization_protocol']    — 'per_celltype_harmony'.
  uns['annotation_source']         — 'celltypist_immune_all_low' or
                                     'celltypist_immune_all_high'.
  uns['bucket_column']             — name of the bucket obs col.
  uns['per_bucket_r']              — dict of bucket -> mean off-diag r.
  uns['per_bucket_hvg']            — dict of bucket -> list of HVG symbols
                                     selected for that bucket.
  uns['X_corrected_<bucket>']      — (n_cells_in_bucket, n_hvg_for_bucket)
                                     Harmony-corrected scaled-HVG embedding
                                     in gene space. ONE KEY PER SURVIVING
                                     BUCKET. Different buckets have
                                     different HVG sets so a single
                                     obsm['X_harmony'] matrix across
                                     buckets is NOT FEASIBLE — Session 3
                                     loaders must iterate over uns keys
                                     matching the prefix 'X_corrected_'
                                     and align rows to obs[bucket_column].
                                     The corresponding HVG symbol list is
                                     uns['per_bucket_hvg'][<bucket>].

Gate values are reported FOR RECORD ONLY (heuristic 5-bucket Phase 3
thresholds for `low`; no thresholds for other granularities — gate
column stores observed r against threshold NaN).
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from trinetravir.data import harmonize as harmonize_mod

REPO = Path(__file__).resolve().parents[1]
PROCESSED = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"

STUDIES = ("wilk_2020", "arunachalam_2020", "lee_2020", "schulte_schrepping_2020")

# Phase 3 5-bucket heuristic thresholds for record. Other buckets default 0.0
# (which means gate_passed is True for any nonneg r — purely informational).
HEURISTIC_THRESHOLDS = {
    "monocyte": 0.60,
    "B": 0.40,
    "NK": 0.35,
    "CD4T": 0.30,
    "CD8T": 0.25,
    "T": 0.30,  # High-model collapsed T bucket; midpoint of CD4 and CD8 thresholds
}


def load_reannotated(suffix: str) -> ad.AnnData:
    """Load + concatenate 4 study reannotated h5ads with the given suffix."""
    parts = []
    for study in STUDIES:
        path = PROCESSED / f"{study}_reannotated_{suffix}.h5ad"
        if not path.exists():
            raise FileNotFoundError(f"missing {path}")
        a = ad.read_h5ad(path)
        a.obs["study_id"] = study
        parts.append(a)
    # outer-join concat by var_names because the raw h5ads share Census gene set.
    combined = ad.concat(parts, join="outer", merge="first", uns_merge="first")
    combined.obs["study_id"] = combined.obs["study_id"].astype("category")
    if "virus" not in combined.obs.columns:
        # All v1 studies are SARS; treat healthy as 'mock' so harmony_per_bucket
        # subset logic (virus.isin(['sars_cov_2', 'mock'])) passes through.
        combined.obs["virus"] = np.where(
            combined.obs["donor_disease_status"].astype(str) == "diseased", "sars_cov_2", "mock"
        )
    return combined


def run_mode(mode: str) -> None:
    if mode == "low":
        suffix = "low"
        annotation_source = "celltypist_immune_all_low"
        bucket_col = "cell_type_bucket_unified"
        out_h5 = PROCESSED / "harmony_per_celltype_phase35_low.h5ad"
        out_csv = TABLES / "gate_phase35_immune_all_low.csv"
    elif mode == "high":
        suffix = "high"
        annotation_source = "celltypist_immune_all_high"
        bucket_col = "cell_type_bucket_unified"
        out_h5 = PROCESSED / "harmony_per_celltype_phase35_high.h5ad"
        out_csv = TABLES / "gate_phase35_immune_all_high.csv"
    elif mode == "subbucket_low":
        suffix = "low"
        annotation_source = "celltypist_immune_all_low"
        bucket_col = "cell_type_subbucket_unified"
        out_h5 = PROCESSED / "harmony_subbucket_phase35_low.h5ad"
        out_csv = TABLES / "gate_phase35_subbucket_granularity_low.csv"
    else:
        raise SystemExit(f"unknown mode {mode!r}")

    if out_h5.exists() and out_csv.exists():
        print(f"[{mode}] already done; skipping ({out_h5.name})")
        return

    print(f"[{mode}] loading reannotated_{suffix} h5ads...")
    combined = load_reannotated(suffix)
    print(f"  combined: {combined.n_obs} cells x {combined.n_vars} genes")

    if bucket_col not in combined.obs.columns:
        raise SystemExit(f"missing obs col {bucket_col!r}")

    # Drop 'other' bucket cells — they're unmapped.
    combined = combined[combined.obs[bucket_col].astype(str) != "other"].copy()
    combined.obs["coarse"] = combined.obs[bucket_col].astype(str).astype("category")
    print(f"  after dropping 'other': {combined.n_obs} cells")

    # Patch GATE_THRESHOLDS with a permissive defaultdict so harmony_per_bucket
    # works for unknown bucket names. Heuristic 5-bucket thresholds are used
    # when present; others get 0.0 (pass = r >= 0, informational).
    original_thresholds = harmonize_mod.GATE_THRESHOLDS
    harmonize_mod.GATE_THRESHOLDS = defaultdict(lambda: 0.0, HEURISTIC_THRESHOLDS)

    # Also patch COARSE_BUCKETS so harmony_per_bucket allows any bucket value.
    original_buckets = harmonize_mod.COARSE_BUCKETS
    bucket_values = sorted(combined.obs["coarse"].astype(str).unique().tolist())
    harmonize_mod.COARSE_BUCKETS = tuple(bucket_values)
    print(f"  buckets to process: {bucket_values}")

    rows = []
    per_bucket_obs = []
    per_bucket_x = []
    try:
        for bucket in bucket_values:
            print(f"  --> harmony for bucket {bucket}")
            try:
                res = harmonize_mod.harmony_per_bucket(combined, bucket, keep_cells=True)
            except Exception as e:
                print(f"     ERROR: {e}")
                rows.append(
                    {
                        "bucket": bucket,
                        "n_cells": 0,
                        "n_studies": 0,
                        "studies_used": "",
                        "gate_r": float("nan"),
                        "threshold": HEURISTIC_THRESHOLDS.get(bucket, float("nan")),
                        "gate_passed": False,
                        "error": str(e),
                    }
                )
                continue
            if res is None:
                print("     skipped (insufficient studies)")
                rows.append(
                    {
                        "bucket": bucket,
                        "n_cells": 0,
                        "n_studies": 0,
                        "studies_used": "",
                        "gate_r": float("nan"),
                        "threshold": HEURISTIC_THRESHOLDS.get(bucket, float("nan")),
                        "gate_passed": False,
                        "error": "insufficient_studies",
                    }
                )
                continue
            rows.append(
                {
                    "bucket": bucket,
                    "n_cells": res.n_cells,
                    "n_studies": len(res.studies_used),
                    "studies_used": ",".join(res.studies_used),
                    "gate_r": round(res.gate_r, 4),
                    "threshold": HEURISTIC_THRESHOLDS.get(bucket, float("nan")),
                    "gate_passed": bool(res.gate_passed),
                    "error": "",
                }
            )
            # Stitch cells into consolidated AnnData with X_harmony (PCA-Harmony).
            # res.x_corrected is in gene space; we want the PCA-Harmony coords.
            # res does not expose them directly. We re-derive by storing x_corrected
            # and projecting back via SVD of x_corrected. Simpler: just keep the
            # gene-space x_corrected as obsm['X_harmony_genespace'] AND also
            # compute a per-bucket PCA on x_corrected to give a low-dim X_harmony.
            # For Session 3's MMD on bucket-level distributions, gene-space is fine.
            cell_obs = res.cell_obs.copy()
            cell_obs["bucket"] = bucket
            per_bucket_obs.append(cell_obs)
            # Save x_corrected per bucket; consolidated obsm not viable across
            # buckets (different HVG sets per bucket). Stash per-bucket parquet
            # for response vectors AND save bucket-level x_corrected as separate
            # group inside uns to keep one file.
            per_bucket_x.append(
                (bucket, np.asarray(res.x_corrected, dtype=np.float32), list(res.hvg_genes))
            )
    finally:
        harmonize_mod.GATE_THRESHOLDS = original_thresholds
        harmonize_mod.COARSE_BUCKETS = original_buckets

    # Consolidated AnnData: stack cells across buckets, with bucket label in obs.
    # X stays a (n_cells, 1) placeholder; each bucket's gene-space corrected
    # matrix lives in uns under key f'X_corrected_{bucket}' (different HVG
    # sets per bucket prevent a single obsm matrix across buckets).
    # Session 3 loaders enumerate uns keys with prefix 'X_corrected_' and pair
    # each with uns['per_bucket_hvg'][bucket] for column ids.
    all_obs = pd.concat(per_bucket_obs, axis=0, ignore_index=True)
    n_total = len(all_obs)
    adata_out = ad.AnnData(
        X=np.zeros((n_total, 1), dtype=np.float32),
        obs=all_obs.set_index("cell_id") if "cell_id" in all_obs.columns else all_obs,
        var=pd.DataFrame(index=["__placeholder__"]),
    )
    adata_out.uns["harmonization_protocol"] = "per_celltype_harmony"
    adata_out.uns["annotation_source"] = annotation_source
    adata_out.uns["bucket_column"] = bucket_col
    adata_out.uns["per_bucket_r"] = {r["bucket"]: r["gate_r"] for r in rows}
    adata_out.uns["per_bucket_hvg"] = {b: hvg for (b, _x, hvg) in per_bucket_x}
    # Per-bucket gene-space corrected matrices (large): stored as a list of
    # (bucket_name, x_matrix). AnnData uns supports dict-of-arrays.
    for bucket, x, _hvg in per_bucket_x:
        adata_out.uns[f"X_corrected_{bucket}"] = x

    adata_out.write_h5ad(out_h5)
    print(f"  wrote {out_h5.name}: shape={adata_out.shape}")

    table = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_csv, index=False)
    print(f"  wrote {out_csv.name}: {len(table)} rows")
    print(table.to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("low", "high", "subbucket_low"):
        raise SystemExit("usage: phase3_5_harmony.py {low|high|subbucket_low}")
    run_mode(sys.argv[1])
