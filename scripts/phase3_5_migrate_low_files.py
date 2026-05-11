"""Migrate existing _reannotated.h5ad (single-model Phase 3.5 outputs) to the
dual-model naming convention: <study>_reannotated_low.h5ad.

Also stamps:
  uns['annotation_source'] = 'celltypist_immune_all_low'
  uns['celltypist_model'] = 'Immune_All_Low.pkl'
  uns['celltypist_majority_voting'] = True
  obs['cell_type_subbucket_unified']  -- finer-granularity assignment for Issue 2.

Idempotent: if <study>_reannotated_low.h5ad already exists with annotation_source
set, the script does nothing for that study.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

from trinetravir.data.bucket_map import map_label_to_subbucket_low

REPO = Path(__file__).resolve().parents[1]
PROCESSED = REPO / "data" / "processed"

STUDIES = ("wilk_2020", "arunachalam_2020", "lee_2020", "schulte_schrepping_2020")


def migrate_one(study: str) -> None:
    src = PROCESSED / f"{study}_reannotated.h5ad"
    dst = PROCESSED / f"{study}_reannotated_low.h5ad"
    if dst.exists():
        # Re-stamp uns if missing; do not rewrite obs.
        a = ad.read_h5ad(dst)
        changed = False
        if a.uns.get("annotation_source") != "celltypist_immune_all_low":
            a.uns["annotation_source"] = "celltypist_immune_all_low"
            changed = True
        if a.uns.get("celltypist_model") != "Immune_All_Low.pkl":
            a.uns["celltypist_model"] = "Immune_All_Low.pkl"
            changed = True
        if "celltypist_majority_voting" not in a.uns:
            a.uns["celltypist_majority_voting"] = True
            changed = True
        if (
            "cell_type_subbucket_unified" not in a.obs.columns
            and "cell_type_unified" in a.obs.columns
        ):
            a.obs["cell_type_subbucket_unified"] = pd.Categorical(
                a.obs["cell_type_unified"].astype(str).map(map_label_to_subbucket_low)
            )
            changed = True
        if changed:
            a.write_h5ad(dst)
            print(f"  {study}: {dst.name} updated in place")
        else:
            print(f"  {study}: {dst.name} already complete; skip")
        return
    if not src.exists():
        print(f"  {study}: NEITHER {src.name} NOR {dst.name} on disk; skipping")
        return
    print(f"  {study}: loading {src.name}")
    a = ad.read_h5ad(src)
    a.uns["annotation_source"] = "celltypist_immune_all_low"
    a.uns["celltypist_model"] = "Immune_All_Low.pkl"
    a.uns["celltypist_majority_voting"] = True
    if "cell_type_subbucket_unified" not in a.obs.columns:
        a.obs["cell_type_subbucket_unified"] = pd.Categorical(
            a.obs["cell_type_unified"].astype(str).map(map_label_to_subbucket_low)
        )
    a.write_h5ad(dst)
    print(f"  {study}: wrote {dst.name} ({a.shape})")


def main() -> int:
    for s in STUDIES:
        migrate_one(s)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
