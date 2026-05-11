"""Session 6A Part C harmonization for Yoshida 2022 cohort.

Inputs:
  data/raw/yoshida_2022/yoshida_2022_pbmc.h5ad  (cellxgene PBMC h5ad)

Outputs:
  data/processed/yoshida_2022_processed_v6.h5ad
  results/tables/yoshida_2022_bucket_mismatch_log.csv

Pipeline (per SESSION_6A_CHECKLIST C-pre.5 + C-pre.6):
  1. Load cellxgene h5ad (422,220 cells x 32,344 genes).
  2. Apply Yoshida age stratification (Issue 28 amendment aff75d4):
     - pediatric = Age_group in {Young child, Child, Adolescent}
     - adult = Age_group = Adult
     - drop Neonate + Infant + Elderly
     - drop post-COVID-19 disorder (different disease state)
  3. Map cell types: use Yoshida's annotation_broad column directly to v1's
     5-bucket scheme. CellTypist consistency check via mismatch rate:
       monocyte | CD4T | CD8T | B | NK | other
     If 'other' fraction < 5% of total cells, use direct mapping; else flag
     for re-run (CellTypist Immune_All_Low).
  4. Apply schema_v6_migration:
     - donor_response_design = cross_sectional
     - exposure_type = natural_infection
     - exposure_duration_hours = NaN
     - age_years = subject.ageAtFirstDraw OR development_stage parsing
     - age_group_category = derived from Age_group
     - infection_state = acute (COVID-19) | naive (normal)
     - donor_serostatus = unknown
   5. Save data/processed/yoshida_2022_processed_v6.h5ad.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from trinetravir.data.schema_v6_migration import V6_CATEGORIES

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw" / "yoshida_2022" / "yoshida_2022_pbmc.h5ad"
OUT = REPO / "data" / "processed" / "yoshida_2022_processed_v6.h5ad"
MISMATCH_LOG = REPO / "results" / "tables" / "yoshida_2022_bucket_mismatch_log.csv"

# Yoshida annotation_broad -> v1 5-bucket scheme.
YOSHIDA_BROAD_TO_BUCKET = {
    "T CD4+": "CD4T",
    "T CD8+": "CD8T",
    "Monocyte": "monocyte",
    "B": "B",
    "NK": "NK",
    "T reg": "CD4T",
    "T g/d": "other",
    "DC": "other",
    "MAIT": "CD8T",
    "Cycling": "other",
    "Platelets": "other",
    "Plasma": "B",
    "HPC": "other",
    "ILC": "other",
    "RBC": "other",
    "Baso/Eos": "other",
}

PEDIATRIC_AGE_GROUPS = {"Young child", "Child", "Adolescent"}
ADULT_AGE_GROUPS = {"Adult"}


def derive_age_group_category(age_group: str) -> str:
    if age_group in PEDIATRIC_AGE_GROUPS:
        # Use v1 schema's pediatric categories. Map Yoshida's broad pediatric
        # labels to the v1 categorical: child_1to12yr by default; adolescent_12to18yr
        # for "Adolescent" specifically.
        if age_group == "Adolescent":
            return "adolescent_12to18yr"
        return "child_1to12yr"
    if age_group in ADULT_AGE_GROUPS:
        return "adult"
    if age_group in {"Neonate", "Infant"}:
        return "infant_<1yr"
    if age_group == "Elderly":
        return "older_adult_>65yr"
    return "adult"  # fallback (shouldn't happen)


def main() -> int:
    logger.info("loading %s", RAW)
    a = ad.read_h5ad(RAW)
    n_raw = a.n_obs
    logger.info("raw shape: %s", a.shape)

    # ---- Primary stratification filter (Issue 28 amendment aff75d4) ----
    keep_age = a.obs["Age_group"].astype(str).isin(PEDIATRIC_AGE_GROUPS | ADULT_AGE_GROUPS)
    keep_disease = a.obs["disease"].astype(str).isin(["normal", "COVID-19"])
    keep = keep_age & keep_disease
    logger.info(
        "dropping %d cells: age=%d, disease=%d",
        int((~keep).sum()),
        int((~keep_age).sum()),
        int((~keep_disease).sum()),
    )
    a = a[keep].copy()
    logger.info("after primary stratification filter: %d cells", a.n_obs)

    # ---- Bucket assignment via Yoshida annotation_broad ----
    bk = a.obs["annotation_broad"].astype(str).map(YOSHIDA_BROAD_TO_BUCKET)
    if bk.isna().any():
        unmapped = a.obs.loc[bk.isna(), "annotation_broad"].astype(str).value_counts()
        logger.warning("unmapped annotation_broad labels: %s", dict(unmapped))
        bk = bk.fillna("other")
    a.obs["cell_type_bucket_unified"] = pd.Categorical(
        bk.values, categories=["monocyte", "CD4T", "CD8T", "B", "NK", "other"]
    )

    # ---- Cross-corpus consistency check (C-pre.5) ----
    bucket_counts = a.obs["cell_type_bucket_unified"].value_counts()
    other_count = int(bucket_counts.get("other", 0))
    total = int(a.n_obs)
    other_pct = 100.0 * other_count / total if total else 0.0
    logger.info(
        "bucket counts: %s | 'other' fraction: %.2f%%",
        dict(bucket_counts),
        other_pct,
    )
    consistency_pass = other_pct < 5.0

    # Mismatch log: each annotation_broad value -> bucket, count, pct
    mismatch_rows = []
    for broad_label, count in a.obs["annotation_broad"].astype(str).value_counts().items():
        bucket = YOSHIDA_BROAD_TO_BUCKET.get(broad_label, "other")
        mismatch_rows.append(
            {
                "yoshida_annotation_broad": broad_label,
                "v1_bucket": bucket,
                "n_cells": int(count),
                "pct_of_total": round(100.0 * count / total, 3),
                "routed_to_other": bucket == "other",
            }
        )
    mismatch_df = pd.DataFrame(mismatch_rows).sort_values("n_cells", ascending=False)
    MISMATCH_LOG.parent.mkdir(parents=True, exist_ok=True)
    mismatch_df.to_csv(MISMATCH_LOG, index=False)
    logger.info("wrote %s", MISMATCH_LOG)
    logger.info(
        "C-pre.5 consistency check: other_pct=%.2f%% (threshold 5%%); PASS=%s",
        other_pct,
        consistency_pass,
    )

    # ---- Apply schema v6 ----
    n = a.n_obs

    # Build v6 obs columns. donor_disease_status from disease.
    a.obs["donor_disease_status"] = pd.Categorical(
        np.where(a.obs["disease"].astype(str) == "COVID-19", "diseased", "healthy_control"),
        categories=["diseased", "healthy_control"],
    )

    # study_id is required for harmonize / calibration framework
    a.obs["study_id"] = "yoshida_2022"
    a.obs["study_id"] = a.obs["study_id"].astype("category")

    a.obs["donor_response_design"] = pd.Categorical(
        ["cross_sectional"] * n, categories=V6_CATEGORIES["donor_response_design"]
    )
    a.obs["exposure_pair_id"] = pd.array([""] * n, dtype="string")
    a.obs["exposure_type"] = pd.Categorical(
        ["natural_infection"] * n, categories=V6_CATEGORIES["exposure_type"]
    )
    a.obs["exposure_duration_hours"] = np.full(n, np.nan, dtype=np.float64)

    # age_years: parse from development_stage where possible (e.g. "25-year-old stage")
    ds = a.obs["development_stage"].astype(str)

    def _parse_age(s: str) -> float:
        import re

        m = re.match(r"^(\d+)-year-old stage$", s)
        if m:
            return float(m.group(1))
        return float("nan")

    a.obs["age_years"] = ds.map(_parse_age).astype(np.float64).values

    a.obs["age_group_category"] = pd.Categorical(
        a.obs["Age_group"].astype(str).map(derive_age_group_category).values,
        categories=V6_CATEGORIES["age_group_category"],
    )

    a.obs["infection_state"] = pd.Categorical(
        np.where(a.obs["disease"].astype(str) == "COVID-19", "acute", "naive"),
        categories=V6_CATEGORIES["infection_state"],
    )

    a.obs["donor_serostatus"] = pd.Categorical(
        ["unknown"] * n, categories=V6_CATEGORIES["donor_serostatus"]
    )

    # Stamp annotation provenance
    a.uns["annotation_source"] = "yoshida_2022_annotation_broad_direct_mapping"
    a.uns["consistency_check_other_pct"] = other_pct
    a.uns["consistency_check_pass"] = bool(consistency_pass)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    logger.info("writing %s", OUT)
    a.write_h5ad(OUT)
    logger.info(
        "done. processed shape %s. cells dropped from raw: %d (-> %d)",
        a.shape,
        n_raw - a.n_obs,
        a.n_obs,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
