"""Session 6A Part C harmonization for Allen Institute Immune Health Atlas.

Inputs:
  data/raw/allen_atlas/human_immune_health_atlas_mono.h5ad  (monocyte bucket only)

Outputs:
  data/processed/allen_atlas_monocyte_processed_v6.h5ad

Pipeline (per SESSION_6A_CHECKLIST C-pre.6):
  1. Load Allen Atlas monocyte h5ad (327,919 cells x 33,538 genes).
  2. AIFI_L1 = "Monocyte" for all cells. cell_type_bucket_unified = "monocyte" uniformly.
  3. Apply schema_v6_migration:
     - donor_response_design = cross_sectional
     - exposure_type = chronic_carriage
     - infection_state = chronic_latent (CMV+) | naive (CMV-)
     - donor_serostatus = positive (CMV+) | negative (CMV-)
     - age_group_category from subject.ageGroup
     - age_years from subject.ageAtFirstDraw
  4. donor_disease_status = diseased (CMV+) | healthy_control (CMV-)
     Per Issue 29 amendment: CMV+ chronic latent carriers are the "diseased"
     side of the discrimination test; CMV- are the naive controls.
  5. Save data/processed/allen_atlas_monocyte_processed_v6.h5ad
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
RAW = REPO / "data" / "raw" / "allen_atlas" / "human_immune_health_atlas_mono.h5ad"
OUT = REPO / "data" / "processed" / "allen_atlas_monocyte_processed_v6.h5ad"


def derive_age_group_category(age_group_label: str) -> str:
    """Allen Atlas subject.ageGroup -> v1 schema age_group_category."""
    if age_group_label == "Young Adult":
        return "adult"
    if age_group_label == "Older Adult":
        return "older_adult_>65yr"
    if age_group_label == "Children":
        return "child_1to12yr"
    return "adult"  # fallback


def main() -> int:
    logger.info("loading %s", RAW)
    a = ad.read_h5ad(RAW)
    n_raw = a.n_obs
    logger.info("raw shape: %s", a.shape)

    # ---- Filter: keep adults only for Issue 29 primary analysis ----
    # Children stratum has only 3 CMV+ subjects (below Issue 4 ≥4 minimum).
    # Per Issue 29 amendment: Young Adult + Older Adult combined for primary analysis.
    keep = a.obs["subject.ageGroup"].astype(str).isin(["Young Adult", "Older Adult"])
    dropped = int((~keep).sum())
    logger.info("dropping %d Children cells (below Issue 4 for CMV+ stratum)", dropped)
    a = a[keep].copy()
    logger.info("after adult-only filter: %d cells", a.n_obs)

    # ---- Bucket assignment: all monocyte ----
    n = a.n_obs
    a.obs["cell_type_bucket_unified"] = pd.Categorical(
        ["monocyte"] * n,
        categories=["monocyte", "CD4T", "CD8T", "B", "NK", "other"],
    )

    # ---- donor_disease_status from subject.cmv (Issue 29 mapping) ----
    cmv = a.obs["subject.cmv"].astype(str)
    a.obs["donor_disease_status"] = pd.Categorical(
        np.where(cmv == "Positive", "diseased", "healthy_control"),
        categories=["diseased", "healthy_control"],
    )

    # ---- study_id + donor_id ----
    a.obs["study_id"] = "allen_atlas"
    a.obs["study_id"] = a.obs["study_id"].astype("category")
    a.obs["donor_id"] = a.obs["subject.subjectGuid"].astype(str)
    a.obs["donor_id"] = a.obs["donor_id"].astype("category")

    # ---- Schema v6 obs ----
    a.obs["donor_response_design"] = pd.Categorical(
        ["cross_sectional"] * n, categories=V6_CATEGORIES["donor_response_design"]
    )
    a.obs["exposure_pair_id"] = pd.array([""] * n, dtype="string")
    a.obs["exposure_type"] = pd.Categorical(
        ["chronic_carriage"] * n, categories=V6_CATEGORIES["exposure_type"]
    )
    a.obs["exposure_duration_hours"] = np.full(n, np.nan, dtype=np.float64)

    a.obs["age_years"] = a.obs["subject.ageAtFirstDraw"].astype(np.float64).values

    a.obs["age_group_category"] = pd.Categorical(
        a.obs["subject.ageGroup"].astype(str).map(derive_age_group_category).values,
        categories=V6_CATEGORIES["age_group_category"],
    )

    a.obs["infection_state"] = pd.Categorical(
        np.where(cmv == "Positive", "chronic_latent", "naive"),
        categories=V6_CATEGORIES["infection_state"],
    )

    a.obs["donor_serostatus"] = pd.Categorical(
        np.where(
            cmv == "Positive",
            "positive",
            np.where(cmv == "Negative", "negative", "unknown"),
        ),
        categories=V6_CATEGORIES["donor_serostatus"],
    )

    # Stamp provenance
    a.uns["annotation_source"] = "allen_atlas_AIFI_L1_monocyte_bucket_file"
    a.uns["bucket_column"] = "cell_type_bucket_unified"

    # Inventory
    donor_cmv = a.obs[["donor_id", "donor_serostatus"]].drop_duplicates()
    pos_n = int((donor_cmv["donor_serostatus"].astype(str) == "positive").sum())
    neg_n = int((donor_cmv["donor_serostatus"].astype(str) == "negative").sum())
    logger.info("donor CMV serostatus: %d positive / %d negative (Issue 4 >=4/>=4)", pos_n, neg_n)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    logger.info("writing %s", OUT)
    a.write_h5ad(OUT)
    logger.info(
        "done. processed shape %s. cells dropped from raw: %d -> %d",
        a.shape,
        n_raw - a.n_obs,
        a.n_obs,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
