"""Schema v6 migration — add obs columns to support held-out cohort designs.

Session 6A Part B. Atomic schema change per Issue 17.

New obs columns introduced in v6 (all optional with explicit defaults so v1
corpus migrates non-destructively):

  donor_response_design : Categorical
      "cross_sectional" (default; v1 corpus + most natural-infection cohorts)
      "paired_within_donor" (Randolph 2021 — mock + IAV per donor)

  exposure_pair_id : str
      Unique identifier linking mock + infected samples from the same donor.
      Empty string for cross_sectional designs.

  exposure_type : Categorical
      "natural_infection" (default; v1 corpus + GSE283744 RSV/SARS + Wang CMV
        + Lee 2025 HIV — all are natural exposures even though biologically
        distinct)
      "ex_vivo_challenge" (Randolph)
      "chronic_carriage" (Wang CMV)
      "retroviral_infection" (Lee 2025 HIV — distinct from acute respiratory)

  exposure_duration_hours : float
      Hours post-exposure when the sample was collected. 6.0 for Randolph
      ex vivo. NaN for natural infection where time-from-infection is
      variable / unknown.

  age_years : float
      Continuous donor age. NaN for cohorts without detailed age metadata.

  age_group_category : Categorical
      "adult" (default; v1 corpus + Randolph + Wang CMV older adult)
      "infant_<1yr"
      "child_1to12yr"
      "adolescent_12to18yr"
      "older_adult_>65yr"

  infection_state : Categorical
      "acute" (default for symptomatic disease, including v1 corpus)
      "chronic_latent" (Wang CMV+ carriage)
      "convalescent" (post-acute recovery; reserved)
      "naive" (healthy unexposed; v1 corpus healthy donors)
      Note: donor_disease_status (diseased/healthy_control) is the primary
      label; infection_state adds biological context for non-acute designs.

  donor_serostatus : Categorical
      "positive", "negative", "unknown" (default)

Default values for v1 corpus (Wilk, Lee, Arunachalam, Schulte-Schrepping —
all acute COVID PBMC):
  - donor_response_design = "cross_sectional"
  - exposure_pair_id = ""  (empty)
  - exposure_type = "natural_infection"
  - exposure_duration_hours = NaN
  - age_years = NaN  (study metadata varies; not retroactively populated)
  - age_group_category = "adult"  (all v1 cohorts are adult by inclusion)
  - infection_state = "acute" for diseased donors, "naive" for healthy_control
  - donor_serostatus = "unknown"

Atomic per Issue 17: this migration script + the corresponding tests in
src/tests/test_schema_v6_migration.py + the METHODS_CHOICES Issue 27-30
pre-specifications land in a single commit.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# Canonical v6 schema column names + dtypes.
V6_COLUMNS: dict[str, str] = {
    "donor_response_design": "category",
    "exposure_pair_id": "string",
    "exposure_type": "category",
    "exposure_duration_hours": "float64",
    "age_years": "float64",
    "age_group_category": "category",
    "infection_state": "category",
    "donor_serostatus": "category",
}

# Allowed values for categorical columns.
V6_CATEGORIES: dict[str, tuple[str, ...]] = {
    "donor_response_design": ("cross_sectional", "paired_within_donor"),
    "exposure_type": (
        "natural_infection",
        "ex_vivo_challenge",
        "chronic_carriage",
        "retroviral_infection",
    ),
    "age_group_category": (
        "adult",
        "infant_<1yr",
        "child_1to12yr",
        "adolescent_12to18yr",
        "older_adult_>65yr",
    ),
    "infection_state": ("acute", "chronic_latent", "convalescent", "naive"),
    "donor_serostatus": ("positive", "negative", "unknown"),
}


def v1_corpus_defaults(disease_status: pd.Series) -> dict[str, pd.Series]:
    """Default values for v1 corpus (Wilk/Lee/Arunachalam/Schulte acute COVID).

    Returns a dict mapping v6 column name -> Series of values matching
    ``disease_status``'s index. The infection_state column is acute for
    diseased donors and naive for healthy_control donors.
    """
    n = len(disease_status)
    idx = disease_status.index
    diseased_mask = disease_status.astype(str) == "diseased"
    infection_state = np.where(diseased_mask, "acute", "naive")
    return {
        "donor_response_design": pd.Categorical(
            ["cross_sectional"] * n,
            categories=V6_CATEGORIES["donor_response_design"],
        ),
        "exposure_pair_id": pd.array([""] * n, dtype="string"),
        "exposure_type": pd.Categorical(
            ["natural_infection"] * n,
            categories=V6_CATEGORIES["exposure_type"],
        ),
        "exposure_duration_hours": pd.Series([np.nan] * n, index=idx, dtype="float64"),
        "age_years": pd.Series([np.nan] * n, index=idx, dtype="float64"),
        "age_group_category": pd.Categorical(
            ["adult"] * n, categories=V6_CATEGORIES["age_group_category"]
        ),
        "infection_state": pd.Categorical(
            infection_state, categories=V6_CATEGORIES["infection_state"]
        ),
        "donor_serostatus": pd.Categorical(
            ["unknown"] * n, categories=V6_CATEGORIES["donor_serostatus"]
        ),
    }


def migrate_v1_to_v6(adata: ad.AnnData, *, in_place: bool = False) -> ad.AnnData:
    """Apply v6 schema migration to a v1 AnnData (idempotent + non-destructive).

    Adds v6 obs columns with v1-corpus defaults. Existing v6 columns are
    preserved (idempotent re-runs are no-ops). All existing obs columns
    untouched. Pre-existing donor_disease_status, study_id, donor_id are
    required for default derivation.
    """
    if not in_place:
        adata = adata.copy()
    if "donor_disease_status" not in adata.obs.columns:
        raise ValueError("schema v6 migration requires donor_disease_status obs column")
    defaults = v1_corpus_defaults(adata.obs["donor_disease_status"])
    n_added = 0
    for col, val in defaults.items():
        if col in adata.obs.columns:
            continue
        adata.obs[col] = val
        n_added += 1
    if n_added:
        logger.info(
            "v6 migration added %d obs columns: %s",
            n_added,
            [c for c in defaults if c not in adata.obs.columns or c in defaults],
        )
    return adata


def has_v6_schema(adata: ad.AnnData) -> bool:
    """Return True iff all v6 schema columns are present in adata.obs."""
    return all(c in adata.obs.columns for c in V6_COLUMNS)


def migrate_file(src: str | Path, dst: str | Path | None = None) -> Path:
    """Read an h5ad, migrate to v6 schema, write to ``dst`` (default in-place)."""
    src = Path(src)
    dst = Path(dst) if dst is not None else src
    logger.info("loading %s", src)
    adata = ad.read_h5ad(src)
    migrated = migrate_v1_to_v6(adata, in_place=True)
    logger.info("writing %s", dst)
    migrated.write_h5ad(dst)
    return dst
