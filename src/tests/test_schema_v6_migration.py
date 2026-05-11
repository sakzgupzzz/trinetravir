"""Tests for src/trinetravir/data/schema_v6_migration.py (Session 6A Part B5)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from trinetravir.data.schema_v6_migration import (
    V6_CATEGORIES,
    V6_COLUMNS,
    has_v6_schema,
    migrate_v1_to_v6,
    v1_corpus_defaults,
)


def _make_v1_adata(n_cells: int = 50, n_diseased: int = 30) -> ad.AnnData:
    """Build a tiny synthetic v1 AnnData with the minimum required obs columns."""
    rng = np.random.default_rng(0)
    obs = pd.DataFrame(
        {
            "study_id": pd.Categorical(["lee_2020"] * n_cells),
            "donor_id": pd.Categorical([f"donor_{i // 5}" for i in range(n_cells)]),
            "donor_disease_status": pd.Categorical(
                ["diseased"] * n_diseased + ["healthy_control"] * (n_cells - n_diseased)
            ),
        }
    )
    return ad.AnnData(X=rng.normal(size=(n_cells, 10)).astype(np.float32), obs=obs)


def test_migration_adds_all_v6_columns() -> None:
    a = _make_v1_adata()
    assert not has_v6_schema(a)
    migrated = migrate_v1_to_v6(a)
    assert has_v6_schema(migrated)
    for col in V6_COLUMNS:
        assert col in migrated.obs.columns, f"missing v6 column {col}"


def test_migration_preserves_original_obs() -> None:
    a = _make_v1_adata()
    n_obs_cols_before = len(a.obs.columns)
    pre = a.obs.copy()
    migrated = migrate_v1_to_v6(a)
    # Original columns untouched + only the v6 columns added.
    for col in pre.columns:
        assert col in migrated.obs.columns
        # Categorical equality
        assert (migrated.obs[col].astype(str) == pre[col].astype(str)).all()
    assert len(migrated.obs.columns) == n_obs_cols_before + len(V6_COLUMNS)


def test_migration_idempotent() -> None:
    a = _make_v1_adata()
    m1 = migrate_v1_to_v6(a)
    m2 = migrate_v1_to_v6(m1)
    # Idempotent: no new columns or shape changes.
    assert m1.shape == m2.shape
    assert list(m1.obs.columns) == list(m2.obs.columns)


def test_v1_corpus_defaults_infection_state_split() -> None:
    """diseased -> acute, healthy_control -> naive."""
    status = pd.Series(["diseased", "diseased", "healthy_control", "healthy_control"])
    d = v1_corpus_defaults(status)
    inf = list(d["infection_state"])
    assert inf == ["acute", "acute", "naive", "naive"]


def test_migrated_columns_are_categorical_with_correct_categories() -> None:
    a = _make_v1_adata()
    m = migrate_v1_to_v6(a)
    for col, cats in V6_CATEGORIES.items():
        series = m.obs[col]
        # pandas Categorical
        assert hasattr(series, "cat"), f"{col} should be Categorical"
        assert set(series.cat.categories) == set(cats), (
            f"{col} categories mismatch: got {sorted(series.cat.categories)} expected {sorted(cats)}"
        )


def test_migration_requires_donor_disease_status() -> None:
    """Migration on an AnnData missing donor_disease_status should raise."""
    rng = np.random.default_rng(0)
    obs = pd.DataFrame({"study_id": ["a"] * 10})
    a = ad.AnnData(X=rng.normal(size=(10, 5)).astype(np.float32), obs=obs)
    try:
        migrate_v1_to_v6(a)
    except ValueError:
        return
    raise AssertionError("expected ValueError for missing donor_disease_status")


def test_exposure_duration_hours_is_nan_for_v1_defaults() -> None:
    a = _make_v1_adata()
    m = migrate_v1_to_v6(a)
    assert m.obs["exposure_duration_hours"].isna().all(), (
        "v1 corpus should have all-NaN exposure_duration_hours (natural infection, no defined exposure time)"
    )


def test_age_group_category_default_is_adult() -> None:
    a = _make_v1_adata()
    m = migrate_v1_to_v6(a)
    assert (m.obs["age_group_category"].astype(str) == "adult").all()


def test_donor_serostatus_default_is_unknown() -> None:
    a = _make_v1_adata()
    m = migrate_v1_to_v6(a)
    assert (m.obs["donor_serostatus"].astype(str) == "unknown").all()


def test_donor_response_design_default_is_cross_sectional() -> None:
    a = _make_v1_adata()
    m = migrate_v1_to_v6(a)
    assert (m.obs["donor_response_design"].astype(str) == "cross_sectional").all()


def test_exposure_type_default_is_natural_infection() -> None:
    a = _make_v1_adata()
    m = migrate_v1_to_v6(a)
    assert (m.obs["exposure_type"].astype(str) == "natural_infection").all()
