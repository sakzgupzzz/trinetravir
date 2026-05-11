"""Unit tests for trinetravir.data.download.

Network-free. Builds tiny synthetic AnnData fixtures that mimic the
cellxgene Census schema (var has ``feature_name``, obs has ``disease`` and
``tissue_general``) and exercises every transformation step.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from trinetravir.data.download import (
    REQUIRED_OBS_COLUMNS,
    StudyConfig,
    apply_infection_status,
    build_manifest,
    filter_pbmc,
    load_dataset_config,
)


def _make_census_like_adata(
    n_cells: int = 20,
    diseases: list[str] | None = None,
    tissues: list[str] | None = None,
) -> ad.AnnData:
    rng = np.random.default_rng(0)
    diseases = diseases or ["COVID-19"] * 10 + ["normal"] * 10
    tissues = tissues or ["blood"] * n_cells
    assert len(diseases) == n_cells
    assert len(tissues) == n_cells

    X = rng.poisson(1.0, size=(n_cells, 6)).astype(np.float32)
    # Mimic Census var: integer-string soma_joinid as index, HGNC in feature_name
    var = pd.DataFrame(
        {"feature_name": ["IFI6", "IFIT1", "OAS1", "ISG15", "ACE2", "GAPDH"]},
        index=[str(i) for i in range(6)],
    )
    obs = pd.DataFrame(
        {
            "disease": diseases,
            "tissue_general": tissues,
            "cell_type": ["monocyte"] * n_cells,
            "donor_id": [f"d{i // 5}" for i in range(n_cells)],
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    return ad.AnnData(X=X, obs=obs, var=var)


def _lee_like_study() -> StudyConfig:
    return StudyConfig(
        study_id="lee_2020",
        source="cellxgene",
        accession="de2c780c-1747-40bd-9ccf-9588ec186cee",
        virus_map={"COVID-19": "sars_cov_2", "influenza": "iav"},
        infection_status_rule="disease_proxy",
        citation="Lee et al. 2020",
        notes="anchor study",
    )


# ---- filter_pbmc ----------------------------------------------------------


def test_filter_pbmc_keeps_blood_only() -> None:
    a = _make_census_like_adata(
        n_cells=4,
        diseases=["COVID-19", "COVID-19", "normal", "normal"],
        tissues=["blood", "lung", "blood", "lung"],
    )
    out = filter_pbmc(a)
    assert out.n_obs == 2
    assert (out.obs["tissue_general"] == "blood").all()


def test_filter_pbmc_raises_when_no_blood() -> None:
    a = _make_census_like_adata(
        n_cells=2, diseases=["COVID-19", "normal"], tissues=["lung", "lung"]
    )
    with pytest.raises(ValueError, match="No PBMC cells"):
        filter_pbmc(a)


def test_filter_pbmc_requires_tissue_general() -> None:
    a = _make_census_like_adata()
    del a.obs["tissue_general"]
    with pytest.raises(KeyError, match="tissue_general"):
        filter_pbmc(a)


# ---- apply_infection_status ----------------------------------------------


def test_apply_infection_status_lee_dual_virus() -> None:
    a = _make_census_like_adata(
        n_cells=6,
        diseases=["COVID-19", "COVID-19", "influenza", "influenza", "normal", "normal"],
        tissues=["blood"] * 6,
    )
    out = apply_infection_status(a, study=_lee_like_study())
    assert out.n_obs == 6
    assert out.obs["virus"].tolist() == ["sars_cov_2", "sars_cov_2", "iav", "iav", "mock", "mock"]
    assert out.obs["donor_disease_status"].tolist() == [
        "diseased",
        "diseased",
        "diseased",
        "diseased",
        "healthy",
        "healthy",
    ]
    assert (out.obs["label_source"] == "disease_proxy").all()
    assert (out.obs["study_id"] == "lee_2020").all()


def test_apply_infection_status_drops_unmapped_diseases() -> None:
    # Lee's Census record contains many unrelated disease tokens; we must drop them.
    a = _make_census_like_adata(
        n_cells=4,
        diseases=["COVID-19", "Alzheimer disease", "normal", "luminal B breast carcinoma"],
        tissues=["blood"] * 4,
    )
    out = apply_infection_status(a, study=_lee_like_study())
    assert out.n_obs == 2
    assert set(out.obs["virus"]) == {"sars_cov_2", "mock"}


def test_apply_infection_status_healthy_token_case_insensitive() -> None:
    a = _make_census_like_adata(
        n_cells=2, diseases=["NORMAL", "Healthy"], tissues=["blood", "blood"]
    )
    out = apply_infection_status(a, study=_lee_like_study())
    assert out.n_obs == 2
    assert (out.obs["donor_disease_status"] == "healthy").all()


def test_apply_infection_status_required_columns_present() -> None:
    a = _make_census_like_adata()
    out = apply_infection_status(a, study=_lee_like_study())
    for col in REQUIRED_OBS_COLUMNS:
        if col == "cell_type":
            continue  # supplied by Census, not by apply_infection_status
        assert col in out.obs.columns, f"{col} missing after annotation"


def test_apply_infection_status_empty_match_raises() -> None:
    a = _make_census_like_adata(
        n_cells=2, diseases=["Crohn disease", "ulcerative colitis"], tissues=["blood", "blood"]
    )
    with pytest.raises(ValueError, match="No cells matched"):
        apply_infection_status(a, study=_lee_like_study())


def test_apply_infection_status_unknown_rule_raises() -> None:
    bad = StudyConfig(
        study_id="x",
        source="cellxgene",
        accession="x",
        virus_map={"COVID-19": "sars_cov_2"},
        infection_status_rule="viral_read_threshold",  # not implemented in v1.1
        citation="",
        notes="",
    )
    a = _make_census_like_adata()
    with pytest.raises(NotImplementedError, match="viral_read_threshold"):
        apply_infection_status(a, study=bad)


# ---- load_dataset_config -------------------------------------------------


def test_load_dataset_config_parses_real_yaml() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    defaults, studies = load_dataset_config(repo_root / "configs" / "datasets.yaml")
    assert defaults["census_version"] == "2025-11-08"
    assert "lee_2020" in studies
    assert studies["lee_2020"].virus_map == {"COVID-19": "sars_cov_2", "influenza": "iav"}
    assert studies["wilk_2020"].virus_map == {"COVID-19": "sars_cov_2"}


# ---- build_manifest ------------------------------------------------------


def test_build_manifest_roundtrip(tmp_path: Path) -> None:
    a = _make_census_like_adata()
    out = apply_infection_status(filter_pbmc(a), study=_lee_like_study())
    out.write_h5ad(tmp_path / "lee_2020.h5ad")
    manifest = build_manifest(tmp_path)
    assert len(manifest) == 1
    row = manifest.iloc[0]
    assert row["study_id"] == "lee_2020"
    assert row["n_cells"] == out.n_obs
    assert "sars_cov_2" in row["viruses"]
    assert "mock" in row["viruses"]
    assert row["n_donors"] >= 1
