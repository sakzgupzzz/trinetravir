"""Dataset download orchestration for cellxgene Census and GEO sources.

Phase 2 entry point. Loads PBMC scRNA-seq datasets defined in
``configs/datasets.yaml`` into a uniform AnnData layout with HGNC gene symbols
and study-specific ``virus`` + ``donor_disease_status`` annotations.

Schema reality check (v1.1, PBMC scope): ``donor_disease_status`` carries
values ``diseased`` / ``healthy_control``. It is a *donor-level* label
inferred from the disease ontology, not a per-cell viral-read measurement.
``mock_control`` is reserved for future in-vitro mock-infected studies
(no PBMC study in v1 produces that value). The ``infected`` /
``bystander`` / ``mock`` per-cell semantics from PLAN §2 belong to
airway-epithelium studies and are out of scope until v2.

See METHODS_CHOICES.md Issue 1 for the rationale behind the label
vocabulary (donor-level proxy, not cell-autonomous infection state).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import cellxgene_census
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

# Census version is pinned in the dataset config; this module never reads a
# floating "stable" pointer to keep the loader reproducible. The pin lives in
# configs/datasets.yaml under defaults.census_version.

REQUIRED_OBS_COLUMNS = (
    "virus",
    "donor_disease_status",
    "label_source",
    "cell_type",
    "donor_id",
    "study_id",
)


@dataclass(frozen=True)
class StudyConfig:
    """Resolved per-study configuration drawn from ``configs/datasets.yaml``.

    Attributes
    ----------
    study_id
        Short identifier used as filename and ``study_id`` obs column value.
    source
        ``cellxgene`` or ``geo``. Drives which loader path is taken.
    accession
        Census ``dataset_id`` for cellxgene; GEO accession (e.g. GSE158055) for
        GEO sources.
    virus_map
        Map from disease-ontology label (or other obs token) to the canonical
        virus label written to ``adata.obs['virus']``.
    infection_status_rule
        How to derive cell labels. Only ``disease_proxy`` is implemented in
        v1.1: cells from a virally-infected donor are labelled
        ``donor_disease_status='diseased'``, cells from a healthy donor are
        ``donor_disease_status='healthy_control'``. The rule name retains the
        historical ``infection_status_rule`` key for forward-compat with the
        v2 cell-level rule (``viral_read_threshold``) which will populate a
        separate ``infection_status`` column for airway-epithelium studies.
    citation
        Free-text human-readable citation.
    notes
        Free-text caveats (e.g. severity-strata, longitudinal timepoints).
    excluded
        If True, ``--all`` skips this study. Direct ``--study-id`` requests
        still load it (override path for ablation / re-evaluation).
    exclusion_reason
        Free-text reason recorded when ``excluded=True``. Required when
        excluded; checked at config-load time.
    """

    study_id: str
    source: str
    accession: str
    virus_map: dict[str, str]
    infection_status_rule: str
    citation: str
    notes: str
    excluded: bool = False
    exclusion_reason: str = ""


def load_dataset_config(config_path: str | Path) -> tuple[dict[str, Any], dict[str, StudyConfig]]:
    """Parse ``configs/datasets.yaml`` into defaults dict + study map.

    Parameters
    ----------
    config_path
        Path to ``configs/datasets.yaml``.

    Returns
    -------
    defaults
        Top-level ``defaults`` block (QC thresholds, census_version, etc).
    studies
        Mapping ``study_id -> StudyConfig``.
    """
    with open(config_path) as fh:
        raw = yaml.safe_load(fh)

    defaults = raw.get("defaults", {})
    studies_raw = raw.get("studies", {}) or {}

    studies: dict[str, StudyConfig] = {}
    for study_id, payload in studies_raw.items():
        excluded = bool(payload.get("excluded", False))
        exclusion_reason = str(payload.get("exclusion_reason", "")).strip()
        if excluded and not exclusion_reason:
            raise ValueError(
                f"{study_id}: excluded=true requires exclusion_reason in datasets.yaml"
            )
        studies[study_id] = StudyConfig(
            study_id=study_id,
            source=payload["source"],
            accession=payload["accession"],
            virus_map=payload.get("virus_map", {}),
            infection_status_rule=payload.get("infection_status_rule", "disease_proxy"),
            citation=payload.get("citation", ""),
            notes=payload.get("notes", ""),
            excluded=excluded,
            exclusion_reason=exclusion_reason,
        )
    return defaults, studies


def load_cellxgene_dataset(
    dataset_id: str,
    census_version: str,
    *,
    organism: str = "Homo sapiens",
    obs_value_filter: str | None = None,
) -> ad.AnnData:
    """Pull a single Census dataset and rebind ``var_names`` to HGNC symbols.

    cellxgene Census returns AnnData with ``var_names`` set to soma_joinid
    integers, which is correct for the storage layer but unusable downstream
    where every model expects HGNC gene symbols. This loader is the single
    point in the codebase that knows about that quirk.

    Parameters
    ----------
    dataset_id
        Census ``dataset_id`` (UUID string).
    census_version
        Pinned Census release (e.g. ``"2025-11-08"``). Never pass ``"stable"``
        or ``"latest"`` — that breaks reproducibility.
    organism
        Census organism label. Defaults to human.
    obs_value_filter
        Optional extra filter pushed into the Census query. The dataset_id
        filter is always added on top of whatever is passed here.

    Returns
    -------
    adata
        AnnData with ``var_names`` set to HGNC ``feature_name`` strings.

    Raises
    ------
    AssertionError
        If any feature_name is null after the rebind. That would silently
        corrupt downstream gene-set lookups (ISG masks, viral receptor masks),
        so we fail loudly instead.
    """
    base_filter = f"dataset_id == '{dataset_id}'"
    full_filter = f"({base_filter}) and ({obs_value_filter})" if obs_value_filter else base_filter

    logger.info("Opening Census %s for dataset %s", census_version, dataset_id)
    with cellxgene_census.open_soma(census_version=census_version) as census:
        adata = cellxgene_census.get_anndata(
            census=census,
            organism=organism,
            obs_value_filter=full_filter,
        )

    if "feature_name" not in adata.var.columns:
        raise RuntimeError(
            f"Census schema regression: 'feature_name' missing from var for {dataset_id}"
        )
    feature_names = adata.var["feature_name"].astype(str)
    assert feature_names.notna().all(), (
        f"NaN feature_name found in {dataset_id}; would corrupt gene lookups"
    )
    adata.var_names = pd.Index(feature_names.values)
    adata.var_names_make_unique()  # rare duplicates exist in Census; suffix them rather than crash

    logger.info("Loaded %s: %d cells x %d genes", dataset_id, adata.n_obs, adata.n_vars)
    return adata


def filter_pbmc(adata: ad.AnnData) -> ad.AnnData:
    """Subset to peripheral-blood cells.

    Uses the ``tissue_general`` ontology column from cellxgene, which collapses
    fine-grained tissue terms ("peripheral blood mononuclear cell",
    "venous blood", etc) into a single ``blood`` token.

    Parameters
    ----------
    adata
        Input AnnData with cellxgene-style obs columns.

    Returns
    -------
    adata_pbmc
        View / copy filtered to blood cells only.
    """
    if "tissue_general" not in adata.obs.columns:
        raise KeyError("tissue_general missing from obs; not a Census-loaded AnnData")
    mask = adata.obs["tissue_general"].astype(str).str.lower().eq("blood")
    n_before, n_after = adata.n_obs, int(mask.sum())
    if n_after == 0:
        raise ValueError("No PBMC cells after tissue_general='blood' filter")
    if n_after < n_before:
        logger.info(
            "PBMC filter: %d -> %d cells (dropped %d non-blood)",
            n_before,
            n_after,
            n_before - n_after,
        )
    return adata[mask].copy()


def apply_infection_status(
    adata: ad.AnnData,
    *,
    study: StudyConfig,
) -> ad.AnnData:
    """Annotate ``adata.obs`` with canonical ``virus`` and ``donor_disease_status``.

    For v1.1 the only implemented rule is ``disease_proxy``. A cell's virus
    is derived from ``adata.obs['disease']`` via ``study.virus_map``;
    ``donor_disease_status`` is ``diseased`` for cells whose disease maps to
    a virus and ``healthy_control`` for cells whose disease is ``normal`` /
    ``healthy``. Cells with diseases not in ``virus_map`` and not in the
    healthy token set are dropped (Census records often contain unrelated
    diseases like Alzheimer or breast carcinoma cells from atlas merges).

    The ``donor_disease_status`` label is *donor-level*, not cell-level.
    PBMCs from a COVID donor are mostly NOT directly infected by SARS-CoV-2;
    the signal we are modelling is systemic cytokine / interferon response,
    not cell-autonomous infection. The framing matters for the eventual
    paper: this is cross-virus transfer of donor-level disease response in
    peripheral blood, not direct cellular infection. Per-cell
    ``infection_status`` (infected / bystander / mock) requires
    airway-epithelium datasets where viral reads can be assigned per cell —
    out of scope until v2.

    The ``label_source`` obs column records that the labels were derived
    via the disease-proxy rule, so downstream code can branch on label
    semantics without re-inferring from study metadata.

    Parameters
    ----------
    adata
        AnnData with cellxgene ``disease`` obs column.
    study
        Resolved StudyConfig.

    Returns
    -------
    adata
        Copy with ``virus``, ``donor_disease_status``, ``label_source``, and
        ``study_id`` obs columns added.
    """
    if study.infection_status_rule != "disease_proxy":
        raise NotImplementedError(
            f"infection_status_rule={study.infection_status_rule!r} not implemented in v1.1"
        )
    if "disease" not in adata.obs.columns:
        raise KeyError("disease column missing from obs; cannot apply disease_proxy rule")

    disease = adata.obs["disease"].astype(str)

    healthy_tokens = {"normal", "healthy"}
    virus_assignment = pd.Series(pd.NA, index=adata.obs.index, dtype="object")
    status_assignment = pd.Series(pd.NA, index=adata.obs.index, dtype="object")

    for disease_token, virus_label in study.virus_map.items():
        match = disease.str.lower().eq(disease_token.lower())
        virus_assignment = virus_assignment.mask(match, virus_label)
        status_assignment = status_assignment.mask(match, "diseased")

    healthy_match = disease.str.lower().isin(healthy_tokens)
    virus_assignment = virus_assignment.mask(healthy_match, "mock")
    status_assignment = status_assignment.mask(healthy_match, "healthy_control")

    keep = virus_assignment.notna()
    n_keep = int(keep.sum())
    if n_keep == 0:
        raise ValueError(
            f"No cells matched virus_map={study.virus_map} or healthy tokens for {study.study_id}"
        )

    out = adata[keep].copy()
    out.obs["virus"] = virus_assignment[keep].astype("category")
    out.obs["donor_disease_status"] = status_assignment[keep].astype("category")
    out.obs["label_source"] = pd.Categorical([study.infection_status_rule] * n_keep)
    out.obs["study_id"] = pd.Categorical([study.study_id] * n_keep)
    return out


def download_study(
    study: StudyConfig,
    *,
    census_version: str,
    out_dir: str | Path,
    overwrite: bool = False,
) -> Path:
    """End-to-end download for one study: Census fetch -> PBMC filter -> annotate -> persist.

    Parameters
    ----------
    study
        Resolved StudyConfig.
    census_version
        Pinned Census release.
    out_dir
        Directory to write ``{study_id}.h5ad`` into.
    overwrite
        If False (default), skip studies already present on disk.

    Returns
    -------
    out_path
        Path to the persisted ``.h5ad``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{study.study_id}.h5ad"

    if out_path.exists() and not overwrite:
        logger.info("%s already exists; skipping (overwrite=False)", out_path)
        return out_path

    if study.source != "cellxgene":
        raise NotImplementedError(
            f"source={study.source!r} not implemented; v1.1 only supports cellxgene Census"
        )

    adata = load_cellxgene_dataset(study.accession, census_version=census_version)
    adata = filter_pbmc(adata)
    adata = apply_infection_status(adata, study=study)

    missing = [c for c in REQUIRED_OBS_COLUMNS if c not in adata.obs.columns]
    if missing:
        raise RuntimeError(f"Required obs columns missing after annotation: {missing}")

    adata.write_h5ad(out_path, compression="gzip")
    logger.info(
        "Wrote %s (%d cells, %.1f MB)", out_path, adata.n_obs, out_path.stat().st_size / 1e6
    )
    return out_path


def build_manifest(out_dir: str | Path) -> pd.DataFrame:
    """Scan ``data/raw/`` and produce a provenance manifest of downloaded studies.

    Parameters
    ----------
    out_dir
        Directory containing ``{study_id}.h5ad`` files.

    Returns
    -------
    manifest
        DataFrame with one row per study: study_id, n_cells, viruses, file_size_mb.
    """
    out_dir = Path(out_dir)
    rows: list[dict[str, Any]] = []
    for h5ad in sorted(out_dir.glob("*.h5ad")):
        adata = ad.read_h5ad(h5ad, backed="r")
        viruses = (
            sorted(set(adata.obs["virus"].astype(str))) if "virus" in adata.obs.columns else []
        )
        n_donors = int(adata.obs["donor_id"].nunique()) if "donor_id" in adata.obs.columns else 0
        rows.append(
            {
                "study_id": h5ad.stem,
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "n_donors": n_donors,
                "viruses": ",".join(viruses),
                "file_size_mb": round(h5ad.stat().st_size / 1e6, 1),
            }
        )
    return pd.DataFrame(rows)
