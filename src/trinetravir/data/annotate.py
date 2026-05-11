"""Phase 3.5 unified re-annotation of PBMC studies with CellTypist.

Addresses the annotation divergence surfaced by the Phase 3 gate
diagnostic (Lee + wilk lacked memory/naive subdivisions for lymphoid
lineages, breaking cross-study B + CD8T comparability).

Strategy:
- Run CellTypist with the ``Immune_All_Low.pkl`` pretrained model.
- Use ``majority_voting=True`` so labels are assigned per cluster
  rather than per cell (reduces noise for lymphoid sublineages where
  per-cell predictions are unstable).
- CellTypist expects log-normalized counts at target_sum=1e4. Our
  cellxgene Census-derived h5ads carry raw integer counts in ``.X``;
  this module normalizes a copy in place, runs prediction, then writes
  three new ``obs`` columns to the returned AnnData:
    cell_type_unified         predicted label (string from Immune_All_Low)
    cell_type_bucket_unified  one of 5 coarse buckets or 'other'
    cell_type_original        snapshot of the pre-existing cell_type col
  ``.X`` is *not* mutated on the returned AnnData; the normalization is
  done on a working copy passed to CellTypist.

Output is persisted to ``data/processed/<study>_reannotated.h5ad``.
``data/raw/`` stays immutable per project rules.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

from trinetravir.data.bucket_map import (
    is_unexpected_label,
    is_unexpected_label_high,
    map_label_to_bucket,
    map_label_to_bucket_high,
    map_label_to_subbucket_low,
)

_MODEL_TAG = {
    "Immune_All_Low.pkl": "celltypist_immune_all_low",
    "Immune_All_High.pkl": "celltypist_immune_all_high",
}

logger = logging.getLogger(__name__)


def annotate_unified(
    adata: ad.AnnData,
    model_name: str = "Immune_All_Low.pkl",
    majority_voting: bool = True,
) -> ad.AnnData:
    """Run CellTypist on ``adata`` and attach unified-label obs columns.

    Parameters
    ----------
    adata
        Per-study AnnData with raw integer counts in ``.X`` and existing
        ``cell_type`` obs column (cellxgene ontology label).
    model_name
        CellTypist model filename. Default ``Immune_All_Low.pkl`` — coarse
        immune labels with consistent vocabulary across studies. Switch
        to ``Immune_All_High.pkl`` only if granularity is insufficient.
    majority_voting
        Pass through to ``celltypist.annotate``. True = over-cluster +
        per-cluster majority vote; recommended for lymphoid sublineages.

    Returns
    -------
    adata_out
        Shallow copy of input with three new obs columns:
        ``cell_type_unified``, ``cell_type_bucket_unified``,
        ``cell_type_original``. ``.X`` unchanged.
    """
    # Local import — celltypist pulls scanpy.__version__ deprecation
    # warnings at module load; we want the warning emitted from the
    # function call rather than at module import time so callers can
    # silence it locally if they want.
    import celltypist

    if "cell_type" not in adata.obs.columns:
        raise ValueError("adata must have an existing 'cell_type' obs column")

    # Working copy for normalization. CellTypist mutates .X.
    work = adata.copy()
    work.X = work.X.astype(np.float32)
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)

    logger.info(
        "running celltypist (%s, majority_voting=%s) on %d cells",
        model_name,
        majority_voting,
        work.n_obs,
    )
    predictions = celltypist.annotate(
        work,
        model=model_name,
        majority_voting=majority_voting,
    )
    pred_adata = predictions.to_adata()

    # majority_voting=True adds 'majority_voting' obs col with cluster-level
    # majority label. predicted_labels is the raw per-cell label. We use
    # the majority-voted label as cell_type_unified when available.
    if "majority_voting" in pred_adata.obs.columns and majority_voting:
        unified = pred_adata.obs["majority_voting"].astype(str)
    else:
        unified = pred_adata.obs["predicted_labels"].astype(str)

    # Sanity: index alignment between work and adata.
    if not unified.index.equals(adata.obs.index):
        raise RuntimeError("celltypist output index does not match input adata index")

    if model_name == "Immune_All_High.pkl":
        bucket_fn = map_label_to_bucket_high
        unexpected_fn = is_unexpected_label_high
        subbucket_fn = None
    else:
        bucket_fn = map_label_to_bucket
        unexpected_fn = is_unexpected_label
        subbucket_fn = map_label_to_subbucket_low

    adata_out = adata.copy()
    adata_out.obs["cell_type_original"] = adata.obs["cell_type"].astype(str)
    adata_out.obs["cell_type_unified"] = unified.astype("category")
    adata_out.obs["cell_type_bucket_unified"] = unified.map(bucket_fn).astype("category")
    if subbucket_fn is not None:
        adata_out.obs["cell_type_subbucket_unified"] = unified.map(subbucket_fn).astype("category")

    # Log unexpected labels so we can audit and extend the bucket map.
    unexpected = sorted({lbl for lbl in unified.unique() if unexpected_fn(lbl)})
    if unexpected:
        logger.warning(
            "unexpected celltypist labels (mapped to 'other'): %s",
            unexpected,
        )

    # Stamp uns so downstream loaders know which model produced these labels.
    adata_out.uns["annotation_source"] = _MODEL_TAG.get(model_name, model_name)
    adata_out.uns["celltypist_model"] = model_name
    adata_out.uns["celltypist_majority_voting"] = bool(majority_voting)

    # Log per-bucket counts for quick sanity.
    counts = adata_out.obs["cell_type_bucket_unified"].value_counts()
    logger.info("unified bucket counts: %s", dict(counts))
    return adata_out


def annotate_and_save(
    raw_path: str | Path,
    out_path: str | Path,
    model_name: str = "Immune_All_Low.pkl",
    majority_voting: bool = True,
) -> ad.AnnData:
    """Convenience wrapper: load -> annotate_unified -> write h5ad.

    Parameters
    ----------
    raw_path
        Path to input ``<study>.h5ad`` (in ``data/raw/``).
    out_path
        Path to output ``<study>_reannotated.h5ad`` (in ``data/processed/``).
    model_name
        Passed to ``annotate_unified``.
    majority_voting
        Passed to ``annotate_unified``.

    Returns
    -------
    adata_out
        Annotated AnnData (also written to disk).
    """
    raw_path = Path(raw_path)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("loading %s", raw_path)
    adata = ad.read_h5ad(raw_path)
    adata_out = annotate_unified(adata, model_name=model_name, majority_voting=majority_voting)
    logger.info("writing %s", out_path)
    adata_out.write_h5ad(out_path)
    return adata_out
