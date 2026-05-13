"""Predict-mean baseline (Issue 24 Category 1; pre-spec 2026-05-11; Phase 5 impl 2026-05-13).

Strongest predict-mean baseline: per-cell-type per-virus mean across training
cells in HVG space. Establishes the trivial lower bound that the factorized
model must beat to demonstrate value beyond per-stratum averaging.

API operates at the response-vector level (per-donor mean expression in a
(cell_type_bucket, virus, condition) stratum) since v1 calibration framework
is response-vector based per Issues 3, 8-11, 26.

For per-cell prediction (per PLAN.md spec literal reading), the per-cell
prediction = the per-(cell_type, virus) mean. Same value broadcasts to all
cells in the stratum. Cell-level eval reduces to response-vector eval.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np


class PredictMeanBaseline:
    """Predict per-(cell_type_bucket, virus) mean post-infection response from training cells.

    Fit:
      X_post: (n_train_cells, n_genes) — post-infection expression vectors
      cell_type: (n_train_cells,) — bucket assignment per cell
      virus_id: (n_train_cells,) — virus per cell

    Stores per-(bucket, virus) mean post-infection expression vectors.

    Predict:
      For each test cell, look up mean[bucket, virus] and return that vector.
    """

    def __init__(self) -> None:
        self.means_: dict[tuple[str, str], np.ndarray] = {}
        self.global_mean_: np.ndarray | None = None  # fallback for unseen (bucket, virus)

    def fit(self, X_post, cell_type, virus_id) -> PredictMeanBaseline:
        X_post = np.asarray(X_post, dtype=np.float32)
        cell_type = np.asarray(cell_type)
        virus_id = np.asarray(virus_id)

        groups: dict[tuple[str, str], list[int]] = defaultdict(list)
        for i, (ct, vi) in enumerate(zip(cell_type, virus_id, strict=False)):
            groups[(str(ct), str(vi))].append(i)

        for key, idx in groups.items():
            self.means_[key] = X_post[idx].mean(axis=0)
        self.global_mean_ = X_post.mean(axis=0)
        return self

    def predict(self, X_baseline, cell_type, virus_id) -> np.ndarray:
        """Return predicted post-infection expression per cell.

        X_baseline is unused (predict-mean ignores input cell per PLAN.md spec).
        Shape: (n_cells, n_genes).
        """
        n = len(cell_type)
        cell_type = np.asarray(cell_type)
        virus_id = np.asarray(virus_id)
        if self.global_mean_ is None:
            raise RuntimeError("Must call fit() before predict().")
        n_genes = self.global_mean_.shape[0]
        out = np.empty((n, n_genes), dtype=np.float32)
        for i, (ct, vi) in enumerate(zip(cell_type, virus_id, strict=False)):
            key = (str(ct), str(vi))
            out[i] = self.means_.get(key, self.global_mean_)
        return out
