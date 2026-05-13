"""KNN baseline (Issue 24 Category 1; pre-spec 2026-05-11; Phase 5 impl 2026-05-13).

Per PLAN.md Phase 5 spec: "for each test baseline cell, find K nearest training
baseline cells (K = 25, 50, 100), return their post-infection mean."

Issue 24 spec: cosine distance on log-normalized HVG, k=10 default with distance-
weighted predictions. Within-virus training neighbors only for cross-virus
evaluation per Issue 15 protocol.
"""

from __future__ import annotations

import numpy as np
from sklearn.neighbors import NearestNeighbors


class KNNBaseline:
    """Cosine-distance KNN with distance-weighted post-infection prediction."""

    def __init__(self, k: int = 10, within_virus_only: bool = True) -> None:
        self.k = k
        self.within_virus_only = within_virus_only
        self.X_baseline_: np.ndarray | None = None
        self.X_post_: np.ndarray | None = None
        self.virus_id_: np.ndarray | None = None
        self._knn_per_virus: dict[str, NearestNeighbors] = {}
        self._global_knn: NearestNeighbors | None = None

    def fit(self, X_baseline, X_post, virus_id) -> KNNBaseline:
        self.X_baseline_ = np.asarray(X_baseline, dtype=np.float32)
        self.X_post_ = np.asarray(X_post, dtype=np.float32)
        self.virus_id_ = np.asarray(virus_id)

        if self.within_virus_only:
            for v in np.unique(self.virus_id_):
                mask = (self.virus_id_ == v).nonzero()[0]
                if len(mask) < self.k:
                    continue
                self._knn_per_virus[str(v)] = NearestNeighbors(
                    n_neighbors=self.k, metric="cosine"
                ).fit(self.X_baseline_[mask])
        else:
            self._global_knn = NearestNeighbors(n_neighbors=self.k, metric="cosine").fit(
                self.X_baseline_
            )
        return self

    def predict(self, X_baseline, virus_id) -> np.ndarray:
        """For each test cell, distance-weighted mean of K nearest training cells' post."""
        if self.X_baseline_ is None or self.X_post_ is None or self.virus_id_ is None:
            raise RuntimeError("Must call fit() before predict().")
        X_baseline = np.asarray(X_baseline, dtype=np.float32)
        virus_id = np.asarray(virus_id)
        n = len(X_baseline)
        n_genes = self.X_post_.shape[1]
        out = np.empty((n, n_genes), dtype=np.float32)

        for i, (xb, v) in enumerate(zip(X_baseline, virus_id, strict=False)):
            v_str = str(v)
            if self.within_virus_only:
                if v_str not in self._knn_per_virus:
                    if self._global_knn is None:
                        out[i] = self.X_post_.mean(axis=0)
                        continue
                    knn = self._global_knn
                    pool_idx = np.arange(len(self.X_post_))
                else:
                    knn = self._knn_per_virus[v_str]
                    pool_idx = (self.virus_id_ == v).nonzero()[0]
            else:
                knn = self._global_knn
                pool_idx = np.arange(len(self.X_post_))

            distances, neighbor_idx = knn.kneighbors(xb.reshape(1, -1))
            distances = distances.flatten()
            neighbor_idx = neighbor_idx.flatten()
            weights = 1.0 / (distances + 1e-8)
            weights /= weights.sum()
            real_idx = pool_idx[neighbor_idx]
            out[i] = (self.X_post_[real_idx] * weights[:, None]).sum(axis=0)
        return out
