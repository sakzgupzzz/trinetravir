"""Linear-delta baseline (Issue 24 Category 1; pre-spec 2026-05-11; Phase 5 impl 2026-05-13).

Per PLAN.md Phase 5 spec: "compute (post − pre) shift vector in PCA space from
training; add to test cell baseline."

Implementation: ridge regression on baseline expression in HVG space with
cell-type one-hot + virus one-hot per Issue 24 spec. Predicts (post - pre)
response vector. Final prediction = baseline + predicted_response.

Alpha tuned via Issue 14 held-out validation by external Phase 5 harness; this
class exposes single-alpha fit API.
"""

from __future__ import annotations

import numpy as np
from sklearn.linear_model import Ridge


class LinearDeltaBaseline:
    """Ridge regression predicting (post - pre) shift from baseline + one-hots."""

    def __init__(self, alpha: float = 1.0) -> None:
        self.alpha = alpha
        self.model_: Ridge | None = None
        self.cell_types_: list[str] = []
        self.viruses_: list[str] = []

    def _onehot(self, cell_type, virus_id) -> np.ndarray:
        ct = np.asarray(cell_type)
        vi = np.asarray(virus_id)
        ct_oh = np.zeros((len(ct), len(self.cell_types_)), dtype=np.float32)
        vi_oh = np.zeros((len(vi), len(self.viruses_)), dtype=np.float32)
        ct_lookup = {c: i for i, c in enumerate(self.cell_types_)}
        vi_lookup = {v: i for i, v in enumerate(self.viruses_)}
        for i, c in enumerate(ct):
            if str(c) in ct_lookup:
                ct_oh[i, ct_lookup[str(c)]] = 1.0
        for i, v in enumerate(vi):
            if str(v) in vi_lookup:
                vi_oh[i, vi_lookup[str(v)]] = 1.0
        return np.concatenate([ct_oh, vi_oh], axis=1)

    def fit(self, X_baseline, X_response, cell_type, virus_id) -> LinearDeltaBaseline:
        """X_baseline + X_response shapes (n_cells, n_genes); response = post - pre."""
        X_baseline = np.asarray(X_baseline, dtype=np.float32)
        X_response = np.asarray(X_response, dtype=np.float32)
        self.cell_types_ = sorted({str(c) for c in cell_type})
        self.viruses_ = sorted({str(v) for v in virus_id})
        onehots = self._onehot(cell_type, virus_id)
        X = np.concatenate([X_baseline, onehots], axis=1)
        self.model_ = Ridge(alpha=self.alpha, solver="auto")
        self.model_.fit(X, X_response)
        return self

    def predict(self, X_baseline, cell_type, virus_id) -> np.ndarray:
        """Return predicted post-infection expression = baseline + predicted_response."""
        if self.model_ is None:
            raise RuntimeError("Must call fit() before predict().")
        X_baseline = np.asarray(X_baseline, dtype=np.float32)
        onehots = self._onehot(cell_type, virus_id)
        X = np.concatenate([X_baseline, onehots], axis=1)
        delta = self.model_.predict(X)
        return X_baseline + delta
