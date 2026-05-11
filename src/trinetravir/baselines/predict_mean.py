"""Predict-mean baseline (Issue 24 Category 1; pre-spec 2026-05-11).

Strongest predict-mean baseline: per-cell-type per-virus mean across training
cells in HVG space. Establishes the trivial lower bound that the factorized
model must beat to demonstrate value beyond per-stratum averaging.

Implementation deferred to Phase 5. See METHODS_CHOICES.md Issue 24.
"""

from __future__ import annotations


class PredictMeanBaseline:
    """Predict per-cell-type per-virus mean across training cells in HVG space."""

    def __init__(self) -> None:
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def fit(self, X, y, cell_type, virus_id):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def predict(self, X, cell_type, virus_id):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")
