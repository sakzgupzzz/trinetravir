"""Linear-delta baseline (Issue 24 Category 1; pre-spec 2026-05-11).

Ridge regression on baseline expression in HVG space + cell-type one-hot +
virus one-hot for within-virus training. Alpha cross-validated within Issue 14
20-config budget (held-out donor split).

Implementation deferred to Phase 5. See METHODS_CHOICES.md Issue 24.
"""

from __future__ import annotations


class LinearDeltaBaseline:
    """Ridge regression predicting response from baseline + cell-type + virus features."""

    def __init__(self, alpha: float = 1.0) -> None:
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def fit(self, X_baseline, X_response, cell_type, virus_id):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def predict(self, X_baseline, cell_type, virus_id):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")
