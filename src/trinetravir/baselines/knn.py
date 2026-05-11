"""KNN baseline (Issue 24 Category 1; pre-spec 2026-05-11).

Cosine distance on log-normalized HVG expression, k=10 with distance-weighted
predictions. Within-virus training neighbors only for cross-virus evaluation
per Issue 15 protocol. Sensitivity at k=5 and k=20 in Phase 5 supplementary.

Implementation deferred to Phase 5. See METHODS_CHOICES.md Issue 24.
"""

from __future__ import annotations


class KNNBaseline:
    """Cosine-distance KNN with distance-weighted predictions, within-virus neighbors only."""

    def __init__(self, k: int = 10) -> None:
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def fit(self, X, y, virus_id):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def predict(self, X, virus_id):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")
