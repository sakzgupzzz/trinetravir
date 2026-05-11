"""ISG-score regression baseline (Issue 24 Category 2; pre-spec 2026-05-11).

Per-donor ISG signature score (mean expression of Khatri MVS genes per
Issue 18) as single feature; linear regression to predict virus identity +
cell-type. The simplest possible "ISG explains everything" baseline.

Critique-document concern 2 (deep learning necessity): if the factorized
model fails to beat this baseline by ≥0.05 cross-study Pearson r averaged
across buckets, the paper acknowledges that the model architecture does
not add value beyond gene-set restriction.

Implementation deferred to Phase 5. See METHODS_CHOICES.md Issue 24.
"""

from __future__ import annotations


class ISGScoreRegressionBaseline:
    """Linear regression on per-donor mean Khatri MVS score predicting virus + cell-type."""

    def __init__(self) -> None:
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def fit(self, X, y, virus_id, cell_type):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def predict(self, X, virus_id, cell_type):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")
