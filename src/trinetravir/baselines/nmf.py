"""NMF baseline (Issue 24 Category 2; pre-spec 2026-05-11).

scikit-learn NMF on log-normalized HVG counts (max(0, log1p(x))).
k components matching factorized model's shared latent dim search
(Issue 21, k ∈ {16, 32, 64}). init='nndsvd' for reproducibility.

Key comparison: tests whether nonlinear factorization in the v1 factorized
model adds value over non-negative linear factorization.

Implementation deferred to Phase 5. See METHODS_CHOICES.md Issue 24.
"""

from __future__ import annotations


class NMFBaseline:
    """Non-negative matrix factorization baseline; k matches factorized model latent dim search."""

    def __init__(self, n_components: int = 32) -> None:
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def fit(self, X, y, virus_id):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def predict(self, X, virus_id):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")
