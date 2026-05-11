"""Sparse PCA baseline (Issue 24 Category 2; pre-spec 2026-05-11).

scikit-learn SparsePCA with α tuned via Issue 14 held-out validation.
k components ∈ {16, 32, 64} matching factorized model's shared latent
dim search (Issue 21).

Key comparison: tests whether nonlinear factorization in the v1 factorized
model adds value over linear sparse factorization.

Implementation deferred to Phase 5. See METHODS_CHOICES.md Issue 24.
"""

from __future__ import annotations


class SparsePCABaseline:
    """Sparse PCA factorization baseline; k matches factorized model latent dim search."""

    def __init__(self, n_components: int = 32, alpha: float = 1.0) -> None:
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def fit(self, X, y, virus_id):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")

    def predict(self, X, virus_id):
        raise NotImplementedError("Phase 5 implementation. See METHODS_CHOICES.md Issue 24.")
