"""Baseline implementations per METHODS_CHOICES.md Issue 24 (pre-spec 2026-05-11).

Two categories:
  Category 1 — Trivial baselines (lower bound):
    - predict_mean: per-cell-type per-virus mean across training cells.
    - linear_delta: ridge regression on baseline + cell-type + virus one-hots.
    - knn: cosine distance on log-normalized HVG, k=10, within-virus neighbors only.
  Category 2 — Simpler factorization baselines (key comparison):
    - sparse_pca: scikit-learn SparsePCA, alpha tuned via Issue 14 held-out validation.
    - nmf: scikit-learn NMF on log1p counts, init='nndsvd'.
    - isg_score_regression: simplest "ISG explains everything" baseline.

All stubs raise NotImplementedError; full implementation at Phase 5.
Implementation must pass Issue 14 hyperparameter policy (held-out donor split,
20-config budget) and Issue 15 cross-virus evaluation protocol.
"""
