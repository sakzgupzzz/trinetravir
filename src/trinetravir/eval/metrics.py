"""Cross-study response-vector coherence metrics for the Phase 3 gate.

Per METHODS_CHOICES Issue 3, the headline gate metric is mean off-diagonal
Pearson r across per-study response vectors. This module exposes that
metric plus alternative metrics for sensitivity analysis: Spearman r
(rank-based, robust to outliers), top-100 differentially-expressed gene
Jaccard overlap (set-based, interpretable), and an MMD-RBF placeholder
that requires per-cell distributions (deferred to v1.5; not implemented).

Each summary function takes a ``response_vectors`` dict (study_id ->
ndarray of shape (n_genes,)) and returns a single scalar — the mean
off-diagonal pairwise value across studies. This shape is intentionally
the same as ``calibration._pairwise_mean_off_diag`` so the calibration
permutation-null and split-half infrastructure can be re-used per metric
when x_corrected is cached at harmonize time.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

MetricFn = Callable[[dict[str, np.ndarray]], float]


def pearson_off_diag(response_vectors: dict[str, np.ndarray]) -> float:
    """Mean off-diagonal Pearson r across the per-study response vectors."""
    keys = sorted(response_vectors.keys())
    if len(keys) < 2:
        return float("nan")
    mat = np.stack([response_vectors[k] for k in keys], axis=1)
    r = np.corrcoef(mat, rowvar=False)
    off = r[~np.eye(len(keys), dtype=bool)]
    return float(off.mean())


def spearman_off_diag(response_vectors: dict[str, np.ndarray]) -> float:
    """Mean off-diagonal Spearman r across per-study response vectors.

    Rank-based; insensitive to outlier genes and to monotone scale changes.
    """
    keys = sorted(response_vectors.keys())
    if len(keys) < 2:
        return float("nan")
    mat = np.stack([response_vectors[k] for k in keys], axis=1)
    # scipy.spearmanr on 2 columns returns a scalar; on 3+ returns matrix.
    rho_raw, _ = spearmanr(mat, axis=0)
    if np.isscalar(rho_raw) or (hasattr(rho_raw, "ndim") and rho_raw.ndim == 0):
        return float(rho_raw)
    rho_mat = np.asarray(rho_raw)
    off = rho_mat[~np.eye(len(keys), dtype=bool)]
    return float(off.mean())


def de_jaccard_off_diag(
    response_vectors: dict[str, np.ndarray],
    *,
    top_k: int = 100,
    direction: str = "abs",
) -> float:
    """Mean off-diagonal Jaccard overlap of top-K DE genes per study.

    Each study contributes a set of gene indices ranked by the magnitude
    of its response vector. Off-diagonal pairs compute Jaccard:
    |A intersect B| / |A union B|.

    Parameters
    ----------
    response_vectors
        ``study_id -> (n_genes,) ndarray``.
    top_k
        Number of top-ranked genes per study to include in each set.
    direction
        ``"abs"`` (default) ranks by ``|x|`` (captures both up and down
        regulation); ``"up"`` ranks by ``x`` ascending-then-reverse (only
        induced genes); ``"down"`` ranks by ``x`` ascending (only
        suppressed genes).
    """
    keys = sorted(response_vectors.keys())
    if len(keys) < 2:
        return float("nan")
    n_genes = len(next(iter(response_vectors.values())))
    if top_k > n_genes:
        top_k = n_genes
    sets: dict[str, set[int]] = {}
    for k in keys:
        v = response_vectors[k]
        if direction == "abs":
            ranked = np.argsort(-np.abs(v))[:top_k]
        elif direction == "up":
            ranked = np.argsort(-v)[:top_k]
        elif direction == "down":
            ranked = np.argsort(v)[:top_k]
        else:
            raise ValueError(f"unknown direction {direction!r}")
        sets[k] = set(ranked.tolist())
    pairs: list[float] = []
    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            inter = len(sets[a] & sets[b])
            union = len(sets[a] | sets[b])
            pairs.append(inter / union if union else 0.0)
    return float(np.mean(pairs)) if pairs else float("nan")


def mmd_rbf_off_diag(
    embeddings_per_study: dict[str, np.ndarray],
    *,
    gamma: float | None = None,
) -> float:
    """Mean off-diagonal MMD with RBF kernel across per-study cell sets.

    Unlike the response-vector metrics, MMD operates on per-cell distributions:
    each study contributes an (n_cells, n_dims) matrix and pairwise MMD is
    computed via the empirical biased estimator:

        MMD^2 = E[k(x,x')] + E[k(y,y')] - 2*E[k(x,y)]

    with k(u,v) = exp(-gamma * ||u-v||^2). Lower MMD = more similar
    distributions. We negate so the convention matches the other metrics
    (higher = better cross-study coherence): score = -MMD^2.

    Parameters
    ----------
    embeddings_per_study
        Mapping study_id -> (n_cells, n_dims) Harmony-corrected embedding
        (typically obsm['X_harmony'] sliced to the bucket).
    gamma
        RBF bandwidth parameter. If None, set via median heuristic on the
        pooled embedding: gamma = 1 / (2 * sigma^2) where sigma is the
        median pairwise Euclidean distance in a random subsample.

    Returns
    -------
    score : float
        Negative mean off-diagonal pairwise MMD^2. Higher = better.
    """
    keys = sorted(embeddings_per_study.keys())
    if len(keys) < 2:
        return float("nan")

    # Median heuristic gamma: subsample at most 1000 cells from pooled and
    # compute pairwise distances. Cheap + reproducible by construction since
    # the subsample is the first 1000 by row order.
    if gamma is None:
        pooled = np.concatenate(
            [
                embeddings_per_study[k][: max(1, min(len(embeddings_per_study[k]), 250))]
                for k in keys
            ],
            axis=0,
        )
        # Distance via cdist; clip to a manageable size to avoid n^2 blow up.
        n = min(len(pooled), 1000)
        sub = pooled[:n]
        # Pairwise distances; flat upper triangle to estimate median.
        sq = ((sub[:, None, :] - sub[None, :, :]) ** 2).sum(axis=-1)
        sigma_sq = float(np.median(sq[np.triu_indices(n, k=1)]))
        if sigma_sq <= 0:
            sigma_sq = 1.0
        gamma = 1.0 / (2.0 * sigma_sq)

    def mmd2(x: np.ndarray, y: np.ndarray) -> float:
        # Biased empirical MMD^2.
        n_x, n_y = len(x), len(y)
        if n_x == 0 or n_y == 0:
            return float("nan")
        # Subsample to cap compute: at most 500 cells per study.
        rng = np.random.default_rng(0)
        if n_x > 500:
            x = x[rng.choice(n_x, 500, replace=False)]
        if n_y > 500:
            y = y[rng.choice(n_y, 500, replace=False)]
        dxx = ((x[:, None, :] - x[None, :, :]) ** 2).sum(axis=-1)
        dyy = ((y[:, None, :] - y[None, :, :]) ** 2).sum(axis=-1)
        dxy = ((x[:, None, :] - y[None, :, :]) ** 2).sum(axis=-1)
        kxx = np.exp(-gamma * dxx).mean()
        kyy = np.exp(-gamma * dyy).mean()
        kxy = np.exp(-gamma * dxy).mean()
        return float(kxx + kyy - 2.0 * kxy)

    pairs: list[float] = []
    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            m = mmd2(embeddings_per_study[a], embeddings_per_study[b])
            if not np.isnan(m):
                pairs.append(m)
    if not pairs:
        return float("nan")
    return float(-np.mean(pairs))


# Registry of metric_fn used by the sensitivity tooling.
METRICS: dict[str, MetricFn] = {
    "pearson": pearson_off_diag,
    "spearman": spearman_off_diag,
    "de_jaccard_top100": de_jaccard_off_diag,
}

# MMD operates on per-cell embeddings, not on response vectors, so it lives
# outside the response-vector METRICS registry. Calibration runners pass it
# as a separate code path.


def load_phase3_response_vectors(
    parquet_dir, buckets: list[str]
) -> dict[str, dict[str, np.ndarray]]:
    """Load cached Phase 3 response vectors per bucket from parquet.

    Files are produced by ``notebooks/04_phase3_harmonization.ipynb`` and
    named ``phase3_response_vectors_<bucket>.parquet`` (index = HVG genes,
    columns = study_ids). Returns ``{bucket: {study_id: ndarray}}``.
    """
    from pathlib import Path

    parquet_dir = Path(parquet_dir)
    out: dict[str, dict[str, np.ndarray]] = {}
    for bucket in buckets:
        path = parquet_dir / f"phase3_response_vectors_{bucket}.parquet"
        if not path.exists():
            continue
        df = pd.read_parquet(path)
        out[bucket] = {sid: df[sid].to_numpy() for sid in df.columns}
    return out


def metric_sensitivity_table(
    response_vectors: dict[str, dict[str, np.ndarray]],
    *,
    metrics: dict[str, MetricFn] | None = None,
    heuristic_thresholds: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Build a long-format sensitivity table: bucket x metric.

    Columns: bucket, metric, observed, heuristic_threshold, exceeds_threshold,
             passes_pearson_threshold_scaled (for non-Pearson metrics, applies
             Pearson threshold * 0.5 as a soft floor since metric units differ).

    The Pearson row per bucket reuses the Phase 3 calibrated thresholds in
    ``GATE_THRESHOLDS``. Other metrics are reported alongside; the verdict
    consistency check is "does the bucket pass under at least N of the
    metrics" (computed downstream).
    """
    from trinetravir.data.harmonize import GATE_THRESHOLDS

    if metrics is None:
        metrics = METRICS
    if heuristic_thresholds is None:
        heuristic_thresholds = GATE_THRESHOLDS
    rows = []
    for bucket, rvs in response_vectors.items():
        for metric_name, fn in metrics.items():
            obs = fn(rvs)
            if metric_name == "pearson":
                th = heuristic_thresholds.get(bucket, float("nan"))
                exceeds = obs >= th
            else:
                # For non-Pearson metrics, use a soft scaling: the metric
                # passes if observed >= 0.5 * Pearson threshold for the same
                # bucket. This is an interim rule documented in
                # METHODS_CHOICES Issue 3 sensitivity resolution; full
                # per-metric thresholds require permutation + split-half
                # calibration on x_corrected (v1.5).
                th = 0.5 * heuristic_thresholds.get(bucket, float("nan"))
                exceeds = obs >= th
            rows.append(
                {
                    "bucket": bucket,
                    "metric": metric_name,
                    "observed": round(float(obs), 4) if not np.isnan(obs) else float("nan"),
                    "heuristic_threshold": round(float(th), 3)
                    if not np.isnan(th)
                    else float("nan"),
                    "exceeds_threshold": bool(exceeds) if not np.isnan(obs) else False,
                }
            )
    return pd.DataFrame(rows)
