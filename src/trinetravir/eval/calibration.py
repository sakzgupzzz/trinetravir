"""Phase 3 / 3.5 threshold calibration for per-bucket cross-study Pearson r.

Three complementary calibration approaches:

1. **Donor-level permutation null.** For each bucket, generate N random
   reassignments of disease/healthy labels at the donor level (preserving
   per-study marginal counts). For each permutation, recompute the per-study
   response vector on the *same* Harmony-corrected scaled-HVG embedding and
   take the mean off-diagonal cross-study Pearson r. The empirical null is
   the distribution of these N values. Report null mean / SD / 95th / 99th
   percentile and empirical p-value = fraction of null >= observed.

2. **Donor-level within-study split-half reliability ceiling.** For each
   bucket and each study, repeatedly split the study's donors into two
   stratified halves (preserving disease/healthy ratio). Compute the
   response vector for each half independently. The mean of the per-split
   Pearson r between halves is the within-study reliability. The
   cross-study ceiling = mean of within-study reliabilities across studies.
   Any cross-study correlation above the ceiling would imply better signal
   replication across cohorts than within a single cohort, which is
   biologically implausible — the ceiling is the noise floor on cross-study
   r given finite per-study sampling.

3. **Combined pass criterion.** Bucket passes the calibrated gate iff
   (a) permutation p-value < 0.01 (signal exceeds chance) AND
   (b) observed r >= 0.5 * split-half ceiling (signal captures at least
   half the within-study reliability).

Donor-level (not cell-level) shuffling/splitting is mandatory. Cell-level
permutation breaks within-donor correlation structure, producing tight,
over-confident null distributions. Cell-level split-half similarly inflates
reliability because half-1's cells from donor X are still highly correlated
with half-2's cells from donor X.

All randomness is seeded from a single configuration parameter so the
calibration is reproducible. Permutation results are cached to disk per
(bucket, seed, n_perm) so re-runs are cheap.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from trinetravir.data.harmonize import BucketResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Common helpers
# ---------------------------------------------------------------------------


def _per_study_donor_status(cell_obs: pd.DataFrame) -> dict[str, pd.Series]:
    """Per-study mapping from donor_id to its (unique) donor_disease_status.

    Asserts each donor appears under one status in the study (true for the
    4 Phase-3 studies by construction).
    """
    out: dict[str, pd.Series] = {}
    for sid, grp in cell_obs.groupby("study_id", observed=True):
        donor_status = (
            grp[["donor_id", "donor_disease_status"]]
            .drop_duplicates()
            .set_index("donor_id")["donor_disease_status"]
        )
        n_donors = grp["donor_id"].nunique()
        if len(donor_status) != n_donors:
            # A donor with mixed status across cells. Should not happen, but
            # if it does we degrade gracefully by taking the majority label.
            logger.warning("study %s has donors with non-unique status; taking mode", sid)
            donor_status = grp.groupby("donor_id")["donor_disease_status"].agg(
                lambda s: s.mode().iat[0]
            )
        out[str(sid)] = donor_status
    return out


def _response_vector(
    x_corrected: np.ndarray,
    cell_donor: np.ndarray,
    donor_to_label: dict[str, str],
) -> np.ndarray | None:
    """Mean(diseased) - mean(healthy) over the cells whose donor is labeled.

    Cells whose donor maps to neither status are dropped. Returns None if
    either pool is empty.
    """
    labels = np.array([donor_to_label.get(d) for d in cell_donor])
    mask_d = labels == "diseased"
    mask_h = labels == "healthy_control"
    if not mask_d.any() or not mask_h.any():
        return None
    return x_corrected[mask_d].mean(axis=0) - x_corrected[mask_h].mean(axis=0)


def _pairwise_mean_off_diag(rvs: dict[str, np.ndarray]) -> float | None:
    """Pearson r matrix of stacked response vectors -> mean off-diagonal r.

    Returns None if fewer than 2 valid response vectors are present.
    """
    keys = [k for k, v in rvs.items() if v is not None]
    if len(keys) < 2:
        return None
    mat = np.stack([rvs[k] for k in keys], axis=1)
    r = np.corrcoef(mat, rowvar=False)
    off = r[~np.eye(len(keys), dtype=bool)]
    return float(off.mean())


# ---------------------------------------------------------------------------
# Approach 1: donor-level permutation null
# ---------------------------------------------------------------------------


@dataclass
class PermutationNullResult:
    bucket: str
    n_perm_requested: int
    n_perm_actual: int
    observed_r: float
    null_mean: float
    null_sd: float
    null_p95: float
    null_p99: float
    p_value: float
    null_values: np.ndarray


def donor_permutation_null(
    bucket_result: BucketResult,
    *,
    n_perm: int = 1000,
    seed: int = 42,
    cache_dir: str | Path | None = None,
    label: str = "phase3",
) -> PermutationNullResult:
    """Empirical null distribution of mean off-diag cross-study r.

    Permutes donor disease/healthy labels per study (preserves marginals).
    Each permutation recomputes per-study response vectors using the same
    Harmony-corrected ``x_corrected`` embedding and computes mean off-diag
    Pearson r across studies.

    Parameters
    ----------
    bucket_result
        Output of ``harmony_per_bucket(..., keep_cells=True)``. Must have
        non-None ``x_corrected`` and ``cell_obs``.
    n_perm
        Target number of permutations. If the per-study combinatorial limit
        is lower (e.g. small donor counts), at most C(n_d, n_h) per study
        and we draw uniformly from the product space.
    seed
        Single integer seeding all per-permutation randomness.
    cache_dir
        Optional directory. When set, the null array is cached to
        ``<cache_dir>/perm_null_<label>_<bucket>_<n_perm>_<seed>.npz`` so
        re-runs of the same configuration are instant.
    label
        Identifier (``phase3`` / ``phase35``) embedded in the cache filename
        to separate calibration runs on different harmonized inputs.

    Returns
    -------
    PermutationNullResult
    """
    if bucket_result.x_corrected is None or bucket_result.cell_obs is None:
        raise ValueError("bucket_result needs keep_cells=True at harmonize time")

    x = bucket_result.x_corrected
    obs = bucket_result.cell_obs
    bucket = bucket_result.bucket
    observed_r = bucket_result.gate_r

    cache_path: Path | None = None
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir / f"perm_null_{label}_{bucket}_{n_perm}_{seed}.npz"
        if cache_path.exists():
            logger.info("loading cached permutation null %s", cache_path)
            data = np.load(cache_path)
            null = data["null"]
            return _pack_perm_result(bucket, n_perm, int(data["n_actual"]), observed_r, null)

    rng = np.random.default_rng(seed)

    per_study_donors = _per_study_donor_status(obs)
    studies = sorted(per_study_donors.keys())

    # Per-study list of (donor_ids array, original n_diseased)
    study_specs: dict[str, tuple[np.ndarray, int]] = {}
    for sid in studies:
        ds = per_study_donors[sid]
        donors = ds.index.to_numpy()
        n_d = int((ds == "diseased").sum())
        if n_d == 0 or n_d == len(donors):
            logger.warning(
                "study %s has all-one-class donors (%d diseased of %d); "
                "skipping its contribution to permutation null",
                sid,
                n_d,
                len(donors),
            )
            continue
        study_specs[sid] = (donors, n_d)

    if len(study_specs) < 2:
        raise RuntimeError(f"bucket {bucket}: fewer than 2 studies have both donor classes")

    # Pre-compute per-cell (study_id, donor_id) once
    cell_study = obs["study_id"].astype(str).to_numpy()
    cell_donor = obs["donor_id"].astype(str).to_numpy()
    study_to_mask = {sid: (cell_study == sid) for sid in study_specs}
    study_to_donors = {sid: cell_donor[m] for sid, m in study_to_mask.items()}
    study_to_x = {sid: x[m] for sid, m in study_to_mask.items()}

    null_values: list[float] = []
    n_failed = 0
    for _ in range(n_perm):
        rvs: dict[str, np.ndarray] = {}
        for sid, (donors, n_d) in study_specs.items():
            permuted_donors = rng.permutation(donors)
            diseased_set = set(permuted_donors[:n_d].tolist())
            donor_to_label = {
                d: ("diseased" if d in diseased_set else "healthy_control") for d in donors
            }
            rv = _response_vector(study_to_x[sid], study_to_donors[sid], donor_to_label)
            if rv is not None:
                rvs[sid] = rv
        mean_off = _pairwise_mean_off_diag(rvs)
        if mean_off is None:
            n_failed += 1
            continue
        null_values.append(mean_off)

    null = np.asarray(null_values, dtype=np.float64)
    n_actual = len(null)
    if n_failed:
        logger.warning(
            "bucket %s: %d/%d permutations failed (insufficient cells per pool)",
            bucket,
            n_failed,
            n_perm,
        )
    if cache_path is not None:
        np.savez(cache_path, null=null, n_actual=np.int64(n_actual), seed=np.int64(seed))
        logger.info("cached permutation null -> %s", cache_path)
    return _pack_perm_result(bucket, n_perm, n_actual, observed_r, null)


def _pack_perm_result(bucket, n_perm, n_actual, observed_r, null) -> PermutationNullResult:
    if len(null) == 0:
        return PermutationNullResult(
            bucket=bucket,
            n_perm_requested=n_perm,
            n_perm_actual=0,
            observed_r=float(observed_r),
            null_mean=float("nan"),
            null_sd=float("nan"),
            null_p95=float("nan"),
            null_p99=float("nan"),
            p_value=float("nan"),
            null_values=null,
        )
    p = float(((null >= observed_r).sum() + 1) / (len(null) + 1))
    return PermutationNullResult(
        bucket=bucket,
        n_perm_requested=n_perm,
        n_perm_actual=n_actual,
        observed_r=float(observed_r),
        null_mean=float(null.mean()),
        null_sd=float(null.std(ddof=1)) if len(null) > 1 else float("nan"),
        null_p95=float(np.percentile(null, 95)),
        null_p99=float(np.percentile(null, 99)),
        p_value=p,
        null_values=null,
    )


# ---------------------------------------------------------------------------
# Approach 2: literature anchoring (data class only — literature notes live
# in references/notes/calibration_*.md and the consolidated summary table)
# ---------------------------------------------------------------------------


@dataclass
class LiteratureAnchor:
    """Literature reference value used as a soft anchor for one bucket."""

    bucket: str
    reference_value: float | None  # comparable cross-cohort metric, if any
    source_paper: str
    metric_description: str
    comparability: str  # "direct", "indirect", "field-context-only"


# ---------------------------------------------------------------------------
# Approach 3: donor-level split-half reliability ceiling
# ---------------------------------------------------------------------------


@dataclass
class SplitHalfStudyResult:
    study: str
    n_splits_attempted: int
    n_splits_completed: int
    n_donors: int
    n_donors_per_half: tuple[int, int]
    mean_r: float
    ci_low: float
    ci_high: float


@dataclass
class SplitHalfCeilingResult:
    bucket: str
    per_study: list[SplitHalfStudyResult]
    ceiling: float  # mean across studies of mean within-study split-half r
    ceiling_sd: float


def split_half_ceiling(
    bucket_result: BucketResult,
    *,
    n_splits: int = 50,
    seed: int = 42,
) -> SplitHalfCeilingResult:
    """Per-study within-study donor-level split-half reliability.

    For each study and each of ``n_splits`` repeats:
      - Stratify-split the study's donors into two halves preserving the
        disease/healthy ratio.
      - Compute response vector for each half in ``x_corrected`` space.
      - Pearson r between the two halves.

    Per-study reliability = mean over splits. Bucket ceiling = mean across
    studies of per-study reliability.
    """
    if bucket_result.x_corrected is None or bucket_result.cell_obs is None:
        raise ValueError("bucket_result needs keep_cells=True at harmonize time")

    x = bucket_result.x_corrected
    obs = bucket_result.cell_obs
    bucket = bucket_result.bucket
    per_study_donors = _per_study_donor_status(obs)
    rng = np.random.default_rng(seed)

    cell_study = obs["study_id"].astype(str).to_numpy()
    cell_donor = obs["donor_id"].astype(str).to_numpy()

    per_study_results: list[SplitHalfStudyResult] = []
    for sid, ds in per_study_donors.items():
        m_study = cell_study == sid
        study_donors_x = x[m_study]
        study_donors_arr = cell_donor[m_study]
        diseased = ds[ds == "diseased"].index.to_numpy()
        healthy = ds[ds == "healthy_control"].index.to_numpy()
        n_d = len(diseased)
        n_h = len(healthy)
        # Need at least 2 donors per class to split into two halves each
        # containing at least one of each class.
        if n_d < 2 or n_h < 2:
            logger.warning(
                "study %s bucket %s: insufficient donors per class (%d diseased / %d healthy); skipping split-half",
                sid,
                bucket,
                n_d,
                n_h,
            )
            continue
        n_d1 = n_d // 2
        n_d2 = n_d - n_d1
        n_h1 = n_h // 2
        n_h2 = n_h - n_h1

        rs: list[float] = []
        for _ in range(n_splits):
            rng.shuffle(diseased)
            rng.shuffle(healthy)
            half1_d = diseased[:n_d1]
            half2_d = diseased[n_d1:]
            half1_h = healthy[:n_h1]
            half2_h = healthy[n_h1:]
            donor_to_label1 = {
                **dict.fromkeys(half1_d, "diseased"),
                **dict.fromkeys(half1_h, "healthy_control"),
            }
            donor_to_label2 = {
                **dict.fromkeys(half2_d, "diseased"),
                **dict.fromkeys(half2_h, "healthy_control"),
            }
            rv1 = _response_vector(study_donors_x, study_donors_arr, donor_to_label1)
            rv2 = _response_vector(study_donors_x, study_donors_arr, donor_to_label2)
            if rv1 is None or rv2 is None:
                continue
            r = float(np.corrcoef(rv1, rv2)[0, 1])
            if not np.isnan(r):
                rs.append(r)

        if not rs:
            logger.warning("study %s bucket %s: no successful splits", sid, bucket)
            continue
        rs_arr = np.asarray(rs)
        ci_low, ci_high = np.percentile(rs_arr, [2.5, 97.5])
        per_study_results.append(
            SplitHalfStudyResult(
                study=sid,
                n_splits_attempted=n_splits,
                n_splits_completed=len(rs_arr),
                n_donors=n_d + n_h,
                n_donors_per_half=(n_d1 + n_h1, n_d2 + n_h2),
                mean_r=float(rs_arr.mean()),
                ci_low=float(ci_low),
                ci_high=float(ci_high),
            )
        )

    if not per_study_results:
        return SplitHalfCeilingResult(
            bucket=bucket, per_study=[], ceiling=float("nan"), ceiling_sd=float("nan")
        )

    per_study_means = np.array([s.mean_r for s in per_study_results])
    return SplitHalfCeilingResult(
        bucket=bucket,
        per_study=per_study_results,
        ceiling=float(per_study_means.mean()),
        ceiling_sd=float(per_study_means.std(ddof=1)) if len(per_study_means) > 1 else float("nan"),
    )


# ---------------------------------------------------------------------------
# Combined calibration row + table
# ---------------------------------------------------------------------------


@dataclass
class CalibrationRow:
    bucket: str
    observed_r: float
    threshold_heuristic: float
    permutation: PermutationNullResult
    split_half: SplitHalfCeilingResult
    literature: LiteratureAnchor | None
    p_value_threshold: float
    ceiling_fraction_threshold: float

    @property
    def passes_permutation(self) -> bool:
        return self.permutation.p_value < self.p_value_threshold

    @property
    def passes_split_half(self) -> bool:
        if np.isnan(self.split_half.ceiling):
            return False
        return self.observed_r >= self.ceiling_fraction_threshold * self.split_half.ceiling

    @property
    def passes_heuristic(self) -> bool:
        return self.observed_r >= self.threshold_heuristic

    @property
    def passes_calibrated(self) -> bool:
        return self.passes_permutation and self.passes_split_half

    def to_dict(self) -> dict:
        return {
            "bucket": self.bucket,
            "observed_r": round(self.observed_r, 3),
            "heuristic_threshold": self.threshold_heuristic,
            "heuristic_pass": self.passes_heuristic,
            "perm_null_mean": round(self.permutation.null_mean, 3),
            "perm_null_p95": round(self.permutation.null_p95, 3),
            "perm_null_p99": round(self.permutation.null_p99, 3),
            "perm_p_value": round(self.permutation.p_value, 4),
            "perm_n_actual": self.permutation.n_perm_actual,
            "perm_pass": self.passes_permutation,
            "split_half_ceiling": round(self.split_half.ceiling, 3)
            if not np.isnan(self.split_half.ceiling)
            else float("nan"),
            "ceiling_fraction": round(self.observed_r / self.split_half.ceiling, 3)
            if self.split_half.ceiling and not np.isnan(self.split_half.ceiling)
            else float("nan"),
            "split_half_pass": self.passes_split_half,
            "literature_value": (self.literature.reference_value if self.literature else None),
            "literature_source": (self.literature.source_paper if self.literature else None),
            "literature_comparability": (
                self.literature.comparability if self.literature else None
            ),
            "calibrated_pass": self.passes_calibrated,
        }


def build_calibration_row(
    bucket_result: BucketResult,
    *,
    n_perm: int = 1000,
    n_splits: int = 50,
    seed: int = 42,
    cache_dir: str | Path | None = None,
    label: str = "phase3",
    literature: LiteratureAnchor | None = None,
    threshold_heuristic: float | None = None,
    p_value_threshold: float = 0.01,
    ceiling_fraction_threshold: float = 0.5,
) -> CalibrationRow:
    from trinetravir.data.harmonize import GATE_THRESHOLDS

    th = (
        GATE_THRESHOLDS[bucket_result.bucket]
        if threshold_heuristic is None
        else threshold_heuristic
    )
    perm = donor_permutation_null(
        bucket_result,
        n_perm=n_perm,
        seed=seed,
        cache_dir=cache_dir,
        label=label,
    )
    sh = split_half_ceiling(bucket_result, n_splits=n_splits, seed=seed)
    return CalibrationRow(
        bucket=bucket_result.bucket,
        observed_r=bucket_result.gate_r,
        threshold_heuristic=th,
        permutation=perm,
        split_half=sh,
        literature=literature,
        p_value_threshold=p_value_threshold,
        ceiling_fraction_threshold=ceiling_fraction_threshold,
    )


def build_calibration_table(rows: Iterable[CalibrationRow]) -> pd.DataFrame:
    """Tabulate per-bucket calibration rows into a single DataFrame."""
    return pd.DataFrame([r.to_dict() for r in rows]).set_index("bucket")


# ---------------------------------------------------------------------------
# Metric-fn aware permutation + split-half (Session 3 extension)
#
# The original donor_permutation_null + split_half_ceiling are Pearson-only
# (via _pairwise_mean_off_diag). The functions below accept a generic
# metric_fn(dict[study_id -> response_vector]) -> float and run the same
# donor-level shuffle / stratified-split protocol against any metric.
# They take raw (x_corrected, cell_obs, bucket_name) rather than a
# BucketResult so they work uniformly for Phase 3 per-bucket files and the
# Phase 3.5 consolidated uns['X_corrected_<bucket>'] entries.
# ---------------------------------------------------------------------------


MetricFn = Callable[[dict[str, np.ndarray]], float]


def permutation_null_with_metric(
    x_corrected: np.ndarray,
    cell_obs: pd.DataFrame,
    bucket: str,
    *,
    metric_fn: MetricFn,
    n_permutations: int = 1000,
    seed: int = 42,
    cache_dir: str | Path | None = None,
    cache_label: str = "metric",
) -> dict:
    """Generic permutation null. Returns {'null': ndarray, 'n_actual': int}.

    Same donor-level shuffle protocol as donor_permutation_null but the metric
    is pluggable.  Caches per (cache_label, bucket, n_permutations, seed).
    """
    cache_path: Path | None = None
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir / (f"perm_null_{cache_label}_{bucket}_{n_permutations}_{seed}.npz")
        if cache_path.exists():
            data = np.load(cache_path)
            return {"null": data["null"], "n_actual": int(data["n_actual"])}

    rng = np.random.default_rng(seed)
    per_study_donors = _per_study_donor_status(cell_obs)
    studies = sorted(per_study_donors.keys())
    study_specs: dict[str, tuple[np.ndarray, int]] = {}
    for sid in studies:
        ds = per_study_donors[sid]
        donors = ds.index.to_numpy()
        n_d = int((ds == "diseased").sum())
        if n_d == 0 or n_d == len(donors):
            continue
        study_specs[sid] = (donors, n_d)
    if len(study_specs) < 2:
        raise RuntimeError(f"bucket {bucket}: fewer than 2 valid studies")

    cell_study = cell_obs["study_id"].astype(str).to_numpy()
    cell_donor = cell_obs["donor_id"].astype(str).to_numpy()
    study_to_mask = {sid: (cell_study == sid) for sid in study_specs}
    study_to_donors = {sid: cell_donor[m] for sid, m in study_to_mask.items()}
    study_to_x = {sid: x_corrected[m] for sid, m in study_to_mask.items()}

    null_values: list[float] = []
    for _ in range(n_permutations):
        rvs: dict[str, np.ndarray] = {}
        for sid, (donors, n_d) in study_specs.items():
            permuted = rng.permutation(donors)
            diseased_set = set(permuted[:n_d].tolist())
            donor_to_label = {
                d: ("diseased" if d in diseased_set else "healthy_control") for d in donors
            }
            rv = _response_vector(study_to_x[sid], study_to_donors[sid], donor_to_label)
            if rv is not None:
                rvs[sid] = rv
        if len(rvs) < 2:
            continue
        v = metric_fn(rvs)
        if not np.isnan(v):
            null_values.append(v)

    null = np.asarray(null_values, dtype=np.float64)
    if cache_path is not None:
        np.savez(cache_path, null=null, n_actual=np.int64(len(null)), seed=np.int64(seed))
    return {"null": null, "n_actual": len(null)}


def split_half_with_metric(
    x_corrected: np.ndarray,
    cell_obs: pd.DataFrame,
    bucket: str,
    *,
    metric_fn: MetricFn,
    n_splits: int = 50,
    seed: int = 42,
) -> dict:
    """Per-study donor-stratified split-half reliability under metric_fn.

    Returns {'per_study_rs': list of arrays per study, 'ceiling': float,
             'split_half_distribution': ndarray of all split values}.
    """
    rng = np.random.default_rng(seed)
    per_study_donors = _per_study_donor_status(cell_obs)
    cell_study = cell_obs["study_id"].astype(str).to_numpy()
    cell_donor = cell_obs["donor_id"].astype(str).to_numpy()

    per_study_rs: dict[str, np.ndarray] = {}
    all_rs: list[float] = []
    for sid, ds in per_study_donors.items():
        m_study = cell_study == sid
        study_x = x_corrected[m_study]
        study_donors_arr = cell_donor[m_study]
        diseased = ds[ds == "diseased"].index.to_numpy()
        healthy = ds[ds == "healthy_control"].index.to_numpy()
        if len(diseased) < 2 or len(healthy) < 2:
            continue
        n_d1 = len(diseased) // 2
        n_h1 = len(healthy) // 2

        rs: list[float] = []
        for _ in range(n_splits):
            rng.shuffle(diseased)
            rng.shuffle(healthy)
            half1_d, half2_d = diseased[:n_d1], diseased[n_d1:]
            half1_h, half2_h = healthy[:n_h1], healthy[n_h1:]
            l1 = {**dict.fromkeys(half1_d, "diseased"), **dict.fromkeys(half1_h, "healthy_control")}
            l2 = {**dict.fromkeys(half2_d, "diseased"), **dict.fromkeys(half2_h, "healthy_control")}
            rv1 = _response_vector(study_x, study_donors_arr, l1)
            rv2 = _response_vector(study_x, study_donors_arr, l2)
            if rv1 is None or rv2 is None:
                continue
            # Apply metric on a 2-element "study set" {half1, half2}.
            v = metric_fn({"half1": rv1, "half2": rv2})
            if not np.isnan(v):
                rs.append(v)
        if rs:
            per_study_rs[sid] = np.asarray(rs)
            all_rs.extend(rs)
    if not per_study_rs:
        return {
            "per_study_rs": {},
            "ceiling": float("nan"),
            "split_half_distribution": np.asarray([]),
        }
    ceiling = float(np.mean([rs.mean() for rs in per_study_rs.values()]))
    return {
        "per_study_rs": per_study_rs,
        "ceiling": ceiling,
        "split_half_distribution": np.asarray(all_rs, dtype=np.float64),
    }


def bootstrap_ci_overlap(
    observed: float,
    split_half_distribution: np.ndarray,
    *,
    alpha: float = 0.05,
) -> dict:
    """Test whether ``observed`` is AT-OR-ABOVE lower bound of split-half CI.

    Session 5 Part A1 correction: the prior implementation treated observed
    *within* the (1-alpha) split-half CI as PASS and *above* the upper CI
    as FAIL. That was wrong — observed r > split-half upper CI represents
    cross-study coherence that EXCEEDS within-study reliability, which is
    the strongest possible signal (interesting, not failure).

    The corrected criterion: observed ≥ lower bound of split-half (1-alpha)
    CI. Above the upper bound is still PASS. Returns the CI bounds plus the
    ``at_or_above_ci_low`` flag (replaces the old ``in_ci`` flag).

    Returns {'at_or_above_ci_low': bool, 'in_ci': bool (legacy, retained
    for audit trail in v1 tables), 'ci_low': float, 'ci_high': float}.
    """
    if len(split_half_distribution) == 0:
        return {
            "at_or_above_ci_low": False,
            "in_ci": False,
            "ci_low": float("nan"),
            "ci_high": float("nan"),
        }
    lo = float(np.percentile(split_half_distribution, 100 * alpha / 2))
    hi = float(np.percentile(split_half_distribution, 100 * (1 - alpha / 2)))
    return {
        "at_or_above_ci_low": bool(observed >= lo),
        "in_ci": bool(lo <= observed <= hi),
        "ci_low": lo,
        "ci_high": hi,
    }


def calibrated_gate_verdict(
    observed: float,
    null_dist: np.ndarray,
    split_half_dist: np.ndarray,
    *,
    percentile: float = 99,
    alpha: float = 0.05,
    use_corrected_ci: bool = True,
) -> dict:
    """Two-criterion calibrated gate verdict.

    Pass iff:
      (1) observed exceeds (percentile)th percentile of null distribution
          (equivalently, permutation p < (1-percentile/100)).
      (2) (corrected) observed ≥ lower bound of (1-alpha) split-half CI.
          (legacy) observed within (1-alpha) split-half CI.

    The corrected criterion (Session 5 Part A1) is the default. The legacy
    "in CI" criterion remains accessible via use_corrected_ci=False for
    audit-trail v1 verdict reproduction.

    Returns dict with keys: ``pass``, ``null_threshold``, ``p_value``,
    ``at_or_above_ci_low``, ``in_split_half_ci`` (legacy), ``ci_low``,
    ``ci_high``.
    """
    if len(null_dist) == 0:
        return {
            "pass": False,
            "null_threshold": float("nan"),
            "p_value": float("nan"),
            "at_or_above_ci_low": False,
            "in_split_half_ci": False,
            "ci_low": float("nan"),
            "ci_high": float("nan"),
        }
    thr = float(np.percentile(null_dist, percentile))
    p = float(((null_dist >= observed).sum() + 1) / (len(null_dist) + 1))
    ci = bootstrap_ci_overlap(observed, split_half_dist, alpha=alpha)
    ci_pass = ci["at_or_above_ci_low"] if use_corrected_ci else ci["in_ci"]
    return {
        "pass": bool(observed >= thr and ci_pass),
        "null_threshold": thr,
        "p_value": p,
        "at_or_above_ci_low": ci["at_or_above_ci_low"],
        "in_split_half_ci": ci["in_ci"],
        "ci_low": ci["ci_low"],
        "ci_high": ci["ci_high"],
    }


def bootstrap_observed_r(
    x_corrected: np.ndarray,
    cell_obs: pd.DataFrame,
    *,
    metric_fn: Callable[[dict[str, np.ndarray]], float],
    n_bootstrap: int = 1000,
    seed: int = 42,
) -> dict:
    """Donor-level bootstrap CI on observed cross-study r.

    Session 5 Part A2. Each iteration resamples donors WITH REPLACEMENT
    within each study (preserving per-study donor count) and computes
    bootstrap mean(diseased donors) - mean(healthy donors) by averaging
    pre-computed per-donor cell means. Per-iteration cost is O(n_donors *
    n_genes), independent of n_cells. Returns dict with observed_ci_low,
    observed_ci_high, n_bootstrap_completed, bootstrap_values.
    """
    rng = np.random.default_rng(seed)
    per_study = _per_study_donor_status(cell_obs)
    cell_study = cell_obs["study_id"].astype(str).to_numpy()
    cell_donor = cell_obs["donor_id"].astype(str).to_numpy()

    studies = sorted(per_study.keys())
    # Pre-compute per-donor mean for each (study, donor). Bootstrap then
    # resamples donor-level means rather than re-aggregating cells.
    per_study_donor_means: dict[str, dict[str, np.ndarray]] = {}
    per_study_donor_status: dict[str, dict[str, str]] = {}
    for sid in studies:
        m = cell_study == sid
        sub_x = x_corrected[m]
        sub_donor = cell_donor[m]
        ds = per_study[sid].to_dict()
        donor_means: dict[str, np.ndarray] = {}
        for d in set(sub_donor.tolist()):
            mask_d = sub_donor == d
            donor_means[d] = sub_x[mask_d].mean(axis=0)
        per_study_donor_means[sid] = donor_means
        per_study_donor_status[sid] = ds

    boot_values: list[float] = []
    for _ in range(n_bootstrap):
        rvs: dict[str, np.ndarray] = {}
        for sid in studies:
            donor_means = per_study_donor_means[sid]
            ds = per_study_donor_status[sid]
            donors = np.array(list(donor_means.keys()))
            if len(donors) < 2:
                continue
            sampled = rng.choice(donors, size=len(donors), replace=True)
            d_means, h_means = [], []
            for d in sampled:
                if ds.get(d) == "diseased":
                    d_means.append(donor_means[d])
                elif ds.get(d) == "healthy_control":
                    h_means.append(donor_means[d])
            if not d_means or not h_means:
                continue
            rv = np.mean(np.stack(d_means, axis=0), axis=0) - np.mean(
                np.stack(h_means, axis=0), axis=0
            )
            rvs[sid] = rv
        if len(rvs) < 2:
            continue
        v = metric_fn(rvs)
        if not np.isnan(v):
            boot_values.append(v)
    arr = np.asarray(boot_values, dtype=np.float64)
    if len(arr) == 0:
        return {
            "observed_ci_low": float("nan"),
            "observed_ci_high": float("nan"),
            "n_bootstrap_completed": 0,
            "bootstrap_values": arr,
        }
    return {
        "observed_ci_low": float(np.percentile(arr, 2.5)),
        "observed_ci_high": float(np.percentile(arr, 97.5)),
        "n_bootstrap_completed": len(arr),
        "bootstrap_values": arr,
    }


def fdr_bh(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction.

    Session 5 Part A3. Standard step-up BH on the input p-value array.
    NaN inputs propagate to NaN outputs. Returns the FDR-corrected
    (adjusted) p-values.
    """
    p = np.asarray(p_values, dtype=np.float64)
    out = np.full_like(p, np.nan)
    mask = ~np.isnan(p)
    pv = p[mask]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(pv)
    ranks = np.arange(1, n + 1, dtype=np.float64)
    sorted_p = pv[order]
    bh = sorted_p * n / ranks
    # Step-up: enforce monotonicity (running min from right).
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.clip(bh, 0.0, 1.0)
    inv = np.empty_like(order)
    inv[order] = np.arange(n)
    out[mask] = bh[inv]
    return out


def summary_stats_off_diag(rvs: dict[str, np.ndarray]) -> dict:
    """Mean / median / min off-diagonal Pearson r (Issue 11)."""
    keys = sorted(rvs.keys())
    if len(keys) < 2:
        return {"mean": float("nan"), "median": float("nan"), "min": float("nan")}
    mat = np.stack([rvs[k] for k in keys], axis=1)
    r = np.corrcoef(mat, rowvar=False)
    off = r[~np.eye(len(keys), dtype=bool)]
    return {"mean": float(off.mean()), "median": float(np.median(off)), "min": float(off.min())}


def to_paper_markdown(rows: Iterable[CalibrationRow]) -> str:
    """Render a paper-ready markdown table from the calibration rows."""
    header = (
        "| Bucket | Observed r | Heuristic threshold | Perm null p99 | Perm p-value | "
        "Split-half ceiling | Obs / ceiling | Literature ref | Calibrated verdict |"
    )
    sep = "|" + "|".join(["---"] * 9) + "|"
    lines = [header, sep]
    for r in rows:
        sh_ceiling = f"{r.split_half.ceiling:.3f}" if not np.isnan(r.split_half.ceiling) else "n/a"
        fraction = (
            f"{r.observed_r / r.split_half.ceiling:.2f}"
            if r.split_half.ceiling and not np.isnan(r.split_half.ceiling)
            else "n/a"
        )
        lit = (
            f"{r.literature.reference_value:.2f} ({r.literature.source_paper}; {r.literature.comparability})"
            if r.literature and r.literature.reference_value is not None
            else (r.literature.source_paper + "; " + r.literature.comparability)
            if r.literature
            else "—"
        )
        verdict = "PASS" if r.passes_calibrated else "FAIL"
        lines.append(
            f"| {r.bucket} | {r.observed_r:.3f} | {r.threshold_heuristic:.2f} | "
            f"{r.permutation.null_p99:.3f} | {r.permutation.p_value:.4f} | "
            f"{sh_ceiling} | {fraction} | {lit} | {verdict} |"
        )
    return "\n".join(lines)
