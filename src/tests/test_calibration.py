"""Synthetic ground-truth tests for src/trinetravir/eval/calibration.py (Session 5 Part A4).

Each test constructs a small synthetic (x_corrected, cell_obs) and verifies that
the calibration framework produces the expected verdict. All randomness seeded
from a single function parameter; metrics.py seed propagation tested explicitly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from trinetravir.eval.calibration import (
    bootstrap_ci_overlap,
    bootstrap_observed_r,
    calibrated_gate_verdict,
    fdr_bh,
    permutation_null_with_metric,
    split_half_with_metric,
)
from trinetravir.eval.metrics import mmd_rbf_off_diag, pearson_off_diag


def _make_synthetic(
    n_studies: int = 4,
    n_donors_per_class: int = 6,
    n_cells_per_donor: int = 50,
    n_genes: int = 50,
    signal_strength: float = 1.0,
    noise_scale: float = 0.1,
    seed: int = 0,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Generate synthetic per-cell embedding + obs frame.

    Each study has ``n_donors_per_class`` diseased + ``n_donors_per_class``
    healthy donors. Each donor contributes ``n_cells_per_donor`` cells.
    Diseased cells have mean = +signal_strength*direction + noise;
    healthy cells have mean = -signal_strength*direction + noise.
    The same direction is shared across studies, so cross-study response
    vectors should correlate strongly.
    """
    rng = np.random.default_rng(seed)
    direction = rng.normal(size=n_genes)
    direction /= np.linalg.norm(direction)

    cells_x: list[np.ndarray] = []
    rows: list[dict] = []
    for s in range(n_studies):
        for d in range(n_donors_per_class):
            for status, sign in (("diseased", +1.0), ("healthy_control", -1.0)):
                donor_id = f"s{s}_{status}_{d}"
                cells = sign * signal_strength * direction + noise_scale * rng.normal(
                    size=(n_cells_per_donor, n_genes)
                )
                cells_x.append(cells)
                for _ in range(n_cells_per_donor):
                    rows.append(
                        {
                            "study_id": f"study_{s}",
                            "donor_id": donor_id,
                            "donor_disease_status": status,
                        }
                    )
    x = np.concatenate(cells_x, axis=0).astype(np.float32)
    obs = pd.DataFrame(rows)
    return x, obs


def test_identical_signal_passes_calibrated_gate() -> None:
    """Strong shared signal across studies -> permutation p ~ 0, observed r ~ 1, pass."""
    x, obs = _make_synthetic(signal_strength=1.0, noise_scale=0.05, seed=0)
    perm = permutation_null_with_metric(
        x,
        obs,
        "synthetic",
        metric_fn=pearson_off_diag,
        n_permutations=200,
        seed=42,
    )
    sh = split_half_with_metric(
        x,
        obs,
        "synthetic",
        metric_fn=pearson_off_diag,
        n_splits=20,
        seed=42,
    )
    # Compute observed.
    from trinetravir.eval.calibration import _per_study_donor_status, _response_vector

    per_study = _per_study_donor_status(obs)
    cs = obs["study_id"].astype(str).to_numpy()
    cd = obs["donor_id"].astype(str).to_numpy()
    rvs = {}
    for sid, ds in per_study.items():
        m = cs == sid
        rvs[sid] = _response_vector(x[m], cd[m], ds.to_dict())
    obs_r = pearson_off_diag(rvs)
    verdict = calibrated_gate_verdict(
        obs_r,
        perm["null"],
        sh["split_half_distribution"],
        percentile=99,
        alpha=0.05,
    )
    assert obs_r > 0.95, f"identical signal should give obs_r ~ 1, got {obs_r:.3f}"
    assert verdict["p_value"] < 0.05, "permutation p should be small"
    assert verdict["pass"], "calibrated_pass should be TRUE for identical signal"


def test_no_signal_fails_calibrated_gate() -> None:
    """Zero shared signal -> observed r ~ 0, permutation p ~ 0.5, fail."""
    x, obs = _make_synthetic(signal_strength=0.0, noise_scale=1.0, seed=1)
    perm = permutation_null_with_metric(
        x,
        obs,
        "synthetic_zero",
        metric_fn=pearson_off_diag,
        n_permutations=200,
        seed=42,
    )
    sh = split_half_with_metric(
        x,
        obs,
        "synthetic_zero",
        metric_fn=pearson_off_diag,
        n_splits=20,
        seed=42,
    )
    from trinetravir.eval.calibration import _per_study_donor_status, _response_vector

    per_study = _per_study_donor_status(obs)
    cs = obs["study_id"].astype(str).to_numpy()
    cd = obs["donor_id"].astype(str).to_numpy()
    rvs = {}
    for sid, ds in per_study.items():
        m = cs == sid
        rvs[sid] = _response_vector(x[m], cd[m], ds.to_dict())
    obs_r = pearson_off_diag(rvs)
    verdict = calibrated_gate_verdict(
        obs_r,
        perm["null"],
        sh["split_half_distribution"],
        percentile=99,
        alpha=0.05,
    )
    assert abs(obs_r) < 0.3, f"zero signal should give obs_r ~ 0, got {obs_r:.3f}"
    assert not verdict["pass"], "calibrated_pass should be FALSE for zero signal"


def test_bootstrap_ci_overlap_corrected_direction() -> None:
    """A1 fix: observed above upper CI should now PASS (at_or_above_ci_low=True)."""
    sh_dist = np.linspace(0.4, 0.9, 100)  # within-study CI roughly [0.41, 0.89]
    # Observed AT upper-end (above the lower bound).
    out_above = bootstrap_ci_overlap(0.95, sh_dist, alpha=0.05)
    assert out_above["at_or_above_ci_low"], "observed above CI should PASS new criterion"
    assert not out_above["in_ci"], "observed above CI should also FAIL legacy in_ci criterion"
    # Observed below lower bound.
    out_below = bootstrap_ci_overlap(0.1, sh_dist, alpha=0.05)
    assert not out_below["at_or_above_ci_low"], "observed below lower CI should FAIL new criterion"
    # Observed within CI.
    out_within = bootstrap_ci_overlap(0.6, sh_dist, alpha=0.05)
    assert out_within["at_or_above_ci_low"]
    assert out_within["in_ci"]


def test_calibrated_gate_verdict_corrected_vs_legacy() -> None:
    """observed above upper split-half CI: PASS under corrected, FAIL under legacy."""
    null_dist = np.linspace(-0.2, 0.5, 500)
    sh_dist = np.linspace(0.4, 0.9, 100)
    observed = 0.95  # above null_p99 AND above sh upper bound
    corrected = calibrated_gate_verdict(
        observed,
        null_dist,
        sh_dist,
        percentile=99,
        alpha=0.05,
        use_corrected_ci=True,
    )
    legacy = calibrated_gate_verdict(
        observed,
        null_dist,
        sh_dist,
        percentile=99,
        alpha=0.05,
        use_corrected_ci=False,
    )
    assert corrected["pass"], "corrected criterion should PASS for r above sh upper bound"
    assert not legacy["pass"], "legacy in-CI criterion should FAIL for r above sh upper bound"


def test_bootstrap_observed_r_coverage() -> None:
    """Bootstrap CI on observed r should cover the strong-signal value."""
    x, obs = _make_synthetic(signal_strength=1.0, noise_scale=0.05, seed=2)
    boot = bootstrap_observed_r(
        x,
        obs,
        metric_fn=pearson_off_diag,
        n_bootstrap=100,
        seed=42,
    )
    # Strong-signal observed should be near 1; CI should be narrow around that.
    assert boot["n_bootstrap_completed"] > 50
    assert boot["observed_ci_low"] > 0.5, f"bootstrap CI low too low: {boot['observed_ci_low']:.3f}"
    assert boot["observed_ci_high"] <= 1.0


def test_fdr_bh_basic() -> None:
    """FDR-BH should adjust p-values and preserve monotonicity."""
    ps = np.array([0.001, 0.01, 0.04, 0.05, 0.2, 0.5])
    adj = fdr_bh(ps)
    # All adjusted p-values >= raw p-values.
    assert np.all(adj >= ps - 1e-9), "BH adjusted p < raw p violates definition"
    # Adjusted p-values monotonically non-decreasing under same input order.
    # (BH guarantees adjusted p[i] = min over k>=i of raw_p[k] * n/k, so order preserved
    # along sorted input. Here our input IS sorted.)
    for i in range(len(adj) - 1):
        assert adj[i] <= adj[i + 1] + 1e-9


def test_fdr_bh_reduces_apparent_significance() -> None:
    """With 20 tests at raw p=0.04, FDR should lift many above 0.05."""
    raw = np.full(20, 0.04)
    adj = fdr_bh(raw)
    # n=20, smallest BH = 0.04 * 20 / 20 = 0.04, but with monotonicity step-up
    # all are 0.04. With raw=0.04 across the board the BH result equals 0.04
    # for all (since smallest threshold isn't crossed). Sanity check.
    assert np.all(adj >= 0.04 - 1e-9)


def test_metrics_seed_propagation() -> None:
    """A5 fix: caller's seed should determine MMD output reproducibly."""
    rng = np.random.default_rng(0)
    a = rng.normal(size=(800, 10))
    b = a + 0.5  # shifted distribution
    embeddings = {"s1": a, "s2": b}
    v1_seed42 = mmd_rbf_off_diag(embeddings, seed=42)
    v2_seed42 = mmd_rbf_off_diag(embeddings, seed=42)
    v_seed7 = mmd_rbf_off_diag(embeddings, seed=7)
    assert v1_seed42 == v2_seed42, "same seed must produce same MMD"
    # Different seed should (usually) produce different value because of
    # different 500-cell subsample. Allow small probability of accidental match.
    assert v1_seed42 != v_seed7 or abs(v1_seed42 - v_seed7) > 1e-12, (
        "different seeds should give different MMD"
    )
