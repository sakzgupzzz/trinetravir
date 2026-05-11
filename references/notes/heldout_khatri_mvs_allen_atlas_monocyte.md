# Allen Atlas (CMV serostatus) — Khatri MVS external validation

Date: 2026-05-11
Cohort: Allen Institute / OneK1K subset. Healthy donor PBMC with CMV+/− serostatus stratification. Chronic-latent herpesvirus context.
Design: `cross_sectional`; `donor_disease_status` derived from `infection_state` (CMV+ = chronic_latent diseased; CMV− = naive healthy_control). Adult-only stratum (Children excluded per Issue 29 amendment, applied at C-pre.6).

## Bucket coverage (monocyte only — other buckets failed n_cells ≥ 50 gate)

| Bucket | n_HVG_common | n_MVS_common | r_full | r_MVS | Δ (MVS − full) |
|---|---|---|---|---|---|
| monocyte | 3460 | 57 | 0.1523 | **−0.0102** | **−0.1625** |

## Khatri MVS interpretation

- **No ISG lift; near-zero / mildly anti-correlated MVS signal.** r_MVS = −0.01 vs r_full = 0.15. Inverted Δ.
- Pattern indicates Khatri ISG subset does NOT capture shared biology between chronic-latent CMV monocyte response and v1 corpus (which is acute-RNA-virus-dominant: SARS-CoV-2, IAV, RSV).
- Chronic CMV serostatus → low-level innate stimulation pattern is biologically distinct from acute paracrine IFN-α/β cascade. No reason to expect canonical type-I-IFN ISG induction at chronic-latent steady state.

## Issue 29 verdict (mechanically applied)

Pre-committed rule: `r_mvs in [0.10, 0.40] APPROPRIATE; > 0.50 OVER_PREDICTION; < 0.05 NO_SHARED_BIOLOGY`.
Observed r_mvs = −0.0102 → **CONCERNING_NO_SHARED_BIOLOGY**.

## Reinterpretation as scope-limitation finding

CONCERNING verdict was pre-specified before knowing chronic-vs-acute distinction would map cleanly to r_MVS≈0. The verdict is mechanically correct under the rule, but the underlying biology is the *expected* outcome for an acute-disease-specific framework probed against chronic-latent context. Reframed (per METHODS_CHOICES.md Issue 29 resolution):

- v1 framework's domain of validity = **acute respiratory viral infection** (paracrine IFN-α/β cascade <72h post-infection).
- Chronic-latent herpesvirus → not in domain. r_MVS≈0 is appropriate boundary marker, not framework failure.
- v1.5 + v2 should formalize "domain of validity" as an explicit deliverable.

## Calibration (N=200 perm)

- monocyte: perm_p_value_mvs = 0.4925, FDR-corrected p = 0.5304. CI = [−0.4469, 0.4197] (very wide, centered near zero).
- Null cannot be rejected. Consistent with no shared biology.

## Scope caveat

Single-bucket coverage (monocyte only) limits this cohort to a one-data-point test of chronic-latent boundary. v1.5 may relax the n_cells ≥ 50 gate to recover B/NK/T buckets and confirm whether the chronic-naive pattern holds across cell types.
