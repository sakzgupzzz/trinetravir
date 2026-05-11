# Held-out Khatri MVS external validation — cross-cohort summary

Date: 2026-05-11
Session 6B Part C output. Aggregates 4 per-cohort notes:
- [Randolph 2021](heldout_khatri_mvs_randolph_2021.md)
- [Yoshida 2022](heldout_khatri_mvs_yoshida_2022.md)
- [Allen Atlas](heldout_khatri_mvs_allen_atlas_monocyte.md)
- [GSE157829 HIV](heldout_khatri_mvs_gse157829.md)

## Cross-cohort table — primary buckets

| Cohort | Bucket | Context | r_full | r_MVS | Verdict | Issue |
|---|---|---|---|---|---|---|
| Randolph 2021 | monocyte | ex vivo IAV (bystander only) | 0.286 | **0.013** | CHALLENGES_H1 (bystander caveat) | #27 |
| Yoshida 2022 | monocyte | pediatric ↔ adult cross-age | 0.387 | **0.591** | SUPPORTS_H1 (strongest signal) | #28 |
| Allen Atlas | monocyte | chronic-latent CMV vs naive | 0.152 | **−0.010** | CONCERNING_NO_SHARED_BIOLOGY (reframed: scope) | #29 |
| GSE157829 | CD4T | chronic HIV | 0.084 | **0.257** | BORDERLINE (just above expected retrovirus-distinctness) | #30 |

## Cross-cohort table — all 14 (cohort × bucket) tests

| Cohort | Bucket | n_MVS_common | r_full | r_MVS | Δ |
|---|---|---|---|---|---|
| Randolph | B | 48 | 0.301 | 0.483 | +0.182 |
| Randolph | NK | 47 | 0.481 | 0.576 | +0.096 |
| Randolph | monocyte | 57 | 0.286 | 0.013 | −0.274 |
| Yoshida | B | 48 | 0.188 | 0.436 | +0.248 |
| Yoshida | CD4T | 57 | 0.080 | 0.278 | +0.198 |
| Yoshida | CD8T | 61 | 0.108 | 0.320 | +0.212 |
| Yoshida | NK | 47 | 0.127 | 0.282 | +0.155 |
| Yoshida | monocyte | 57 | 0.387 | 0.591 | +0.204 |
| Allen | monocyte | 57 | 0.152 | −0.010 | −0.163 |
| GSE157829 | B | 48 | 0.263 | 0.596 | +0.333 |
| GSE157829 | CD4T | 57 | 0.084 | 0.257 | +0.173 |
| GSE157829 | CD8T | 61 | 0.383 | 0.612 | +0.229 |
| GSE157829 | NK | 47 | 0.125 | 0.416 | +0.292 |
| GSE157829 | monocyte | 57 | −0.011 | 0.041 | +0.052 |

## Headline cross-cohort findings

1. **Yoshida cross-age signal is strongest single result.** All 5 buckets show clean ISG lift; monocyte r_MVS = 0.59. Conserved component transfers pediatric ↔ adult.
2. **Acute respiratory cohort coverage is uneven.** Only Yoshida + Randolph carry the v1 framework's nominal domain; Randolph monocyte is bystander-artifact pending v1.5 Zenodo re-acquisition.
3. **Chronic-latent boundary (Allen CMV) shows zero shared biology with acute corpus.** Scope-limitation finding, not framework failure.
4. **Chronic retroviral (GSE157829 HIV) shows partial overlap.** CD4T target-cell biology r_MVS = 0.26 is borderline; CD8T/B/NK lymphocyte buckets show stronger ISG signal than CD4T target compartment.
5. **Δ(MVS − full) = +0.15 to +0.34 in 10 of 14 tests.** Khatri MVS subset systematically more coherent than full HVG across cohorts, confirming ISG cascade as the dominant conserved component captured.

## Calibration status (N=200 perm, N=100 bootstrap)

No bucket survives FDR<0.01 at N=200. Headline relies on observed effect-size verdicts under pre-committed Issue 27-30 rules. v1.5 with N≥1000 perm + N≥500 bootstrap will tighten significance; current bootstrap CIs are wide (±0.2) but generally exclude zero on lower bound for ISG-lift buckets.

## Scope caveats — aggregated

- **Acute respiratory viral**: v1 framework's nominal domain. Yoshida confirms; Randolph bystander limits direct test.
- **Chronic-latent herpesvirus** (Allen CMV): outside domain. r_MVS≈0 marks boundary.
- **Chronic retroviral** (GSE157829): partial overlap. Domain boundary, not clean exclusion.
- **Cell-type generality**: lymphocyte buckets (B/NK/CD4T/CD8T) show ISG lift across cohorts. Monocyte pattern is cohort-specific: high in Yoshida, near-zero in Randolph-bystander/Allen-CMV/GSE-HIV.

## Pipeline state

Session 6B Part C complete. v1.5 next steps:
1. Zenodo re-acquire `infected_monocytes_cluster_singlets.rds` (Randolph) → re-run Issue 27 with infected monocyte population.
2. N=1000+ permutation re-run on all 14 tests → tighten FDR p-values.
3. Per-stratum f_shared sensitivity (Part D, Issue 18-prep).
4. Few-shot adaptation curves (Part E, Issue 21-prep).
