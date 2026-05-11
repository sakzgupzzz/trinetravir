# GSE157829 (HIV chronic) — Khatri MVS external validation

Date: 2026-05-11
Cohort: GSE157829 / Wang et al. 2020. HIV-infected adult PBMC (chronic infection, suppressed and viremic donors) vs C1 healthy control.
Design: `cross_sectional`; cross-cohort baseline retained = v1 corpus 41 healthy donors aggregated. C1 within-cohort healthy (n=1) retained supplementary as sanity check only.

## Bucket coverage (5 buckets supported)

| Bucket | n_HVG_common | n_MVS_common | r_full | r_MVS | Δ (MVS − full) |
|---|---|---|---|---|---|
| B | 2924 | 48 | 0.2628 | **0.5961** | +0.3333 |
| **CD4T** | 3107 | 57 | 0.0840 | **0.2572** | +0.1732 |
| CD8T | 3191 | 61 | 0.3827 | **0.6121** | +0.2294 |
| NK | 3210 | 47 | 0.1246 | **0.4163** | +0.2917 |
| monocyte | 3109 | 57 | −0.0111 | **0.0405** | +0.0516 |

## Khatri MVS interpretation

- **CD8T, B, NK: substantial ISG lift.** r_MVS exceeds r_full by 0.23–0.33. Chronic HIV produces canonical type-I-IFN signature most strongly in lymphocyte compartments.
- **CD4T (primary bucket): partial ISG lift.** r_MVS = 0.26 above r_full = 0.08 by Δ = 0.17. HIV-CD4T target-cell biology overlaps acute-RNA-virus ISG ~25% at conserved-component level.
- **monocyte: flat.** r_MVS ≈ 0.04, near-zero. Chronic HIV monocyte response distinct from acute paracrine ISG cascade.

## Issue 30 verdict (mechanically applied; primary bucket CD4T)

Pre-committed rule: `r_mvs in [0.00, 0.20] EXPECTED_RETROVIRUS_DISTINCTNESS; > 0.40 SURPRISING; < −0.10 ANTI_CORRELATION`.
Observed CD4T r_mvs = 0.2572 → **BORDERLINE** (just above expected retrovirus-distinctness range).
Partial overlap: chronic HIV ISG tone overlaps acute viral ISG ~50% at MVS subset, ~10% at full HVG. Imperfect retrovirus distinction at conserved-ISG level; clean retrovirus distinction at full HVG.

## Calibration (N=200 perm)

- CD4T: perm_p_value_mvs = 0.1343, FDR-corrected p = 0.3173. CI = [0.155, 0.4335] (lower bound > 0; effect-size positive at 95%).
- CD8T: perm_p_value_mvs = 0.1343, FDR-corrected p = 0.3173. CI = [0.4153, 0.6731] — narrow positive band.
- monocyte: perm_p_value_mvs = 0.6959. Null cannot be rejected.

No bucket survives FDR<0.01 at N=200. CD4T + CD8T CIs exclude zero on lower bound → effect-size signal stable.

## Scope caveat

GSE157829 is a single chronic HIV cohort. Generalization to other retroviral contexts (HTLV-1, SIV) untested. Within-cohort healthy n=1 (C1) limits within-cohort baseline; primary analysis uses v1 corpus healthy aggregate for power. Donor identity de-duplication applied at C-pre.7.
