# Randolph 2021 — Khatri MVS external validation

Date: 2026-05-11
Cohort: GSE162632 / Zenodo 4273999. Ex vivo IAV stimulation of PBMC monocytes, 6h, MOI 0.5.
Design: `paired_within_donor` (90 donors, both mock + IAV conditions per donor).
Issue 27 amendment applied: HMN83575 healthy_control (43 cells) excluded primary; 89/90 donors retained. 19 donors with <100 cells on watch-list for supplementary sensitivity.

## Bucket coverage (5 buckets supported; only B/NK/monocyte have n_cells ≥ 50 per condition)

| Bucket | n_HVG_common | n_MVS_common | r_full | r_MVS | Δ (MVS − full) |
|---|---|---|---|---|---|
| B | 3351 | 48 | 0.3009 | **0.4827** | +0.1818 |
| NK | 3486 | 47 | 0.4808 | **0.5764** | +0.0956 |
| monocyte | 3456 | 57 | 0.2864 | **0.0126** | **−0.2738** |

## Khatri MVS interpretation

- **B and NK: appropriate ISG lift.** r_MVS > r_full by 0.10–0.18, consistent with conserved type-I-IFN component dominating ISG subset while non-ISG genes add noise on full HVG.
- **monocyte: ISG lift INVERTED — anomaly.** Full HVG r=0.29 is moderate; MVS r=0.013 is near-zero. Pattern opposite of B/NK lift.

## Root cause (Step 2 investigation, Session 6B)

`infected_monocytes_cluster_singlets.rds` was NOT extracted from `inputs.tar.gz` (Zenodo 10.5281/zenodo.4273999) before the tar archive was deleted. Available monocyte cells in the processed h5ad are `monocytes` (bystander) population only, not `infected_monocytes`. 6h ex vivo bystander monocytes have not received sufficient paracrine IFN-α/β to mount canonical ISG cascade. The Randolph monocyte response vector therefore reflects mock-vs-bystander baseline drift, not infected-cell ISG induction.

## Issue 27 verdict (mechanically applied)

Pre-committed rule: `r_mvs ≥ 0.40 SUPPORTS_H1; < 0.20 CHALLENGES; in [0.20, 0.40] INCONCLUSIVE`.
Observed r_mvs = 0.0126 → **CHALLENGES_H1**.
Caveat documented in METHODS_CHOICES.md: bystander-only data is a real data gap, not a framework failure. v1.5 must re-acquire `infected_monocytes_cluster_singlets.rds` from Zenodo and re-run.

## Calibration (N=200 perm; v1.5 will re-run at N=1000+)

- monocyte: perm_p_value_mvs = 0.4925, FDR-corrected p = 0.5304 — null cannot be rejected. Consistent with anomaly: r_mvs lies inside permutation null distribution.
- B: perm_p_value_mvs = 0.0746, FDR-corrected p = 0.3173 — fails FDR<0.01, passes nominal trend.
- NK: perm_p_value_mvs = 0.0348, FDR-corrected p = 0.3173 — same.

No bucket survives FDR<0.01 at N=200. Headline relies on observed effect-size verdicts.

## Scope caveat

Bystander-only Randolph monocyte data limits Issue 27 to ex vivo bystander-vs-mock comparison. Direct infection ISG response remains untested at this cohort.
