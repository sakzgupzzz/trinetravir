# Session 6B — held-out calibrated evaluation + per-stratum sensitivity + few-shot

Multi-session, 3-5 weeks total. Per Block #5 of SHORT_TERM_PLAN.md.

## Status snapshot (2026-05-11)

- Part A+B initial sweep: ✅ DONE (commit `eb64886`). 4 cohorts response-vector comparison vs training corpus on gene-space, no permutation null yet.
- Parts A+B follow-up + C + D + E + F + G: pending.

## Next-session sequence (executed mechanically; pre-committed rules)

### Step 1 — Fix Yoshida gene naming
- Apply `var['gene_symbol']` remap before MVS intersection in `scripts/heldout_calibrated_evaluation.py`. Yoshida cellxgene h5ad uses Ensembl IDs as var index; the `var.gene_symbol` column has the symbols.
- Recompute r_MVS for all 5 Yoshida buckets. **Issue 28 verdict blocked on this fix.**
- Expected: pediatric-stratum monocyte r_MVS in [0.20, 0.50] (per-pre-spec decision rule range); adult-stratum monocyte r_MVS similar.

### Step 2 — Investigate Randolph monocyte r_MVS anomaly
- Current: r_full = 0.287, r_mvs = 0.013 (substantially LOWER, opposite of expected ISG lift pattern).
- Check: gene-symbol dedup in `harmonize_randolph_2021.py` `load_geo_pool()` — `pd.DataFrame.var.duplicated()` may have dropped MVS genes preferentially. Inspect.
- Check: which Randolph monocyte cells contribute? If `infected_monocytes` ≠ `monocytes` are pooled together in bucket, may dilute ISG signal.
- If technical fix restores pattern: rerun + report normal r_MVS.
- If anomaly persists after technical fixes: this is a REAL biological finding (Randolph 6h ex vivo IAV monocyte at MOI 0.5 may have *direct cell-autonomous* infection signature dominating, not paracrine ISG). Document as finding rather than bug.

### Step 3 — Run permutation null + bootstrap CI + FDR-BH (Part B full calibration)
- Extend `scripts/heldout_calibrated_evaluation.py` to use Session 5 v2 framework:
  - `permutation_null_with_metric` on held-out cohort donor labels (donor-level shuffle).
  - `bootstrap_observed_r` for observed-r CI (Session 5 v2 fix).
  - `bootstrap_ci_overlap` for ≥ lower CI bound criterion (Session 5 v2 fix, NOT in-CI).
  - `fdr_bh` across all (cohort × bucket × metric) tests.
- Per-cohort design overrides (Issues 27-30):
  - **Randolph (paired_within_donor)**: permutation shuffles condition WITHIN donor, not across donors. Apply Issue 27 amendment exclusion (HMN83575 healthy excluded primary; 89/90 donors). Sensitivity rows for ≥50, ≥100, no-exclusion thresholds.
  - **GSE157829**: cross-cohort baseline = v1 corpus 41 healthy donors aggregated (not within-cohort C1 healthy). C1 supplementary sanity check only.
  - **Allen Atlas**: CMV+ vs CMV- (chronic_latent vs naive); Children stratum already filtered.
  - **Yoshida**: two strata (pediatric primary + adult sensitivity). Each computes its own held-out vector + permutation null.
- Output: extend `heldout_calibration_<cohort>.csv` columns: `perm_p99`, `perm_p_value`, `observed_ci_low`, `observed_ci_high`, `fdr_corrected_p`, `calibrated_pass_p99_fdr`.

### Step 4 — Mechanical application of pre-committed decision rules
- Apply Issues 27/28/29/30 decision rules to v2 calibrated verdicts. No interpretation latitude.
  - Issue 27 (Randolph): monocyte cross-context MVS r ≥ 0.40 supports H1; r < 0.20 challenges.
  - Issue 28 (Yoshida): pediatric monocyte cross-age MVS r ≥ 0.30 supports H1; r < 0.10 challenges.
  - Issue 29 (Allen CMV): monocyte chronic-latent-vs-naive MVS r ∈ [0.10, 0.40] = appropriate; >0.50 over-prediction concerning; <0.05 no shared biology concerning.
  - Issue 30 (GSE157829 HIV): CD4T MVS r ∈ [0.00, 0.20] expected retrovirus distinctness; >0.40 surprising; <-0.10 anti-correlation interpretable.

### Step 5 — Report
- `results/tables/heldout_calibration_<cohort>.csv` × 4 with full v2 columns.
- `results/tables/heldout_vs_training_comparison.csv` cross-table with calibrated PASS/FAIL/PARTIAL verdicts.
- For each Issue 27/28/29/30: side-by-side write-up of pre-committed rule + observed value + verdict, in METHODS_CHOICES.md resolution section.

## Remaining Parts (post Steps 1-5)

- **Part C MVS per-cohort write-up**: largely covered by Steps 1-3 output; needs formal `references/notes/heldout_khatri_mvs_<cohort>.md` per cohort summarizing ISG-lift findings.
- **Part D Per-stratum f_shared**: separate ML task (train per-cohort factorized model; compare f_shared cosine similarity). Heavy compute (4-5 model trainings). Stop and assess before launching.
- **Part E Few-shot adaptation**: requires placeholder v1 factorized model (Issue 21 pre-spec architecture). Heavy compute + framework build. Multi-session.
- **Part F Issue resolutions**: 27/28/29/30 resolved with calibrated evidence + biological caveats.
- **Part G Final 6B report**.

## Constraints

- Pre-committed numerical decision rules in Issues 27-30 must be applied mechanically. No post-hoc adjustment.
- Atomic schema-change rule (Issue 17) applies to all calibration output additions.
- If projection (Part A gene-space) fails decisively for any cohort (held-out r near zero on full HVG across all buckets after Yoshida gene-name fix), STOP and document — flag for Session 4 scVI sensitivity.
- If per-stratum f_shared similarity < 0.3 across cohorts in Part D, STOP for human decision.
- If few-shot fails to converge with N=1000 for multiple cohorts, document as domain-of-validity finding.
