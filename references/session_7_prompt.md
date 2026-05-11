# Session 7 — Pre-Modeling Sensitivity Audit

**Status:** Pre-Phase-4 audit gate. Drafted 2026-05-11 after Session 6B partial completion.
**Trigger:** Critique document concerns 4 (Harmony preserving only conserved axes) and partial concern 3 (cross-study integration assumption) need empirical answers before Phase 4 modeling builds on the harmonized embedding.
**Pattern:** Parallel to Session 5's calibration audit response. Pre-specified decision rules set BEFORE running analyses. Results disclosed regardless of outcome.
**Estimated wall-time:** 1-2 Claude Code sessions + 1 audit confirmation chat. ~2-3 weeks calendar time.

---

## Motivation

The current paper framing rests on cross-study response-vector coherence under per-cell-type Harmony correction. The critique document raised a substantive concern: Harmony may be preserving only dominant conserved axes (ISGs) and artificially inflating cross-study coherence, in which case the ISG-conservation finding could be partly an integration artifact rather than biological signal.

This concern is not addressable via re-running the existing calibration framework — the framework operates on post-Harmony embeddings by design. Two new analyses are needed:

1. **Pre/post-harmonization response vector comparison.** Quantify how much of the observed cross-study Pearson r comes from Harmony alignment versus underlying biological coherence.

2. **Within-cohort-only sensitivity.** Verify that the response vector findings (especially ISG conservation in lymphoid populations) hold within individual cohorts without any cross-study integration.

Both analyses are run-once, report-results. The pre-committed decision rules below classify the outcomes BEFORE any computation runs.

---

## Pre-conditions (must be complete before Session 7 launches)

- [ ] Session 6B fully closed: Issue 27 corrected verdict (post-Issue-31 monocyte_infected re-run), Issue 28-30 verdicts all locked at N=1000
- [ ] Issue 31 (cross-bucket healthy reference for cluster-defined subsets) committed to METHODS_CHOICES.md
- [ ] MANUSCRIPT_DRAFT.md Sections 3-5 updated with N=1000 verdicts
- [ ] Yoshida CI caveat language inserted per Saksham's verbatim wording
- [ ] All Session 6B commits pushed to origin/main with audit-confirmation
- [ ] This Session 7 prompt itself committed to repo as `references/session_7_prompt.md` (atomic Issue 17)

If any pre-condition is missing, **Session 7 does not launch.** This is enforced per Issue 17 atomic discipline.

---

## Scope overview

Two analyses, both empirical, both pre-Phase-4 critical:

| Part | Analysis | Purpose | Wall-time estimate |
|---|---|---|---|
| A | Pre/post-harmonization response vector comparison | Quantify Harmony's contribution to cross-study coherence | ~2-3h compute + 1h reporting |
| B | Within-cohort-only sensitivity | Verify findings hold without cross-study integration | ~2-3h compute + 1h reporting |
| C | Audit gate | Verify deliverables, atomic commits, manuscript updates | ~30 min chat |

Out of scope (deferred to Phase 5/6 supplementary work):

- Alternative MVS gene set (Khatri vs Interferome 2.0)
- HMN83575 cell-count sensitivity (rule already pre-specified)
- Per-batch effect quantification beyond study-level
- Alternative HVG selection methods

---

## Part A — Pre/post-harmonization response vector comparison

### Scope

Compute per-bucket cross-study Pearson r on response vectors derived from **pre-Harmony** counts (raw normalized log1p) and compare to the equivalent post-Harmony result from Session 5's calibration_phase3_v2.csv.

**Per-cohort method:**

1. Load v1 corpus harmony h5ads
2. For each bucket (monocyte, CD4T, CD8T, NK, B):
   - Compute per-donor pseudo-bulk response vectors using `X_raw` (pre-Harmony normalized counts) — diseased mean minus healthy mean per donor, then average across donors per study
   - Compute equivalent vectors using `X_harmony` (post-Harmony embedding) — already cached from Session 5
   - Compute cross-study Pearson r on each
3. Report Δr = r_post − r_pre per bucket
4. Repeat for full-HVG and MVS-restricted analyses

**Deliverable:** `results/tables/sensitivity_pre_post_harmony.csv` with columns: bucket, gene_set (full / MVS), r_pre, r_post, delta_r, n_studies, n_donors_total.

### Pre-committed decision rule (locked before observation)

| Δr (post minus pre) | Interpretation | Implication for manuscript |
|---|---|---|
| Δr ≤ 0.10 across most buckets | Harmony adds minor smoothing; biological coherence is dominant | ISG-conservation finding holds as biology; manuscript framing unchanged |
| Δr ∈ (0.10, 0.30] across most buckets | Harmony adds substantial smoothing; mixed biological + integration contribution | ISG-conservation finding holds qualitatively but framing must explicitly acknowledge Harmony contribution; sensitivity analysis becomes a load-bearing limitation rather than reassurance |
| Δr > 0.30 across most buckets | Harmony does most of the work | ISG-conservation finding requires significant revision; manuscript must reframe as "post-harmonization coherence" rather than "biological conservation"; reconsider whether the finding is publishable in current form |

**"Most buckets" = ≥3 of 5 buckets.** Mixed patterns (some buckets ≤0.10, others >0.30) trigger per-bucket disclosure rather than aggregate verdict.

**Restricted to MVS gene set:** Δr_MVS interpretation has higher stakes because the ISG-restriction finding is the load-bearing contribution. If Δr_MVS > 0.30 specifically, the methodology contribution claim weakens substantially.

### What this Part DOES NOT test

- It does NOT test whether Harmony is "correct" or "good." Harmony is doing its job of removing batch effects; the question is how much of the observed cross-study coherence comes from that vs from real biology.
- It does NOT test the held-out cohort transfer findings. Those were computed using the same Harmony pipeline on held-out data; their interpretation is separate.

---

## Part B — Within-cohort-only sensitivity

### Scope

Run the full calibration framework on each v1 cohort independently. Compute per-cohort cross-bucket response vector comparisons. Verify that the per-cohort effects align with the cross-study harmonized findings.

**Per-cohort method:**

1. For each v1 cohort (Wilk, Lee, Arunachalam, Schulte-Schrepping):
   - Load raw normalized counts (no Harmony)
   - For each bucket: compute per-donor pseudo-bulk response vectors (diseased minus healthy)
   - Apply v2 calibration framework: permutation null (N=500 for time efficiency), bootstrap CI (N=200), observed Pearson r between buckets within-cohort
   - Generate per-cohort within-study coherence statistics
2. Aggregate: do the within-cohort coherence patterns mirror the cross-study harmonized patterns?

**Deliverable:** `results/tables/sensitivity_within_cohort.csv` with columns: cohort, bucket_pair, observed_r, perm_p_raw, bootstrap_ci_low, bootstrap_ci_high, gene_set (full / MVS).

Plus aggregate comparison: `results/tables/sensitivity_within_vs_cross.csv` with columns: bucket_pair, mean_within_cohort_r, cross_study_harmonized_r, sign_concordance, magnitude_alignment.

### Pre-committed decision rule (locked before observation)

| Pattern | Interpretation | Implication for manuscript |
|---|---|---|
| Within-cohort effects align with cross-study (sign concordance ≥80%, mean within-r within 0.20 of cross-study r) | Biology is consistent within and across studies; cross-study findings reflect real signal | ISG-conservation finding holds; cross-study integration is a useful tool but not creating the signal |
| Within-cohort effects partially align (sign concordance 50-80%, magnitude divergence 0.20-0.50) | Biology shows within-cohort but cross-study integration changes magnitudes | Findings hold qualitatively; manuscript discusses cross-study integration as amplifying rather than creating the signal |
| Within-cohort effects disappear or reverse (sign concordance <50%, or systematic magnitude reversal) | Cross-study findings are artifactual or dependent on integration | Major reframing required; the finding becomes "post-harmonization analysis surfaces coherence not visible in raw within-cohort data" — much weaker contribution |

**Sign concordance** = fraction of bucket pairs where within-cohort mean r and cross-study harmonized r have the same sign.

**Magnitude alignment** = mean absolute difference |r_within - r_cross| across bucket pairs.

### What this Part DOES NOT test

- It does NOT replicate Phase 3 calibration. The within-cohort comparison is between buckets WITHIN a cohort, not across studies. This is a different statistic than the cross-study analysis.
- It does NOT use Harmony embeddings. The point is to verify findings without harmonization.

---

## Part C — Audit gate

### Verification checklist

Before Session 7 closes and Session 3.5 unblocks:

- [ ] `sensitivity_pre_post_harmony.csv` exists with rows for all bucket × gene_set combinations
- [ ] `sensitivity_within_cohort.csv` exists with rows for all cohort × bucket_pair × gene_set combinations
- [ ] `sensitivity_within_vs_cross.csv` aggregate table exists
- [ ] Pre-committed decision rules from this prompt verbatim quoted in METHODS_CHOICES.md (new Issue 32: pre/post-harmonization sensitivity; new Issue 33: within-cohort sensitivity)
- [ ] Verdicts mechanically applied per decision rules; no post-hoc reframing
- [ ] MANUSCRIPT_DRAFT.md Section "Methods supplementary" updated with sensitivity analysis results
- [ ] MANUSCRIPT_DRAFT.md "Limitations" section updated to reflect what the sensitivity analyses revealed
- [ ] If Δr_MVS > 0.30 across most buckets OR within-cohort effects disappear: MANUSCRIPT_DRAFT.md gets a substantive reframing pass before Session 3.5 launches
- [ ] All commits atomic per Issue 17
- [ ] Pushed to origin/main

### Human-audit confirmation

Saksham reviews the verdict and confirms:
- The verdicts match the pre-committed rules
- The manuscript updates reflect honest disclosure of what the sensitivity analyses found
- No post-hoc rationalization snuck in
- Session 3.5 is genuinely unblocked (no remaining sensitivity questions)

If any of the above fail, Session 7 stays open and Session 3.5 does not launch.

---

## Issues to be opened in METHODS_CHOICES.md

### Issue 32: Pre/post-harmonization sensitivity analysis design

**Choice:** Compute response vectors on pre-Harmony (raw normalized log1p) counts and post-Harmony embeddings; quantify Δr per bucket per gene set. Pre-committed decision rule classifies Harmony's contribution at thresholds Δr ≤ 0.10 (biology dominant), Δr ∈ (0.10, 0.30] (mixed), Δr > 0.30 (Harmony dominant).

**Rationale:** Cross-study integration could artificially inflate coherence by preserving only dominant conserved axes. Quantifying Harmony's contribution separates biological signal from integration smoothing.

**Validation:** Test_calibration.py extended with synthetic ground-truth case: known-correlated synthetic data with study-batch noise; verify pre-Harmony r < post-Harmony r as expected; verify Δr magnitude scales with noise level.

**Date opened:** [Session 7 launch]
**Date resolved:** [Session 7 close]

### Issue 33: Within-cohort-only sensitivity analysis design

**Choice:** Run v2 calibration framework on each v1 cohort independently (no cross-study integration). Compute per-cohort cross-bucket response vector comparisons. Aggregate to sign concordance and magnitude alignment metrics. Pre-committed decision rule classifies alignment at thresholds sign concordance ≥80% with magnitude divergence ≤0.20 (biology consistent), partial alignment, and disappear/reverse patterns.

**Rationale:** If within-cohort effects don't replicate cross-study findings, the cross-study coherence may be an integration artifact. Within-cohort effects are a more conservative baseline.

**Validation:** Per-cohort calibration uses same v2 framework as cross-study (permutation null, bootstrap CI), just restricted to within-cohort data. Framework-level validation already complete (test_calibration.py 8/8 passing).

**Date opened:** [Session 7 launch]
**Date resolved:** [Session 7 close]

---

## Atomic commit discipline

Per Issue 17, each commit groups: code + tests + result tables in a single atomic commit. Suggested commit sequence:

1. **Commit:** Open Issues 32 + 33 in METHODS_CHOICES.md (pre-spec before run). Include the pre-committed decision rule tables verbatim. Add `references/session_7_prompt.md` to repo.

2. **Commit:** Part A code + sensitivity_pre_post_harmony.csv + test extension for pre-Harmony response vector computation.

3. **Commit:** Part B code + sensitivity_within_cohort.csv + sensitivity_within_vs_cross.csv + test extension for within-cohort framework usage.

4. **Commit:** Apply pre-committed decision rules; update Issue 32 + 33 resolution in METHODS_CHOICES.md with verdicts. Update MANUSCRIPT_DRAFT.md Methods supplementary + Limitations sections.

5. **Commit (conditional):** If Δr > 0.30 OR within-cohort disappears, substantive MANUSCRIPT_DRAFT.md reframing in separate atomic commit before Session 3.5 launches.

6. **Push** all commits to origin/main with audit-confirmation in chat.

---

## Deliverables checklist

After Session 7 closes, the following exist and are committed:

- [ ] `results/tables/sensitivity_pre_post_harmony.csv`
- [ ] `results/tables/sensitivity_within_cohort.csv`
- [ ] `results/tables/sensitivity_within_vs_cross.csv`
- [ ] `METHODS_CHOICES.md` Issues 32 and 33 with pre-committed rules and applied verdicts
- [ ] `MANUSCRIPT_DRAFT.md` Methods supplementary section with sensitivity analysis results
- [ ] `MANUSCRIPT_DRAFT.md` Limitations section updated to reflect findings
- [ ] `src/tests/` extensions for pre-Harmony response vector and within-cohort framework usage
- [ ] `references/session_7_prompt.md` (this document)
- [ ] All commits atomic per Issue 17, pushed to origin/main

---

## Timeline estimate

| Step | Wall-time | Blocking |
|---|---|---|
| Pre-conditions verified | 30 min | Session 6B closure |
| Issues 32 + 33 pre-spec committed | 30 min | Pre-conditions |
| Part A code + run | 2-3h | Pre-spec |
| Part B code + run | 2-3h | Part A (sequential to avoid resource contention) |
| Decision rule application | 30 min | Part A + B |
| MANUSCRIPT_DRAFT.md updates | 1-2h | Decision rules |
| Audit gate confirmation | 30 min chat | All above |
| Session 3.5 unblocks | - | After audit confirmation |

Total: ~8-10h Claude Code time across 1-2 sessions, plus chat audit. Calendar time ~1-2 weeks including audit confirmation latency.

---

## Conditional outcomes

### If Part A shows Δr ≤ 0.10 across most buckets AND Part B shows sign concordance ≥80%

**Outcome:** Manuscript framing holds. Sensitivity analyses become a reassurance in the Methods supplementary. Limitations section can confidently state that the ISG-conservation finding reflects biology rather than integration artifact.

**Time impact:** Minimal. Session 3.5 launches on schedule.

### If Part A shows Δr ∈ (0.10, 0.30] OR Part B shows partial alignment

**Outcome:** Manuscript framing holds qualitatively but requires acknowledgment that Harmony contributes meaningfully to cross-study coherence. The ISG-conservation finding becomes "the cross-study integration framework, including per-cell-type Harmony correction, produces lymphoid response vectors that align with the canonical ISG signature; raw within-cohort analysis shows partial alignment, indicating biology is present but integration amplifies the signal."

**Time impact:** ~1 week additional for MANUSCRIPT_DRAFT.md substantive revision. Session 3.5 launches after revision committed.

### If Part A shows Δr > 0.30 OR Part B shows within-cohort effects disappear

**Outcome:** Major reframing required. The finding becomes "post-harmonization analysis reveals coherence patterns not visible in raw within-cohort or pre-Harmony cross-study data." This is a substantially weaker contribution. Three options:

1. **Reframe v1 around methodology contribution only.** The calibration framework remains a contribution; the ISG-conservation biology framing gets demoted.

2. **Investigate alternative harmonization methods.** Try scVI (Issue 6 sensitivity already pre-specified) or scgen for v1 alongside Harmony. If alternative harmonizations show smaller Δr, the finding may be Harmony-specific rather than general.

3. **Acknowledge as critical limitation and proceed.** Honest disclosure in limitations; let reviewers decide whether the contribution is sufficient.

**Time impact:** ~2-4 weeks additional. Session 3.5 launches only after a decision among options 1/2/3 is made and committed.

---

## How this session relates to Session 5 audit pattern

Session 5 was the audit response to hostile-reviewer concerns about calibration framework correctness. Session 7 follows the same pattern for harmonization-induced bias concerns. Both:

- Triggered by external/critical input
- Run before the next pipeline block (Session 5 before Session 6; Session 7 before Session 3.5)
- Use pre-committed decision rules to prevent post-hoc rationalization
- Result in MANUSCRIPT_DRAFT updates that disclose the analyses regardless of outcome
- Close with audit gate before unlocking the next block

The pattern is intentional and protects the manuscript's credibility against predictable reviewer concerns.

---

## Notes on what Session 7 deliberately does not address

- **Factorized model necessity (critique concern 2):** Tested in Phase 5/6 via baseline comparisons. Already scheduled.
- **Held-out cohort biological independence (critique concern 3):** Empirically tested in Session 6B (Allen CMV scope-limitation demonstrates real boundary).
- **PBMC tissue limitations (critique concern 5):** Scope statement, not a testable hypothesis at v1 level. Disclosed in limitations.
- **Alternative gene sets / cell-count thresholds:** Phase 5/6 supplementary work.
- **scVI sensitivity (Issue 6):** Already scheduled for Session 4.

Session 7 is specifically about concern 4 (Harmony preserving only conserved axes) plus the within-cohort baseline. Scope is intentionally tight.

---

## How to use this prompt

When ready to launch Session 7:

1. Verify all pre-conditions in the checklist above are met
2. Commit this prompt to the repo as `references/session_7_prompt.md` (if not already)
3. Open a fresh Claude Code session with: "Launch Session 7 per references/session_7_prompt.md. Open Issues 32 and 33 with the pre-committed decision rules verbatim before running any analyses."
4. Audit-confirm in chat after each Part completes
5. Audit-confirm final Session 7 close before launching Session 3.5

Apply Issue 17 atomic commit discipline throughout.
