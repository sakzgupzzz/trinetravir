# SESSION 5 — audit response (resolves audit findings, opens Issue 25, revises Issues 3, 7, 9)

**Ordering**: runs AFTER Session 3 closeout (subbucket + Issue 2 committed) and BEFORE Session 3.5. CPU-only. Decision-documentation + calibration framework correctness fixes + external validation. 6-10 hours focused work.

## Prerequisites

- [ ] Session 3 complete (commit 233c068)
- [ ] Session 3 closeout complete (Issue 2 resolved, phase35_subbucket.csv committed)
- [ ] Hostile-reviewer audit findings reviewed and acknowledged

## Constraint

Audit response only. No Phase 4 work. No Phase 5 model implementation. No scVI (Session 4). Stop and report. The strategic v1 framing decision (Issue 25) is human-only — Session 5 opens Issue 25, presents the options, and stops.

---

## PART A — Critical calibration framework correctness fixes

### A1. Bootstrap CI direction fix

Current `calibrated_gate_verdict` in `src/trinetravir/eval/calibration.py:715-753` treats observed r **within** the split-half 95% CI as PASS and outside CI as FAIL. This is wrong — above the upper CI represents BETTER-than-within-study coherence (interesting, not failure). The NK per-cell-type FAIL for r=0.3845 < sh_ci_low=0.4422 demonstrates the bug.

**Fix**: change criterion from "within CI" to "≥ lower CI bound". Two criteria:
1. observed r ≥ 99th percentile of permutation null (unchanged)
2. observed r ≥ lower bound of split-half ceiling 95% CI (revised)

Re-run all calibrated verdicts. Save to `results/tables/calibration_*_v2.csv`. Old tables preserved for audit trail; new tables marked `_v2`.

### A2. Bootstrap CI on observed r

Add bootstrap CI on observed r via donor-level resampling (B=1000 iterations). Each iteration: resample donors with replacement, recompute per-study response vectors, recompute cross-study mean off-diagonal r.

**Implementation**: extend `permutation_null_with_metric` or add `bootstrap_observed_r` to calibration.py. Add CI columns to all calibration output tables.

### A3. Multiple-testing correction

Current calibration reports 5 buckets × 4 metrics = 20 tests as independent p-values. Add FDR-BH correction on permutation null p-values per dataset.

**Implementation**: add `fdr_corrected_p` column to all calibration output tables. Update calibrated_pass criteria to use FDR-corrected p < 0.01 for p99 equivalent.

### A4. test_calibration.py

Write `src/tests/test_calibration.py` with synthetic ground-truth examples:
- Test 1: Two identical response vectors → permutation null p ≈ 1, observed r ≈ 1, calibrated_pass = TRUE
- Test 2: Two orthogonal response vectors → permutation null p ≈ 0, observed r ≈ 0, calibrated_pass = FALSE
- Test 3: Controlled noise → permutation null behaves as expected
- Test 4: Split-half ceiling on known-replicable signal recovers high r
- Test 5: Bootstrap CI on observed r covers true value at expected coverage rate
- Test 6: FDR correction reduces apparent significance proportional to number of tests

All tests use seeded numpy RNG from function parameter (not metrics.py's hardcoded `default_rng(0)` — fixed in A5).

### A5. metrics.py seed leak fix

`src/trinetravir/eval/metrics.py:172` hardcodes `rng=default_rng(0)` ignoring caller seed. Fix to thread caller's seed through. Add test in test_calibration.py demonstrating caller-seed propagation.

---

## PART B — External validation via Khatri MVS reanalysis

Khatri MVS is currently cited as theoretical anchor but never empirically tested. Audit identifies as critical gap.

### B1. Acquire Khatri MVS

- MVS module gene list (~400 genes) → `data/reference/khatri_mvs_module_genes.txt`
- ≥1 independent external PBMC bulk RNA-seq cohort: viral vs healthy. Influenza preferred (matches Lee IAV anchor); SARS-CoV-2 acceptable; other RNA viruses fallback.

Source: Andres-Terre et al. 2015 Immunity + follow-up. Datasets typically GEO.

### B2. External validation analysis

Limited because v1 factorized model not yet implemented (Phase 5 post-Session 4). External validation scope for Session 5:
- (a) Show v1's harmonized cross-study response vectors correlate with Khatri MVS module score on external cohort.
- (b) Show calibration framework's permutation null is not artifactually tight on external data.

Outputs:
- `results/tables/external_validation_khatri.csv`
- `references/notes/external_validation_summary.md`

### B3. Document the result honestly

**Success** (r ≥ 0.5 on external cohort): document as supporting evidence for calibration framework.

**Failure** (r < 0.3 or unexpected): document explicitly. **Do NOT proceed to Session 3.5/4/Phase 4 until failure mode is understood.** Reviewer ammunition is worse than reviewer-finding-it-themselves.

---

## PART C — Reframe Phase 3 as exploratory/discovery

Audit identifies Phase 3 gate thresholds as fit-to-data (HARKing-light). The PLAN.md annotation showing each threshold sits above the pre-Harmony r value is post-hoc rationalization.

### C1. Update METHODS_CHOICES.md

For Issues 3 and the Phase 3 results in general, add to each relevant issue's "Validation strategy" noting Phase 3 thresholds were set after observing Harmony output and are exploratory/discovery, NOT confirmatory.

Add new **Issue 26** (Phase 3 threshold provenance — RESOLVED at acknowledgment level):

```
The choice as it stands: Phase 3 buckets were declared PASS/FAIL using
thresholds annotated post-Harmony as "above the pre-Harmony r." This
is fit-to-data, not pre-specification.

Acknowledgment: Phase 3 results are reframed as exploratory/discovery
evidence. The Phase 3 PASS/FAIL verdicts indicate which buckets have
signal worth pursuing in downstream phases, not which buckets have
confirmed cross-study coherence at pre-specified thresholds.

Forward commitment: Phase 5 thresholds will be set from external
literature (Khatri MVS r=0.45 for module preservation; other cited
literature anchors) BEFORE running Phase 5. The Phase 5 pre-
registration commits to thresholds and to the v1 paper's primary
claims before any Phase 5 calibration runs.

Validation: methods section reports Phase 3 results as exploratory
and Phase 5 results as confirmatory at pre-registered thresholds.

Date opened: <today>
Date resolved: <today> (at acknowledgment level; full validation at
Phase 5 launch)
```

### C2. Update PLAN.md

Add §1.8 (or similar) titled "Exploratory vs confirmatory distinction" explicitly noting:
- Phases 1-3 produced exploratory evidence (heuristic thresholds, data-driven decisions).
- Phase 4 onwards produces confirmatory evidence with pre-registered protocols.
- The v1 paper distinguishes these clearly in methods and results.

---

## PART D — Revise post-hoc rationalizations in resolved issues

### D1. Issue 7 (per-cell-type vs global Harmony)

Current resolution: "the pre-specified rule would favour Global" followed by methodological-preference override. This is post-hoc.

**Revision options** (Claude Code picks one):

**Option 1**: Accept the pre-specified rule's verdict. Switch to global Harmony as v1 primary. Per-cell-type reported as sensitivity. Honest but reverses several Session 3 decisions.

**Option 2**: Acknowledge override is post-hoc. Rewrite Issue 7 resolution to lead with "the pre-specified rule favored global; we override to per-cell-type for methodological reasons; this is a post-hoc decision and we report it as such."

**Option 3**: Run both protocols in parallel through Phase 5+ and report both in the paper. Slightly more work but fully defensible.

**Recommendation**: Option 3 if compute allows; Option 2 otherwise. Document choice with rationale.

### D2. Issue 3 (DE-Jaccard dismissal)

Current resolution frames DE-Jaccard as a "different question" (top-100 ranking vs full vector). Audit calls this cherry-pick.

**Revision**: rewrite DE-Jaccard discussion to lead with the **degeneracy explanation** (top-100 of 50-dim PCA = all; metric is mathematically degenerate on this embedding) rather than "different question" framing. The degeneracy is a real methodological problem, not a post-hoc dismissal.

---

## PART E — Open Issue 25: v1 paper framing decision

Audit identifies v1's "first systematic cross-virus benchmark" framing as unsupportable with current corpus (4 SARS-CoV-2 studies + 1 IAV study within Lee).

Open **Issue 25** in METHODS_CHOICES.md with two options for human decision. Session 5 does NOT make this decision; opens issue and presents analysis.

```
Issue 25 (v1 paper framing — OPEN, requires human decision)

The choice as it stands: PLAN.md frames v1 as "cross-virus
generalization for single-cell host response prediction" with H1-H5
about cross-virus transfer learning. Current corpus (Wilk, Lee,
Arunachalam, Schulte-Schrepping) contains 4 SARS-CoV-2 studies and 1
IAV study (Lee). RSV and HSV/CMV planned but not acquired/harmonized.

Why this matters: audit identifies that "cross-virus transfer
learning" claims require multiple non-SARS-CoV-2 studies AND multiple
IAV studies. With n=1 IAV study, demonstrated cross-virus result (Lee
within-study SARS-vs-IAV monocyte r=0.651) is single within-study
data point, not benchmark evidence.

Resolution required (human decision):

Option A: Reframe v1 honestly. v1 becomes "PBMC SARS-CoV-2 cross-
study harmonization benchmark with Lee within-study cross-virus
exploration." Factorized model demonstrates methodology on SARS-CoV-2
cross-study task with Lee IAV exploration as cross-virus data point.
v1.5 then becomes the proper cross-virus paper after acquiring
additional viral data.

Pro: defensible at peer review. Honest about what data supports.
Allows v1 to ship in roughly planned timeline (Phase 4-7 over 8-12
weeks).
Con: smaller-claim paper. Cross-virus framing was the project's
novelty hook. Reframing loses some.

Option B: Acquire additional viral data before v1 ships. Add ≥1 more
IAV study (or RSV, or HSV/CMV) meeting Issue 4 criteria. Re-run
harmonization, Phase 3.5 re-annotation, Phase 3 calibration. Update
PLAN.md scope. Ship v1 with original framing.

Pro: preserves original cross-virus framing. Stronger paper.
Con: 2-4 weeks additional data acquisition + harmonization before
Phase 4. Pushes v1 timeline out. Reintroduces scope expansion.

Decision authority: human only. Session 5 stops at Issue 25 open;
human reviews and decides; subsequent sessions (3.5 revised, 4) are
re-scoped based on decision.

Validation: chosen option's methodology pre-registered before Phase 5
launch.

Date opened: <today>
Date resolved: <pending human decision>
```

---

## PART F — STOP AND REPORT

Report:

1. Calibration framework v2 — `calibration_phase3_v2.csv`, `calibration_phase35_low_v2.csv`, `calibration_phase35_high_v2.csv`, `calibration_phase3_global_harmony_v2.csv`. Compare to v1 tables; document where corrected bootstrap CI direction changes any verdict. Document where FDR-BH correction changes any verdict. Document per-bucket observed-r CI ranges.
2. test_calibration.py passing with synthetic ground truth tests. metrics.py seed leak fixed.
3. Khatri MVS external validation outputs. Verdict: succeeded / failed / unexpected. Specific evidence cited.
4. Issue 26 added (Phase 3 threshold provenance acknowledgment). PLAN.md updated with §1.8 exploratory-vs-confirmatory section.
5. Issue 7 resolution revised per Option 1, 2, or 3 chosen with rationale. Issue 3 DE-Jaccard discussion revised to lead with degeneracy explanation.
6. Issue 25 opened with both framing options in full and analysis.
7. State after Session 5:
   - 16 Sessions 1-3 issues resolved (1, 2, 3, 4, 5, 7-revised, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)
   - 1 Session 5 issue resolved (26)
   - 1 Session 5 issue open (25, awaiting human decision)
   - 1 still open (6, Session 4)
   - Session 3.5 and Session 4 BLOCKED until Issue 25 decision
8. Recommended next step: human reviews Issue 25 options and decides. After decision, Session 3.5 (potentially revised) and Session 4 can be revisited.

## Time budget

6-10 hours focused work.
- Part A (calibration fixes): 2-3 hours
- Part B (Khatri MVS): 2-3 hours incl. data acquisition
- Part C (Phase 3 reframing): 1 hour
- Part D (post-hoc revisions): 1-2 hours
- Part E (Issue 25 open): 30 min
- Part F (stop and report): 30 min

## Constraints

- No Phase 4 / Phase 5 / scVI work.
- Atomic schema-change rule (Issue 17): calibration framework v2 + tests + result tables in same commits.
- Old _v1 calibration tables preserved; new _v2 tables produced alongside.
- If Khatri MVS validation fails decisively (r < 0.3 or unexpected), STOP and report before Part C. Failure may indicate calibration framework problems beyond audit findings.
- Issue 25 is human-only decision. Session 5 opens and stops.
- Session 3.5 and Session 4 BLOCKED until Issue 25 is resolved.
