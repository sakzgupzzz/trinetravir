# Session 4.5 prompt — Audit-response cleanup + manuscript propagation + Part B global scVI

**Status:** Drafted 2026-05-12 during Session 4. Session 4.5 follows Session 4 close (Option Narrow scope) per pipeline rev 3.

This session closes the audit-response work that emerged during Session 4 (Issues 38, 39, Amendment 1), runs the deferred Part B global scVI sensitivity via a precomputed-cache unblock for the 16 GB Modal CLI upload failure, and propagates findings into `MANUSCRIPT_DRAFT.md` before Phase 5 modeling implementation begins.

**Integration choice:** standalone file pattern per Session 4 + 3.5 + 7 precedent. Commit as `references/session_4_5_prompt.md`.

---

## Entry conditions

**At Session 4.5 entry (post-Session 4 close):**

- Session 4 closed under Option Narrow scope: Issue 6 resolved (Tier IV HARMONY_PREFERRED), Issue 35 resolved (INCLUDE foundation models per projection-based assessment), Issue 34 Amendment 1 (two-sided Tier I) landed.
- Pipeline 8/8 done per pipeline rev 3.
- 33 base issues + Amendment 1 + Issues 36/37/38/39 in METHODS_CHOICES.md.
- Wilk depth watchpoint clean (Part A.5b complete in Session 4; no artifact concern).
- `scvi_input_global.h5ad` (118 MB) already in Modal volume from Session 4 (uploaded successfully).
- `harmony_global_embedding.h5ad` (16 GB) Modal CLI upload BLOCKED (two attempts, empty logs). Session 4.5 implements the cache-based unblock.

**Empirical anchors carried in:**

- Sessions 5/6B/7 calibration verdicts at canonical N=1000 (Issue 38 partial recompute confirmed CIs slightly wider at canonical N, revealing N=200 over-confidence in earlier reports — validates canonical N=1000 standard retroactively).
- Tier IV HARMONY_PREFERRED verdict from Session 4 Part A (locked; Amendment 1 audit trail preserved).
- Issue 34 Amendment 1: Tier I corrected to two-sided `max(|Δr|) ≤ 0.05` rule.
- Issue 35 INCLUDE verdict for Geneformer + scGPT in Phase 7 baselines (projection-based per Theodoris 2023 + Cui 2024 throughputs; not measured empirically because pre-spec math >10× under budget).
- Issue 39 Wilk mito sensitivity: status pending Session 4.5 escalation decision.

**Part B infrastructure preserved as untracked files in working tree from Session 4:**

- `scripts/extract_global_counts.py` — works; generated `scvi_input_global.h5ad` (118 MB, in Modal volume)
- `scripts/session4_part_b_global_sweep.py` — complete code, tested locally on small subset
- `scripts/session4_part_b_modal.py` — Modal wrapper ready

These files need modification per Part F protocol (cache-based approach, not 16 GB h5ad upload). They move from untracked to tracked as part of F.2's atomic commit.

---

## SESSION 4.5 — Audit-response cleanup + manuscript propagation + Part B global

```
CPU + ~$5-10 GPU. Decision-pre-specification + cleanup + manuscript work
+ deferred GPU sweep via cache-based unblock. Roughly 5-7 hours focused
work distributed across 1-2 chat sessions.

Prerequisite: Session 4 closed cleanly under Option Narrow (Issue 6 +
Issue 35 + Amendment 1 resolved; commit #13 pipeline closure landed).
Issues 38, 39, Amendment 1 manuscript propagation, and deferred Part B
carried forward per Option Narrow scoping decision.

CONSTRAINT: Audit-response cleanup + manuscript propagation + deferred
Part B global scVI only. No Phase 5 modeling work. No new methodological
pre-specs beyond Issue 39 decision rule. Atomic schema-change rule
(Issue 17) applies. Pre-registration discipline preserved: Issue 39
decision rule committed BEFORE measuring Wilk mito fraction.

═══════════════════════════════════════════════════════════════════

PART A — Randolph Issue 31 canonical N=1000 recompute

  Issue 31 (cross-bucket healthy reference design for cluster-defined
  subsets) was originally committed during Session 6B at heuristic N.
  Issue 38 audit (Session 4 audit-response sweep) flagged that Sessions
  5/6B/7 results should all be at canonical N=1000 for cross-cohort
  verdict-table consistency.

  v1 corpus + Yoshida + Allen + GSE157829 already at canonical N per
  Issues 36/37/38 audit-response work. Randolph remains pending — Part
  A closes the gap.

  Inputs:
    data/processed/randolph_2021_processed_v6.h5ad
    references/khatri_mvs_gene_list.csv (Khatri 86-gene Table S2)
    3 buckets only: monocyte, B, NK (CD4T/CD8T deferred to v1.5 per
      Issue 27 C-pre.6 amendment).

  Output:
    results/tables/session4_5_randolph_issue31_canonical.csv
    Updated verdict for Issue 27 monocyte (Randolph cross-context IAV)
    applies pre-committed rule from Issue 27 (≥0.40 SUPPORTS_H1, <0.20
    CHALLENGES, intermediate INCONCLUSIVE).

  Decision rule: pre-committed in Issue 27. No new pre-spec needed.

  Wall-time: ~30-60 min CPU. Atomic commit.

═══════════════════════════════════════════════════════════════════

PART B — Issue 38 verdict-comparison table completion

  Cross-cohort verdict consistency check at canonical N=1000. Partial
  table (3 of 4 cohorts) already produced in commit 793344e. Part B
  closes the table by adding the Randolph row from Part A.

  Coverage:
    - v1 corpus (Sessions 5/7 reanalysis at canonical N)
    - Yoshida 2022 (Issue 28)
    - Allen Immune Health Atlas (Issue 29)
    - GSE157829 (Issue 30)
    - Randolph 2021 (Part A output)

  Output:
    results/tables/session4_5_verdict_comparison_canonical.csv
    Per-row: cohort, primary bucket, verdict at heuristic N, verdict at
      canonical N=1000, verdict change flag.

  Verdict consistency assessment: verify Sessions 6B headline verdicts
  are robust to canonical N. Document any verdict changes with
  rationale.

  If ≥1 verdict flip: open Issue 40 documenting the N-sensitivity
  finding. If all verdicts stable: confirm Sessions 6B headline
  conclusions hold at canonical N.

  Wall-time: ~30 min CPU. Atomic commit.

═══════════════════════════════════════════════════════════════════

PART C — Issue 39 Wilk mito sensitivity decision

  Issue 39 was opened during Session 4 audit-response sweep based on
  Wilk's distinct sequencing depth profile (158-293 vs 1500-3600 for
  other studies). Wilk depth watchpoint (Part A.5b in Session 4)
  verified depth is non-confounding for Harmony/scVI integration.
  Question remaining: does Wilk's MITOCHONDRIAL gene proportion differ
  from other studies in a way that could affect scVI's NB likelihood
  inference?

  Pre-commit decision rule BEFORE measurement:

    Compute mito-gene fraction per study (using MT- prefix on gene
    symbols, n=37 mt genes in 4000-HVG space if present).

    Tier I — REDUNDANT_DEFER:
      Wilk mito fraction within 1 SD of mean(other studies)
      -> Issue 39 closed as redundant given Part A.5b watchpoint
         outcome. No new sensitivity analysis. Documented as v1
         limitation in manuscript discussion.

    Tier II — RUN_SENSITIVITY:
      Wilk mito fraction > 1 SD from mean(other studies)
      -> Run per-study response vector correlation excluding mito
         genes; verify |delta_r_no_mito| <= 0.05 across all 5
         buckets.
      -> If passes: Issue 39 resolved as "no mito-driven artifact."
      -> If fails: escalate to Phase 5 supplementary with explicit
         manuscript disclosure.

    Tier III — ESCALATE_NOW:
      Wilk mito fraction > 2 SD from mean(other studies)
      -> Run mito sensitivity AND consider whether Tier IV
         HARMONY_PREFERRED verdict requires Wilk-specific caveat
         beyond the existing depth caveat.

  Rationale: Wilk depth watchpoint verified one dimension of QC
  concern. Mito sensitivity is a parallel but distinct dimension. Pre-
  committing the rule prevents post-hoc escalation framing.

  Threshold anchoring note: 1 SD / 2 SD thresholds are conventional
  but pulled from intuition. If Wilk mito fraction is dramatically
  outside expected range and the thresholds feel mis-calibrated,
  amend pre-commit before measurement (data-direction-independent
  fix per Amendment 1 precedent).

  Output:
    results/tables/session4_5_wilk_mito_sensitivity.csv
    Issue 39 resolution committed atomically per tier.

  Wall-time: 30 min measurement + variable depending on tier (max ~2h
  if Tier III). Two atomic commits (decision rule pre-commit, then
  resolution post-measurement).

═══════════════════════════════════════════════════════════════════

PART D — Manuscript propagation (Steps 6a, 6b, 6d, 6e, then 6c)

  Apply Session 4 + 4.5 findings to MANUSCRIPT_DRAFT.md. EXECUTION
  ORDER: 6a -> 6b -> 6d -> 6e -> 6c. Step 6c (smoke-test) runs LAST
  so it can validate all prior changes within the same session.

  Step 6a — Issue 38 canonical N propagation:
    Update reported CIs in manuscript tables where calibration N
    affects the values. Reference results/tables/session4_5_*.csv as
    new source of truth for canonical N values. Methods section notes
    that canonical N=1000 is used throughout for consistency. The N=200
    over-confidence finding (CIs were ~10-20% narrower at heuristic N
    than at canonical) is documented as audit-response justification
    for the canonical N standard.

  Step 6b — Session 7 sensitivity audit propagation:
    Apply Issue 32 BIOLOGY_DOMINANT verdict + Issue 33
    BIOLOGY_CONSISTENT verdict to manuscript discussion section.
    Harmony lift disclosure: pre-Harmony cross-study r already
    substantial (0.13-0.58 across buckets); Harmony adds 0.02-0.25 on
    top (BIOLOGY_DOMINANT for load-bearing monocyte MVS at delta_r=
    0.08; MIXED for full HVG space). Methods section discloses this
    honestly per post-Session-7 framing.

  Step 6d — Issue 34 Amendment 1 propagation:
    Methods section discloses the Tier I two-sided correction. Frame:
    "Pre-spec gap surfaced post-launch; structural fix applied (data-
    direction-independent) before final verdict computation. Amendment
    commit hash 9ba7c25 + verdict commit a92404a available for audit."

  Step 6e — Issue 6 resolution propagation:
    Methods section reflects Tier IV HARMONY_PREFERRED verdict with
    three supporting lines: (1) literature anchor (scIB Luecken 2022 +
    Briefings 2022 finding that scVI ~ Harmony on immune integration
    with this paper landing on the Harmony-favored end of the range),
    (2) Wilk depth watchpoint clean per Part A.5b, (3) Amendment 1
    disclosure per Step 6d.

  Step 6c — Manuscript smoke-test (LAST):
    Verify all cited issue numbers exist in METHODS_CHOICES.md. Verify
    all referenced result tables exist on disk. Verify all empirical
    claims trace to committed result files. If practical, CI workflow
    at .github/workflows/manuscript-smoke-test.yaml that fails on
    dangling references. Smoke-test catches dangling refs introduced
    by 6a/6b/6d/6e propagation in the same session.

  Each step is its own atomic commit per Issue 17, OR 6a+6b+6d+6e
  bundle into one commit + 6c its own commit if changes are small.

═══════════════════════════════════════════════════════════════════

PART F — Part B global scVI (deferred from Session 4)

  Background: Session 4's original Part B path uploaded harmony_
  global_embedding.h5ad (16 GB) to Modal and ran scVI global. The
  16 GB Modal CLI upload exceeded the practical limit and failed
  silently in two attempts (b3am0lm30, brha6vfiv) during Session 4.
  Session 4.5 implements the cache-based unblock.

  F.1 — Build precomputed Harmony response vector cache locally

    Inputs: harmony_global_embedding.h5ad (16 GB, stays local; never
    uploaded). harmony_per_celltype_<bucket>.h5ad files for bucket
    assignment if needed.

    Process: For each cell in harmony_global_embedding, compute the
    per-(bucket, study, donor_disease_status) mean of
    obsm['X_harmony']. Aggregate into a small parquet table:
    ~5 buckets x 4 studies x 2 statuses = 40 rows x n_harmony_dim
    columns. Total ~5 MB.

    Bucket assignment from obs['coarse'] (NOT cell_type_bucket — per
    Issue 7 status note).

    Output:
      data/processed/harmony_global_response_vector_cache.parquet
      Schema: bucket (str), study_id (str), donor_disease_status
      (str), and one column per Harmony PC (n_harmony_dim).

    Atomic commit (cache file + build script).
    Wall-time: ~30 min CPU (single pass over 244K cells).

  F.2 — Modify session4_part_b_global_sweep.py to use cache

    Update the existing untracked script to:
      - Load harmony_global_response_vector_cache.parquet instead of
        full harmony_global_embedding.h5ad as the Harmony comparator.
      - scVI training still uses scvi_input_global.h5ad (118 MB,
        already in Modal volume from Session 4).
      - delta_r comparison computes scVI global response vectors per-
        (bucket, study, status) and compares against cached Harmony
        global response vectors.

    Files moving from untracked to tracked: the 3 Part B infra
    scripts (extract_global_counts.py, session4_part_b_global_sweep
    .py, session4_part_b_modal.py) now commit cleanly as part of
    F.2's atomic unit.

    Atomic commit. Wall-time: ~30 min CPU.

  F.3 — Upload cache + launch Modal sweep

    Upload data/processed/harmony_global_response_vector_cache
    .parquet (~5 MB, well within Modal CLI limit) to /inputs/ on
    Modal volume. Launch session4_part_b_global_sweep.py via Modal.

    Wall-time: ~2-4h GPU. Cost: ~$5-10.

  F.4 — Verdict commit

    Apply four-tier verdict structure per Session 4 spec
    (HARMONY_ADEQUATE / MIXED / SCVI_PREFERRED / HARMONY_PREFERRED).
    Supplementary verdict (doesn't drive Issue 6 headline). Expected
    outcome: similar to Part A (Harmony preferred in global mode
    too) but could surprise.

    If verdict differs substantively from Part A per-bucket result,
    document the asymmetry in manuscript discussion section. If
    similar, manuscript cites Part B as supplementary confirmation
    that Part A's per-bucket result holds in scVI's native global
    mode.

    Atomic commit.

═══════════════════════════════════════════════════════════════════

PART E — STOP AND REPORT (audit gate)

Report:

  1. Randolph Issue 31 canonical N verdict + table committed (Part A).
  2. Issue 38 cross-cohort verdict-comparison table complete (Part B).
     Verdict consistency: stable at canonical N OR documented changes.
     (Issue 40 opened only if verdict flips surfaced.)
  3. Issue 39 decision: REDUNDANT_DEFER / RUN_SENSITIVITY /
     ESCALATE_NOW per pre-committed rule (Part C). Resolution
     committed.
  4. Manuscript propagation steps 6a + 6b + 6d + 6e committed in
     execution order, then 6c (smoke-test) committed last and
     validates all prior changes (Part D).
  5. Part B global scVI verdict + commit (Part F). Supplementary
     tier. Cache-based unblock executed successfully.
  6. State after Session 4.5:
     - All Session-4-emergent issues resolved
     - Part B global scVI completed via cache-based approach
     - 0 open methodological issues
     - MANUSCRIPT_DRAFT.md reflects all empirical findings +
       amendments
     - Audit trail clean: zero placeholders, zero pending markers
     - Ready for Phase 5 (factorized model implementation per PLAN.md
       numbering; "Phase 4" per handoff numbering)
  7. Recommended next step: Phase 5 launch in fresh chat.

CONSTRAINTS

- No Phase 5 modeling work. No model code. No training runs.
- Part B global scVI is supplementary; doesn't drive Issue 6 headline
  verdict (already locked at Tier IV from Session 4 Part A).
- Atomic schema-change rule (Issue 17) applies to each commit.
- Pre-registration discipline preserved: Issue 39 decision rule
  committed BEFORE measuring Wilk mito fraction (Part C).
- Manuscript propagation respects post-Session-7 framing: primary
  calibration framework + secondary factorized model + empirical ISG-
  conservation as replicated finding with Harmony lift disclosed.
- Smoke-test (Step 6c) runs LAST in manuscript propagation to validate
  all prior changes within the same session.
- 16 GB Modal CLI upload path stays dead. Cache-based approach (~5 MB)
  is the only Part B path in this session.
```

---

## Issues opened / resolved in Session 4.5

**Resolved:**

- Issue 31 (Randolph cross-bucket healthy reference) at canonical N — Part A
- Issue 38 (verdict-comparison table) — Part B closes the table at canonical N
- Issue 39 (Wilk mito sensitivity) — per Part C decision rule

**Conditionally opened:**

- Issue 40 (N-sensitivity finding) — only if Part B surfaces verdict flips between heuristic N and canonical N. Otherwise not opened.

**Deferred work completed:**

- Part B global scVI sensitivity (Session 4 → 4.5 deferral resolved via cache-based approach).

---

## Atomic commit sequence (expected 10-12 commits)

1. `references/session_4_5_prompt.md` (this document)
2. Randolph Issue 31 canonical N recompute result + verdict (Part A)
3. Issue 38 verdict-comparison table + consistency assessment (Part B)
4. Issue 39 decision rule pre-commit (Part C pre-measurement)
5. Issue 39 resolution per measured tier (Part C post-measurement)
6. F.1 — Harmony global response vector cache (`data/processed/harmony_global_response_vector_cache.parquet` + build script)
7. F.2 — Part B sweep script + Modal wrapper modified to use cache (lifts 3 untracked files into tracked)
8. F.3 / F.4 — Part B sweep results + verdict commit
9. Manuscript steps 6a + 6b + 6d + 6e (canonical N + Session 7 + Amendment 1 + Issue 6 — bundled or atomic per file)
10. Manuscript step 6c (smoke-test runs last; validates all prior changes)
11. (Optional) Issue 40 N-sensitivity finding if verdict flips surfaced
12. Session 4.5 pipeline closure commit

---

## Timeline estimate

- Spec drafting: this commit
- Part A: ~30–60 min CPU
- Part B: ~30 min CPU
- Part C: 30 min measurement + variable depending on tier (max ~2h if Tier III)
- Part D: ~1–2h focused manuscript work
- Part F (Part B global scVI): ~30 min cache build + ~30 min script modification + ~2–4h GPU async + ~30 min verdict
- Audit gate: ~1h chat session
- **Total: 1–2 chat sessions, ~5–7h focused work + 2–4h GPU async**

---

## Pattern relation to prior sessions

Session 4.5 follows the cleanup-plus-deferred-work pattern: Sessions 1–3 base methodology → Session 5 audit → Session 6A/B held-out → Session 7 sensitivity → Session 3.5 pre-specs → Session 4 GPU sensitivity (Option Narrow close) → **Session 4.5 audit-response cleanup + manuscript propagation + deferred Part B** → Phase 5 modeling implementation.

Distinct from earlier sessions in two ways: (a) absorbs work that emerged during Session 4 audit-response sweep (Issues 38/39/Amendment 1) under the Option Narrow scope decision, and (b) executes the cache-based Part B unblock that the original 16 GB Modal CLI upload path couldn't deliver.

---

## How to start Session 4.5

In fresh Claude.ai chat (recommended — Session 4 chat will be bloated by close):

1. Paste `trinetravir_chat_handoff.md` (post-Session-4 version — needs refresh after Session 4 closure to reflect 8/8 pipeline complete + Session 4.5 as next).
2. Paste this Session 4.5 prompt document.
3. Add opening line:

```
I'm continuing Trinetravir from Session 4 close. Pipeline 8/8 done.
Ready to launch Session 4.5: audit-response cleanup + manuscript
propagation + deferred Part B global scVI.

Let's start with Part A (Randolph Issue 31 canonical N=1000 recompute).
```

The untracked Part B infra files from Session 4 (`scripts/extract_global_counts.py`, `scripts/session4_part_b_global_sweep.py`, `scripts/session4_part_b_modal.py`) should be visible in `git status` when fresh chat starts — they get modified and lifted into tracked state during Part F.2.
