# Short-term plan (pipeline enforcement, read FIRST every session)

**Read this file at the start of every Claude Code session before any other work.** This file defines the strict pipeline order for the remaining work blocks. Do NOT start a block until its prerequisites are checked off.

## Pipeline order (LOCKED — no reordering without human approval)

| # | Block | Spec | Prereqs | Status |
|---|-------|------|---------|--------|
| 1 | **Session 3 closeout** | Inline (this file §1) | Phase35_subbucket bg job `bfs5zrjht` completion | **DONE 2026-05-11** |
| 2 | **Session 5 (audit response)** | `SESSION_5_SPEC.md` | Block #1 fully committed | **DONE 2026-05-11** |
| 3 | **Human decision on Issue 25** | METHODS_CHOICES.md Issue 25 entry | Block #2 complete + Issue 25 opened | **DONE 2026-05-11** (Option B hybrid) |
| 4 | **Session 6A (held-out cohort acquisition + harmonization + schema v6)** | `SESSION_6A_CHECKLIST.md` | Block #3 (Issue 25 Option B hybrid) | **DONE 2026-05-11** (Parts A/B/C/D/E all complete; 4/4 cohorts harmonized; audit 14/0 PASS) |
| 5 | **Session 6B (held-out calibration + per-stratum sensitivity + few-shot)** | inline spec (chat) | Block #4 audit gate passes | **DONE 2026-05-11 (Parts A/B/C/F/G; Parts D/E deferred to Phase 5+ — require v1 factorized model)** |
| 6 | **Session 7 (pre-modeling sensitivity audit)** | `references/session_7_prompt.md` | Block #5 closure + Session 7 prompt committed | **DONE 2026-05-11** — Issue 32 MIXED (Δr 0.02-0.25; monocyte MVS BIOLOGY_DOMINANT) + Issue 33 BIOLOGY_CONSISTENT (100% within-cohort sign concordance); ISG-conservation framework empirically defended; human audit-confirmed |
| 7 | **Session 3.5 (pre-specs)** | `SESSION_3_5_CHECKLIST.md` | Block #6 audit-confirm | **DONE 2026-05-11** — Issues 18-24 pre-specified; Khatri MVS gene-list canonicalized as `references/khatri_mvs_gene_list.csv`; baseline stubs at `src/trinetravir/baselines/`; PLAN.md v1.3 with §1.1 expansion + §1.5 + §1.6 + §1.7 |
| 8 | **Session 4 (GPU/scVI)** | Not yet drafted | Block #7 complete | **NEXT** — Issue 6 scVI sensitivity analysis + GPU environment setup |

**Block #6 (Session 7) is the new gate** between Session 6B closure and Session 3.5. Triggered by critique-document concern 4 (Harmony preserving only conserved axes). Pattern mirrors Session 5 audit response: pre-committed decision rules before computation; results disclosed regardless of outcome.

**Block #8 is the LAST block**, not the first. Session 4 (GPU/scVI work) does NOT run earlier in the pipeline.

## Hard rules

- **NEVER skip ahead.** If a block is BLOCKED, do not start it even if it looks easier. The dependency exists for a reason (Session 5 fixes calibration framework bugs before Session 3.5 makes pre-spec decisions; Issue 25 framing decides what Session 3.5 even pre-specifies).
- **Atomic commits per block.** Each block's work commits as one (or a small number of clearly-themed) commit(s). Do NOT mix work from different blocks in a single commit.
- **Update this file before starting any block.** Mark the previous block as DONE, mark the new block as IN PROGRESS, and check the prereqs are met. If prereqs are missing, list what's missing and STOP.
- **Update this file at end of every session** with current status of in-progress block. If a session ends mid-block, note exactly what's done and what remains.
- **NEVER delete this file.** When all blocks complete, replace contents with new plan; never remove the file itself.

## Current status (2026-05-11)

**Pipeline rev 2 (2026-05-11)**: Issue 25 resolved as **Option B hybrid** — acquire 4 held-out cohorts (Randolph 2021 IAV, GSE283744 pediatric RSV+SARS, Wang 2025 CMV, Lee 2025 HIV) for cross-context validation, NOT to retrain v1 on. Original Wilk+Lee+Arunachalam+Schulte training corpus intact. New blocks #4 (Session 6A: acquisition + harmonization + schema v6) and #5 (Session 6B: held-out calibration) inserted before Session 3.5 + Session 4. Sessions 6A/6B span 4-6 weeks of work, multi-Claude-session.

- **Block #1 (Session 3 closeout)**: **DONE 2026-05-11**
  - Phase35_subbucket calibration bg job `bfs5zrjht` completed (~80 min CPU; 12 sub-buckets × 3 metrics + 12 split-half + 12 MMD).
  - Issue 2 resolution written citing within-Immune_All_Low granularity sweep. Headline: 5-bucket primary; sub-bucket sensitivity reveals additional B-cell signal (B_naive + B_memory both PASS calibrated where 5-bucket B FAILS).
  - Committed: calibration_phase35_subbucket.csv + Issue 2 resolution + SHORT_TERM_PLAN.md + SESSION_5_SPEC.md + SESSION_3_5_CHECKLIST.md (commit hash to be filled).
  - State at end of Block #1: 16 issues resolved (1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17); 1 open (Issue 6 for Session 4).

- **Block #2 (Session 5, audit response)**: BLOCKED on #1
  - Spec: `SESSION_5_SPEC.md`
  - Tasks: #15-24 already created
  - Key gate inside Block #2: if Khatri MVS external validation (Part B) FAILS decisively (r<0.3 or unexpected), STOP and report before proceeding to Part C+. Failure may indicate calibration framework problems beyond audit findings.
  - End-state: Issue 25 OPEN (v1 paper framing decision), Issue 26 RESOLVED (Phase 3 threshold provenance acknowledgment), Issues 3+7 revised for post-hoc rationale, calibration framework v2 produced.

- **Block #3 (Issue 25 human decision)**: BLOCKED on #2
  - Human reviews Issue 25's two options (A: reframe v1 honestly as SARS-CoV-2 cross-study + Lee within-study exploration; B: acquire additional viral data before v1 ships).
  - No Claude Code work in this block — wait for human.
  - If Option B: 2-4 weeks additional data acquisition + harmonization before Block #4 begins.

- **Block #4 (Session 3.5, pre-specs)**: BLOCKED on #3
  - Spec: `SESSION_3_5_CHECKLIST.md`
  - Pre-specifies Issues 18-24 (Phase 5/7/9 prerequisites) + METHODS_CHOICES hygiene cleanup + PLAN.md v1.3 integration + stub baseline files at `src/trinetravir/baselines/`.
  - End-state: 7 issues at pre-specification level (18-24). Compute begins in Block #5.

- **Block #5 (Session 4, GPU/scVI)**: BLOCKED on #4
  - Spec: not yet drafted. Likely covers scVI sensitivity analysis (Issue 6), GPU environment setup, possibly initial Phase 4 work.
  - End-state: Issue 6 resolved. Project ready for Phase 4 implementation.

- **Block #6 (Session 7, pre-modeling sensitivity audit)**: **DONE 2026-05-11**
  - Spec: `references/session_7_prompt.md` (committed 2026-05-11)
  - Atomic commits: `9d363d8` (Issues 32+33 pre-spec gate) → `6abc7bd` (Part A) → `1382811` (Part B) → `5796130` (mechanical rules + manuscript update).
  - Issue 32 (pre/post-Harmony Δr): MIXED both gene sets. NO bucket × gene_set crosses Δr > 0.30 HARMONY_DOMINANT threshold. Monocyte MVS Δr=0.08 → BIOLOGY_DOMINANT at the load-bearing grain.
  - Issue 33 (within-cohort sensitivity): BIOLOGY_CONSISTENT both gene sets. Sign concordance = 100% across 20 bucket-pair × gene_set aggregate tests; mean magnitude alignment 0.077 (full HVG) / 0.136 (MVS).
  - Conditional atomic commit #5 (substantive reframing) NOT triggered.
  - Human audit-confirmed: verdicts match pre-committed rules; manuscript updates honest disclosure; no post-hoc rationalization in commit sequence.
  - ISG-conservation framework empirically defended against critique-document concern 4 (Harmony preserving only conserved axes). Biology with Harmony amplification, NOT integration artifact.

- **Block #7 (Session 3.5, pre-specs)**: **DONE 2026-05-11**
  - Spec: `SESSION_3_5_CHECKLIST.md` + session-prompt-as-pasted (decision-documentation work).
  - Atomic commits (10 total this session):
    - `d86585c` Issue 18 (ISG gene set: Khatri MVS primary)
    - `f73dfda` Issue 19 (pathway: REACTOME R-HSA-913531)
    - `932b984` Issue 20 (reconstruction loss: MSE primary, NB sensitivity)
    - `fed4d19` Issue 21 (factorized model hyperparameter search space)
    - `af95d94` Issue 22 (few-shot adaptation protocol)
    - `2c0373d` Issue 23 (comparison method versions + reproducibility)
    - `eb8564f` Issue 24 + 6 baseline stub files
    - `1e07453` METHODS_CHOICES.md hygiene cleanup (Issues 8-11 `<fill in>`)
    - `c44d1ff` Khatri MVS gene-list canonical CSV + Issue 18 amendment
    - `332f023` PLAN.md v1.3 (§1.1 expansion + §1.5 + §1.6 + §1.7)
  - State at close:
    - 25 issues resolved at entry (1-5, 7-17, 25-33)
    - 7 issues resolved at pre-specification level (18-24)
    - 1 open: Issue 6 (Session 4 GPU work)
    - 33 issues total in METHODS_CHOICES.md
  - All Phase 5/7/9 modeling decisions are now pre-registered against the Session 5+6B+7 audit trail. Implementation in Phases 5+ cannot drift from these pre-specs without an atomic Issue 17 audit-trail change.

- **Block #8 (Session 4, GPU/scVI)**: NEXT, after Session 3.5 close
  - Spec: not yet drafted.
  - Scope: Issue 6 (scVI sensitivity analysis vs Harmony), GPU environment setup, possibly initial Phase 4+ implementation work.
  - Compute envelope: foundation model baselines (Issue 23: Geneformer, scGPT) gated on Session 4 GPU access.
  - End-state: Issue 6 resolved. Pipeline cleared for Phase 4+ modeling implementation.

## Enforcement protocol

At the start of every Claude Code session that touches this project:

1. **Read this file** (the auto-memory system makes this happen automatically via the entry pointing here).
2. **Identify the current IN PROGRESS block** from the table above.
3. **Check that block's prereqs are met.** If any prereq is unchecked, STOP and either complete the prereq or escalate to the user.
4. **Work only within the in-progress block's scope.** Do not start a future block.
5. **Update this file at end of session** with current state of in-progress block before terminating.

If a user prompt during a session requests work from a future block, respond: "Block #X is BLOCKED on Block #Y completion. Current in-progress is Block #Z. Please confirm whether to deviate from the pipeline order or finish Block #Z first."

## Why this file exists

The pipeline now has 8 blocks (rev 3, 2026-05-11), each with non-trivial dependencies. Skipping a block (e.g., starting Session 3.5 before Session 7's harmonization-bias audit) would commit pre-specifications against a known-uncertain harmonization signal contribution, locking in errors that Session 7 was designed to catch. The user explicitly requires pipeline order be followed "to the dot." This file is the enforcement mechanism.

Pipeline rev 3 (2026-05-11): Session 7 inserted between Session 6B and Session 3.5. Triggered by critique-document concern 4 (Harmony preserving only conserved axes). Pattern matches Session 5 audit-response: pre-committed decision rules before run; results disclosed regardless of outcome; audit gate before unlocking next block. See `references/session_7_prompt.md` for full spec.
