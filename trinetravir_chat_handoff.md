# Trinetravir chat handoff — post-Session 3.5 state

**Updated:** 2026-05-11 (post-Session-3.5-closure, post-Session-4-spec-commit).
**Last commit at handoff write:** `d3b2a6b` (`references/session_4_prompt.md`).
**Purpose:** drop into a fresh Claude Code / Claude.ai chat to bring the new session up to current state without re-litigating settled work.

---

## Pipeline state (8-block, rev 3, 2026-05-11)

| # | Block | Status |
|---|-------|--------|
| 1 | Session 3 closeout | DONE |
| 2 | Session 5 (calibration audit response) | DONE |
| 3 | Issue 25 human decision (Option B hybrid) | DONE |
| 4 | Session 6A (4 held-out cohort acquisition + harmonization + Schema v6) | DONE |
| 5 | Session 6B (held-out calibration, Issues 27-31 verdicts at N=1000) | DONE (Parts A/B/C/F/G; Parts D/E deferred to Phase 5+) |
| 6 | Session 7 (pre-modeling sensitivity audit, Issues 32+33) | DONE |
| 7 | Session 3.5 (pre-specs Issues 18-24, gene-list canonicalization, PLAN.md v1.3) | DONE |
| 8 | **Session 4 (Issue 6 scVI sensitivity + GPU + Issue 23 compute-envelope)** | **NEXT** — spec at `references/session_4_prompt.md` |

After Session 4 closes: 0 open methodological issues; Phase 4 modeling implementation unblocked.

---

## Issue counts (METHODS_CHOICES.md, 33 total)

| Status | Count | Issues |
|--------|-------|--------|
| Resolved (full validation) | 25 | 1–5, 7–17, 25–33 |
| Resolved at pre-specification level (Phase 5/7/9 validation pending) | 7 | 18–24 |
| Open (compute-gated, Session 4 deliverable) | 1 | 6 |
| **Total** | **33** | |

Zero `<fill in>` / `pending Session N` / `<pending human decision>` markers in METHODS_CHOICES.md (hygiene pass `7f75218`).

**Issue 6 is the only methodological open issue.** Closes at Session 4 (scVI sensitivity vs Harmony).

Issues 34 + 35 will be opened in Session 4:
- **Issue 34** (NEW): scVI comparison design pre-spec (hyperparameter grid, output normalization, four-tier verdict thresholds).
- **Issue 35** (NEW): foundation model compute-envelope decision (INFERENCE-only scope; fine-tuning deferred to v1.5).

---

## Key file pointers

| File | Purpose | State |
|------|---------|-------|
| `SHORT_TERM_PLAN.md` | Pipeline enforcement (read first every session) | rev 3 with Block #8 as NEXT |
| `PLAN.md` | Project plan | **v1.3** (§1.1 expanded; §1.5 v1.5 pointer; §1.6 related work; §1.7 factorized model motivation-tier; §1.8 exploratory-vs-confirmatory) |
| `METHODS_CHOICES.md` | All methodological decisions with five-field structure | 33 issues, all `<fill in>` / `pending` markers cleared |
| `MANUSCRIPT_DRAFT.md` | Working manuscript document | Sections 4 + 5 updated with N=1000 verdicts; Limitations includes FDR<0.05 disclosure + Harmony Δr disclosure + Yoshida CI caveat + Issue 30 CI mirror language |
| `v1_5_PLAN.MD` | v1.5 forward planning (Ad-vector extension) | Pointed-to from PLAN.md §1.5 |
| `references/session_4_prompt.md` | Session 4 spec | NEW (commit `d3b2a6b`) |
| `references/session_7_prompt.md` | Session 7 spec | Closed |
| `references/khatri_mvs_gene_list.csv` | **CANONICAL** Khatri MVS 86-gene list with provenance header | Commit `c44d1ff`. Phase 5+ model code MUST load from this single source. |
| `references/notes/heldout_khatri_mvs_*.md` | Per-cohort MVS write-up notes | 5 cohort notes from Session 6B Part C |
| `src/trinetravir/baselines/*.py` | 6 baseline stubs (Issue 24 Category 1+2) | All raise NotImplementedError; Phase 5 implementation |
| `src/trinetravir/eval/calibration.py` | Calibration framework v2 | 8/8 unit tests passing; used by Sessions 5+6B+7 |
| `data/processed/harmony_per_celltype_<bucket>.h5ad` | Per-bucket Harmony embeddings | 5 buckets. **NO raw counts** — must extract from per-study reannotated h5ads in Session 4 Part A.1. |
| `data/processed/<study_id>_reannotated.h5ad` | Per-study raw counts | 4 v1 studies + 4 held-out cohorts processed_v6. Used as raw-counts source. |

---

## Empirical anchors (cite from this list rather than re-running)

**Session 5 (calibration audit response)**:
- Calibration framework v2 fixes: bootstrap CI direction, FDR-BH, observed-r CI lower bound.
- Khatri MVS Part B: ISG-restricted r lift +0.06 to +0.23 across 4 of 5 v1 buckets.

**Session 6B (held-out cohort calibration at N=1000)**:
- Issue 27 PRIMARY (Randolph monocyte_infected, post-Issue-31): r_full = 0.129 (p=0.001, FDR p_full = 0.027 — only test crossing FDR<0.05 across 15-test panel; signal in non-ISG genes); r_MVS = −0.011 → **CHALLENGES_H1**. Bystander sensitivity row same direction.
- Issue 28 (Yoshida monocyte cross-age): r_MVS = 0.591, CI [0.017, 0.684] → **SUPPORTS_H1**. CI lower bound below 0.30 threshold disclosed honestly.
- Issue 29 (Allen Atlas chronic-latent CMV): r_MVS = −0.010 → CONCERNING reframed as scope-limitation finding.
- Issue 30 (GSE157829 chronic HIV CD4T): r_MVS = 0.257, CI [0.157, 0.513] → **BORDERLINE** (CI straddles EXPECTED + SURPRISING_HIGH bands).

**Session 7 (pre-modeling sensitivity audit)**:
- Issue 32 (pre/post-Harmony Δr): MIXED both gene sets. Δr range 0.02–0.25 across 5 buckets. **Monocyte MVS Δr = 0.08 → BIOLOGY_DOMINANT** at the load-bearing grain. NO bucket crosses Δr > 0.30 HARMONY_DOMINANT threshold.
- Issue 33 (within-cohort sensitivity): BIOLOGY_CONSISTENT both gene sets. **Sign concordance = 100% across 20 aggregate bucket-pair × gene_set tests**. Cross-study integration amplifies, not creates, the signal.

**Session 3.5 (pre-specs Issues 18-24 + canonicalization)**:
- Issue 18: Khatri MVS Andres-Terre 2015 Table S2 high-confidence core subset (86 genes) as primary. Canonical at `references/khatri_mvs_gene_list.csv`.
- Issue 19: REACTOME R-HSA-913531 pathway for pathway-aware regularization.
- Issue 20: MSE on response vectors primary; NB sensitivity at Phase 5.
- Issue 21: 16-config hyperparameter grid (within Issue 14 20-config budget).
- Issue 22: Few-shot adaptation protocol (Phase 9 deliverable; absorbed deferred Session 6B Part E work).
- Issue 23: Method versions pinned at Phase 7 launch. Foundation model baselines (Geneformer, scGPT) gated on Session 4 compute envelope (Issue 35).
- Issue 24: 6 baselines (3 trivial Category 1 + 3 factorization Category 2). Bar-to-beat: factorized model must exceed Category 2 by Δr ≥ 0.05 cross-study Pearson averaged across buckets.

---

## ISG-conservation framework defense (post-Session-7 framing)

The cross-study coherence finding (lymphoid +0.06 to +0.23 ISG lift across cohorts; held-out replication across Yoshida + GSE157829) is **biology with Harmony amplification, not an integration artifact**.

Defense chain:
1. **Pre-Harmony r already substantial** (0.13–0.58 across 5 buckets, Session 7 Issue 32). The signal exists before Harmony.
2. **Harmony adds Δr 0.02–0.25** (mixed, never > 0.30). Harmony amplifies; doesn't manufacture.
3. **Monocyte MVS Δr = 0.08 = BIOLOGY_DOMINANT** at the load-bearing grain. Where ISG-conservation is most claimed (monocyte canonical-ISG response), Harmony's contribution is below the 0.10 mixed-band ceiling.
4. **Within-cohort sign concordance 100%** (Session 7 Issue 33). Every bucket pair in every cohort matches the sign of the cross-study harmonized result. Removing integration entirely doesn't collapse the pattern.

These four points constitute the manuscript's response to critique concern 4 (Harmony preserving only conserved axes). Don't re-litigate this in new chats.

---

## Manuscript framing (post-Session-7)

**Primary contribution**: calibration framework + methodology audit trail (pre-registered held-out validation; sensitivity audit; documented field-gap responses).

**Secondary contribution**: factorized model architecture (Phase 5/6 empirical question whether it adds value over Issue 24 Category 2 baselines).

Working title (top choice): "Calibrated cross-study integration of PBMC SARS-CoV-2 single-cell transcriptomics reveals ISG-anchored response coherence and cross-context generalization".

Title finalizes after Session 4 closes Issue 6.

---

## What NOT to re-litigate in new chats

- **Bucket granularity**: 5-bucket primary + sub-bucket sensitivity (Issue 2 resolved).
- **Cross-study metric**: mean off-diagonal Pearson r primary + Spearman/DE-Jaccard/MMD supplementary (Issue 3 resolved; DE-Jaccard degeneracy documented per Issue 3 revision).
- **Per-cell-type vs global Harmony**: per-cell-type primary, global supplementary (Issue 7 resolved + revision).
- **Khatri MVS gene set**: 86-gene Table S2 core from Andres-Terre 2015. Canonical CSV at `references/khatri_mvs_gene_list.csv`. Phase 5+ MUST load from this single source.
- **Issue 25 v1 framing**: Option B hybrid resolved 2026-05-11. 4 held-out cohorts acquired for cross-context validation; original v1 corpus preserved as training.
- **Issue 27 bystander caveat**: closed via Zenodo re-acquisition of `infected_monocytes_cluster_singlets.rds` + Issue 31 cross-bucket healthy reference. Both bystander + infected verdicts CHALLENGES_H1. Documented as kinetic boundary, not data gap.
- **FDR<0.05 panel-level significance**: structural power limit at N=1000 × 15 tests. Disclosed in Limitations. Verdicts apply pre-committed effect-size rules, not FDR p-values.
- **Yoshida CI [0.017, 0.684]**: disclosed honestly in Section 4 + Limitations. Mechanical SUPPORTS_H1 verdict holds at point estimate; CI lower below 0.30 threshold acknowledged.
- **Session 7 audit**: ISG-conservation defended. See above section.

---

## Session 4 launch context

**Trigger:** "Launch Session 4 per `references/session_4_prompt.md`. Open Issues 34 and 35 with the pre-committed decision rules verbatim before running any compute. Then write `scripts/extract_per_bucket_counts.py` and run extraction. Then begin Part A sweep."

**Pre-conditions verified at Session 4 entry:**
- ✅ Session 3.5 closed (commit `ec3dfe9` + hygiene pass `7f75218`)
- ✅ Session 4 spec committed (`d3b2a6b`)
- ⚠ GPU environment NOT yet provisioned (user TODO before launching Part A)
- ⚠ scvi-tools NOT yet pinned in `configs/methods_versions.yaml` (Session 4 atomic commit #4)

**Atomic commit sequence (10 expected in Session 4):**
1. ✅ `references/session_4_prompt.md` (DONE, `d3b2a6b`)
2. Issue 34 pre-spec
3. Issue 35 pre-spec
4. `configs/methods_versions.yaml` + `configs/gpu_environment.yaml`
5. `scripts/extract_per_bucket_counts.py` + 5 `scvi_input_<bucket>.h5ad`
6. Part A scVI per-bucket sweep results
7. Part B scVI global supplementary results
8. Issue 6 resolution
9. Part C compute-envelope + Issue 35 resolution
10. Session 4 pipeline closure

**Estimated wall-time:** ~7-14h GPU + 4-6h focused work, 2-3 chat sessions over 1-2 calendar days.

**After Session 4 closes:** 0 open methodological issues; Phase 4 modeling unblocked. Manuscript title locks. Phase 5 hyperparameter sweep begins.

---

## Memory pointer

Auto-memory at `/Users/sakshamgupta/.claude/projects/-Users-sakshamgupta-Documents-coding-projects-trinetravir/memory/MEMORY.md` indexes:
- `short_term_plan_pointer.md` — pipeline rev 3 at 7/8 (Session 4 NEXT)
- Other entries: user role, naming, schema decisions, Phase 2/3 scope decisions

Auto-memory is loaded on every Claude Code session via the memory system. Fresh chats inherit this automatically.

---

## Quick orientation for fresh chat

1. Read `SHORT_TERM_PLAN.md` (auto-loaded). Confirm Block #8 NEXT.
2. Read `references/session_4_prompt.md`. Note 10-commit atomic sequence.
3. Read this handoff doc for state snapshot + pointers.
4. Issue 34 + 35 pre-spec drafting first. No compute until pre-specs commit.
5. After pre-specs: GPU env setup, scvi-tools pin, extraction script, Part A sweep.
