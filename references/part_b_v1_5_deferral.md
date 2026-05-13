---
title: Part B (scVI global mode) — v1.5 deferral
status: Deferred to v1.5
date: 2026-05-13
cross-references:
  - Issue 6 (METHODS_CHOICES.md)
  - Block #8 + Block #9 (SHORT_TERM_PLAN.md)
  - Session 4 spec (references/session_4_prompt.md §B)
  - Session 4.5 spec (references/session_4_5_prompt.md Part B / Part F)
---

# Part B (scVI global mode) — v1.5 deferral

## Decision

Part B (scVI global supplementary sweep against `harmony_global` reference) is **deferred to v1 → v1.5 follow-up**, not part of v1 release. Task #24 closed against this rationale.

## Rationale

**Issue 6 verdict does not depend on Part B.** Session 4 Part A delivered the per-bucket scVI sweep, Δr_mvs computation across 5 buckets, and calibration framework v2 application. Under the two-sided structural rule (Issue 34 Amendment 1), 4/5 buckets returned Δr_mvs < −0.10 with CD8T marginal at −0.046, yielding **Tier IV HARMONY_PREFERRED** as the Issue 6 resolution. This verdict is mechanically locked at Session 4 close.

Part B was scoped as a **supplementary reviewer-concern hedge**: anticipating reviewer questions of the form "did you check scVI in global integration mode, not just per-bucket?" — a different parameter setting of the same comparison, not a different scientific question. Its outcome cannot overturn the Part A verdict because Issue 6 was specified per-bucket per the pre-registered comparison; global mode is a robustness check, not a re-arbitration.

**Resource economics rejected Part B at Session 4.5.** Modal A100 run started 2026-05-13 00:44 EDT, killed at 08:44 EDT 8h timeout mid-config 1, zero outputs landed, total spend ~$22.24. Re-running requires either re-architecting the Modal pipeline to chunk under the timeout ceiling, or moving to a different GPU substrate. Neither is justified ahead of v1 release given that the verdict isn't load-bearing.

## Infrastructure preserved

The following artifacts remain in-tree for v1.5 resumption — no code is removed:

- `scripts/extract_global_counts.py` — builds `scvi_input_global.h5ad` (244,389 cells × 4000 HVG aligned to `harmony_global_embedding.h5ad`). Run once to regenerate input.
- `scripts/build_harmony_global_response_vector_cache.py` — replaces the 16 GB `harmony_global_embedding.h5ad` Modal upload (CLI couldn't handle) with a 3.4 MB precomputed Harmony response-vector parquet cache.
- `scripts/session4_part_b_global_sweep.py` — local scVI sweep driver against global mode.
- `scripts/session4_part_b_modal.py` — Modal GPU wrapper (A100 / 8h timeout window).
- `references/session_4_5_prompt.md` Part B + Part F specs — full pre-registered protocol.

## v1.5 resumption plan

When v1.5 work opens:

1. **Restructure for Modal timeout ceiling.** Chunk the global sweep into per-config or per-shard tasks landing partial outputs ≤6h each, with explicit checkpointing. The 8h ceiling killed cleanly because there were no intermediate writes; v1.5 work must persist after each config.
2. **Compare against fixed Harmony cache.** The precomputed `harmony_global_response_vector_cache.parquet` from `build_harmony_global_response_vector_cache.py` is the reference. No re-running Harmony.
3. **Report as supplementary table only.** Outcome enters MANUSCRIPT_DRAFT.md as a v1.5 supplement, not a primary result. If global-mode scVI dramatically outperforms Harmony in global mode (Δr_mvs > +0.10 across multiple buckets), reconsider Issue 6 in a v1.5 revision letter; otherwise document as robustness check.

## Manuscript treatment in v1

`MANUSCRIPT_DRAFT.md` Methods Harmonization paragraph notes: "scVI global-mode comparison deferred to v1.5; per-bucket scVI vs Harmony comparison (Tier IV HARMONY_PREFERRED) reported in §Results." No further claim about global mode is made.

## Audit-trail timeline

- 2026-05-12 22:59 EDT — `extract_global_counts.py` built input file.
- 2026-05-13 00:37 EDT — `build_harmony_global_response_vector_cache.py` replaced 16 GB upload bottleneck.
- 2026-05-13 00:43 EDT — Modal A100 launched (`session4_part_b_modal.py`).
- 2026-05-13 08:44 EDT — Modal A100 8h timeout fired, killed mid-config 1, zero outputs.
- 2026-05-13 — User decision: defer to v1.5. Task #24 closed.
