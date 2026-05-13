# Session 5.5 prompt — Phase 6 metric pre-spec (cell-level evaluation framework)

**Status:** Drafted 2026-05-13 during Phase 5 closeout. Session 5.5 follows Phase 5 closeout per pipeline rev 4 (Phase 5 verdict landed in commit b32fe37 + protocol §3.1 amendment 67f5875).

This session produces the Phase 6 metric pre-spec — the cell-level evaluation framework Phase 6 methods (scGen, CPA, scVI conditional, foundation-model projections) will be scored against. Phase 5 demonstrated that the aggregate grain in `references/phase5_protocol.md` §3 cannot distinguish baselines from sophisticated methods because `predict_mean ≡ linear_delta` byte-identical on all six aggregate metrics (commit 67f5875 §3.1 amendment). Phase 6 requires metrics computed without the §3 aggregation step.

**Integration choice:** standalone file pattern per Sessions 4 / 4.5 / 7 precedent. Commit as `references/session_5_5_prompt.md`. Session output: `references/phase6_protocol.md`.

---

## Entry conditions

**At Session 5.5 entry (post-Phase 5 close):**

- Phase 5 v6 baseline eval committed (b32fe37): 5 baselines × 20 conditions × 6 metrics, FAISS-optimized KNN, N=1000 calibration framework v2.
- Phase 5 protocol §3.1 amendment committed (67f5875): aggregation-grain identity + scope condition + KNN reading note + Phase 6 metric deferral statement.
- Task #24 Part B v1.5 deferral committed (f3c9cad): Issue 6 verdict mechanically locked at Tier IV HARMONY_PREFERRED.
- Pipeline 8/8 blocks done at Session 4 close (Option Narrow); Session 4.5 audit-response cascade closed; Phase 5 baseline floor established.
- Phase 6 modeling work BLOCKED until this session produces `references/phase6_protocol.md` with explicit metric selection + win conditions.

**Empirical anchors carried in:**

- `predict_mean ≡ linear_delta` byte-identical across 20 conditions × 6 metrics under §3 aggregate grain (Phase 5 v6 commit b32fe37; algebraic identity documented in protocol §3.1).
- KNN baselines differ from `predict_mean` but often score worse on aggregate metrics — evidence that aggregation discards per-cell signal, not that KNN is deficient.
- Phase 5 baseline ceiling under aggregate grain: within-virus Pearson r ≈ 0.998, cross-virus r ≈ 0.957. Any aggregate-grain metric for Phase 6 sits on this same ceiling regardless of method sophistication.

**Field-level empirical anchors:**

- Ahlmann-Eltze et al. 2025 (*Nat Methods* 22:1657–1661, DOI 10.1038/s41592-025-02772-6): five foundation models + two deep-learning methods failed to outperform "no change" / "additive" baselines on perturbation benchmarks. GitHub: https://github.com/const-ae/linear_perturbation_prediction-Paper. Per-cell L1 + pseudobulk computation lives in their analysis notebooks.
- bioRxiv 2025.10.20.683304 (October 2025): argues Ahlmann-Eltze tightness is a metric-calibration artifact; proposes "interpolated duplicate" positive control + "dynamic range fraction" calibration measure.
- Bunne et al. 2023 (*Nat Methods*, CellOT): single-cell perturbation as optimal transport problem; uses Wasserstein-2 distance between predicted and observed cell distributions as primary metric. Distribution-matching grain rather than point-prediction.
- Lotfollahi et al. 2019 (*Nat Methods*, scGen): predictive log-likelihood of held-out cells under latent-space generative model. Method-specific; only applies to methods with explicit likelihoods.
- Lotfollahi et al. 2023 (*Mol Syst Biol*, CPA): factorized embedding evaluation; per-cell reconstruction + held-out perturbation prediction with both L1 and Wasserstein components.

---

## SESSION 5.5 — Phase 6 metric pre-spec

```
CPU only. Decision-pre-specification + literature anchoring + smoke
test. Roughly 3-5 hours focused work. No Phase 6 implementation in
this session.

Prerequisite: Phase 5 closeout (commits b32fe37 + 67f5875 + f3c9cad)
landed, origin/main current.

CONSTRAINT: metric pre-spec only. No Phase 6 method implementation.
No new methodological pre-specs beyond Phase 6 metric framework. Atomic
schema-change rule (Issue 17) applies. Pre-registration discipline:
metric choice + win conditions + decision rules committed BEFORE Phase
6 implementation begins.

═══════════════════════════════════════════════════════════════════

PART A — Primary metric pre-spec: per-cell L1

Anchor: Ahlmann-Eltze et al. 2025 GitHub analysis notebooks.

Define:
  per-cell L1 error_i = ||y_pred_i - y_obs_i||_1 / n_genes
  Aggregate: mean over test cells in (bucket, virus, within/cross) condition.

Rationale (must be stated explicitly in phase6_protocol.md §2):
  - Broadest method coverage: applies to any method producing per-cell
    predictions, including foundation-model projection methods that lack
    explicit likelihood.
  - Field-canonical: directly comparable to Ahlmann-Eltze 2025 leaderboard.
  - Cell-level grain: does not collapse under predict_mean = linear_delta
    aggregation identity (per Phase 5 §3.1 amendment).
  - Simplest interpretation: mean absolute error per gene per cell, units
    are log-normalized counts.

Win condition (Issue 24 discipline):
  - Phase 6 method M beats baseline floor B in condition C iff
    per-cell L1(M, C) < per-cell L1(B, C) at canonical N=1000
    permutation null + bootstrap CI + FDR-BH across all (method × bucket
    × virus × within/cross) tests.
  - "Floor" defined as predict_mean (which ≡ linear_delta at aggregate;
    at cell-level grain they differ — pre-spec must commit to which
    cell-level baseline expression of predict_mean to use as floor).

Output: phase6_protocol.md §2 with metric definition, rationale, win
condition, and explicit anchoring to Ahlmann-Eltze 2025 implementation.

═══════════════════════════════════════════════════════════════════

PART B — Supporting metric pre-spec: Wasserstein-2 / energy distance

Anchor: Bunne et al. 2023 CellOT.

Define:
  W2(P_pred, P_obs) where P_pred = {y_pred_i}, P_obs = {y_obs_i}
  Use either Sinkhorn-regularized W2 (entropic OT, faster) or exact W2 via
  POT library; document choice with rationale.

Rationale:
  - Distribution-matching grain: scores whether the predicted cell ensemble
    matches the observed cell ensemble in shape, not just first moment.
  - Cross-check for per-cell L1: a method can have low L1 (cells close to
    targets) but high W2 (wrong distribution shape), or vice versa.
  - Anchors to CellOT 2023 paradigm, complements Ahlmann-Eltze L1.

Win condition:
  - Same structure as Part A. Method beats baseline iff W2(M, C) <
    W2(B, C) at calibration framework v2 thresholds.

Constraint: Wasserstein computation cost. Pre-spec must commit to either
exact W2 (full pairwise cost matrix) or Sinkhorn-regularized W2 with
explicit regularization parameter ε. Document compute envelope.

Output: phase6_protocol.md §3.

═══════════════════════════════════════════════════════════════════

PART C — Calibration metric pre-spec: dynamic-range fraction

Anchor: bioRxiv 2025.10.20.683304 (October 2025 rebuttal).

Define:
  dynamic_range_fraction = (metric(M) - metric(baseline_uninformative)) /
                           (metric(positive_control) - metric(baseline_uninformative))
  where baseline_uninformative = predict_mean of training cells,
  positive_control = "interpolated duplicate" per bioRxiv rebuttal §2.

Rationale:
  - Addresses Ahlmann-Eltze 2025 tightness as metric-calibration artifact.
  - Complements rather than replaces per-cell L1: same primary metric,
    but normalized so a "competitive baseline" reads as ~0 and a perfect
    predictor reads as ~1.
  - Pre-emptively addresses reviewer concerns about Phase 6 method ranking
    being driven by aggregate-vs-cell-grain choice.

Win condition:
  - Method M demonstrates value over baseline iff dynamic_range_fraction(M)
    > 0.2 (threshold to be set in this session; current placeholder).
  - This metric is reported alongside per-cell L1, not as primary.

Output: phase6_protocol.md §4 with positive-control construction protocol
(interpolated duplicate from bioRxiv rebuttal §2) and threshold rationale.

═══════════════════════════════════════════════════════════════════

PART D — Method-specific supplementary: predictive log-likelihood

Anchor: Lotfollahi 2019 scGen, Lotfollahi 2023 CPA.

Define:
  Per held-out cell, log p(y_obs_i | x_i, perturbation; θ) under model M's
  generative density. Sum over test cells.

Rationale:
  - Only available for methods with explicit likelihood: scGen, CPA,
    scVI conditional, optionally Geneformer/scGPT under masked-token
    pseudo-likelihood reframing.
  - Cannot compare across methods with different output spaces (e.g., scGen
    latent vs scGPT token logits). Report within-method as model-fit
    diagnostic, not cross-method ranking.
  - Explicitly excluded from primary leaderboard.

Constraint: report as supplementary table per method. Do not enter primary
comparison.

Output: phase6_protocol.md §5 with method coverage matrix (which methods
support which metric).

═══════════════════════════════════════════════════════════════════

PART E — Aggregation strategy + win conditions + calibration framework v2

Define:
  Per (bucket, virus, within/cross) condition (20 conditions per Phase 5):
    - Compute per-cell L1, W2, dynamic-range-fraction for each method.
    - N=1000 donor-level permutation null (calibration framework v2).
    - N=1000 cell-level bootstrap CI on observed metric.
    - FDR-BH across all (method × bucket × virus × within/cross × metric)
      tests.
  Report:
    - Primary table: per-cell L1 per (method × condition) with CI + q-value.
    - Supporting table: W2 same shape.
    - Calibration table: dynamic-range-fraction same shape.
    - Method-specific: predictive log-likelihood for methods that support
      it, supplementary appendix.

Win condition aggregation rule (Issue 24 + 36 discipline):
  - Method demonstrates value over baseline floor iff:
    (a) per-cell L1 < baseline floor in ≥4 of 5 buckets within-virus AND
        ≥3 of 5 buckets cross-virus, at FDR-BH q < 0.05.
    OR
    (b) Tier system parallel to Issue 34 four-tier verdict — to be
        committed in this session.

Output: phase6_protocol.md §6.

═══════════════════════════════════════════════════════════════════

PART F — Atomic commit sequence

Commit 1: phase6_protocol.md initial draft with §1-5 (metric pre-specs)
Commit 2: phase6_protocol.md §6 (aggregation + win conditions + tier
          system parallel to Issue 34)
Commit 3: METHODS_CHOICES.md Issue 40 entry — "Phase 6 metric framework"
          referencing phase6_protocol.md as authoritative spec
Commit 4: SHORT_TERM_PLAN.md update marking Session 5.5 DONE, Phase 6
          UNBLOCKED, pipeline rev 5

═══════════════════════════════════════════════════════════════════

PART G — Smoke test before Phase 6 implementation begins

Verify in this order before closing Session 5.5:

1. phase6_protocol.md self-consistent: every metric defined has a win
   condition, a rationale, and a literature anchor.
2. Method coverage matrix (Part D) lists every Phase 6 method planned
   (scGen, CPA, scVI conditional, Geneformer projection, scGPT projection)
   and marks which metrics each supports.
3. Pre-spec discipline: no Phase 6 implementation work begun in this
   session. Only protocol authoring + literature anchoring.
4. Calibration framework v2 application consistent with Issue 38 canonical
   N=1000 standard.
5. Tier system (if adopted) symmetric across method vs baseline directions
   per Amendment 1 lesson from Issue 34.
6. METHODS_CHOICES.md Issue 40 entry references protocol file (not
   inline duplication).

Output: smoke test PASS confirmation in commit message of final Session
5.5 commit (Part F commit 4).
```

---

## Out-of-scope explicitly

- **No Phase 6 method implementation.** scGen / CPA / scVI conditional / Geneformer / scGPT projection code lives in subsequent sessions per Session 5.5 protocol output.
- **No Phase 7 foundation-model probe revisions.** Issue 35 INCLUDE verdict stands; foundation-model evaluation happens after Phase 6 metric framework is fixed.
- **No new methodological issues beyond Issue 40 (Phase 6 metric framework reference).**
- **No manuscript Discussion paragraph for Phase 5 pm ≡ ld finding.** Staged until Phase 6 outcome per closing-checklist instruction from Phase 5 closeout chat (2026-05-13).

## Session 5.5 outputs (committed)

1. `references/phase6_protocol.md` — Phase 6 metric framework, full pre-spec.
2. `METHODS_CHOICES.md` Issue 40 — Phase 6 metric framework reference entry.
3. `SHORT_TERM_PLAN.md` updates — Session 5.5 closure, Phase 6 unblocked, pipeline rev 5.

## Decision rule for entering Phase 6 implementation

Phase 6 work (model implementation) begins only after Session 5.5 commits all four atomic commits per Part F. Pre-registration discipline: implementation cannot precede pre-spec.

---

## References

- Ahlmann-Eltze, C., Huber, W., Anders, S. (2025). Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nat Methods* 22, 1657–1661. DOI: 10.1038/s41592-025-02772-6. GitHub: https://github.com/const-ae/linear_perturbation_prediction-Paper
- bioRxiv 2025.10.20.683304 (October 2025). Deep Learning-Based Genetic Perturbation Models Do Outperform Uninformative Baselines on Well-Calibrated Metrics.
- Bunne, C., et al. (2023). Learning single-cell perturbation responses using neural optimal transport (CellOT). *Nat Methods*.
- Lotfollahi, M., Wolf, F.A., Theis, F.J. (2019). scGen predicts single-cell perturbation responses. *Nat Methods* 16, 715–721.
- Lotfollahi, M., et al. (2023). Predicting cellular responses to complex perturbations in high-throughput screens (CPA). *Mol Syst Biol*.

**Project cross-references:**

- `references/phase5_protocol.md` §3.1 — aggregation-grain identity + Phase 6 metric deferral statement.
- `METHODS_CHOICES.md` Issue 24 — baseline win-condition framework.
- `METHODS_CHOICES.md` Issue 36 — calibration framework v2 (canonical N=1000).
- `METHODS_CHOICES.md` Issue 38 — N=1000 reconciliation.
- `SHORT_TERM_PLAN.md` — pipeline state tracker.
