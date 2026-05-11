# Methods Choices Log

This document tracks every non-trivial methodological choice made in the Trinetravir project. Every entry requires:

1. **The choice** — what was decided
2. **The alternatives considered** — what else could have been done
3. **The scientific rationale** — why this choice was made, with citations where possible
4. **The validation strategy** — how this choice will be defended (sensitivity analysis, calibration, citation, or pre-specified rule)
5. **The status** — open / in-progress / validated / paper-ready

No methodological choice is "arbitrary" or "by intuition." If a justification can't be written, the choice must be revisited until one can be. This document IS the methods section of the eventual paper, written incrementally as decisions are made.

When a choice is added to this log, it must include the decision date and the phase of the project when it was made. When a choice is modified, the prior version is preserved with a strikethrough or change note — this is the audit trail.

---

## Open issues requiring immediate resolution

The following choices were made earlier in the project with insufficient justification. Each must be resolved before proceeding to Phase 4. Resolution means: justification written, alternatives evaluated, validation strategy specified, and the resulting evidence collected.

### Issue 1: `infection_status` label semantics (LOAD-BEARING)

**Status**: open. Resolution required before Phase 4.

**The choice as it stands**: cells in the harmonized AnnData are labeled `infected` or `mock` in the `infection_status` obs column, based on whether they came from a diseased or healthy donor. For PBMC studies, this is a *donor-level* disease state proxy, not a *cell-level* viral infection state. PBMCs rarely contain directly virally infected cells.

**Why this is a problem**: the label name is misleading. A reviewer reading "infected cells" in the methods will reasonably assume per-cell viral RNA detection, which is not what the column represents. The mismatch between label and meaning will be flagged as imprecise at best, misleading at worst.

**Resolution required**:
- Rename the obs column from `infection_status` to `donor_disease_status`.
- Allowed values become `diseased`, `healthy_control`, `mock_control` (the last for in-vitro mock-infected controls if any).
- Update the loader, harmonization, and downstream code to use the new name.
- Document in this file that the project measures *systemic immune response to viral disease in PBMCs*, not *cell-autonomous response to direct viral infection*. The factorized model still works on this signal but the framing of the eventual paper must be honest about it.

**Alternatives considered**:
- Keep `infection_status` but document the proxy clearly: rejected because the label name itself is the source of confusion, and documentation cannot fully fix that.
- Use per-cell viral read detection to define `infected` cells properly: not possible at scale in our data because PBMC viral read counts are extremely sparse and most studies did not align reads to viral genomes.

**Validation strategy**: the rename is a refactoring task, not a scientific claim, so no sensitivity analysis required. The validation is that the methods section explicitly defines what the label means.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 2: cell-type bucket granularity (LOAD-BEARING)

**Status**: open. Resolution required before Phase 6.

**The choice as it stands**: cells are bucketed into five coarse-grained categories — monocyte, B, NK, CD4T, CD8T — for all downstream analyses including the Phase 3 gate check and the eventual cross-virus benchmark.

**Why this is a problem**: the bucket granularity is itself a methodological decision that affects results. Monocytes have meaningful subtype heterogeneity (classical, non-classical, intermediate) with different antiviral responses. B cells span naive, memory, plasma. The choice of "five buckets at this granularity" was driven by what mapped cleanly across the original study labels, not by a pre-specified scientific criterion.

**Resolution required**:
- Pre-specify the criterion for bucket choice in this file before running Phase 6.
- Run a sensitivity analysis at finer granularity for at least the headline cell type (monocyte) using CellTypist's Immune_All_High labels or equivalent, where unified labels permit. Show that the cross-virus result is qualitatively preserved across granularity choices.
- Document the trade-off: finer granularity yields more biological precision but lower cell counts per bucket per donor, increasing response vector noise.

**Alternatives considered**:
- Fine-grained labels throughout (15-20 buckets at the cell subtype level): rejected for v1 because per-bucket per-donor cell counts become too small in some studies to compute stable response vectors. Documented as a v2 extension.
- Bulk PBMC response without bucketing: rejected because the Phase 3 stratified diagnostic showed that bulk cross-study correlation (r=0.054) is dominated by composition drift, not biological signal. Bucketing is necessary to recover meaningful signal.

**Validation strategy**: sensitivity analysis at finer granularity for monocyte in Phase 6 or Phase 7. Supplementary figure showing cross-virus results are robust to granularity choice.

**Pre-specified criterion to be written**: "Bucket granularity is the coarsest level at which (a) all four studies' annotations can be reliably mapped to the same vocabulary, (b) each bucket contains at least N=200 cells per donor per study on average, and (c) prior PBMC integration literature (Khatri MVS, scIB benchmark) reports comparable groupings."

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 3: Pearson r as the primary cross-study coherence metric (LOAD-BEARING)

**Status**: open. Resolution required before Phase 6.

**The choice as it stands**: cross-study response-vector coherence is measured by mean off-diagonal Pearson r across the per-study response vectors. The Phase 3 gate uses this metric. The Phase 4 cross-virus benchmark will likely use it too.

**Why this is a problem**: Pearson r is sensitive to outlier genes, scale-dependent, and does not capture distribution-shape differences. The Feb 2026 metrics-failure literature (Evaluating Single-Cell Perturbation Response Models Is Far from Straightforward, bioRxiv 2026.02.14.705879) explicitly shows that metric choice substantially affects apparent method rankings in this domain. A reviewer will ask why Pearson and not Spearman, MMD, Wasserstein, or differential expression overlap.

**Resolution required**:
- Pre-specify the primary metric and the rationale for using it.
- Run a sensitivity analysis using at least three additional metrics: Spearman r (rank-based, robust to outliers), differential-expression-set overlap (interpretable, captures top-N gene agreement), and MMD with multiple kernel choices (distribution-aware, addresses the Feb 2026 paper's critiques of Wasserstein and Energy distance).
- Report all metrics in the eventual paper. The headline result must be qualitatively consistent across at least three of them.

**Alternatives considered**:
- MMD as primary: rejected as primary because it requires kernel choice, which adds another arbitrary parameter. Use as secondary.
- Wasserstein distance: explicitly avoided per the Feb 2026 paper's finding that Wasserstein fails in high-dimensional gene expression spaces under variance scaling.
- Energy distance: avoided for the same reason; the Feb 2026 paper shows Energy distance can overlook disruptions in gene-gene relationships.

**Validation strategy**: sensitivity analysis with at least three additional metrics; cite the Feb 2026 paper in the methods section as the rationale for choosing Pearson over distribution-distance metrics.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 4: study exclusion criteria for Guo and MGH (MODERATE)

**Status**: open. Resolution required before Phase 4.

**The choice as it stands**: Guo 2020 was excluded because it has 2 diseased donors and 0 healthy controls. MGH was excluded because it has 14 diseased and 1 healthy donor. The exclusions were made on a case-by-case basis after seeing the data.

**Why this is a problem**: the exclusion criterion was not pre-specified. A reviewer could ask whether the exclusions are post-hoc rationalizations or principled.

**Resolution required**:
- Pre-specify the exclusion rule and document it in this file.
- Apply the rule retroactively to confirm Guo and MGH exclusions are consistent with it.

**Pre-specified rule (to be confirmed by the user)**: "A study is included in the harmonization corpus if and only if it has at least 4 healthy donors and at least 4 diseased donors. The threshold of 4 is set by the requirement that within-study donor-level split-half analysis (for the threshold calibration in Issue 9) produces at least 2 donors per split per class."

**Alternatives considered**:
- Lower threshold (≥3 donors): rejected because donor-level statistical analysis becomes unstable.
- Higher threshold (≥6 donors): rejected because it would unnecessarily exclude studies like Wilk that have 6 healthy and 7 diseased donors.
- No threshold, include all studies: rejected because Guo (0 healthy) cannot produce a response vector and MGH (1 healthy donor) produces donor-confounded response vectors.

**Validation strategy**: the criterion is pre-specified before the rest of the analysis; the rule justifies both exclusions cleanly.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 5: Gate 1 sanity-check threshold of r < 0.7 (MINOR but principle-bearing)

**Status**: open. Resolution required for paper.

**The choice as it stands**: at the start of the project, the Gate 1 sanity check threshold was set to "if SARS-vs-IAV response vector r < 0.7, the project has signal; if > 0.9, stop and reconsider." The observed value was 0.46, well below 0.7.

**Why this is a problem**: the 0.7 threshold was hand-picked. The decision (proceed) was unambiguous given the observed 0.46, but the threshold itself is undefended.

**Resolution required**:
- Re-evaluate the Gate 1 result against a calibrated threshold derived from the permutation null distribution (see Issue 9). If observed r exceeds the 99th percentile of the permutation null for cross-virus comparison, the gate is passed on calibrated grounds.
- Document the recalibration in the paper.

**Validation strategy**: subsumed by the permutation null calibration in Issue 9.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 6: Harmony as the batch correction method (MODERATE)

**Status**: open. Resolution required before Phase 6.

**The choice as it stands**: Harmony was used for batch correction with study_id as the batch key, per-cell-type. The choice was driven by CPU-friendliness and speed.

**Why this is a problem**: scVI / scANVI / BBKNN / FastMNN / LIGER are also standard methods in this domain, with different inductive biases. A reviewer will ask why Harmony was chosen.

**Resolution required**:
- Document the principled reasons for choosing Harmony over alternatives.
- Run a sensitivity analysis with at least scVI as an alternative. Show that the cross-study response-vector correlation results are qualitatively preserved.

**Justification (to be refined)**: Harmony was chosen because (a) it runs on CPU within the project's compute constraints, (b) it does not require per-batch reference data unlike scANVI, (c) the Theis lab and Korsunsky et al. 2019 benchmarks show Harmony performs competitively with deep-learning alternatives for PBMC integration tasks, and (d) the per-cell-type harmonization protocol avoids Harmony's known weakness of mixing cell types when used globally with cell-type-aware covariates.

**Alternatives considered**:
- scVI: GPU-required for reasonable wall time, but produces a learned latent space that may be more useful for downstream modeling. Should be run as sensitivity analysis.
- scANVI: requires cell-type labels as semi-supervision; available but introduces another methodological dependency.
- BBKNN, MNN: less commonly used as the primary integration method in recent PBMC work.

**Validation strategy**: scVI sensitivity analysis in Phase 6 or 7; one supplementary figure showing qualitative consistency.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 7: per-cell-type harmonization vs joint harmonization (MODERATE)

**Status**: open. Resolution required before Phase 6.

**The choice as it stands**: Harmony is run separately on each cell-type bucket (monocyte, B, NK, CD4T, CD8T), with study_id as the batch key within each. This was chosen over global Harmony with study_id as the batch key on all cells together.

**Why this is a problem**: the choice was made for clean statistical interpretation per bucket and to avoid the risk of accidentally mixing cell types when using global Harmony. But the choice is itself a methodological decision with implications.

**Resolution required**:
- Document the rationale for per-cell-type over global.
- Run global Harmony as a sensitivity analysis. Confirm that response-vector recovery is qualitatively similar.

**Justification (to be refined)**: per-cell-type harmonization (a) avoids the risk of Harmony mixing cells across cell type boundaries when batch effects are larger than cell-type effects in some regions of the embedding, (b) allows per-bucket Harmony parameters to be tuned to each cell type's batch effect magnitude, and (c) produces cleaner statistical interpretation because response vectors are computed within the same harmonized space they're evaluated in.

**Validation strategy**: global Harmony sensitivity analysis; supplementary figure.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 8: permutation N=1000 and split-half N=50 (MINOR)

**Status**: open. Resolution: one sentence each in methods.

**The choice as it stands**: the calibration prompt specifies N=1000 permutations for the null distribution and N=50 random splits for split-half reliability.

**Why this is acceptable**: 1000 is conventional for permutation testing in genomics (it produces stable p-values down to ~0.001 resolution, sufficient for our 99th percentile threshold). 50 random splits produces stable mean estimates for split-half reliability at our cell counts.

**Resolution**: cite convention in methods. If a reviewer pushes, increase to 10,000 permutations and 100 splits as a sensitivity check.

**Validation strategy**: one sentence in methods citing the convention.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 9: combined pass criterion threshold (LOAD-BEARING)

**Status**: open. Resolution required before Phase 4.

**The choice as it stands**: in the calibration prompt, a bucket passes if (a) observed cross-study r exceeds the 99th percentile of the permutation null AND (b) observed r is at least 50% of the within-study split-half ceiling.

**Why this is a problem**: the 50% threshold for "fraction of ceiling" is itself a hand-picked number. Why not 40% or 60%?

**Resolution required**:
- Replace the fraction-of-ceiling criterion with a more principled test: observed cross-study r should be statistically indistinguishable from within-study split-half r at some confidence level (e.g., observed r within the 95% CI of split-half r), OR formally test whether observed r is significantly lower than split-half r using a bootstrap of both quantities.
- If the principled test is too strict (most buckets fail because cross-study can rarely match within-study), retain the fraction-of-ceiling criterion but justify the fraction with a citation or empirical argument.

**Validation strategy**: replace with bootstrap CI overlap test if practical; otherwise document the fraction choice with rationale.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 10: 99th percentile in permutation null (MINOR)

**Status**: open. Resolution: report both 95th and 99th percentile and pick one.

**The choice as it stands**: the permutation null criterion uses the 99th percentile (equivalent to p < 0.01).

**Why this is acceptable**: 99th percentile is more conservative than 95th and reduces false-positive rate. Multiple testing across 5 buckets and 6 cross-study pairs (10 comparisons per bucket) increases family-wise error rate, which 99th percentile partially accounts for without explicit Bonferroni correction.

**Resolution**: state in methods that 99th percentile was chosen for conservatism given the multiple comparison structure, and report both thresholds in supplementary.

**Validation strategy**: methods sentence + supplementary table.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 11: mean off-diagonal r as the cross-study summary statistic (MINOR)

**Status**: open. Resolution: report mean, median, and minimum, justify mean.

**The choice as it stands**: cross-study r matrices are summarized by their mean off-diagonal value (across all study pairs).

**Why this is a problem**: mean, median, and minimum are all reasonable choices. Minimum is most conservative (worst study pair sets the bound). Median is robust to outliers.

**Resolution**: report all three in supplementary. Use mean in headline figures because it's conventional and matches what most PBMC integration literature reports.

**Validation strategy**: supplementary table.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 12: CellTypist model choice (Immune_All_Low) (MODERATE)

**Status**: open. Resolution required after Phase 3.5 completes.

**The choice as it stands**: Phase 3.5 uses CellTypist's Immune_All_Low pretrained model for unified cell-type re-annotation.

**Why this is a problem**: Immune_All_High has finer granularity, and other PBMC-specific CellTypist models exist. Azimuth is also a defensible alternative.

**Resolution required**:
- Justify Immune_All_Low based on bucket granularity matching (Issue 2) and per-study cell counts.
- Run Immune_All_High as a sensitivity analysis if Phase 3.5 surfaces unified labels that are too coarse to distinguish biologically meaningful subtypes.
- Document the choice and cite the CellTypist paper (Domínguez Conde et al. 2022, Science).

**Validation strategy**: methods justification + optional Immune_All_High sensitivity analysis.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 13: cellxgene Census as the data source (MINOR)

**Status**: open. Resolution: one paragraph in methods.

**The choice as it stands**: scRNA-seq datasets are downloaded from cellxgene Census, which provides processed and harmonized AnnData files with standardized metadata.

**Why this is acceptable**: Census is the field's standard repository for harmonized scRNA-seq data, maintained by CZI Biohub. Using it ensures reproducibility of data ingestion. Reprocessing from raw FASTQs would be possible but is outside the scope of v1.

**Resolution**: document the use of Census, cite the Census paper / resource, acknowledge that downstream results inherit Census's processing choices (specifically alignment to a fixed reference genome version, gene annotation version, and QC defaults applied by the original study authors before submission).

**Validation strategy**: methods paragraph; cite Census as data source.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 14: hyperparameter defaults for benchmark methods (LOAD-BEARING)

**Status**: open. Resolution required before Phase 6.

**The choice as it stands**: to be determined. When scGen, scCausalVI, CoupleVAE, and foundation models are run in Phase 6, hyperparameters must be chosen.

**Why this is a problem**: tuning hyperparameters favors the methods we tune the hardest. Using published defaults treats all methods symmetrically but may not reflect each method's best achievable performance.

**Resolution required**: pre-specify the hyperparameter policy before Phase 6 begins. Choose one of:
- *Published defaults policy*: every method is run at its published default hyperparameters from the original paper or GitHub repo. No tuning. Justification: this protocol is symmetric across methods and reflects out-of-the-box performance, which is the relevant comparison for a new user evaluating which method to adopt.
- *Light tuning policy*: each method gets a small hyperparameter sweep (≤20 configurations) on a held-out validation set. The best configuration per method is used for evaluation. Justification: this gives each method a fair shot at its best performance.
- *Held-out validation policy*: each method tunes hyperparameters using a within-virus held-out validation split, and the tuned hyperparameters are evaluated on the cross-virus test split. Justification: this is the standard ML protocol and most defensible scientifically.

The held-out validation policy is the most defensible but adds compute cost. Decide before Phase 6 begins and document the chosen policy here with rationale.

**Validation strategy**: pre-specified policy; methods section explicitly states the protocol.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

### Issue 15: cross-virus training/test split protocol (LOAD-BEARING)

**Status**: open. Resolution required before Phase 4.

**The choice as it stands**: to be determined. The cross-virus benchmark requires choosing which virus(es) to train on and which to hold out.

**Why this is a problem**: train-on-SARS-test-on-IAV and train-on-IAV-test-on-SARS will produce different results because of asymmetric sample sizes, biological differences, and severity distributions. A single direction is arbitrary; choosing both with equal weight implies a symmetric protocol.

**Resolution required**:
- Pre-specify the cross-virus split protocol. Recommended: leave-one-virus-out cross-validation, where each virus in the benchmark serves as the held-out target in turn, and the model is trained on all other viruses. Report mean and per-virus performance.
- For v1 with two viruses (SARS, IAV), this becomes train-SARS-test-IAV and train-IAV-test-SARS. Report both and the mean.
- When RSV or other viruses are added later, the protocol extends naturally.

**Validation strategy**: pre-specified protocol; methods section describes the leave-one-virus-out scheme.

**Date opened**: 2026-05-10
**Date resolved**: <fill in>

---

## Process rules for future methodological choices

The following process applies to any methodological choice made after this document is initialized.

**Rule 1: every choice gets logged before it's used.** When a new threshold, parameter, exclusion rule, metric, or modeling decision is made, an entry is added to this file with the five required fields (choice, alternatives, rationale, validation strategy, status). The choice does not enter the codebase or the analysis until the entry is written.

**Rule 2: no choice is "intuitive."** If the rationale field reads "based on intuition" or "seems reasonable," the choice is incomplete. Find a citation, a pre-specified rule, an empirical anchor, or a sensitivity analysis. If no rationale can be found, revisit the choice.

**Rule 3: every load-bearing choice gets a validation strategy.** Load-bearing = "if this choice were different, the headline result might change." For load-bearing choices, the validation strategy must include either (a) a sensitivity analysis run as part of the project, or (b) an explicit principled justification with citation, or (c) both.

**Rule 4: the audit trail is preserved.** When a choice is modified, the prior version is kept in the document with a change note explaining the modification. Reviewers will sometimes ask "did you change methodology after seeing results?" and the audit trail answers that question honestly.

**Rule 5: this document IS the methods section.** When the paper is drafted, the methods section is built by lifting entries from this file. If an entry doesn't read like methods-section prose, refine it now, not at submission time.

**Rule 6: read this file before starting any new phase.** Phase transitions are when methodological choices are typically made. Reading this file at the start of each phase ensures that open issues are addressed before new choices accumulate on top of unresolved ones.

---

## Choice log (chronological)

This section records resolved choices in the order they were made. As issues above are resolved, their content is moved here with the resolution date.

*(initially empty; populated as choices are resolved)*

---

## Pending revisions

This section tracks choices that have been resolved but may need revisiting based on later findings.

*(initially empty)*
