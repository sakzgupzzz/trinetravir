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

### Resolved Issue 1: `donor_disease_status` label semantics (LOAD-BEARING) — 2026-05-10

**Final choice**: cells in every harmonized AnnData carry a `donor_disease_status` obs column with allowed values `diseased` and `healthy_control`. `mock_control` is reserved for future in-vitro mock-infected studies (not used in v1; no PBMC study in the v1 corpus produces it). The column was historically named `infection_status` with values `{infected, mock}` in early PLAN drafts; the rename to `donor_disease_status` with values `{diseased, healthy_control}` reflects what the column actually measures: a *donor-level* disease-state proxy derived from the cellxgene `disease` ontology, not a per-cell viral-read detection.

The conceptual rename (`infection_status` -> `donor_disease_status`) was applied in an earlier session (logged in `memory/schema_decisions.md`). The value rename (`healthy` -> `healthy_control`) was applied 2026-05-10 in this session:
- Source updates: `src/trinetravir/data/download.py` (label writer), `src/trinetravir/data/harmonize.py` (label readers in `harmony_per_bucket`), `src/trinetravir/eval/calibration.py` (donor-level masks in permutation null and split-half ceiling), `src/tests/test_download.py` + `src/tests/test_harmonize.py` (assertions), `scripts/phase3_lee_diagnostic.py`.
- On-disk migration: `scripts/migrate_donor_disease_status_value.py` re-writes the value in existing h5ads under `data/raw/` and `data/processed/` (idempotent).
- The Census remote source files are not modified.

**Why this matters scientifically**: the project models *systemic immune response to viral disease in PBMCs* (interferon and cytokine programs in circulating immune cells whose donor has the virus), not *cell-autonomous response to direct viral infection* (a separate phenomenon that requires per-cell viral-RNA detection and is achievable mainly in airway-epithelium studies). The factorized cross-virus model still works on the systemic-response signal but the framing of the eventual paper must be honest about it. The `donor_disease_status` name is precise; the prior `infection_status` name was misleading and would have been flagged by reviewers.

**Alternatives considered and rejected**:
- Keep `infection_status` and document the proxy in methods: rejected because the label *name* is the source of reviewer confusion; documentation alone cannot fix that.
- Define `infected` cells via per-cell viral read detection: not feasible in v1. PBMC viral read counts are extremely sparse and most v1-corpus studies (Lee, Wilk, Arunachalam, Schulte-Schrepping) did not align reads to viral genomes. v2 may revisit this for airway-epithelium studies.

**Validation strategy**: refactor with full test coverage (`uv run pytest src/tests/` = 39 passed at resolution time). No scientific sensitivity analysis required — this is a naming and vocabulary clarification, not a methods change. The methods section of the eventual paper will define both the column and the allowed values explicitly so a reviewer cannot mistake the donor-level proxy for cell-level infection state.

**Date opened**: 2026-05-10
**Date resolved**: 2026-05-10

---

### Resolved Issue 15: cross-virus train/test split protocol (LOAD-BEARING) — 2026-05-10

**Final choice**: cross-virus evaluation uses leave-one-virus-out cross-validation. Each virus in the benchmark serves as the held-out target in turn; the model is trained on all other viruses and evaluated zero-shot on the held-out virus. Both directions plus the mean are reported.

For v1 (two viruses: SARS-CoV-2 and IAV), the protocol reduces to two directions: train on SARS-CoV-2 / test on IAV, and train on IAV / test on SARS-CoV-2. When RSV / 2nd IAV / DNA-virus control studies are added in v1.5+, the protocol extends naturally to N held-out directions.

The protocol is encoded in `configs/evaluation.yaml` under `cross_virus_protocol:` with `protocol: leave_one_virus_out`, `v1_directions:` explicit list, `report_directions: all`, and an explicit rationale block.

**Why leave-one-virus-out rather than a single direction or random holdout**: SARS and IAV have asymmetric sample sizes in the v1 corpus (SARS has more donors and cells). Picking a single direction (e.g. only train-SARS-test-IAV) would inflate headline performance because that is the easier direction; reviewers would (rightly) ask why only one direction was tested. Random holdout of (donor × virus) pairs would let a method pass by memorizing within-virus structure, which is not the cross-virus generalization question.

**Alternatives considered and rejected**:
- Single direction (the higher-N one): rejected for inflating headlines and reviewer-flagging.
- Random 80/20 holdout: rejected because random holdout does not test cross-virus generalization; it tests within-virus interpolation.

**Validation strategy**: protocol pre-specified in `configs/evaluation.yaml` before Phase 4 begins. Methods section quotes the protocol verbatim. Headline results in the paper report each direction separately AND the mean; if the two directions diverge substantially (>0.1 absolute Pearson r difference), that divergence is itself reported as a finding about cross-virus asymmetry, not concealed by averaging.

**Date opened**: 2026-05-10
**Date resolved**: 2026-05-10

---

### Resolved Issue 14: hyperparameter policy for benchmark methods (LOAD-BEARING) — 2026-05-10

**Final choice**: every benchmark method tunes hyperparameters using a within-virus held-out validation split, then evaluates on the cross-virus test split. Each method receives the same compute budget — at most 20 hyperparameter configurations evaluated on the validation set. The configuration that minimizes the primary cross-study coherence metric on within-virus held-out donors is selected; the cross-virus test split then evaluates that single configuration. The validation split is donor-level (not cell-level), seeded, and holds out at least 2 donors per disease class.

The protocol is encoded in `configs/evaluation.yaml` under `hyperparameter_policy:` with `protocol: held_out_validation`, `validation_split_fraction: 0.2`, `tuning_budget_per_method: 20`, and an explicit rationale block.

**Why held-out validation rather than published defaults or light tuning**: published defaults are typically tuned for the original paper's dataset (not ours), which biases against methods whose original-paper data differs most from PBMC cross-virus — that is a *measurement* bias, not a methodological strength. Light tuning without enforced donor-level isolation lets the same donor appear in tuning and test pools, producing optimistically biased results. Held-out validation with donor-level isolation is the standard ML protocol and the most defensible of the three.

**Why a 20-config budget**: a finite, equal budget per method is what makes the comparison symmetric. The number itself is empirical — large enough for typical method-specific hyperparameter spaces (learning rate, latent dim, batch size, regularization weight) to be reasonably explored, small enough to fit within the project's compute envelope when summed over 6-8 methods. If a reviewer pushes, the budget can be doubled in a sensitivity run; the comparison framework is unchanged.

**Alternatives considered and rejected**:
- Published defaults: rejected because of original-paper-data bias described above.
- Light tuning (no enforced donor isolation): rejected because donor leakage produces optimistic bias.

**Validation strategy**: the policy is pre-specified in `configs/evaluation.yaml` before Phase 6 begins. The methods section of the eventual paper will quote the protocol verbatim from the config file. If any single method's reported best configuration is suspiciously close to its published default, that is reported as an integrity check passing; if it diverges substantially, the divergence and the validation-set performance gap is reported.

**Date opened**: 2026-05-10
**Date resolved**: 2026-05-10

---

### Resolved Issue 4: study exclusion criteria (MODERATE) — 2026-05-10

**Final choice**: a study is included in the harmonization corpus if and only if it has at least 4 healthy donors and at least 4 diseased donors. The criterion is encoded in `configs/datasets.yaml` under the top-level `inclusion_criteria:` section with `min_healthy_donors: 4` and `min_diseased_donors: 4`. A retroactive application table records the donor counts and the resulting include/exclude verdict for every study currently in the registry.

Guo 2020 (0 healthy / 2 diseased) and MGH acute COVID (1 healthy / 14 diseased) fail the rule and remain excluded. Lee 2020 (4/13), Wilk 2020 (6/7), Arunachalam 2020 (5/7), and Schulte-Schrepping 2020 (21/18) satisfy the rule and remain included.

**Why this threshold**: donor-level statistical analyses — within-study split-half reliability (Issue 9 calibration) and donor-level permutation null (Issue 9) — require at least 2 donors per split per class for stable estimates. A 4/4 minimum is the smallest donor count compatible with that requirement. A higher threshold (e.g. 6/6) would unnecessarily exclude Lee (anchor cross-virus study) and Wilk. A lower threshold (3/3) would admit Lee but produce split-half estimates with only 1-2 donors per half — statistically unstable.

**Alternatives considered and rejected**:
- Threshold of 3 donors per class: rejected because split-half analysis on 3 donors produces 1-2 donors per half per class, too small for stable reliability estimates.
- Threshold of 6 donors per class: rejected because it would exclude Lee (4 healthy donors) — Lee is the only cross-virus anchor (SARS + IAV + healthy in one study); excluding it would kill v1.
- No threshold (include everything): rejected because Guo (0 healthy) cannot compute a within-study response vector at all, and MGH (1 healthy donor) produces donor-confounded response vectors with strongly negative cross-study r (the Phase 3 prep notebook surfaced this empirically; see memory/phase3_decisions.md).

**Validation strategy**: pre-specified in the registry; enforced by `src/tests/test_datasets.py` (4 tests, all passing). The test fails informatively if a future study is added that violates the rule, so the criterion cannot silently drift. The methods section of the paper will quote the rule verbatim from `inclusion_criteria.rationale`.

**Date opened**: 2026-05-10
**Date resolved**: 2026-05-10

---

## Pending revisions

This section tracks choices that have been resolved but may need revisiting based on later findings.

*(initially empty)*
