# Phase 5 Protocol: Simple Baselines Evaluation

**File:** `references/phase5_protocol.md`
**Status:** Pre-specification for Phase 5 implementation. Per Issue 24 + PLAN.md §4 Phase 5.
**Date:** 2026-05-13

This document specifies the evaluation protocol for Phase 5 Category 1 baselines (`predict_mean`, `linear_delta`, `knn`) per Issue 24. The protocol follows canonical scGen 2019 + Ahlmann-Eltze 2025 evaluation conventions for unpaired-donor single-cell perturbation prediction.

Referenced from PLAN.md §4 Phase 5 as the authoritative specification for Phase 5 implementation.

---

## 1. Evaluation Protocol

Phase 5 follows the canonical cell-level prediction protocol established in Lotfollahi et al. 2019 (*Nat Methods*): per-cell predictions aggregated to per-(cell_type, virus) mean response vectors for evaluation.

Unpaired-donor designs (Lee 2020 PBMC; equivalent to Kang 2018 PBMC IFN-β) do NOT require donor pairing under this protocol — cell-level grain provides baseline feature variance naturally. Random donor pairing is non-standard and adds noise; not used here.

**Reference implementations:**

- `theislab/scgen` GitHub — scGen baseline functional forms
- Ahlmann-Eltze 2025 GitHub — additive / no-change baseline conventions

## 2. Baseline Functional Forms

For test cell $i$ with baseline expression vector $x_i$, target (cell_type, virus) bucket $b$, training cells $\mathcal{T}_b$ partitioned into:

- $\mathcal{T}_b^{baseline}$: mock donor cells in bucket
- $\mathcal{T}_b^{post}$: diseased donor cells in bucket

### 2.1 `predict_mean`

Cohort-level post-perturbation mean, ignoring cell-level input:

$$\hat{y}_i = \frac{1}{|\mathcal{T}_b^{post}|} \sum_{j \in \mathcal{T}_b^{post}} y_j$$

Constant per (cell_type, virus); cell-level grain produces identical per-cell predictions within a group. Equivalent to Ahlmann-Eltze 2025 "additive" baseline at population mean.

### 2.2 `linear_delta`

Adds training-set shift to per-cell baseline:

$$\hat{y}_i = x_i + \delta_b$$

where:

$$\delta_b = \frac{1}{|\mathcal{T}_b^{post}|} \sum_{j \in \mathcal{T}_b^{post}} y_j - \frac{1}{|\mathcal{T}_b^{baseline}|} \sum_{k \in \mathcal{T}_b^{baseline}} x_k$$

Per-cell predictions vary because $x_i$ varies per cell. The shift $\delta_b$ is bucket-level (cell_type × virus), computed once from training cells. Equivalent to Ahlmann-Eltze 2025 "no change" baseline with population-shift adjustment.

### 2.3 `knn` ($K \in \{25, 50, 100\}$)

For each test cell, find $K$ nearest training baseline cells by cosine distance on HVG-restricted gene expression; return the mean post-perturbation expression from their bucket:

$$\hat{y}_i = \frac{1}{|\mathcal{T}_{b,K}^{post}|} \sum_{n \in \mathcal{T}_{b,K}^{post}} y_n$$

where:

$$\mathcal{T}_{b,K}^{post} = \{\text{post cells of donors of } \text{kNN}(x_i, \mathcal{T}_b^{baseline}; K)\}$$

Under unpaired-donor design, "donors of kNN" reduces to a cell-type-bucket-conditional mean weighted by baseline neighborhood structure of $x_i$. Per-cell predictions vary based on local baseline neighborhood.

**Distance choice:** Cosine distance preferred over Euclidean for sparse single-cell data per scGen convention.

## 3. Evaluation Metrics

Per cell-type bucket, aggregate per-cell predictions to per-virus response vector:

$$\bar{\hat{y}}_b^v = \frac{1}{|\text{cells}_b^v|} \sum_i \hat{y}_i$$

Compare to observed aggregate $\bar{y}_b^v$ via:

1. **Per-gene Pearson r** between $\bar{\hat{y}}_b^v$ and $\bar{y}_b^v$ — primary headline metric per scGen convention.
2. **Per-gene R²** with separate computation on Khatri MVS subset vs non-MVS (Issue 18 anchor) — tests whether baseline captures conserved ISG component.
3. **Top-100 DE-gene Jaccard** between predicted and observed gene rankings on $\text{post} - \text{pre}$ direction.
4. **Direction-of-change accuracy** per gene (sign of $\text{post} - \text{pre}$).

## 4. Calibration Framework v2 Application

Apply per Issue 36 framework to all reported metrics:

- **Donor-level permutation null:** ~5,000 permutations of donor-disease-status labels within bucket-virus
- **Bootstrap CI on observed r:** ~5,000 cell-level bootstrap resamples
- **FDR-BH correction:** across (baseline × bucket × virus × within/cross) tests

Output per Issue 24 win-condition: confirmatory verdict per baseline-bucket-virus tuple.

## 5. Expected Results — Honest Framing

The "mean-baseline-is-competitive" phenomenon is a robust field-level finding:

- **Ahlmann-Eltze 2025 (*Nat Methods*):** Five foundation models + two deep learning models all failed to outperform simple "no change" / "additive" baselines on standard perturbation benchmarks.
- **bioRxiv Oct 2025 rebuttal (`2025.10.20.683304`):** Argues the competitiveness is a metric-calibration artifact. Proposes "interpolated duplicate" positive control + "dynamic range fraction" calibration measure.

**For v1:** Report standard Pearson r headline. Flag in Discussion that the metric-calibration debate is open and may affect interpretation of small $\Delta r$ between baselines and downstream methods. This framing pre-empts reviewer concerns; do not over-interpret tight baseline clustering as model-architecture-irrelevant.

## 6. Implementation Tasks

1. Implement baselines per functional forms above:
   - `src/trinetravir/baselines/predict_mean.py`
   - `src/trinetravir/baselines/linear_delta.py`
   - `src/trinetravir/baselines/knn.py`
   - Per-cell prediction interface: `predict(x_test, train_baselines, train_posts) → y_pred`
2. Cell-level evaluation harness in `scripts/phase5_baselines_eval.py`:
   - Per (cell_type, virus, within/cross) condition
   - Predict per-test-cell
   - Aggregate to per-virus mean response vector
   - Compute metrics + calibration framework v2 outputs
3. Output table `results/tables/phase5_baselines_eval.csv` with columns:
   - `bucket, virus, baseline, metric, value, lower_ci, upper_ci, perm_null_p_value, fdr_q_value, within_or_cross`
4. Cross-reference Ahlmann-Eltze 2025 GitHub + scGen GitHub baseline implementations during code review.

All three baselines run on CPU in seconds-to-minutes per bucket. No GPU needed.

## 7. Deliverable

Per-bucket-virus comparison table with calibration-framework-v2-corrected metrics, framed against the mean-baseline competitiveness pattern documented in field literature.

Phase 5 verdict per Issue 24: baseline performance becomes the empirical floor for Phase 6 comparison methods + Phase 8 factorized model.

## 8. References

- Lotfollahi, M., Wolf, F.A., Theis, F.J. (2019). scGen predicts single-cell perturbation responses. *Nat Methods* 16, 715–721.
- Ahlmann-Eltze, C., Huber, W., Anders, S. (2025). Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nat Methods* 22, 1657–1661.
- bioRxiv 2025.10.20.683304 (October 2025). Deep Learning-Based Genetic Perturbation Models Do Outperform Uninformative Baselines on Well-Calibrated Metrics.
- Lotfollahi, M., et al. (2023). Predicting cellular responses to complex perturbations in high-throughput screens (CPA). *Mol Syst Biol*.
- Bunne, C., et al. (2023). Learning single-cell perturbation responses using neural optimal transport (CellOT). *Nat Methods*.

**Project cross-references:**

- `PLAN.md` §4 Phase 5 — phase definition
- `METHODS_CHOICES.md` Issue 18 — Khatri MVS gene set
- `METHODS_CHOICES.md` Issue 24 — baseline win-condition
- `METHODS_CHOICES.md` Issue 36 — calibration framework v2
