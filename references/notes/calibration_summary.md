---
name: Calibration literature summary
description: Supplementary table mapping batch-correction and meta-analysis literature values to Trinetravir's cross-study response-vector Pearson r metric.
type: reference
---

# Calibration literature summary

Our metric: per-cell-type response vector = mean(diseased) - mean(healthy) in scaled-HVG gene space, computed per study; pairwise Pearson r across studies; mean off-diagonal r is the gate. Per-bucket thresholds: 0.60 monocytes, 0.25 CD8T.

| Paper | Metric | PBMC value | Mapping (Direct / Indirect / Context) |
|---|---|---|---|
| Luecken et al. 2022 *Nat Methods* (scIB) | Composite of 14 embedding metrics (ASW, ARI, NMI, kBET, iLISI, graph connectivity, etc.) | Harmony among top-4 on Immune Cell Hum (10 PBMC+BM batches); exact decimals only as figure panels | Indirect — embedding mixing, not signature replication |
| Tran et al. 2020 *Genome Biol* | kBET, iLISI, cLISI, ASW, ARI on 10 datasets | Harmony tied for first overall on PBMC Dataset 5 (with Seurat 3); ARI_batch > 0.97; top-tier cLISI; 3rd on iLISI and kBET | Indirect — method-choice justification |
| Korsunsky et al. 2019 *Nat Methods* (Harmony original) | iLISI / cLISI on the post-correction embedding | 3 PBMC chemistries: iLISI 1.00 -> 1.96 (95% [1.36, 2.56]); cLISI 1.00 retained; competitors fail to exceed iLISI 1.1 | Indirect — algorithmic citation, embedding-mixing precedent |
| Hao et al. 2021 *Cell* (Seurat v4 / Azimuth) | Reference-mapping prediction score; cross-platform fraction concordance | MAIT fraction R = 0.911 vs CyTOF on COVID-19 PBMC; ~80% cells at confidence >= 0.75 (community report) | Context — upstream cell-type-assignment reproducibility, 1-dim fraction (not response vector); soft ceiling hint |
| Zheng et al. 2021 *Immunity* (Khatri lab MVS) | Per-cohort Pearson R of MVS score vs severity; per-cell-type MVS-correlation at single-cell | Bulk discovery: 0.43 <= R <= 0.93 per cohort; validation: 0.33 <= R <= 0.78; aggregate R = 0.75. scRNA: CD14+ monocyte R = 0.45 (p = 2.7e-14), CD16+ monocyte R = 0.25, myeloid R = 0.28. Severity AUROC >= 0.929 discovery, > 0.98 validation | **Direct** — load-bearing citation; closest published cross-cohort PBMC scRNA response-replication numbers |

## Suggested calibration table interpretation
- Trinetravir's monocyte gate (0.60) sits **above** Khatri's CD14+ monocyte single-cell R (0.45). Either our gate is more conservative (good, since we measure pairwise study agreement directly while Khatri measures score-vs-severity), or we should soften toward 0.45-0.50 if data fails to clear 0.60.
- Trinetravir's CD8T gate (0.25) matches Khatri's CD16+ monocyte R (0.25) and generic myeloid R (0.28) — calibration on rare/heterogeneous cell types is approximately field-consistent.
- scIB / Tran / Korsunsky values cannot be used to derive a numerical r threshold; cite them only as evidence that Harmony is a defensible PBMC integration choice.
- Seurat v4's R = 0.911 on a 1-dim cell-fraction metric is a soft *ceiling* — our 2000-dim response-vector r should not exceed this, and on the most abundant cell types is expected to plateau well below it.
