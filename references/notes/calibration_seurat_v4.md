---
name: Integrated analysis of multimodal single-cell data (Seurat v4 / Azimuth)
description: Multimodal PBMC CITE-seq reference (~162k cells, 228 antibodies); reports MAIT fraction R=0.911 vs CyTOF on COVID query but does not publish per-cell-type cross-study label-transfer F1.
type: reference
---
**Citation**: Hao Y, Hao S, Andersen-Nissen E, Mauck WM 3rd, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, Hoffman P, Stoeckius M, Papalexi E, Mimitou EP, Jain J, Srivastava A, Stuart T, Fleming LM, Yeung B, Rogers AJ, McElrath JM, Blish CA, Gottardo R, Smibert P, Satija R. 2021. *Cell* 184:3573-3587.e29.
**DOI**: 10.1016/j.cell.2021.04.048
**URL**: https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3 ; preprint: https://www.biorxiv.org/content/10.1101/2020.10.12.335331

## Metric reported
Three things relevant to our calibration:
1. **Reference mapping prediction score** (per-cell, continuous in [0, 1]) returned by `MapQuery` based on weighted nearest-neighbor distances to the multimodal reference. Confidence threshold is typically 0.75.
2. **Cross-modality concordance**: protein modality is held out, then predicted protein expression from RNA-only mapped cells is compared to measured CITE-seq protein. Reported as correlations.
3. **External-cohort validation**: COVID-19 PBMC query mapped to healthy reference, predictions validated against orthogonal CyTOF measurements on the same donors.

## PBMC values observed
- Multimodal reference: ~162,000 PBMCs, 228-antibody CITE-seq panel, eight donors (HIV vaccine trial), integrated via WNN.
- **Cross-modality validation**: MAIT cell fraction predicted from scRNA-seq Azimuth annotation correlated R = 0.911 with MAIT fraction measured independently by CyTOF on COVID-19 PBMC samples.
- **Inter-method agreement vs scArches**: "In 73.8% of cases" Seurat v4 annotations had stronger protein-expression support than scArches when the two disagreed on cell type calls.
- Aggregate cross-study cell-type label-transfer F1 / accuracy is **not** reported as a single headline number in the paper text I could extract. Community reports (Azimuth documentation, third-party benchmarks) report ~80% of cells confidently annotated (prediction.score >= 0.75) on typical PBMC queries — but that is a confidence threshold, not an accuracy.
- Per-cell-type accuracy on monocytes / CD4T / CD8T / B / NK is shown only in confusion-matrix figure panels and not as decimals in the abstract/results text I could access.

## Mapping to our metric
INDIRECT but closer than the integration benchmarks above. Reference mapping concordance (R = 0.911 on MAIT fractions cross-platform) tells you that *cell-type assignments* are reproducible cross-study, which is upstream of our metric — if cell-type calls are unstable, response-vector Pearson r will be noise. The MAIT correlation is computed at the *cohort-fraction* level (one number per donor) which is the same scale as cohort-level pseudobulk, so it is a loosely comparable order of magnitude: R ~0.9 is achievable for *abundance* concordance across platforms on PBMC. Note, however, that our response-vector r is computed on a 2000-dim gene-space delta, not a 1-dim fraction — so we should not expect 0.91 on our metric; the dimensionality and the diseased-minus-healthy subtraction both inject noise.

## Reference value for our calibration table
- Context; cite as upstream cell-type-assignment reproducibility precedent (R ~0.91 cross-platform on a 1-dim fraction). Treat as a soft *ceiling* hint for our metric on the most-abundant cell types, not a target.
