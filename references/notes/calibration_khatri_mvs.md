---
name: Multi-cohort analysis of host immune response identifies conserved protective and detrimental modules associated with severity across viruses (Khatri lab MVS)
description: 34-cohort + 3-scRNA-cohort meta-analysis; conserved viral severity signature replicates cross-cohort with Pearson R 0.43-0.93 (discovery), 0.33-0.78 (validation); scRNA monocyte-MVS correlations R=0.25-0.45.
type: reference
---
**Citation**: Zheng H, Rao AM, Dermadi D, Toh J, Murphy Jones L, Donato M, Liu Y, Su Y, Dai CL, Kornilov SA, Karagiannis M, Marantos T, Hasin-Brumshtein Y, Lu Y-F, Zhang Y, Wang J, Liu W, Wang Y, Zhang B, Mao Y, Beck CG, Ouyang Z, Davis MM, Heath JR, Atreya MR, Hotchkiss RS, Remy KE, Standiford TJ, Giamarellos-Bourboulis EJ, Liu Y, Khatri P. 2021. *Immunity* 54:753-768.e5.
**DOI**: 10.1016/j.immuni.2021.03.002
**URL**: https://www.cell.com/immunity/fulltext/S1074-7613(21)00114-X ; PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC7988739/ ; preprint: https://www.medrxiv.org/content/10.1101/2020.10.02.20205880v1

## Metric reported
This is **the most directly comparable paper** to Trinetravir. The Khatri lab pipeline:
1. Computes a per-sample MVS (Meta-Virus Signature) score = standardized geometric-mean expression of over- vs under-expressed genes, applied per-cohort after within-cohort COCONUT-style standardization.
2. Reports **Pearson R of MVS score vs ordinal severity** within each cohort, then aggregates across cohorts — this is conceptually analogous to "compute a response score per study, then check whether it correlates with the biological gradient consistently across studies."
3. Reports module-level cross-cohort **AUROC** for binary severity classification.
4. Validates at single-cell resolution: 702,970 immune cells / 289 PBMC samples / 3 cohorts; computes per-cell-type MVS-correlations.
5. Uses Hedges' g effect sizes for differential cell-proportion analysis across cohorts (cross-cohort meta-analysis effect-size language, same family as our cross-study agreement).

## PBMC values observed
**Bulk multi-cohort (34 cohorts, 4,780 samples, 16 viruses):**
- Cross-cohort Pearson R of MVS score vs severity, discovery (19 datasets, 1,674 samples): **0.43 <= R <= 0.93** per cohort.
- Validation cohorts (SARS-CoV-2, Ebola, chikungunya, additional): **0.33 <= R <= 0.78, p <= 1.8e-05** per cohort.
- Aggregate severity-MVS correlation across discovery: **R = 0.75, p < 2.2e-16**.
- Module-level severity classification: discovery AUROC **>= 0.929** (3,183 samples); validation AUROC **> 0.98** (1,154 samples, 4 viruses).
- Mild-vs-moderate distinction: SoM score AUROC **> 0.75** (vs MVS AUROC < 0.63 for this finer distinction).

**Single-cell (702,970 cells, 3 cohorts, PBMC):**
- CD14+ monocyte MVS correlation with severity: **R = 0.45, p = 2.7e-14**.
- CD16+ monocyte MVS correlation: **R = 0.25, p = 2.4e-05**.
- Generic myeloid cell MVS correlation: **R = 0.28, p = 2.4e-06**.
- Cross-cohort cell-proportion effect sizes (Hedges' g, non-severe vs HC): total monocytes 1.10, CD14+ 1.12, CD16+ -0.88 (severe vs non-severe), neutrophils 1.24 (severe vs HC), NK -0.85 to -1.03.

**Other Khatri-lab related papers (briefly):**
- Andres-Terre et al. 2015 *Immunity* (PMC4684904) — the original MVS, 396-gene signature, validation AUC range 0.84-1.00 across respiratory virus cohorts (RSV, H1N1, H3N2, influenza); 18 microarray datasets, 2,939 samples. Cohort-Pearson values not reported in same form as Zheng 2021.
- Lewis et al. 2025 *Immunity* (PMID 40532705) — "A conserved immune dysregulation signature" — most recent update extending the framework; not directly extracted here.

## Mapping to our metric
**DIRECT — this is the load-bearing citation.** The Khatri lab's per-cohort severity-MVS Pearson R is computed in exactly the same statistical family as our cross-study response-vector Pearson r: both are correlations between a derived per-cohort biological score and a reference (severity for them, the other study's response vector for us). Key differences to flag honestly:
- Khatri's R is severity-vs-score within each cohort (a 1-d signature scored against an ordinal severity outcome). Ours is response-vector-vs-response-vector across two studies in scaled-HVG gene space (a 2000-d cosine-like Pearson). So Khatri's R is a per-cohort generalization metric, ours is a pairwise cross-study replication metric — different but in the same order of magnitude.
- Khatri uses bulk transcriptome standardization (COCONUT / per-cohort z-scoring); we use Harmony on scRNA. Both are batch-correction-plus-meta-analysis pipelines.
- **The scRNA monocyte numbers (R = 0.45 for CD14+, R = 0.25 for CD16+) are the closest single-cell anchor we have**: a state-of-the-art cross-cohort scRNA harmonization pipeline produces R ~0.45 on the dominant myeloid response and R ~0.25 on the rarer one, against severity as the reference. Our 0.60 gate for monocytes and 0.25 for CD8T are in the right neighborhood — possibly slightly optimistic for monocytes (Khatri reports 0.45, we target 0.60) and well-calibrated for the rarer cell types (Khatri 0.25 for CD16+ monocytes ≈ our 0.25 for CD8T).

## Reference value for our calibration table
- **DIRECT anchor**: use Khatri lab Zheng et al. 2021 R = 0.45 (CD14+ monocyte) as the literature soft-floor for monocyte response replication; R = 0.25 (CD16+ monocyte / generic myeloid) as the literature soft-floor for less-abundant cell types. These are cohort-vs-severity, not study-vs-study, so cite as "comparable order of magnitude" not "identical metric." Our per-bucket thresholds (0.60 monocytes, 0.25 CD8T) should be discussed against these values in the calibration table.
