---
name: A benchmark of batch-effect correction methods for single-cell RNA sequencing data
description: 14-method benchmark with kBET/LISI/ASW/ARI; on the PBMC cross-protocol dataset Harmony tied for first overall with Seurat 3 (ARI_batch > 0.97, top-tier cLISI).
type: reference
---
**Citation**: Tran HTN, Ang KS, Chevrier M, Zhang X, Lee NYS, Goh M, Chen J. 2020. *Genome Biology* 21:12.
**DOI**: 10.1186/s13059-019-1850-9
**URL**: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9 ; PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC6964114/

## Metric reported
Five metric families: (1) kBET acceptance rate per cell type (batch mixing); (2) LISI — iLISI (batch mixing, higher = better) and cLISI (cell-type purity, lower = better); (3) ASW (silhouette) for batch and for cell-type; (4) ARI on batch labels and on cell-type labels post-clustering; (5) DEG concordance pre/post correction. They benchmark 14 methods (Harmony, Seurat 2, Seurat 3, LIGER, MNN Correct, fastMNN, BBKNN, scMerge, ZINB-WaVE, scGen, Combat, limma, Scanorama, MMD-ResNet) on 10 datasets including a PBMC cross-protocol scenario (Dataset 5).

## PBMC values observed
Dataset 5 is two PBMC batches from different 10X protocols, containing CD4/CD8 T cells and CD14/FCGR3A monocytes (the cell-type structure most relevant to our task).
- Harmony **tied for best overall** on Dataset 5 (PBMC) via rank-sum across all metrics, tied with Seurat 3.
- Harmony was **third place** on iLISI (after LIGER and Seurat 2) on Dataset 5.
- Harmony was **tied for best cLISI** with Seurat 3 (cell-type purity preserved).
- Harmony was **third place** on kBET acceptance.
- Harmony and Seurat 3 produced **ARI_batch > 0.97** on Dataset 5 (interpretation in this paper: ARI_batch high means batch-label clustering is destroyed post-correction, i.e. good mixing).
- Exact decimal values for Harmony's kBET, iLISI, cLISI on Dataset 5 are reported only in the supplementary heatmaps/tables; I extracted rankings from the PMC text but not the raw numbers without the supplement.
- Overall lab recommendation: Harmony, LIGER, Seurat 3 are the recommended methods.

## Mapping to our metric
INDIRECT. Same caveat as scIB: these are embedding-mixing and cell-type-purity metrics, not response-vector replication. The closest thing they report to "do biological signals survive correction" is DEG concordance pre/post correction within a single dataset — which is again not cross-study response agreement. Useful as evidence that Harmony is a defensible choice for PBMC cross-protocol integration (which is exactly our setting — we have 4 PBMC studies on 10X v2/v3/5'), but not a numerical target.

## Reference value for our calibration table
- Indirect; cite as method-choice justification — "Harmony tied for best overall on PBMC cross-protocol integration in Tran et al. 2020." Do not derive a numerical Pearson r threshold from this paper.
