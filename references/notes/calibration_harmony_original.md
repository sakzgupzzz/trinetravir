---
name: Fast, sensitive and accurate integration of single-cell data with Harmony
description: Original Harmony paper; on 3 PBMC datasets (10X 3' v1, 3' v2, 5') Harmony achieves iLISI median 1.96 (95% [1.36, 2.56]) up from 1.00 pre-integration, with cLISI 1.00.
type: reference
---
**Citation**: Korsunsky I, Millard N, Fan J, Slowikowski K, Zhang F, Wei K, Baglaenko Y, Brenner M, Loh P-R, Raychaudhuri S. 2019. *Nature Methods* 16:1289-1296.
**DOI**: 10.1038/s41592-019-0619-0
**URL**: https://www.nature.com/articles/s41592-019-0619-0 ; PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC6884693/ ; preprint: https://www.biorxiv.org/content/10.1101/461954v2

## Metric reported
The paper introduces LISI (Local Inverse Simpson Index) computed in two flavors:
- **iLISI** (integration LISI): for each cell, the effective number of distinct *batches* in its local neighborhood. Range [1, N_batches]. Higher = better mixing.
- **cLISI** (cell-type LISI): for each cell, the effective number of distinct *cell types* in its local neighborhood. Range [1, N_celltypes]. Lower (closer to 1) = better cell-type purity.

These are local diversity metrics on the post-correction embedding. The paper also shows accuracy on simulated cell-type labels and compares to MNN Correct, BBKNN, MultiCCA, and Scanorama.

## PBMC values observed
The three-PBMC-dataset experiment (Cell Hashing-style, three 10X chemistries: 3' v1, 3' v2, 5' end) — most directly relevant to Trinetravir's setup:
- **Pre-integration iLISI**: median 1.00, 95% [1.00, 1.00] (i.e. zero mixing — cells clustered by chemistry).
- **Post-Harmony iLISI**: median 1.96, 95% [1.36, 2.56] (out of theoretical max 3.0 for three datasets; ~65% of max).
- **Post-Harmony cLISI**: median 1.00, 95% [1.00, 1.02] (cell-type structure fully preserved).
- Competing methods on the same three PBMC datasets: median iLISI "failed to exceed 1.1" for MNN Correct, BBKNN, MultiCCA, Scanorama.

No kBET, ARI, or pseudobulk-response-Pearson-r values are reported on the PBMC datasets in the main paper.

## Mapping to our metric
INDIRECT. LISI is a per-cell local-neighborhood diversity score on the corrected embedding — it tells you "do batch labels mix locally" not "do disease-response signatures replicate." A high iLISI is necessary but not sufficient for our metric: cells must mix in latent space before any meaningful pseudobulk delta can be computed cross-study, but the converse is not true. This paper is the seminal Harmony citation and is appropriate as the algorithmic citation in our methods section, but its reported PBMC numbers (iLISI 1.96, cLISI 1.00) cannot be mapped to a Pearson r threshold.

## Reference value for our calibration table
- Indirect; cite as method origin and embedding-mixing precedent on PBMC. Note in supplementary that Harmony achieves iLISI ~65% of theoretical max on three PBMC chemistries (Korsunsky et al. 2019), consistent with what we observe pre-response-vector computation.
