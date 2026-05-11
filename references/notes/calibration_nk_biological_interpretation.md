# NK divergence under per-cell-type vs global Harmony — biological interpretation

**Empirical observation (Session 1 Step 8b, heuristic threshold)**:
- Per-cell-type Harmony NK r = 0.384 (PASS at threshold 0.35)
- Global Harmony NK r = 0.308 (FAIL at threshold 0.35)
- Delta: per-cell-type higher by 0.077 — the largest divergence among 5 buckets.

For monocyte / B / CD4T / CD8T, per-cell-type and global Harmony produce r values within ±0.06 of each other. NK is the outlier. The Session 3 calibrated comparison confirms whether this divergence is statistically meaningful relative to null + split-half ceiling distributions; the *biological* interpretation below stands independent of that verdict.

## Why NK is the bucket where harmonization protocol matters most

NK cells in human PBMC are not a homogeneous population. Functional subsets defined by surface marker expression include:

- **CD56-bright NK cells** (~10% of circulating NK): cytokine producers, high IFN-γ, lower cytotoxicity, predominantly tissue-trafficking phenotype.
- **CD56-dim CD16+ NK cells** (~90%): cytotoxic effectors, ADCC-capable.
- **Adaptive / memory NK cells**: long-lived expansions in CMV+ donors, distinct transcriptional program.
- **Cycling NK**: actively proliferating, abundant in acute viral infection.

Single-cell PBMC studies routinely surface 3–5 NK subsets per cohort. Within-bucket transcriptional heterogeneity is therefore *high* — substantially higher than within a single CD4T subset or a B-cell sublineage.

## Why this produces the per-cell-type vs global Harmony divergence

**Global Harmony** runs cluster-aware batch correction on ALL 244k PBMC cells together with study_id as the sole batch key. The clustering step (Harmony's internal soft k-means) groups cells by transcriptional similarity, regardless of cell-type label. NK subsets that share transcriptional features with other cell types (e.g., CD56-bright NK has elevated cytokine gene expression overlapping with activated CD4T; cycling NK shares proliferation markers with cycling CD8T) can be assigned to clusters that mix NK with non-NK cells. Harmony's per-cluster correction then *spreads* the NK heterogeneity across cell-type boundaries — diluting the NK-specific signal at the expense of preserving inter-cell-type structure.

**Per-cell-type Harmony** runs the correction *within* the NK bucket only. The clustering step now groups NK subsets among themselves — CD56-bright clusters with CD56-bright, CD16+ with CD16+, cycling NK with cycling NK. Cross-cohort batch effects within each subset are corrected without leaking signal across cell-type boundaries. The NK-specific disease-induced signal (ISG induction, killer-Ig-like receptor modulation, etc.) is preserved at higher fidelity.

The other 4 buckets show smaller divergence because:
- **Monocyte** subsets (classical / non-classical / intermediate) are transcriptionally clean — sub-bucket boundaries align with major marker genes. Global Harmony's cluster assignments rarely cross monocyte subset boundaries.
- **B cells** have less transcriptional overlap with other lineages.
- **CD4T / CD8T** have within-subset gradients (naive → CM → EM → effector) but those gradients run *within* each cell-type bucket; global Harmony's cluster assignments respect the major T-cell axis.

NK is uniquely positioned to suffer from global Harmony because (a) within-bucket heterogeneity is high, and (b) the heterogeneity overlaps transcriptionally with other lineages (cytokine programs, proliferation, activation signatures shared with myeloid + T cells).

## Implication for the paper

If the Session 3 calibrated comparison confirms per-cell-type Harmony is meaningfully better for NK, the paper's methods section should:
1. Justify per-cell-type Harmony as the v1 primary protocol explicitly via the NK-bucket evidence.
2. Cite the NK subset heterogeneity literature (Marquardt et al. 2017 Nature; Crinier et al. 2018 Immunity; Smith et al. 2020 cell-state continuum work).
3. Report the global-Harmony comparison in the supplementary as a sensitivity, with the NK divergence flagged as the load-bearing bucket.

## References for paper write-up
- Marquardt N et al. 2017 *The human NK cell response to yellow fever vaccination*, Nature Comm.
- Crinier A et al. 2018 *High-dimensional single-cell analysis identifies organ-specific signatures and conserved NK cell subsets in human and mouse*, Immunity.
- Smith SL et al. 2020 *Diversity of peripheral blood human NK cells identified by single-cell RNA sequencing*, Blood Adv.
