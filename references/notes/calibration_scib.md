---
name: Benchmarking atlas-level data integration in single-cell genomics (scIB)
description: 14-metric integration benchmark; reports Harmony as a top performer on PBMC/immune tasks, but uses embedding-mixing metrics (ASW/ARI/NMI/kBET) not response-vector replication.
type: reference
---
**Citation**: Luecken MD, Buttner M, Chaichoompu K, Danese A, Interlandi M, Mueller MF, Strobl DC, Zappia L, Dugas M, Colome-Tatche M, Theis FJ. 2022. *Nature Methods* 19:41-50.
**DOI**: 10.1038/s41592-021-01336-8
**URL**: https://www.nature.com/articles/s41592-021-01336-8 ; https://theislab.github.io/scib-reproducibility/

## Metric reported
scIB defines a composite integration score = 0.6 * bio-conservation + 0.4 * batch-correction. The 14 underlying metrics are split into two families: (a) bio-conservation — NMI, ARI, ASW_label (cell-type silhouette), isolated-label F1, isolated-label silhouette, graph cLISI, cell-cycle conservation, HVG conservation, trajectory conservation; (b) batch removal — ASW_batch, graph iLISI, kBET, graph connectivity, PC regression on batch. All metrics are scaled to [0,1] (higher = better). The benchmark was run on 13 integration tasks, including an "Immune Cell Hum" task built from 10 PBMC + bone-marrow batches across 5 datasets.

## PBMC values observed
- Headline qualitative result: for the human immune cell task (PBMC + BM, 10 batches, 5 datasets), the top-ranked methods were Scanorama (embedding), FastMNN (embedding), scANVI, and Harmony. Harmony is described as performing well on the "less complex" tasks (Immune Cell Hum is one of these).
- Specific numerical per-method per-task scores on Immune Cell Hum are published only as figure panels (heatmaps) on the scIB-reproducibility website (`theislab.github.io/scib-reproducibility`) and in Supplementary Table 5 of the Nature Methods paper. I was unable to extract the exact decimal values from the public landing pages or PubMed abstract through WebFetch (both the Nature article and bioRxiv preprint return 403 Forbidden; the scib-reproducibility site renders numbers only inside PNG figures).
- Order-of-magnitude expectation from the published figure panels (qualitative, NOT extracted as decimals): Harmony's overall integration score on Immune Cell Hum is approximately 0.65-0.75, with bio-conservation slightly higher than batch-removal. scANVI and Scanorama-embed sit slightly above Harmony on this task.

## Mapping to our metric
INDIRECT. scIB measures *embedding integration*: do cells from different batches mix in latent space (iLISI, kBET, ASW_batch) while preserving cell-type structure (NMI, ARI, ASW_label). Our metric measures *response signature replication*: does the diseased - healthy contrast vector in scaled-HVG space agree across studies (Pearson r of per-cell-type pseudobulk delta-vectors). High scIB scores guarantee that monocytes-from-Lee and monocytes-from-Wilk land in the same neighborhood post-Harmony; they do NOT guarantee that COVID-induced monocyte shifts replicate across studies. A study could have scIB ~0.80 yet produce wildly different disease-response vectors if cohort case-mix differs. Useful as field context, not as a numerical anchor for our r threshold.

## Reference value for our calibration table
- Indirect; cite as field context — "Harmony is an established top-tier method on PBMC integration benchmarks (Luecken et al. 2022)" — but do NOT use scIB scores as a numerical target for our cross-study Pearson r gate.
