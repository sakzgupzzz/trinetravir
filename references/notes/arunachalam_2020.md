# Systems biological assessment of immunity to mild versus severe COVID-19 infection in humans
**Citation**: Arunachalam PS, Wimmers F, Mok CKP et al. 2020, *Science* 369(6508):1210-1220
**DOI**: 10.1126/science.abc6261
**URL**: https://pmc.ncbi.nlm.nih.gov/articles/PMC7665312/ (PDF downloaded via EuropePMC)
**Status**: downloaded (main text only — Materials and Methods are in supplementary)
**Tier**: 1

## Core claim
Arunachalam et al. apply a multi-modal systems-immunology pipeline (mass cytometry, plasma cytokine multiplex/Olink, bulk RNA-seq, and CITE-seq scRNA-seq) to PBMC and plasma samples from 76 COVID-19 patients and 69 healthy controls across two geographically distinct cohorts (Hong Kong and Atlanta). They argue COVID-19 induces a paradoxical "spatial dichotomy": peripheral blood myeloid cells (monocytes, mDCs, pDCs) show *suppressed* innate function (low HLA-DR, blunted IFN-α via impaired mTOR signaling, reduced TLR-stimulated cytokines) while plasma carries elevated EN-RAGE (S100A12), TNFSF14, and oncostatin M — inflammatory mediators they hypothesize originate from lung tissue rather than peripheral leukocytes. The CITE-seq sub-study (Atlanta only, 7 COVID + 5 healthy) shows ISG induction is early and transient, peaking in moderate (not severe) cases. The contribution to single-cell COVID-19 atlases is the **CITE-seq sub-study**, which is the slice we get from cellxgene.

## Methodology / approach
- **Two cohorts** (intentional cross-population comparison): Princess Margaret Hospital, Hong Kong (n=36 COVID, n=45 healthy) + Hope Clinic at Emory University, Atlanta (n=40 COVID, n=24 healthy + 16 flu/RSV). Age- and sex-matched controls within each cohort.
- **Severity classification**: 4-tier — Mild/moderate, Severe (no ICU), ICU, Convalescent. Hong Kong cohort skews mild (75% mild/moderate); Atlanta skews severe (60% Severe-no-ICU, 18% ICU). Some patients sampled at multiple time points.
- **Treatment confounds**: Hong Kong cohort received IFN-β1 (20%), corticosteroids (19%), antivirals (61%); Atlanta = NA.
- **Assays applied**:
  - **Phospho-CyTOF** (mass cytometry, 22 surface + 12 intracellular markers): 54 PBMC samples (36 patients) HK + 19 samples (16 patients) ATL.
  - **In-vitro TLR stimulation + intracellular staining**: 17 samples ATL only.
  - **Olink 92-cytokine plasma proteomics**: 36 samples ATL + 19 flu/RSV.
  - **CITE-seq scRNA-seq**: **7 PBMC samples (7 COVID patients) + 5 healthy controls — Atlanta only**.
  - **Bulk RNA-seq**: 17 PBMC samples (15 patients) ATL.
  - **Bacterial product (LPS, bacterial DNA) plasma assays**: 51 samples ATL.
- **CITE-seq protocol details (main text)**: DC magnetic enrichment, then mixed ~1:2 with total PBMCs before loading on 10x. Materials and Methods (specific kit version, sequencer, cellranger version, cell-type annotation method) are in the **supplementary file**, not in the main PDF; the main text says "CITE-seq" + "10x" by implication and confirms data deposited to GEO **GSE155673** (CITE-seq) and **GSE152418** (bulk).
- **Cell counts**: CITE-seq totals not reported in the main-text excerpt; Figure 4 panels reference "n=12" combined samples (5 healthy + 7 COVID).
- **Annotation**: manual UMAP + cluster annotation; ISG signature (33 genes) derived from blood transcription modules (BTMs) M75/M127/M150 enriched in DEGs.

## Key findings
- pDC frequency reduced in COVID-19 PBMCs in both cohorts (Wilcoxon p=0.01 ATL, p=0.0015 HK) — but not correlated with time-from-symptom-onset.
- Reduced pS6 in pDCs (low mTOR activation) → reduced IFN-α and TNF-α response to TLR3/TLR7-8 stimulation in COVID-19 pDCs vs healthy.
- Reduced HLA-DR + CD86 expression on monocytes/mDCs; blunted IL-6/TNF-α/IL-1β response to TLR2/4/5 + viral TLR cocktail stimulation.
- Plasma cytokines: 43/71 detected proteins significantly upregulated in COVID-19 (Olink panel); EN-RAGE (S100A12), TNFSF14, oncostatin M correlate with severity.
- Plasma bacterial DNA + LPS positively correlate with EN-RAGE, TNFSF14, OSM, IL-6 — suggests gut/lung bacterial translocation contributes to systemic inflammation.
- CITE-seq: ISGs strongly induced in **moderate** COVID-19 but not severe/ICU — consistent with Hadjadj 2020 et al. (cited).
- ISG signature in bulk RNA-seq has strong temporal dependence: peaks early, declines by ~day 10 post symptom onset.
- TNFSF14 and OSM gene expression are *down* in PBMCs of severe patients despite plasma protein being up → tissue (likely lung) origin, not blood leukocyte origin.
- EN-RAGE (S100A12) gene expression IS high in blood myeloid cells of severe patients (so this one mediator is cell-of-origin in blood).

## Hyperparameters / evaluation details
N/A — dataset paper. Linear modeling for severity-discriminating features used cohort (HK vs ATL) as covariate to identify cross-cohort consistent signals. ISG signature derived from BTM-enriched DEGs (33 genes).

## Relevance to Trinetravir
**Smallest CITE-seq sub-study of our 4 studies**: only **7 COVID + 5 healthy donors** (Atlanta cohort only) for the single-cell modality — even though the full systems-immunology study has 76 + 69 donors. This is the cellxgene record we'll be using and is by far the lowest-N PBMC contribution. Severity stratification is **4-tier ordinal (Moderate / Severe-no-ICU / ICU / Convalescent)** but the CITE-seq subsample size (n=7) is too small to power within-study severity contrasts; we should treat this study primarily as additional cells for the "COVID-positive" pool with weak severity weighting. Healthy controls (n=5) are matched within-cohort. Time-from-symptom-onset is annotated (their Figure 4E uses days 0-60 axis) and is critical because their core finding is that ISG response is time-dependent. **Sequencing platform**: CITE-seq on 10x — the chemistry version (3′ v2 vs v3 vs v3.1) is in supplementary methods we did not extract; we need to either fetch the supplement or read from the GEO GSE155673 metadata to know definitively.

## Caveats / limitations
- Authors' stated: cross-cohort heterogeneity (HK skews mild, ATL skews severe) makes it hard to disentangle population effects from severity effects (they used cohort as covariate to mitigate); treatment heterogeneity in HK cohort (IFN-β1, steroids, antivirals); enrichment-then-mix CITE-seq design for DCs distorts cell-type proportions and they note this as a limitation for compositional analysis.
- For us additionally: very small CITE-seq sample (n=12 total → may not survive scvi/Harmony integration meaningfully; dominated by other cohorts in joint embedding); supplementary methods retrieval needed for sequencing-platform version; CITE-seq antibody panel adds protein covariates that cellxgene may not preserve (we are RNA-only); Atlanta-only CITE-seq subset is geographically narrow despite the wider study span.
