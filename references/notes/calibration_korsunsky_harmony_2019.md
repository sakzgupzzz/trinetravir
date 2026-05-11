# Korsunsky et al. 2019 (Nature Methods) — Harmony

**Paper**: Korsunsky I et al., *Fast, sensitive and accurate integration of single-cell data with Harmony*, Nature Methods 16, 1289–1296 (2019).

**Metric reported**: LISI (Local Inverse Simpson's Index) — iLISI for batch mixing, cLISI for cell-type retention. Reported on the PBMC stimulation atlas (Kang et al. 2018 stim/ctrl PBMCs).

**PBMC-specific values**: Harmony on Kang et al. PBMC stim/ctrl achieves iLISI ≈ 2.0 (good mixing for 2-batch task; max possible = 2), cLISI ≈ 1.05 (good cell-type separation; close to ideal 1.0). The paper documents per-cell-type Harmony integration as the recommended use when cell types are known.

**Mapping to our metric**: Korsunsky's per-cell-type Harmony is the integration approach we use. The paper's PBMC stim/ctrl benchmark is a *2-batch within-study* setup; our setup is *4-study cross-cohort*. The Korsunsky result establishes the floor of expected integration quality but does not directly anchor our cross-study response-vector r.

**Suitability as soft reference**: **indirect**. Establishes the methodological foundation for our per-cell-type Harmony choice (METHODS_CHOICES Issue 6 + Issue 7).

**Reference value field**: n/a (LISI is not Pearson r).
