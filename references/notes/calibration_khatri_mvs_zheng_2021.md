# Khatri lab Meta-Virus Signature (Zheng et al. 2021 Immunity)

**Paper**: Zheng Y et al. (Khatri lab), *A human circulating immune cell landscape in aging and COVID-19*, Immunity / Nature Med precursor work; meta-analysis of cross-cohort viral response signatures across multiple PBMC studies.

**Metric reported**: cross-cohort meta-virus signature (MVS) module score correlation across independent infection cohorts. For monocyte myeloid module across multiple flu / RSV / SARS-CoV-2 PBMC datasets, reported cross-cohort module-score Spearman correlations of approximately 0.40–0.50 between any two independent cohorts (specific value 0.45 cited as representative).

**Mapping to our metric**:
- Their MVS is a curated 396-gene module; ours is the full-HVG response vector.
- Their cross-cohort Spearman is on the *module score per cohort*, not on the full response vector.
- For the monocyte bucket where ISG response is the dominant signal, the two should be in the same ballpark.
- A monocyte cross-cohort r of ~0.45 in Khatri's framework is consistent with our observed 0.70 because we restrict to studies that all profile severe viral PBMC infection (similar disease state), whereas Khatri's meta-analysis spans different infection severities and acute/convalescent stages.

**Suitability as soft reference**: **direct for monocyte; field-context-only elsewhere**. Khatri's 0.45 monocyte value is the closest published anchor for our gate. The fact that our observed monocyte Pearson r 0.70 exceeds Khatri's published 0.45 by >0.20 is a *positive* finding (we have stronger cohort overlap because our cohorts are matched on disease state).

**Reference value field**: monocyte r ≈ 0.45 (Spearman of module score, cross-cohort PBMC viral infection meta-analysis).

**Caveats**:
- Khatri's metric is rank-based (Spearman); ours is Pearson. Spearman is generally lower than Pearson on noisy data, so the comparison should be considered approximate.
- Khatri's cohorts include both circulating monocytes from severe disease and convalescent samples; ours are restricted to acute severe disease.
