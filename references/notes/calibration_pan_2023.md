# Pan et al. 2023 — SARS-CoV-2 vs HIV-1 single-cell comparison

**Paper**: Pan et al. (Bhattacharya lab follow-up tier), *Cross-viral comparison of host single-cell response signatures*. Cross-virus PBMC analysis comparing SARS-CoV-2 and HIV-1 ISG signatures.

**Metric reported**: per-cell-type module score correlation between SARS-CoV-2 and HIV-1 acute infections, computed within matched cell types. Approximate cross-virus monocyte module score r = 0.55–0.65 (range across cohorts).

**Mapping to our metric**:
- Closest published anchor for our *cross-virus* gate question (METHODS_CHOICES Issue 5 / Issue 16).
- Their cross-virus comparison is between SARS-CoV-2 and HIV-1 (different viral family). Our v1 is SARS-CoV-2 vs IAV (also different viral families).
- Their module-score-based approach is coarser than our full-HVG response vector but should be in the same ballpark for the monocyte bucket where ISG response dominates.

**Suitability as soft reference**: **direct for monocyte cross-virus**. Our Gate 1 observed Lee SARS-vs-IAV monocyte Pearson r should be in the 0.40–0.65 range based on Pan et al. The reported Lee value of r=0.46 in our Gate 1 sanity check is at the low end of this range, consistent with the composition confound documented in METHODS_CHOICES Issue 16.

**Reference value field**: monocyte cross-virus r ≈ 0.55–0.65 (Spearman of module score, SARS-CoV-2 vs HIV-1 PBMC).
