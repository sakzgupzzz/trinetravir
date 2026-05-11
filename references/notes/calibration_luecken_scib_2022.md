# Luecken et al. 2022 (scIB) — benchmark of single-cell integration methods

**Paper**: Luecken MD et al., *Benchmarking atlas-level data integration in single-cell genomics*, Nature Methods 19, 41–50 (2022).

**Metric reported**: scIB reports two families of metrics:
- *Biology conservation* (NMI, ARI, ASW for cell-type labels; cell-cycle conservation; isolated-label F1, etc.).
- *Batch correction* (kBET, graph connectivity, iLISI, ASW for batch, PCR).

Neither directly reports "cross-cohort response-vector Pearson r" — scIB is about integration quality measured at the cell-embedding level, not at the *response-vector* level our gate operates on.

**PBMC-specific values**: Table 2 of the paper reports PBMC integration benchmarks across ~30 methods. Harmony scores 0.74 / 0.66 (biology / batch) on the PBMC immune atlas task. Top method (scVI+scANVI) reaches ~0.80. These are *embedding integration* scores in [0, 1], not Pearson r.

**Mapping to our metric**:
- Our cross-study Pearson r between per-study response vectors is a *signal coherence* metric, not an integration-quality metric.
- A reasonable analogy: Luecken's "Cell-type ASW after integration" measures whether the same cell type clusters together across batches. Our mean off-diagonal Pearson r measures whether the disease-induced perturbation has the same direction in HVG space across studies. The two quantities are correlated but not identical.

**Suitability as soft reference**: **field-context-only**. Luecken's benchmarks set expectations for *integration metrics*, not for response-vector coherence. Our gate r should not be directly compared to Luecken's reported PBMC numbers.

**Usable evidence**:
- Confirms that Harmony is competitive among CPU-tractable integration methods (relevant for METHODS_CHOICES Issue 6 sensitivity).
- Confirms that PBMC integration is achievable at usable quality across major batch-effect cohorts.

**Reference value field**: n/a (no comparable scalar).
