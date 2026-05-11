# Tran et al. 2020 (Genome Biology) — scRNA-seq batch correction benchmark

**Paper**: Tran HTN et al., *A benchmark of batch-effect correction methods for single-cell RNA sequencing data*, Genome Biology 21, 12 (2020).

**Metric reported**: kBET (k-nearest-neighbor batch effect test) + ASW (average silhouette width) for batch and cell-type retention. PBMC-specific values across methods.

**PBMC-specific values**: PBMC immune-system benchmark: Harmony achieves kBET ≈ 0.07 (lower = better mixing), ASW-cell-type ≈ 0.55 (higher = better cell-type separation). Scaled to [0,1] composite scores: ~0.72 PBMC integration quality for Harmony.

**Mapping to our metric**: same caveat as Luecken et al. — Tran's metrics evaluate integration quality at the cell-embedding level, not response-vector coherence. Not a direct anchor.

**Suitability as soft reference**: **field-context-only**. Confirms Harmony performs well on PBMC integration tasks at the embedding level.

**Reference value field**: n/a.
