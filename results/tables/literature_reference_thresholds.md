# Literature reference values for cross-cohort PBMC response-vector coherence

Synthesized from publicly available benchmarks. See `references/notes/calibration_*.md` for per-paper detail.

| Paper | Metric | PBMC value | Comparability | Notes |
|-------|--------|-----------|---------------|-------|
| Khatri lab MVS (Zheng 2021 Immunity) | Cross-cohort Spearman of viral-response module score (monocyte) | ≈ 0.45 | **direct** | Closest published anchor for our monocyte bucket. Their meta-analysis spans severities + viruses; ours restricts to acute severe PBMC. Higher r expected for ours. |
| Pan et al. 2023 | Cross-virus monocyte module Spearman (SARS-CoV-2 vs HIV-1) | 0.55 – 0.65 | **direct** | Most direct anchor for our cross-virus Gate 1 question. Lee SARS-vs-IAV monocyte Pearson r in this range is the expected scale. |
| Luecken et al. 2022 (scIB) | Composite integration score (PBMC immune atlas, Harmony) | 0.74 / 0.66 | field-context-only | Embedding integration quality, not response-vector r. Establishes Harmony is competitive. |
| Tran et al. 2020 | kBET / ASW composite (PBMC, Harmony) | ≈ 0.72 | field-context-only | Same caveat as Luecken — integration metric, not response coherence. |
| Korsunsky et al. 2019 (Harmony) | iLISI / cLISI on Kang 2018 PBMC stim/ctrl | 2.0 / 1.05 | indirect | 2-batch within-study setup; not 4-study cross-cohort. Establishes methodological foundation. |

## Headline reference values per bucket

- **monocyte**: 0.40 – 0.65 (Khatri meta-virus + Pan cross-virus). Our gate threshold 0.60 sits at the upper end.
- **B**, **NK**, **CD4T**, **CD8T**: no published direct anchors. Field convention says lymphoid lineages are noisier than myeloid; thresholds should be lower than monocyte.
- **Cross-virus**: 0.45 – 0.65 (Pan et al. monocyte). Our Gate 1 Lee SARS-vs-IAV per-stratum monocyte r 0.65 sits at the upper end.

## Suitability summary

Only Khatri MVS and Pan et al. 2023 are *direct* anchors for our response-vector Pearson r convention (and even those are Spearman-based module scores, not full-HVG Pearson). The remaining literature provides *field context* (Harmony is the standard PBMC tool; PBMC integration is feasible) but does not anchor specific Pearson r thresholds.

This is why our gate thresholds are *primarily* calibrated via the donor-level permutation null + split-half ceiling (calibrated significance), and *secondarily* checked against the Khatri + Pan reference values. The calibrated approach is the headline; literature is a soft sanity check.
