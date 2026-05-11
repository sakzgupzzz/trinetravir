# Khatri MVS external validation (Session 5 Part B) — summary

Date: 2026-05-11
MVS gene list: 86 canonical type-I-IFN ISGs + Khatri core (see `data/reference/khatri_mvs_module_genes.txt`).

## Per-bucket cross-study Pearson r — full HVG vs MVS subset

| Bucket | n_HVG | n_MVS_in_HVG | r_full | r_MVS_subset | Δ |
|---|---|---|---|---|---|
| monocyte | 4000 | 57 | 0.7012 | 0.6566 | -0.0446 |
| B | 4000 | 48 | 0.2971 | 0.3587 | +0.0616 |
| NK | 4000 | 47 | 0.3845 | 0.4690 | +0.0845 |
| CD4T | 4000 | 57 | 0.3214 | 0.4818 | +0.1604 |
| CD8T | 4000 | 61 | 0.1686 | 0.4000 | +0.2315 |

**Verdict — monocyte MVS check: WEAK PASS.** r_MVS (0.6566) is high in absolute terms but lower than r_full (0.7012). This is consistent with the full HVG response vector containing additional non-ISG signal that adds coherence; the MVS subset captures the canonical ISG-dominated component.

## Scope caveat

This is a within-corpus consistency check on our own harmonized response vectors, not a true external-cohort validation against bulk RNA-seq data. A FULL external validation against an independent PBMC bulk cohort is deferred to v1.5; bulk cohort acquisition + alignment to v1's pipeline is outside Session 5's scope.

## Interpretation

- The check is necessary but not sufficient. Passing the MVS subset coherence check confirms our pipeline captures the canonical IFN-stimulated gene signature.
- A failing check (r_MVS << r_full or r_MVS < 0.3) would indicate our pipeline is finding cross-study coherence in non-MVS genes — possible HVG selection artifact or batch-correction over-fitting.
- A passing check supports proceeding with Session 5 Parts C-E and downstream Phase 4+ work.
