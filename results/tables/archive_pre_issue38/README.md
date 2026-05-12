# Archive: pre-Issue-38 calibration outputs

This directory preserves the calibration table outputs from before the **Issue 38 bootstrap N reconciliation** (commit `793344e`, 2026-05-12). Pre-Issue-38 outputs used non-canonical N for bootstrap CIs and (in Session 7) for permutation p-values.

## Files

| Archived file | Original location | Pre-Issue-38 N values |
|---------------|-------------------|----------------------|
| `heldout_v2_calibration_yoshida_2022_n200_archived.csv` | `../heldout_v2_calibration_yoshida_2022_n1000.csv` | N_PERM=1000, N_BOOTSTRAP=200 |
| `heldout_v2_calibration_allen_atlas_monocyte_n200_archived.csv` | `../heldout_v2_calibration_allen_atlas_monocyte_n1000.csv` | N_PERM=1000, N_BOOTSTRAP=200 |
| `heldout_v2_calibration_gse157829_n200_archived.csv` | `../heldout_v2_calibration_gse157829_n1000.csv` | N_PERM=1000, N_BOOTSTRAP=200 |
| `heldout_v2_calibration_randolph_2021_n200_issue31_archived.csv` | `../heldout_v2_calibration_randolph_2021_n1000_issue31.csv` | N_PERM=1000, N_BOOTSTRAP=200 (Issue 31 cross-bucket healthy ref variant) |
| `heldout_v2_calibration_combined_n200_archived.csv` | `../heldout_v2_calibration_combined_n1000.csv` | mixed N_BOOTSTRAP per row (100/63/71/72/59 — variable due to degenerate cases) |
| `sensitivity_within_cohort_n500_perm_n200_boot_archived.csv` | `../sensitivity_within_cohort.csv` | N_PERM=500, N_BOOTSTRAP=200 |
| `heldout_issue_verdicts_n200_archived.csv` | `../heldout_issue_verdicts_n1000.csv` | Verdicts derived from N_BOOTSTRAP=200 source tables |

## Why archived (vs. relying on git history)

Per Issue 38 audit-trail discipline (`METHODS_CHOICES.md`), pre-Issue-38 outputs are preserved as explicit archived files so reviewers can directly compare old-N vs new-canonical-N CI bounds without git operations. Supplementary supports the "verdicts unchanged under N reconciliation" claim (verdicts assigned mechanically on observed `r_mvs` point estimate alone, not CI bounds).

## What replaces these files

After Task #17 recompute (post-Issue-38, 2026-05-12+):
- `../heldout_v2_calibration_*_n1000.csv` will be regenerated at canonical N_BOOTSTRAP=1000.
- `../sensitivity_within_cohort.csv` will be regenerated at canonical N_PERM=1000, N_BOOTSTRAP=1000.
- `../heldout_issue_verdicts_n1000.csv` will be regenerated from canonical-N sources.

## Reference

- `METHODS_CHOICES.md` Issue 38 — full audit map + reconciliation plan
- Commit `793344e` — Issue 38 documentation + config harmonization + script patches
- Recompute task (Task #17) — produces canonical-N outputs that supersede these
