"""Acquire Lee 2025 — HIV-1 early infection multi-omics cohort.

Session 6A Part A4 (TERTIARY held-out cohort). STUB. MARGINAL SAMPLE SIZE.

Citation: Lee et al. eLife 2025, PMC12370253.
Design: 9 early HIV (<6 months) donors, 5 scRNA-seq + 4 snRNA-seq multiome.
Healthy control N TBD on download.

CRITICAL pre-flight check: verify healthy N ≥ 4. If N < 4, falls back to
qualitative-validation status per Issue 30 fallback protocol.

See references/notes/hiv_biology_note.md for the forward-flag on retrovirus
biology (different from v1 corpus's RNA respiratory viruses).
"""

from __future__ import annotations

import sys


def main() -> int:
    print("[STUB] acquire Lee 2025 HIV cohort")
    print("  9 early HIV donors (<6 months from infection)")
    print("  healthy N TBD — verify on download. Falls back to qualitative-only if <4.")
    print("  modality: 5 scRNA + 4 snRNA multiome")
    print("  storage: data/raw/lee_2025_hiv/")
    print()
    print("ACTUAL DOWNLOAD NOT IMPLEMENTED. This is a Session 6A scaffolding stub.")
    print("On acquisition: verify healthy N + flag retrovirus biology per Issue 30.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
