"""Acquire Wang 2025 — chronic CMV carriage cohort.

Session 6A Part A3 (TERTIARY held-out cohort). STUB.

Citation: bioRxiv 2025.06.24.661167.
Design: 19 CMV(-) + 17 CMV(+) older adults, median age 71. 10x scRNA-seq +
flow cytometry. Discovery cohort.

Acquisition: GEO accession TBD (verify on download via project GitHub linked
in preprint).
"""

from __future__ import annotations

import sys


def main() -> int:
    print("[STUB] acquire Wang 2025 CMV cohort (GEO accession TBD)")
    print("  CMV+ n=17 / CMV- n=19 older adults; median age 71")
    print("  design: chronic latent carriage (asymptomatic)")
    print("  storage: data/raw/wang_2025_cmv/")
    print()
    print("ACTUAL DOWNLOAD NOT IMPLEMENTED. This is a Session 6A scaffolding stub.")
    print("On acquisition: verify chronic-carriage caveat in Issue 29.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
