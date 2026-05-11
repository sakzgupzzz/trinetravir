"""Acquire GSE283744 — pediatric RSV + SARS-CoV-2 + healthy cohort.

Session 6A Part A2 (SECONDARY held-out cohort). STUB.

Citation: Research Square 2025 "Comparative Single-Cell Analyses in Infants..."
GEO: GSE283744.

Design: 19 RSV-infected (mild=5, moderate=7, severe=7) + 30 SARS-CoV-2 + 17
healthy infants. Median age 2.3 months. 66 scRNA-seq + 51 snATAC-seq.

scRNA-seq is primary analysis modality; snATAC-seq is v1.5 multi-modal
extension.

Pre-flight verification: confirm publication state has stabilized (current
status is Research Square preprint; may have moved to final publication
by the time Session 6A runs).
"""

from __future__ import annotations

import sys

GEO_ACCESSION = "GSE283744"


def main() -> int:
    print(f"[STUB] acquire {GEO_ACCESSION}")
    print("  RSV n=19 (mild=5, moderate=7, severe=7) + SARS-CoV-2 n=30 + healthy n=17")
    print("  modality: scRNA primary; snATAC deferred")
    print("  storage: data/raw/gse283744/")
    print()
    print("ACTUAL DOWNLOAD NOT IMPLEMENTED. This is a Session 6A scaffolding stub.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
