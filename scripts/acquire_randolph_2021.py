"""Acquire Randolph 2021 (GSE162632) — ex vivo IAV challenge cohort.

Session 6A Part A1 (PRIMARY held-out cohort).

This script is a STUB. Actual acquisition is multi-hour bg work that spans
multiple Claude Code sessions:

  1. Query GEO for GSE162632 supplementary files (10x .h5 or .mtx + barcodes
     + features per sample).
  2. Download per-sample files. ~90 donors × 2 conditions = 180 samples.
  3. Load each into per-sample AnnData, store under data/raw/randolph_2021/.
  4. Build per-donor donor metadata table (sex, ancestry, paired_id).
  5. Validate cell counts match published 235,161 high-quality / 255,731 raw.

Run as: uv run python scripts/acquire_randolph_2021.py [--sample-only N]
        (--sample-only restricts to N donors for development; default = all 90)

Expected wall time: 4-6 hours for full cohort acquisition.

TODO before launch:
  - Verify GEO supplementary file format (10x .h5 vs .mtx per sample).
  - Determine whether to consolidate all donors into one h5ad or keep
    per-sample files (recommend per-sample for memory safety on laptop).
  - Implement parallel download with retry-on-failure.
  - After acquisition: verify per-donor cell counts in cohort_qc_inventory.csv.
"""

from __future__ import annotations

import argparse
import sys

GEO_ACCESSION = "GSE162632"
CITATION = "Randolph HE et al. Science 374:1127-1133 (2021)"
DESIGN = (
    "90 male donors (EUR+AFR ancestry), paired mock + IAV Cal/04/09 (H1N1), 6h ex vivo, MOI 0.5"
)
EXPECTED_CELLS = 235161


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample-only", type=int, default=0, help="restrict to first N donors (0 = all 90)"
    )
    args = parser.parse_args()
    print(f"[STUB] acquire {GEO_ACCESSION} ({CITATION})")
    print(f"  design: {DESIGN}")
    print(f"  expected cells: {EXPECTED_CELLS} high-quality")
    print(f"  sample_only: {args.sample_only or 'all 90 donors'}")
    print("  storage: data/raw/randolph_2021/")
    print()
    print("ACTUAL DOWNLOAD NOT IMPLEMENTED. This is a Session 6A scaffolding stub.")
    print("Launch acquisition in a future session as a background job with proper")
    print("network monitoring + retry-on-failure. Estimated wall time: 4-6h.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
