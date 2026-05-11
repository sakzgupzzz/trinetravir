"""Run CellTypist Immune_All_High on the 4 v1 PBMC studies for Issue 12.

Produces <study>_reannotated_high.h5ad with coarse top-level labels
(Monocytes / B / T / etc.). The High model cannot resolve CD4 vs CD8
or NK; this is documented in METHODS_CHOICES Issue 12 and in the
bucket map module.

Order: wilk (smallest, validates pipeline), arunachalam, lee
(load-bearing diagnostic), schulte_schrepping (largest).
"""

from __future__ import annotations

import logging
from pathlib import Path

from trinetravir.data.annotate import annotate_and_save

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw"
PROCESSED = REPO / "data" / "processed"

STUDIES = ("wilk_2020", "arunachalam_2020", "lee_2020", "schulte_schrepping_2020")


def main() -> int:
    for study in STUDIES:
        out = PROCESSED / f"{study}_reannotated_high.h5ad"
        if out.exists():
            print(f"[{study}] already at {out.name}; skipping")
            continue
        print(f"\n========== {study} (Immune_All_High) ==========")
        annotate_and_save(
            raw_path=RAW / f"{study}.h5ad",
            out_path=out,
            model_name="Immune_All_High.pkl",
            majority_voting=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
