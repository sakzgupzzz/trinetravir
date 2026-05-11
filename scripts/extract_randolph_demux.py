"""Per-cell-type demux metadata extractor for Randolph 2021.

Reads ONE Seurat .rds file via rdata, dumps just the meta.data subset as CSV,
exits. Designed to be called as a subprocess per cell type so memory doesn't
accumulate across all 5 .rds reads (the CD4_T .rds is 10 GB compressed and
explodes to >30 GB in memory).

Usage:
  uv run python scripts/extract_randolph_demux.py {monocyte|B|NK|CD4T|CD8T}
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import rdata

warnings.filterwarnings("ignore", category=UserWarning, module="rdata")

REPO = Path(__file__).resolve().parents[1]
ZENODO = (
    REPO / "data" / "raw" / "randolph_2021" / "zenodo_inputs" / "inputs" / "1_calculate_pseudobulk"
)
OUT_DIR = REPO / "data" / "raw" / "randolph_2021"

RDS_BY_BUCKET = {
    "monocyte": ZENODO / "monocytes_cluster_singlets.rds",
    "B": ZENODO / "B_cluster_singlets.rds",
    "NK": ZENODO / "NK_cluster_singlets.rds",
    "CD4T": ZENODO / "CD4_T_cluster_singlets.rds",
    "CD8T": ZENODO / "CD8_T_cluster_singlets.rds",
}


def main() -> int:
    if len(sys.argv) != 2 or sys.argv[1] not in RDS_BY_BUCKET:
        print(f"usage: {sys.argv[0]} {{monocyte|B|NK|CD4T|CD8T}}", file=sys.stderr)
        return 2
    bucket = sys.argv[1]
    rds_path = RDS_BY_BUCKET[bucket]
    out_csv = OUT_DIR / f"demux_meta_{bucket}.csv"
    if out_csv.exists():
        print(f"[{bucket}] already exists at {out_csv.name}; skipping")
        return 0
    print(
        f"[{bucket}] reading {rds_path.name} ({rds_path.stat().st_size / 1e9:.1f} GB)", flush=True
    )
    parsed = rdata.parser.parse_file(str(rds_path))
    conv = rdata.conversion.convert(parsed)
    md = getattr(conv, "meta.data")
    md.columns = [str(c) for c in md.columns]
    keep_cols = [
        c
        for c in ["orig.ident", "MS_indiv_ID", "MS_infection_status", "celltype"]
        if c in md.columns
    ]
    sub = md[keep_cols].copy()
    sub["barcode"] = md.index.astype(str)
    sub["bucket"] = bucket
    sub.to_csv(out_csv, index=False)
    print(
        f"[{bucket}] wrote {out_csv.name}: {len(sub)} rows, {sub['MS_indiv_ID'].nunique()} donors",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
