"""One-shot migration: rename donor_disease_status value 'healthy' -> 'healthy_control'.

Resolves METHODS_CHOICES.md Issue 1 by re-writing on-disk h5ads to match the
updated download.py output vocabulary. Idempotent: skips files whose values
are already migrated.

Touches both:
  - data/raw/<study>.h5ad        (output of download_dataset)
  - data/processed/<study>_reannotated.h5ad  (Phase 3.5 reannotated copies)

The original Census source files (remote) are untouched.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

REPO = Path(__file__).resolve().parents[1]
TARGET_DIRS = [REPO / "data" / "raw", REPO / "data" / "processed"]

OLD_VALUE = "healthy"
NEW_VALUE = "healthy_control"
COLUMN = "donor_disease_status"


def migrate_file(path: Path) -> str:
    """Migrate one h5ad. Returns status string for logging."""
    a = ad.read_h5ad(path)
    if COLUMN not in a.obs.columns:
        return f"skip (no {COLUMN} column)"
    series = a.obs[COLUMN].astype(str)
    n_old = int((series == OLD_VALUE).sum())
    n_new = int((series == NEW_VALUE).sum())
    if n_old == 0:
        return f"skip (already migrated; {n_new} {NEW_VALUE!r} cells)"
    new_series = series.replace({OLD_VALUE: NEW_VALUE}).astype("category")
    a.obs[COLUMN] = new_series
    a.write_h5ad(path)
    return f"migrated {n_old} {OLD_VALUE!r} -> {NEW_VALUE!r}"


def main() -> int:
    paths: list[Path] = []
    for d in TARGET_DIRS:
        if not d.exists():
            continue
        paths.extend(sorted(d.glob("*.h5ad")))
    if not paths:
        print("No h5ad files found under data/raw or data/processed.")
        return 0
    for p in paths:
        try:
            status = migrate_file(p)
        except Exception as exc:  # noqa: BLE001
            print(f"[ERROR] {p.name}: {exc}")
            continue
        print(f"{p.relative_to(REPO)}: {status}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
