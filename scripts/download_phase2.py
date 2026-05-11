"""CLI orchestrator for Phase 2 PBMC dataset acquisition.

Reads ``configs/datasets.yaml``, downloads the requested Census studies via
``trinetravir.data.download``, persists each as ``data/raw/{study_id}.h5ad``,
and writes a provenance manifest.

Usage
-----
List configured studies and exit (no network):
    uv run python scripts/download_phase2.py --list

Dry run (resolve config + show plan, no downloads):
    uv run python scripts/download_phase2.py --all --dry-run

Download a single study:
    uv run python scripts/download_phase2.py --study-id wilk_2020

Download all configured Census studies:
    uv run python scripts/download_phase2.py --all
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from trinetravir.data.download import (
    StudyConfig,
    build_manifest,
    download_study,
    load_dataset_config,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = REPO_ROOT / "configs" / "datasets.yaml"
DEFAULT_OUT_DIR = REPO_ROOT / "data" / "raw"
DEFAULT_MANIFEST = REPO_ROOT / "data" / "raw" / "manifest.csv"


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="Path to datasets.yaml")
    p.add_argument(
        "--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Where to write {study_id}.h5ad"
    )
    p.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST, help="Manifest CSV output path"
    )
    p.add_argument(
        "--study-id", action="append", help="Specific study_id(s) to download. Repeatable."
    )
    p.add_argument(
        "--all",
        action="store_true",
        help="Download every Census study in the config (skips excluded).",
    )
    p.add_argument(
        "--include-excluded",
        action="store_true",
        help="When combined with --all, also download studies marked excluded in datasets.yaml.",
    )
    p.add_argument("--list", action="store_true", help="List configured studies and exit.")
    p.add_argument(
        "--dry-run", action="store_true", help="Resolve plan and report without downloading."
    )
    p.add_argument(
        "--overwrite", action="store_true", help="Re-download even if {study_id}.h5ad exists."
    )
    p.add_argument("--verbose", "-v", action="count", default=0, help="Increase log verbosity.")
    return p.parse_args(argv)


def _select_studies(
    studies: dict[str, StudyConfig],
    requested_ids: list[str] | None,
    all_flag: bool,
    include_excluded: bool = False,
) -> list[StudyConfig]:
    census_studies = {sid: s for sid, s in studies.items() if s.source == "cellxgene"}
    if all_flag:
        # --all skips excluded studies. Excluded studies are still available
        # via explicit --study-id (override path for ablation work).
        return [s for s in census_studies.values() if include_excluded or not s.excluded]
    if requested_ids:
        missing = [sid for sid in requested_ids if sid not in census_studies]
        if missing:
            raise SystemExit(f"Unknown or non-Census study_id(s): {missing}")
        return [census_studies[sid] for sid in requested_ids]
    raise SystemExit("Pass --study-id <id> (repeatable) or --all. Use --list to see options.")


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    log_level = {0: logging.WARNING, 1: logging.INFO}.get(args.verbose, logging.DEBUG)
    logging.basicConfig(level=log_level, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    log = logging.getLogger("download_phase2")

    defaults, studies = load_dataset_config(args.config)
    census_version = defaults["census_version"]

    if args.list:
        print(f"{'study_id':<28} {'source':<10} {'excluded':>8}  accession")
        for sid, s in studies.items():
            mark = "EXCL" if s.excluded else ""
            print(f"{sid:<28} {s.source:<10} {mark:>8}  {s.accession}")
        return 0

    selected = _select_studies(
        studies, args.study_id, args.all, include_excluded=args.include_excluded
    )
    log.info("Census version pinned: %s", census_version)
    log.info("Selected %d study(ies) for download", len(selected))
    for s in selected:
        log.info("  - %s  (%s)", s.study_id, s.accession)

    if args.dry_run:
        log.info("Dry run; exiting without downloading.")
        return 0

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for s in selected:
        log.info("=== %s ===", s.study_id)
        download_study(
            s,
            census_version=census_version,
            out_dir=args.out_dir,
            overwrite=args.overwrite,
        )

    manifest = build_manifest(args.out_dir)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.manifest, index=False)
    pd.set_option("display.max_colwidth", 60)
    print()
    print("Manifest written to", args.manifest)
    print(manifest.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
