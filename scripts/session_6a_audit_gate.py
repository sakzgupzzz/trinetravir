"""Session 6A Part E audit gate.

Verifies before Session 6B launch:
  (1) All 4 held-out cohorts processed v6 h5ads exist + load + have v6 obs cols.
  (2) Issue 4 sample sizes verified per cohort (with documented deviations).
  (3) Schema v6 migration test suite passes (62/62 tests).
  (4) v1 corpus still loads under schema v6 migration (regression check).
  (5) Existing v1 calibration outputs untouched.

Outputs audit summary to stdout; PASS/FAIL per check. Stops with non-zero
exit if any check fails.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import anndata as ad

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
RAW = REPO / "data" / "raw"
TABLES = REPO / "results" / "tables"

V6_COLS = (
    "donor_response_design",
    "exposure_pair_id",
    "exposure_type",
    "exposure_duration_hours",
    "age_years",
    "age_group_category",
    "infection_state",
    "donor_serostatus",
)

EXPECTED_COHORTS = {
    "yoshida_2022": {
        "file": PROC / "yoshida_2022_processed_v6.h5ad",
        "min_diseased": 4,
        "min_healthy": 4,
    },
    "allen_atlas_monocyte": {
        "file": PROC / "allen_atlas_monocyte_processed_v6.h5ad",
        "min_diseased": 4,
        "min_healthy": 4,
    },
    "gse157829": {
        "file": PROC / "gse157829_processed_v6.h5ad",
        "min_diseased": 4,
        "min_healthy": 1,  # cross-cohort design — Issue 30 amendment
        "deviation_note": "Issue 30 amendment 104a688: cross-cohort integration design per field precedent; n=1 healthy in GEO + v1 corpus 41 healthy as primary baseline",
    },
    "randolph_2021": {
        "file": PROC / "randolph_2021_processed_v6.h5ad",
        "min_diseased": 4,
        "min_healthy": 4,
    },
}

results = []


def check(name: str, passed: bool, detail: str) -> None:
    results.append((name, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}")


# ---- (1) Processed h5ads exist + load + v6 obs cols ----
print("\n=== (1) Processed v6 h5ads ===")
for cohort, cfg in EXPECTED_COHORTS.items():
    p = cfg["file"]
    if not p.exists():
        check(f"{cohort} file exists", False, f"missing {p}")
        continue
    try:
        a = ad.read_h5ad(p, backed="r")
    except Exception as e:
        check(f"{cohort} loads", False, f"load error: {e}")
        continue
    missing = [c for c in V6_COLS if c not in a.obs.columns]
    if missing:
        check(f"{cohort} v6 obs cols", False, f"missing cols: {missing}")
    else:
        check(f"{cohort} v6 obs cols", True, f"shape={a.shape}, all 8 v6 obs cols present")

# ---- (2) Issue 4 sample sizes ----
print("\n=== (2) Issue 4 sample sizes (≥4 diseased + ≥4 healthy, with documented deviations) ===")
for cohort, cfg in EXPECTED_COHORTS.items():
    p = cfg["file"]
    if not p.exists():
        check(f"{cohort} Issue 4", False, "file missing")
        continue
    a = ad.read_h5ad(p, backed="r")
    if "donor_disease_status" not in a.obs.columns or "donor_id" not in a.obs.columns:
        check(f"{cohort} Issue 4", False, "missing donor_disease_status or donor_id")
        continue
    df = a.obs[["donor_id", "donor_disease_status"]].drop_duplicates()
    dis = int((df["donor_disease_status"].astype(str) == "diseased").sum())
    hc = int((df["donor_disease_status"].astype(str) == "healthy_control").sum())
    ok = dis >= cfg["min_diseased"] and hc >= cfg["min_healthy"]
    detail = f"{dis} diseased + {hc} healthy_control donors"
    if "deviation_note" in cfg:
        detail += f" [{cfg['deviation_note']}]"
    check(f"{cohort} Issue 4", ok, detail)

# ---- (3) Schema v6 migration test suite ----
print("\n=== (3) pytest src/tests/ ===")
try:
    r = subprocess.run(
        ["uv", "run", "pytest", "src/tests/", "-q"],
        cwd=str(REPO),
        capture_output=True,
        text=True,
        timeout=120,
    )
    last_line = (
        r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr.strip().splitlines()[-1]
    )
    check("pytest src/tests/", r.returncode == 0, last_line)
except Exception as e:
    check("pytest src/tests/", False, f"exception: {e}")

# ---- (4) v1 corpus still loads under v6 migration ----
print("\n=== (4) v1 corpus schema_v6_migration regression check ===")
try:
    from trinetravir.data.schema_v6_migration import has_v6_schema, migrate_v1_to_v6

    # Use a v1 corpus h5ad (lee_2020 reannotated_low.h5ad)
    v1_p = PROC / "lee_2020_reannotated_low.h5ad"
    if not v1_p.exists():
        check("v1 corpus migration", False, f"missing {v1_p}")
    else:
        a = ad.read_h5ad(v1_p, backed="r")
        # Copy to in-memory + migrate
        a_mem = a.to_memory()
        a.file.close()
        m = migrate_v1_to_v6(a_mem)
        ok = has_v6_schema(m) and m.n_obs == a_mem.n_obs
        check(
            "v1 corpus migration", ok, f"v1 corpus shape={m.shape}; v6 schema applied; no row loss"
        )
except Exception as e:
    check("v1 corpus migration", False, f"exception: {e}")

# ---- (5) v1 calibration outputs untouched ----
print("\n=== (5) v1 calibration outputs intact ===")
v1_cals = [
    "calibration_phase3.csv",
    "calibration_phase3_v2.csv",
    "calibration_phase35_low.csv",
    "calibration_phase35_low_v2.csv",
]
for f in v1_cals:
    p = TABLES / f
    if p.exists() and p.stat().st_size > 0:
        check(f"{f}", True, f"{p.stat().st_size} bytes")
    else:
        check(f"{f}", False, "missing or empty")

# ---- Summary ----
n_pass = sum(1 for _, ok, _ in results if ok)
n_fail = sum(1 for _, ok, _ in results if not ok)
print(f"\n=== Audit summary: {n_pass} PASS / {n_fail} FAIL ===")
if n_fail:
    print("\nFAILED checks:")
    for name, ok, detail in results:
        if not ok:
            print(f"  - {name}: {detail}")
sys.exit(0 if n_fail == 0 else 1)
