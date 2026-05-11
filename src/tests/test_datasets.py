"""Tests for the dataset registry inclusion rule (METHODS_CHOICES Issue 4).

Validates that:
- Every included study satisfies the pre-specified donor-count minimum.
- Every excluded study fails it (consistent with the documented exclusion reason).
- The yaml retroactive_application table matches the inclusion/exclusion fields.

If a future study is added that violates the rule, this test fails informatively.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
CFG = REPO / "configs" / "datasets.yaml"


@pytest.fixture(scope="module")
def cfg() -> dict:
    with CFG.open() as f:
        return yaml.safe_load(f)


def test_inclusion_criteria_section_present(cfg: dict) -> None:
    """The criterion is documented in the registry, not in a separate hidden file."""
    assert "inclusion_criteria" in cfg, "inclusion_criteria section missing from datasets.yaml"
    ic = cfg["inclusion_criteria"]
    for key in (
        "min_healthy_donors",
        "min_diseased_donors",
        "rationale",
        "retroactive_application",
    ):
        assert key in ic, f"inclusion_criteria.{key} missing"
    assert ic["min_healthy_donors"] == 4
    assert ic["min_diseased_donors"] == 4


def test_retroactive_application_internally_consistent(cfg: dict) -> None:
    """Each entry's meets_criterion flag must match the donor counts."""
    ic = cfg["inclusion_criteria"]
    min_h = ic["min_healthy_donors"]
    min_d = ic["min_diseased_donors"]
    for row in ic["retroactive_application"]:
        expected = row["healthy_donors"] >= min_h and row["diseased_donors"] >= min_d
        assert row["meets_criterion"] == expected, (
            f"retroactive_application row for {row['study']!r} reports "
            f"meets_criterion={row['meets_criterion']!s} but donor counts "
            f"({row['healthy_donors']}H/{row['diseased_donors']}D vs minimum "
            f"{min_h}H/{min_d}D) imply {expected!s}"
        )


def test_excluded_studies_fail_criterion(cfg: dict) -> None:
    """Every excluded study must independently fail the criterion."""
    ic = cfg["inclusion_criteria"]
    retro = {row["study"]: row for row in ic["retroactive_application"]}
    for study_id, study in cfg["studies"].items():
        if not study.get("excluded"):
            continue
        assert study_id in retro, (
            f"excluded study {study_id!r} missing from inclusion_criteria.retroactive_application"
        )
        assert not retro[study_id]["meets_criterion"], (
            f"excluded study {study_id!r} actually satisfies the criterion "
            "per retroactive_application; either un-exclude it or fix the row"
        )


def test_included_studies_satisfy_criterion(cfg: dict) -> None:
    """Every non-excluded study in the registry must satisfy the criterion."""
    ic = cfg["inclusion_criteria"]
    retro = {row["study"]: row for row in ic["retroactive_application"]}
    for study_id, study in cfg["studies"].items():
        if study.get("excluded"):
            continue
        assert study_id in retro, (
            f"included study {study_id!r} missing from "
            f"inclusion_criteria.retroactive_application (add donor counts there "
            "so the inclusion rule is auditable)"
        )
        assert retro[study_id]["meets_criterion"], (
            f"included study {study_id!r} fails the donor-count criterion "
            f"({retro[study_id]['healthy_donors']}H/"
            f"{retro[study_id]['diseased_donors']}D) — either exclude it or "
            "revise the criterion. See METHODS_CHOICES Issue 4."
        )
