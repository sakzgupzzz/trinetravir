"""Session 4.5 Part D Step 6c: manuscript smoke-test.

Validates MANUSCRIPT_DRAFT.md after Steps 6a/6b/6d/6e propagation:
  1. All cited Issue numbers exist in METHODS_CHOICES.md
  2. All referenced result tables exist on disk
  3. All commit hashes cited are reachable
  4. No stale unresolved TBD/FIXME/placeholder markers (TBD entries in change log
     for future propagation steps are documented; flag if unexpected)
  5. No dangling N=200/N=500 references in active claims (legitimate
     historical/audit-trail mentions are expected)

Output: stdout summary. Exit 0 if all checks pass; exit 1 if any dangling refs.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MS = REPO / "MANUSCRIPT_DRAFT.md"
MC = REPO / "METHODS_CHOICES.md"
TABLES = REPO / "results" / "tables"


def main() -> int:
    ms = MS.read_text()
    mc = MC.read_text()

    failures = []

    # 1. Issue citation integrity
    issue_cites = {int(m.group(1)) for m in re.finditer(r"Issue (\d+)", ms)}
    issue_defs = {int(m.group(1)) for m in re.finditer(r"### Issue (\d+)", mc)}
    issue_resolved = {int(m.group(1)) for m in re.finditer(r"### Resolved Issue (\d+)", mc)}
    all_defs = issue_defs | issue_resolved
    dangling_issues = issue_cites - all_defs

    print("=== 1. Issue citation integrity ===")
    print(f"  Cited in manuscript: {sorted(issue_cites)}")
    print(f"  Defined in METHODS_CHOICES: {sorted(all_defs)}")
    if dangling_issues:
        print(f"  DANGLING (cited but not defined): {sorted(dangling_issues)}")
        failures.append(f"Dangling issue refs: {sorted(dangling_issues)}")
    else:
        print("  ALL CITATIONS RESOLVE ✓")

    # 2. Result table reference integrity
    table_refs = set(re.findall(r"results/tables/([a-zA-Z0-9_/]+\.csv)", ms))
    missing_tables = []
    for t in table_refs:
        p = TABLES / t
        if not p.exists():
            missing_tables.append(t)

    print("\n=== 2. Result table reference integrity ===")
    print(f"  Tables referenced: {len(table_refs)}")
    if missing_tables:
        print(f"  MISSING on disk: {missing_tables}")
        failures.append(f"Missing tables: {missing_tables}")
    else:
        print("  ALL TABLE REFS EXIST ✓")

    # 3. Commit hash references
    commit_refs = set(re.findall(r"commit `([0-9a-f]{7,40})`", ms))
    print("\n=== 3. Commit hash references ===")
    print(f"  Commit refs cited: {sorted(commit_refs)}")
    missing_commits = []
    for c in commit_refs:
        result = subprocess.run(
            ["git", "rev-parse", "--verify", c], capture_output=True, text=True, cwd=str(REPO)
        )
        if result.returncode != 0:
            missing_commits.append(c)
    if missing_commits:
        print(f"  UNREACHABLE commits: {missing_commits}")
        failures.append(f"Unreachable commits: {missing_commits}")
    else:
        print("  ALL COMMITS REACHABLE ✓")

    # 4. TBD/FIXME markers
    tbd_lines = [
        (i, line)
        for i, line in enumerate(ms.splitlines(), 1)
        if re.search(r"\bTBD\b|\bFIXME\b|<fill in>|\bXXX\b", line, re.IGNORECASE)
    ]
    print("\n=== 4. TBD/FIXME markers ===")
    print(f"  Found: {len(tbd_lines)}")
    # Expected TBDs are change-log entries for future propagation (Step 6c smoke-test note)
    # + section placeholders (Issue 18-24 methods detail TBD, model results TBD, etc.)
    unexpected_tbds = []
    # Pre-existing structural placeholders that are NOT Session 4.5 propagation defects:
    # authorship / target venue / corpus-cell-count to-fill / future-step / detail-pending markers.
    allowed_tbd_patterns = [
        re.compile(r"\| TBD "),  # change-log table placeholder rows
        re.compile(r"TBD details per Issues"),  # Methods detail placeholder
        re.compile(r"Author.*TBD|TBD.*Author|TBD.*co-author|TBD.*senior author"),
        re.compile(r"venue.*TBD|TBD.*venue|Target venue.*TBD"),
        re.compile(r"Total v1 corpus harmonized cells.*TBD"),
        re.compile(r"TBD.*discuss with Lilly|TBD.*collaborator"),
        re.compile(r"Authorship.*TBD|TBD.*Authorship"),
    ]
    for ln, line in tbd_lines:
        stripped = line.strip()
        if any(p.search(stripped) for p in allowed_tbd_patterns):
            continue
        unexpected_tbds.append((ln, stripped[:140]))
    if unexpected_tbds:
        print("  UNEXPECTED TBDs (not change-log/placeholder):")
        for ln, snippet in unexpected_tbds[:10]:
            print(f"    line {ln}: {snippet}")
        failures.append(f"Unexpected TBDs: {len(unexpected_tbds)}")
    else:
        print("  ALL TBDs are expected change-log/placeholder entries ✓")

    # 5. Stale N references in active claims (vs historical/audit-trail mentions)
    # Allow N=200 / N=500 in change-log entries (historical) and in Issue 38 reconciliation context
    print("\n=== 5. Stale N reference scan ===")
    stale_lines = []
    for ln, line in enumerate(ms.splitlines(), 1):
        if re.search(r"\bN=200\b|\bN=500\b|n_bootstrap=200|N_BOOTSTRAP=200", line):
            # Filter: change log + Issue 38 reconciliation context + Yoshida wide-CI caveat are legitimate
            if (
                "Step 6" in line
                or "change log" in line.lower()
                or "Issue 38" in line
                or "N=200 was over-confident" in line
                or "N=200 reading" in line
                or "N=200 bootstrap" in line  # Within historical context
                or "reconciliation" in line
                or "Session 6B" in line  # historical session reference
                or "2026-05-1" in line  # change-log date
            ):
                continue
            stale_lines.append((ln, line.strip()[:140]))
    if stale_lines:
        print("  STALE N references in active claims:")
        for ln, snippet in stale_lines[:10]:
            print(f"    line {ln}: {snippet}")
        failures.append(f"Stale N refs: {len(stale_lines)}")
    else:
        print("  No stale N=200/N=500 in active claims (historical/audit-trail mentions allowed) ✓")

    # Final verdict
    print("\n=== Smoke-test verdict ===")
    if failures:
        print(f"  FAIL: {len(failures)} issues")
        for f in failures:
            print(f"    - {f}")
        return 1
    print("  PASS ✓ All integrity checks clean.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
