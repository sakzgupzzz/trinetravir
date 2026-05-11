# Session 6A — held-out cohort acquisition + harmonization + schema v6 (Option B hybrid, part 1 of 2)

**Multi-session, 4-6 weeks total.** This checklist tracks state across multiple Claude Code sessions. Update status after each session.

## Pre-flight gate

- [x] Session 3 closeout complete (commit `0bdd2c5`)
- [x] Session 5 audit complete (commit `d91c24e`)
- [x] Issue 25 resolved as Option B hybrid (2026-05-11 human decision)
- [ ] Existing training corpus (Wilk + Lee + Arunachalam + Schulte) intact
- [ ] Audit gate at Part E **MUST PASS** before Session 6B begins

## Constraint

Acquisition + harmonization + schema work only. No Phase 4 / Phase 5 / Session 6B / Session 4 work. Stop at audit gate. Atomic schema-change rule (Issue 17) applies to all schema changes — migration script + tests + docs land in coordinated commits.

---

## Part A — held-out cohort acquisition (1.5-2 weeks; multi-session)

For each cohort: log source URL, accession, download size + timing, per-donor cell counts in `results/tables/cohort_inventory.csv`, donor metadata (age, sex, disease, severity, time-from-infection, treatment), acquisition issues.

### A1. Randolph 2021 (PRIMARY) — IAV ex vivo challenge
- [ ] Citation: Randolph HE et al., Science 374, 1127-1133 (2021).
- [ ] GEO: **GSE162632**
- [ ] Design: 90 male donors, paired mock + IAV Cal/04/09 (H1N1), 6h ex vivo, MOI 0.5, 10x.
- [ ] Cells: 235,161 high-quality / 255,731 raw.
- [ ] Sample size: vastly exceeds Issue 4 threshold.
- [ ] Acquisition: GEO direct download, 4-6h. **bg job needed**.
- [ ] Storage: `data/raw/randolph_2021/` + `data/processed/randolph_2021_*.h5ad`

### A2. GSE283744 (SECONDARY) — pediatric RSV + SARS-CoV-2
- [ ] Citation: Research Square 2025 "Comparative Single-Cell Analyses in Infants..." (verify final publication state)
- [ ] GEO: **GSE283744**
- [ ] Design: 19 RSV-infected (mild=5, moderate=7, severe=7) + 30 SARS-CoV-2 + 17 healthy infants. Median age 2.3 months. 66 scRNA-seq + 51 snATAC-seq.
- [ ] Modality: scRNA-seq primary; snATAC-seq deferred to v1.5.
- [ ] Sample size: 19 RSV + 17 healthy → passes Issue 4; 30 SARS-CoV-2 + 17 healthy → passes Issue 4.
- [ ] Acquisition: GEO direct, 4-6h. **bg job needed**.
- [ ] Storage: `data/raw/gse283744/` + `data/processed/gse283744_*.h5ad`

### A3. Wang 2025 (TERTIARY) — chronic CMV carriage
- [ ] Citation: bioRxiv 2025.06.24.661167 (verify final publication).
- [ ] Accession: GEO accession TBD (verify on download) + project GitHub linked in preprint.
- [ ] Design: 19 CMV(-) + 17 CMV(+) older adults, median age 71. 10x scRNA-seq + flow cytometry.
- [ ] Sample size: 17/19 meets Issue 4 with chronic-carriage caveat (see B3).
- [ ] Acquisition: 4-6h. **bg job needed**.
- [ ] Storage: `data/raw/wang_2025_cmv/` + `data/processed/wang_2025_cmv_*.h5ad`

### A4. Lee 2025 (TERTIARY) — HIV-1 early infection
- [ ] Citation: eLife 2025, PMC12370253.
- [ ] Design: 9 individuals with early HIV (<6 months), 5 scRNA-seq + 4 snRNA-seq multiome. 10x.
- [ ] Sample size: 9 HIV donors; healthy control N **TBD on download** — may not meet Issue 4. Fall back to qualitative validation if healthy N < 4. **Flag explicitly.**
- [ ] Acquisition: 3-4h. **bg job needed**.
- [ ] Storage: `data/raw/lee_2025_hiv/` + `data/processed/lee_2025_hiv_*.h5ad`

---

## Part B — schema v6 extensions (1 week)

Atomic per Issue 17: migration script + tests + docs in single commit.

- [x] **B1** `donor_response_design` + `exposure_pair_id` + `exposure_type` + `exposure_duration_hours` — for paired-within-donor designs (Randolph).
- [x] **B2** `age_years` + `age_group_category` — for pediatric stratification (GSE283744).
- [x] **B3** `infection_state` + `donor_serostatus` — for chronic latent carriage (Wang CMV).
- [x] **B4** retrovirus design — extend `exposure_type` with `retroviral_infection` value (Lee 2025 HIV). + `references/notes/hiv_biology_note.md` forward flag.
- [x] **B5** `src/trinetravir/data/schema_v6_migration.py` + `src/tests/test_schema_v6_migration.py`.
- [x] Existing v1 cohorts default-migrated: `donor_response_design=cross_sectional`, `exposure_type=natural_infection`, `age_group_category=adult`, `infection_state=acute` (Lee/Wilk/Arunachalam/Schulte are acute COVID), `donor_serostatus=unknown`.

---

## Part C — per-cohort harmonization (1.5-2 weeks; multi-session)

For each of the 4 cohorts:

- [ ] **C1** QC per Issue 4 thresholds. Document per-cohort pass rate in `results/tables/cohort_qc_inventory.csv`.
- [ ] **C2** CellTypist Immune_All_Low annotation (per Issue 12). For GSE283744 pediatric: verify CellTypist label accuracy against published pediatric PBMC composition.
- [ ] **C3** Apply `schema_v6_migration.py` to set new obs columns appropriately per cohort.
- [ ] **C4** Save per-cohort processed h5ads: `data/processed/<cohort>_processed_v6.h5ad`. Inventory in `results/tables/cohort_inventory_v6.csv`.

**Note**: cross-corpus harmonization (projecting held-out cohorts into training corpus integration space) is Session 6B work. Session 6A produces per-cohort processed h5ads ready for that projection.

---

## Part D — pre-specify Issues 27-30 (2-3 days)

- [x] **Issue 27** ex vivo vs natural infection comparison protocol (Randolph) — open at pre-specification; resolution at Session 6B.
- [x] **Issue 28** pediatric age stratification protocol (GSE283744) — open at pre-specification; resolution at Session 6B.
- [x] **Issue 29** chronic latent carriage analysis protocol (Wang CMV) — open at pre-specification; resolution at Session 6B.
- [x] **Issue 30** retrovirus context evaluation protocol (Lee 2025 HIV) — open at pre-specification; resolution at Session 6B.

---

## Part E — audit gate (STOP before Session 6B)

- [ ] Four cohorts acquired and inventoried per A1-A4.
- [ ] Schema v6 extensions added + tested per B1-B5.
- [ ] Schema migration script + tests committed atomically.
- [ ] Per-cohort QC + CellTypist + processed h5ads ready (Part C).
- [ ] Issues 27-30 opened in METHODS_CHOICES.md (Part D).
- [ ] Audit:
  - [ ] v1 `calibration_phase3_v2.csv` reproduces from migrated schema (regression test).
  - [ ] All held-out cohorts pass QC per Issue 4.
  - [ ] Sample sizes verified: Randolph ≥4/≥4, GSE283744 RSV ≥4/≥4, GSE283744 SARS ≥4/≥4, Wang CMV ≥4/≥4 (with chronic caveat), Lee 2025 HIV verified.
- [ ] **DO NOT auto-launch Session 6B.** Human audit confirmation required.

---

## Current state (2026-05-11, this session)

**Completed this session:**
- Schema v6 migration code (B1-B5)
- Pre-specify Issues 27-30 (Part D)
- SHORT_TERM_PLAN.md + SESSION_6A_CHECKLIST.md updates
- Acquisition scaffolding (data/raw/<cohort>/ directories + acquisition stub scripts ready to launch in future sessions)

**Deferred to future sessions:**
- Actual data acquisitions (Part A1-A4) — each is 3-6h of network + processing; requires background-job orchestration spanning multiple Claude Code sessions.
- Per-cohort harmonization (Part C1-C4) — depends on Part A completion.
- Audit gate (Part E) — depends on Parts A + C completion.

## Next-session priorities (for future Claude Code session)

1. Launch Randolph 2021 GEO download as bg job (highest priority — primary held-out cohort).
2. While Randolph downloads, launch GSE283744 in parallel if disk + RAM allow.
3. Run schema v6 migration tests against current v1 corpus to verify no regressions.
4. After each cohort downloads, run Part C harmonization steps for that cohort.
5. When all 4 cohorts processed, run audit gate verification and stop.
