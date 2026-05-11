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

**State (2026-05-11)**: original cohorts A2/A3/A4 BLOCKED. Substitutions made per Issues 28/29/30 cohort-substitution commits (2a0f322 / da40558 / 840f3db). New acquisition order below.

For each cohort: log source URL, accession, download size + timing, per-donor cell counts in `results/tables/cohort_inventory.csv`, donor metadata, acquisition issues.

### A1. Randolph 2021 (PRIMARY) — IAV ex vivo challenge — **ACQUIRED**
- [x] Citation: Randolph HE et al., Science 374, 1127-1133 (2021).
- [x] GEO: **GSE162632**
- [x] Design: 90 male donors, paired mock + IAV Cal/04/09 (H1N1), 6h ex vivo, MOI 0.5, 10x.
- [x] Cells: 235,161 high-quality / 255,731 raw.
- [x] Acquisition: 789MB GSE162632_RAW.tar downloaded + extracted to data/raw/randolph_2021/.
- [x] State: 30 multiplexed 10x pools on disk; awaiting Session 6B demultiplexing.

### A2. Yoshida 2022 (SUBSTITUTE for GSE283744) — pediatric + adult SARS-CoV-2 — **PENDING**
- [ ] Citation: Yoshida M et al. Nature 602:321 (2022).
- [ ] Source: **covid19cellatlas.org** direct h5ad (open-access Wellcome Sanger / HCA team deposition).
- [ ] Design: pediatric + adult + healthy acute primary SARS-CoV-2 PBMC. n=93 total. 10x 5' (matches Lee 2020 in v1 corpus). PBMC compartment: 317,854 cells.
- [ ] Sample size: pediatric + adult cohorts each have healthy controls; passes Issue 4 within each age stratum.
- [ ] Acquisition: covid19cellatlas.org direct h5ad download. Expected size: few GB.
- [ ] Storage: `data/raw/yoshida_2022/` + `data/processed/yoshida_2022_*.h5ad`

### A3. GSE213516 (SUBSTITUTE for Wang 2025 CMV) — chronic-latent CMV — **PENDING**
- [ ] Source: **GSE213516** (PBMC aging clocks cohort with CMV serostatus).
- [ ] Design: CMV(+) vs CMV(-) healthy older adults. 10x Genomics. Used by Science Advances aging clocks + npj Aging integrated atlas.
- [ ] Acquisition: GEO standard download.
- [ ] Storage: `data/raw/gse213516/` + `data/processed/gse213516_*.h5ad`
- [ ] Verify CMV serostatus annotation is in sample metadata for clean harmonization.

### A4. GSE157829 (SUBSTITUTE for Lee 2025 HIV) — chronic HIV exhaustion — **PENDING**
- [ ] Citation: Wang J et al. (2020) HIV exhaustion atlas (PMC7646563).
- [ ] Source: **GSE157829** public GEO.
- [ ] Design: 4 healthy + 6 HIV-infected (3 high VL + 3 low VL). 10x Genomics. ~66k PBMCs.
- [ ] Sample size: **meets Issue 4** (≥4 healthy + ≥4 diseased). No fallback to qualitative-only needed.
- [ ] Acquisition: GEO standard download.
- [ ] Storage: `data/raw/gse157829/` + `data/processed/gse157829_*.h5ad`

### A-extra. PBMCpedia investigation (30-min time-box, NON-BLOCKING)
- [ ] Check PBMCpedia (NAR Nov 2025, DOI 10.1093/nar/gkaf1245) for pre-harmonized versions of Yoshida 2022, GSE157829, GSE213516. If included in PBMCpedia's 24-study collection with standardized preprocessing, document alternative access path + note potential Part C time savings. Do NOT switch to PBMCpedia as primary access yet — requires verifying PBMCpedia preprocessing matches v1 methodology (CellTypist Immune_All_Low + per-cell-type Harmony).

### A-extra. Grabauskas/Ucar GEO probe (15-min, NON-BLOCKING)
- [ ] Search GEO for "Ucar" + "Grabauskas" CMV deposits. If Grabauskas 2025 Cohort 1 (19 CMV- + 17 CMV+) publicly available under different accession than the controlled-access pattern, document as POTENTIAL additional CMV cohort. If not findable, skip. GSE213516 remains primary substitute for Issue 29.

### A-blocked. Original cohorts (deferred; replaced)
- ❌ GSE283744 pediatric (controlled access; replaced by Yoshida 2022)
- ❌ Wang 2025 CMV Jackson Lab (likely controlled access; replaced by GSE213516)
- ❌ Lee 2025 HIV KRA (gated 2-4wk review + no healthy controls; replaced by GSE157829)

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

**Pre-conditions for next Claude session executing Part C** (added 2026-05-11):

### C-pre.1 Pre-specify Yoshida 2022 age stratification rule in Issue 28 BEFORE calibration
- [ ] **Primary stratification rule**: pediatric = age ≤ 18; adult = age 18-65 (matching v1 corpus). **Drop neonate (≤2 yr) and elderly (≥65 yr) from primary test.** Yoshida's `Age_group` column maps as: Neonate → drop; Infant (1-2 yr cutoff) → drop if Yoshida defines Infant as ≤2 yr; Young child / Child / Adolescent → pediatric; Adult → adult; Elderly → drop.
- [ ] Verify Yoshida's `Age_group` value boundaries via cellxgene metadata or paper supplement before applying the rule.
- [ ] **Supplementary sensitivity**: report alternative cutoffs (e.g., pediatric = ≤12, adult = 18-65; pediatric = ≤18, adult = ≤80) in supplementary table. Headline = the primary rule above.
- [ ] Commit Issue 28 age-stratification pre-specification BEFORE running any Yoshida calibration. Atomic per Issue 17.

### C-pre.2 GSE157829 C1 demultiplexing investigation
- [ ] Investigate paper supplementary / GitHub for C1 healthy-donor demultiplexing approach (hashtag IDs, scrublet calls, or pre-demuxed cell-donor map). PMC7646563.
- [ ] **If clean per-donor extraction is possible**: apply, get 4 healthy donors. Issue 4 satisfied directly.
- [ ] **If C1 is unsplittable pseudo-pool**: two paths:
  - (a) Substitute alternative HIV cohort (Cohn 2020 backup; verify accession before substituting). Repeat Issue 30 substitution commit pattern.
  - (b) Treat C1 as a single pseudo-donor for purposes of cross-cohort comparison. Document Issue 4 deviation rationale explicitly in METHODS_CHOICES Issue 30 amendment. Risk: pseudo-donor inflates sample-size appearance without true biological replication.
- [ ] Decision recorded in METHODS_CHOICES Issue 30 before any calibration.

### C-pre.3 GSE213516 CMV serostatus map (BEFORE calibration)
- [ ] Download Figshare 23789844 (cell type annotations) — verify whether CMV serostatus is in that annotation table OR in GEO Series Matrix file.
- [ ] Parse Series Matrix for sample-level CMV serostatus mapping. Build `data/raw/gse213516/cmv_serostatus_map.csv` linking each GSM sample to {CMV_positive, CMV_negative, unknown}.
- [ ] Document any sample-ID ambiguity in `cohort_inventory.csv` (rows where serostatus cannot be cleanly assigned).
- [ ] If <4 CMV+ or <4 CMV- donors after mapping, flag as Issue 4 deviation and consider substitution.

### C-pre.4 Apply v2 calibration framework per cohort
- [ ] Use Session 5 v2 framework: bootstrap CI direction fix (≥ lower CI bound, not in-CI), observed-r bootstrap CI, FDR-BH correction, MVS-restricted analysis per Khatri MVS gene set.
- [ ] Output: `results/tables/cohort_calibration_{cohort}_v2.csv` per cohort.

### C-pre.5 CellTypist annotation consistency check (Yoshida-specific)
- [ ] Yoshida 2022 h5ad ships with pre-annotated cell types. **Consistency check**: map Yoshida's annotations to v1's Immune_All_Low scheme. Document mismatches in `cohort_inventory.csv`.
- [ ] Decision: use Yoshida's annotations directly (faster, but may mismatch v1 buckets) OR re-run Immune_All_Low for cross-corpus consistency (slower, guarantees identical annotation pipeline). Default = re-run for cross-corpus consistency unless mismatch is < 5%.
- [ ] For other 3 cohorts (raw mtx/tsv): run CellTypist Immune_All_Low per Issue 12.

### C-pre.6 Standard Part C steps per cohort
- [ ] **C1** QC per Issue 4 thresholds. Document per-cohort pass rate in `results/tables/cohort_qc_inventory.csv`.
- [ ] **C2** CellTypist Immune_All_Low annotation (per C-pre.5 decision).
- [ ] **C3** Apply `schema_v6_migration.py` to set new obs columns appropriately per cohort:
  - Randolph: `donor_response_design = paired_within_donor`, `exposure_pair_id` per donor, `exposure_type = ex_vivo_challenge`, `exposure_duration_hours = 6.0`, `age_group_category = adult`.
  - Yoshida 2022: `age_years` populated from cellxgene metadata, `age_group_category` per C-pre.1 mapping, `exposure_type = natural_infection`, `infection_state = acute` for COVID-19 / `naive` for normal.
  - GSE157829: `exposure_type = retroviral_infection`, `infection_state = chronic_latent` for HIV+ donors, `age_group_category = adult`.
  - GSE213516: `exposure_type = chronic_carriage`, `infection_state = chronic_latent` for CMV+ / `naive` for CMV-, `donor_serostatus = positive` for CMV+ / `negative` for CMV-, `age_group_category = older_adult_>65yr`.
- [ ] **C4** Save per-cohort processed h5ads: `data/processed/<cohort>_processed_v6.h5ad`. Inventory in `results/tables/cohort_inventory_v6.csv`.
- [ ] **STOP** at end of Part C. Report state. Do not auto-launch Part E audit gate (user-confirmed gate).

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
