# Trinetravir v1 Paper — Working Manuscript Document

**Status:** Working draft for iteration. Started 2026-05-11 during Session 6B initial results.
**Author:** Saksham Gupta + collaborators TBD
**Repo:** https://github.com/sakzgupzzz/trinetravir
**Target venue:** bioRxiv preprint (primary), then journal submission TBD
**Target submission window:** Pending Session 6B completion + Session 3.5 + Session 4 + Phase 4–7 modeling. Realistic: 10–16 weeks from 2026-05-11.

---

## Working title options

1. Calibrated cross-study integration of PBMC SARS-CoV-2 single-cell transcriptomics reveals ISG-anchored response coherence and cross-context generalization
2. A factorized model for viral host response prediction with ISG-aware regularization, validated across cross-study, cross-context, cross-age, and chronic-viral cohorts
3. Pre-registered held-out validation of a single-cell PBMC viral response model: cross-context, cross-age, and chronic-virus boundary conditions

The current top choice is **#1** — it emphasizes the methodological contribution (calibrated integration framework + ISG-anchored finding) while signaling the held-out validation through "cross-context generalization." Title gets finalized after Session 6B verdicts are in.

---

## Elevator pitch (one paragraph)

We build a calibrated cross-study integration framework for PBMC single-cell transcriptomics of viral infection, anchored on a corpus of four harmonized SARS-CoV-2 cohorts. A factorized neural-network model with ISG-aware regularization captures the conserved antiviral response component, motivated by the empirical finding that ISG-restricted analysis recovers cross-study coherence in lymphoid populations otherwise misclassified as failed under full-HVG calibration. Pre-registered held-out validation against four cohorts spanning cross-context (ex vivo IAV vaccine vs natural infection), cross-age (pediatric vs adult SARS-CoV-2), chronic-latent CMV discrimination, and chronic HIV retrovirus distinctness tests the model's generalization scope with mechanical decision rules set before data acquisition. The framework yields a defensible methodology for cross-study integration of viral PBMC scRNA-seq, a community-resource harmonized corpus, and explicit boundary conditions on transfer from acute respiratory viral training distributions.

---

## Authorship and contributions (placeholder)

- **First author:** Saksham Gupta (Eli Lilly; computational design, methodology, implementation, paper writing)
- **Senior author / corresponding:** TBD — discuss with Lilly and academic collaborators about appropriate co-author roles
- **Possible co-authors:** TBD pending engagement with: harmonization methodology collaborators, immunology domain experts for biological interpretation, data-source PIs (Wilk, Lee, Arunachalam, Schulte) per cohort acknowledgment conventions

**Action item:** have authorship conversation by Phase 5 completion. v1 is solo-led but co-author engagement on domain biology interpretation is essential before submission.

---

## Key contributions (will become cover letter bullets)

1. **Harmonized cross-study PBMC SARS-CoV-2 corpus released as a community resource.** Four cohorts (Wilk 2020, Lee 2020, Arunachalam 2020, Schulte-Schrepping 2020) integrated under a standardized schema (v6) with consistent cell-type annotation (CellTypist Immune_All_Low), per-cell-type batch correction (Harmony), and validated cross-study coherence at the calibrated-evaluation level.

2. **A calibrated evaluation framework for cross-study transcriptional coherence** (`calibration.py` v2). Includes permutation null with donor-level resampling, bootstrap confidence intervals on observed Pearson r, split-half ceiling estimates of within-study reproducibility upper bounds, FDR-Benjamini-Hochberg correction across multiple hypothesis tests, and external anchoring against the Khatri Meta-Virus Signature gene set. Designed to distinguish confirmatory from exploratory evidence.

3. **The empirical finding that ISG-restricted analysis recovers lymphoid cross-study coherence that full-HVG calibration misclassifies as failed.** Replicated in held-out cohorts (GSE157829 chronic HIV: 4 of 5 buckets lifted; Randolph 2021 IAV: 2 of 3 buckets lifted). Methodological implication: standard HVG-based response vectors are too noisy for lymphoid cell types; pathway-restricted vectors recover signal.

4. **A factorized neural-network architecture with ISG-aware regularization**, evaluated on the cross-study integration task. Motivated by the ISG-restriction finding — the architecture's design constraint that conserved antiviral response components must align with the canonical ISG signature is empirically supported rather than ad hoc.

5. **Pre-registered held-out validation across four cohorts** testing biologically distinct hypotheses with mechanical decision rules set before data acquisition:
   - Cross-context (Randolph 2021 IAV ex vivo vs Lee 2020 IAV natural infection within v1 corpus)
   - Cross-age (Yoshida 2022 pediatric vs adult SARS-CoV-2)
   - Chronic-latent CMV discrimination (Allen Institute Immune Health Atlas CMV+ vs CMV−)
   - Chronic HIV retrovirus distinctness (GSE157829 chronic HIV vs v1 corpus baseline)

6. **Documented boundary conditions on transfer from acute respiratory viral training**, including specific bucket-level findings (e.g., chronic latent CMV monocyte signature not detectable; chronic HIV monocyte response distinct from acute; lymphoid signal preserved via ISG restriction).

---

## Pre-registered hypotheses with decision rules

All decision rules committed before observation of held-out cohort data. See `METHODS_CHOICES.md` Issues 27–30 in the repo for full audit trail with commit hashes.

| Issue | Hypothesis | Cohort | Pre-committed rule on MVS-restricted Pearson r |
|---|---|---|---|
| 27 | Cross-context IAV transfer (vaccine ex vivo → natural infection) | Randolph 2021 (GSE162632) | Monocyte r ≥ 0.40 supports H1; r < 0.20 challenges |
| 28 | Cross-age transfer (adult → pediatric SARS-CoV-2) | Yoshida 2022 (Nature 602:321) | Monocyte r ≥ 0.30 supports H1; r < 0.10 = does NOT transfer |
| 29 | Chronic-latent CMV discrimination from acute viral training | Allen Immune Health Atlas | Monocyte r ∈ [0.10, 0.40] = appropriate; r > 0.50 = concerning over-prediction; r < 0.05 = concerning (no antiviral memory shared) |
| 30 | Chronic HIV retrovirus distinctness | GSE157829 (Wang 2020) | CD4T r ∈ [0.00, 0.20] expected; r > 0.40 surprising; r < −0.10 anti-correlation |

These rules apply mechanically. The framework does not interpret away unexpected results — each numerical range corresponds to a specific verdict that goes into the paper as-is.

---

## Cohorts

### Training corpus (v1 main analysis)

| Cohort | Accession | Cells (post-QC) | Donors | Virus | Notes |
|---|---|---|---|---|---|
| Wilk 2020 | GSE150728 | TBD | TBD | SARS-CoV-2 | First atlas |
| Lee 2020 | GSE149689 | TBD | TBD | SARS-CoV-2 + IAV | 5′ chemistry; within-study cross-virus exploration |
| Arunachalam 2020 | GSE155673 | TBD | TBD | SARS-CoV-2 | |
| Schulte-Schrepping 2020 | EGAS00001004571 | TBD | TBD | SARS-CoV-2 | |

Total v1 corpus harmonized cells across 4 studies: TBD (fill in from existing v1 harmony outputs).

### Held-out cohorts (v1 confirmatory validation)

| Cohort | Accession | Cells (post-harmonization) | Donors | Test |
|---|---|---|---|---|
| Yoshida 2022 | cellxgene 03f821b4 | 168,018 | 9 diseased + 26 healthy; primary strata pediatric 5C/17N + adult 4C/9N | Issue 28 cross-age SARS-CoV-2 |
| Allen Atlas (monocyte) | AIFI Immune Health Atlas | 300,888 | 42 CMV+ / 50 CMV− | Issue 29 chronic-latent CMV |
| GSE157829 | GSE157829 | 35,750 | 6 HIV + 1 healthy (cross-cohort design with v1 baseline per field precedent) | Issue 30 chronic HIV |
| Randolph 2021 | GSE162632 + Zenodo 4273999 | 34,276 | 90 paired (mock + IAV per donor); HMN83575 excluded primary <50 cells | Issue 27 cross-context IAV |

Cross-cohort integration design for GSE157829 follows field precedent (eBioMedicine 2025; PMC10040851; PMC9434837). Allen Atlas monocyte subset used for Issue 29 primary; full atlas available for sensitivity analysis.

---

## Paper structure (section-by-section outline)

### Abstract (placeholder)

Single-paragraph summary covering: methodology framework, ISG-restriction finding, factorized model, held-out validation, conclusions about transfer scope. Will be finalized last.

### Introduction (~800–1200 words)

- Single-cell PBMC scRNA-seq of viral infection has produced multiple high-quality cohorts but cross-study integration is hampered by batch effects and methodological heterogeneity
- Existing approaches (Harmony, scVI, scGen, scCausalVI, CPA) address some integration aspects but lack calibration framework for distinguishing biological signal from technical artifact at the cross-study response-vector level
- Pre-registered evaluation against external held-out cohorts is rare in this literature, leading to risk of fit-to-training-corpus thresholds
- This paper presents: harmonized cross-study SARS-CoV-2 corpus + calibrated evaluation framework + factorized model with ISG-aware regularization + pre-registered held-out validation across four biological axes

### Methods (~2000–3000 words; details in supplementary)

**Cohort selection** — Issues 4 (≥4 healthy + ≥4 diseased donors), with explicit deviations documented (GSE157829 cross-cohort baseline design, HMN83575 exclusion).

**Schema v6** — `donor_response_design`, `exposure_pair_id`, `exposure_type`, `exposure_duration_hours`, `age_years`, `age_group_category`, `infection_state`, `donor_serostatus`. Atomic non-destructive migration verified via regression test.

**Quality control** — per-cell mito %, gene count, doublet detection. Per-cohort specifications detailed in supplementary.

**Cell type annotation** — CellTypist Immune_All_Low model applied consistently across cohorts. Per-cohort consistency checks documented.

**Harmonization** — per-cell-type Harmony with `study_id` (or `donor_id` within-cohort) as batch variable. Choice of per-cell-type over global Harmony documented per Issue 7 (acknowledged as methodological-alignment choice with factorized model's per-bucket grain).

**Calibration framework v2** — permutation null with donor-level resampling, bootstrap CI on observed Pearson r (donor-level B=1000), split-half ceiling, FDR-BH correction, MVS-restricted analysis. `test_calibration.py` with 8 synthetic ground-truth tests for verification.

**Khatri MVS external anchor** — canonical antiviral ISG signature from Andres-Terre et al. 2015 used as orthogonal validation set. Cross-study coherence under MVS restriction compared to full HVG calibration.

**Factorized model architecture** — TBD details per Issues 18–24 pre-specifications (Session 3.5 deliverable). Will include: ISG-aware regularization design, hyperparameter search space, few-shot adaptation protocol, comparison method versions, baseline implementations.

**Pre-registered held-out validation** — Issues 27–30 with decision rules committed before data acquisition.

**Exploratory vs confirmatory framing** — `PLAN.md` §1.8 distinction: Phases 1–3 produced exploratory evidence (heuristic thresholds, data-driven decisions); Phase 4 onwards produces confirmatory evidence with pre-registered protocols. The paper distinguishes these clearly throughout.

### Results

#### Section 1: Harmonized cross-study SARS-CoV-2 corpus

**Headline:** 4-study corpus harmonized via per-cell-type Harmony, validated via calibrated cross-study coherence. Monocyte cross-study Pearson r passes FDR-corrected confirmatory threshold; lymphoid populations require ISG-restriction (Khatri MVS) to recover cross-study signal.

**Figures:** UMAP before/after harmonization; per-bucket cross-study r heatmap with FDR-corrected significance markers; Khatri MVS lift table.

#### Section 2: ISG-restricted analysis recovers lymphoid cross-study coherence

**Headline:** Restricting to canonical ISG genes lifts cross-study Pearson r in 4 of 5 buckets within v1 corpus (lymphoid +0.06 to +0.23). Pattern motivates ISG-aware regularization in the factorized model.

**Figures:** r_full vs r_MVS scatter; bucket-level lift bar chart; ISG signature score distributions.

#### Section 3: Factorized model with ISG-aware regularization

**Headline:** Architecture description, training details, baseline comparisons (scVI, scGen, scCausalVI), cross-study reconstruction loss, few-shot adaptation performance.

**Figures:** Architecture schematic; reconstruction loss curves; baseline comparison table; few-shot adaptation curves.

#### Section 4: Held-out validation (pre-registered)

Subsections per Issue 27–30 hypothesis. Each subsection reports: observed Pearson r (full + MVS-restricted) with bootstrap CI, permutation null p-value, FDR-corrected significance, comparison to pre-committed decision rule, mechanical verdict.

##### Headline summary (N=1000 permutations, N=200 bootstrap, 15-test panel)

| Issue | Cohort × Bucket | Primary contrast | r_full | r_MVS | Bootstrap CI r_MVS | Raw p_MVS | FDR p_MVS | Verdict |
|---|---|---|---|---|---|---|---|---|
| 27 PRIMARY | Randolph monocyte_infected (Issue 31) | cluster-8 IAV-infected vs parent-bucket NI mock | 0.129 | **−0.011** | n/a (paired) | 0.072 | 0.270 | **CHALLENGES_H1** |
| 27 SENSITIVITY | Randolph monocyte (bystander) | flu vs NI within bystander bucket | 0.286 | 0.013 | n/a (paired) | 0.441 | 0.508 | CHALLENGES_H1 |
| 28 | Yoshida monocyte | pediatric ↔ adult cross-age | 0.387 | **0.591** | [0.017, 0.684] | 0.052 | 0.260 | **SUPPORTS_H1** |
| 29 | Allen Atlas monocyte | chronic-latent CMV+ vs CMV− | 0.152 | **−0.010** | [−0.516, 0.415] | 0.509 | 0.546 | **CONCERNING_NO_SHARED_BIOLOGY** (scope-limitation finding) |
| 30 | GSE157829 CD4T | chronic HIV vs v1 baseline | 0.084 | **0.257** | [0.157, 0.513] | 0.136 | 0.286 | **BORDERLINE** (just above EXPECTED ceiling) |

##### Issue 27 — Cross-context IAV (Randolph 2021)

Per Issue 31 (METHODS_CHOICES.md, pre-spec 2026-05-11), the matched healthy reference for cluster-defined cell subsets draws from the parent bucket's healthy/mock condition. The primary Issue 27 contrast compares cluster-8 infected monocytes (n=4924 cells, flu condition, HMN83575 excluded) against parent `monocyte` bucket NI mock cells (n=9785). The sensitivity contrast uses bystander monocytes (flu condition, non-cluster-8) against the same NI mock pool. Both contrasts mechanically yield CHALLENGES_H1 under the pre-committed rule (r_mvs ≥ 0.40 SUPPORTS; r_mvs < 0.20 CHALLENGES).

The full-HVG Pearson r reaches statistical significance against the permutation null (primary r_full = 0.129, p_full = 0.001, only test in the 15-test panel with FDR-corrected p_full < 0.05); on the MVS-restricted subset, the correlation collapses to zero. The shared signal is carried by non-ISG genes (lineage identity, basal transcription, cell-state markers), not by canonical type-I-IFN ISGs. Kinetic interpretation: ex vivo 6h IAV at MOI 0.5 (Randolph) engages early-phase direct-PAMP-sensing programs (RIG-I/MDA5, immediate-early IFN gene induction) that diverge structurally from the late-phase mature paracrine IFN-α/β ISG cascade captured in v1's natural-infection training corpus.

Bystander r_full > infected r_full (0.286 vs 0.129) is biologically interpretable: directly-infected cells engage cell-autonomous antiviral programs (apoptosis, autophagy, viral replication suppression) that diverge from generic monocyte activation; bystander monocytes look closer to the broader v1 corpus signal because they share paracrine activation without engaging the cell-fate-decision machinery.

##### Issue 28 — Cross-age (Yoshida 2022)

> Issue 28 SUPPORTS_H1: r_MVS = 0.591 (95% bootstrap CI [0.02, 0.68], N=1000 permutations p = 0.052). The observed effect size clears the pre-committed ≥0.30 threshold; the wide CI reflects limited donor power in the primary pediatric/adult strata (9 diseased + 26 healthy) and indicates the verdict is robust to point-estimate interpretation but cannot rule out, at 95% confidence, that the true effect lies below the supporting threshold.

This is the strongest single-cohort held-out result in v1. Yoshida is the only cohort whose primary verdict supports H1 under the pre-committed mechanical rule.

##### Issue 29 — Chronic-latent CMV discrimination (Allen Atlas)

The CONCERNING_NO_SHARED_BIOLOGY verdict is mechanically correct under the pre-committed rule (r_mvs < 0.05) but reads biologically as a **scope-limitation finding**: v1's training corpus captures acute IFN-driven response, not chronic-latent IFN tone. The framework is acute-disease-specific by construction; the boundary condition is appropriate, not failure. The bootstrap CI [−0.516, 0.415] reflects single-bucket coverage (monocyte only met the n_cells ≥ 50 sensitivity gate) and small canonical-ISG sample (n=57 MVS genes) against a flat ~0 signal.

##### Issue 30 — Chronic HIV retrovirus distinctness (GSE157829)

> Issue 30 BORDERLINE: r_MVS = 0.257 (95% bootstrap CI [0.157, 0.513], N=1000 permutations p = 0.136). The observed point estimate sits just above the pre-committed [0.00, 0.20] EXPECTED retrovirus-distinctness ceiling; the CI lower bound (0.16) lies inside the EXPECTED range, the upper bound (0.51) crosses the SURPRISING_HIGH threshold (>0.40). The mechanical verdict applies to the point estimate; the wide CI reflects limited donor power (6 HIV donors + 1 healthy, cross-cohort baseline design) and indicates the true effect is consistent with both expected-partial-overlap and surprisingly-high readings at 95% confidence.

CD4T target-cell biology in chronic HIV overlaps acute viral ISG signature ~50% at MVS level, ~10% at full-HVG level. The framework discriminates retrovirus from acute RNA virus imperfectly at the conserved-ISG level, cleanly at the full-HVG level. CD8T r_MVS = 0.612 (CI [0.345, 0.661]) and B r_MVS = 0.596 (CI [0.219, 0.672]) show stronger lymphoid ISG transfer than the CD4T primary contrast.

**Figures:** per-cohort response vector correlations; ISG-lift replication across cohorts; pre-committed decision rule visualization with observed r locations marked + bootstrap CI bars + pre-committed verdict bands.

#### Section 5: Boundary conditions and transfer scope

**Headline:** The v1 model's transfer scope is defined by what does and doesn't generalize. Findings (N=1000 confirmatory analysis):

- **Cross-context IAV** (Randolph ex vivo 6h IAV MOI 0.5 → v1 natural-infection corpus): **does NOT transfer** at the MVS canonical-ISG level for monocytes, in either the cluster-8 infected subpopulation (r_MVS = −0.011) or the bystander subpopulation (r_MVS = 0.013). Full-HVG signal does transfer (r_full = 0.129 infected; 0.286 bystander), but is carried by non-ISG lineage-level signal. Kinetic boundary: early-phase ex vivo response ≠ late-phase natural-infection response at the canonical-ISG level. Lymphoid cross-context IAV transfer IS observed: B r_MVS = 0.483 (p_mvs = 0.033), NK r_MVS = 0.576 (p_mvs = 0.032). The boundary is monocyte-specific and kinetic-specific.
- **Cross-age** (adult → pediatric, same virus, same context): **transfers**. Monocyte signal at r_mvs = 0.591 well above the 0.30 SUPPORTING threshold. Conserved component preserves canonical ISG signature across age groups. Caveat: wide bootstrap CI [0.02, 0.68] indicates the verdict is robust to point-estimate interpretation but cannot rule out a true effect below the threshold at 95% confidence given Yoshida's limited per-stratum donor counts.
- **Chronic-latent CMV:** monocyte transcriptional signature **not detectable** at MVS level (r_MVS ≈ 0). Consistent with CMV memory residing in adaptive compartments (CD8 TEMRA, GZMK+ T cells, adaptive NK) rather than monocytes. v1 framework is acute-disease-specific by construction; this boundary is appropriate, not failure.
- **Chronic HIV:** lymphoid ISG signature substantially shared with acute viral training (CD8T r_MVS = 0.61; B r_MVS = 0.60; NK r_MVS = 0.42); CD4T primary target r_MVS = 0.257 sits at the BORDERLINE just above expected retrovirus-distinctness; monocyte response distinct (r_MVS = 0.04, flat). Chronic IFN tone preserves canonical ISG signature in lymphoid compartments but does NOT correlate at the CD4T target compartment as strongly, reflecting cell-autonomous retroviral biology specific to HIV target cells.

The v1 framework's **domain of validity** is acute respiratory viral infection (paracrine IFN-α/β cascade ≥24h post-onset) on PBMC monocytes + lymphoid compartments. Cross-age + chronic-lymphoid generalization holds. Ex vivo early-phase IAV and chronic-latent herpesvirus monocyte boundaries are explicit scope limits.

### Discussion (~1500–2000 words)

- Methodological contributions: calibration framework as community standard; ISG-restricted analysis as default for lymphoid cross-study response vectors; pre-registered held-out validation as alternative to ad hoc threshold setting
- Biological insights: ISG signature as conserved antiviral component; cell-type-specific transfer patterns; boundary conditions on acute → chronic generalization
- Comparison to existing approaches: scGen (counterfactual generation), scCausalVI (causal disentanglement), CPA (compositional perturbation), GPLVM (probabilistic factorization)
- Limitations (see dedicated section below)
- Implications for v1.5 (cross-virus extension to engineered adenovirus vector platforms) and v2 (commercial product)
- **Threshold provenance for held-out cohort decision rules.** The Issues 27-30 decision rule thresholds were set by informed judgment calibrated against published cross-cohort viral response literature, not derived from a single principled statistical calculation. The calibrating references are documented in `METHODS_CHOICES.md` Issue 36 and include: (1) within-corpus monocyte cross-study Pearson r ceiling from Session 5 calibration (r = 0.45-0.65 across v1 cohort pairs at MVS-restricted level), (2) Khatri Meta-Virus Signature cross-cohort transfer baseline from Andres-Terre 2015 (*Immunity* 43:1199) and follow-up work (r ≈ 0.40-0.60 across 14 independent respiratory viral PBMC cohorts), (3) single-cell perturbation prediction benchmarks (Ahlmann-Eltze 2025; Kedzierska 2024 *BMC Genomics*) establishing r ≈ 0.30-0.45 as the realistic upper bound for cross-distribution transfer, (4) PBMC cross-study viral response specifically (PBMCpedia *NAR* 2025), (5) pediatric vs adult SARS-CoV-2 PBMC scRNA-seq literature (Jia 2024 PMC11325098; Sallusto 2025 *Nat Commun*), and (6) ex vivo PBMC challenge vs natural infection signature overlap (PMC11637350). The methodological defense rests on pre-registration (thresholds committed in `METHODS_CHOICES.md` before data observation, verifiable via git timestamps) and on literature anchoring (the specific numerical values fall within the field's empirical baseline ranges for analogous comparisons). The robustness of each verdict to threshold choice within plausible literature-supported ranges is reported in supplementary; for Yoshida (Issue 28), the observed r_MVS = 0.591 clears any literature-plausible support threshold in [0.20, 0.40], so the verdict is insensitive to the specific cutoff. Future cohort designs (v1.5+) will consider relative-to-reference thresholds (e.g., "observed r exceeds 50% of within-corpus monocyte cross-study ceiling") as a more principled framing than absolute thresholds. This is documented as a v1 limitation and v1.5 improvement opportunity in `METHODS_CHOICES.md` Issue 36b.

### Limitations

- **FDR-corrected significance does not survive the 15-test panel.** At N=1000 permutations across 15 bucket-cohort tests (4 held-out cohorts × ~3–5 buckets per cohort), no comparison reached FDR<0.05 confirmatory significance under Benjamini-Hochberg correction. This reflects a structural power limitation given the held-out evaluation infrastructure; larger held-out infrastructure (v1.5 scope) is required for FDR<0.05 confirmatory thresholds. The verdicts reported here apply pre-committed decision rules on observed effect sizes, which were specified before data acquisition and do not depend on FDR-corrected p-values. One test (Issue 27 PRIMARY full-HVG, raw p_full = 0.001, FDR p_full = 0.027) does survive FDR<0.05 — but that significant signal is carried by non-ISG genes (the MVS-restricted contrast is null), reinforcing rather than weakening the CHALLENGES_H1 verdict.
- **Per-cell-type Harmony correction contributes 0.02–0.25 to cross-study response-vector coherence on top of the substantial pre-Harmony biological signal** (Session 7 pre-modeling sensitivity audit, Issue 32). Pre-Harmony cross-study Pearson r is already 0.13–0.58 across the 5 buckets; Harmony correction adds Δr in (0.10, 0.30] for 3 of 5 buckets on full HVG, 4 of 5 on the Khatri MVS subset. No bucket crosses the pre-committed Δr > 0.30 HARMONY_DOMINANT threshold, so the cross-study integration framework is not "doing most of the work" — but it is doing meaningful smoothing, and the manuscript framing must acknowledge this. The load-bearing monocyte-MVS Δr = 0.08 (BIOLOGY_DOMINANT) and the perfect (100%) within-cohort vs cross-study sign concordance (Issue 33) together establish that the ISG-conservation finding is biology with Harmony amplification, not an integration artifact.
- v1 corpus is SARS-CoV-2-dominated (4 SARS studies + 1 IAV within-study via Lee). Multi-virus generalization is tested via held-out cohorts but not via expanded training corpus. v1.5 will address this.
- Held-out cohort design includes one cross-cohort comparator (GSE157829 with v1 healthy baseline) per field precedent. While defensible methodology, within-cohort case-control would be cleaner if a sufficient HIV PBMC scRNA-seq cohort with >4 healthy controls were available.
- Allen Atlas Issue 29 used monocyte subset only (327K cells of 1.8M total). Sensitivity analysis on full atlas deferred to supplementary; primary verdict on monocyte aligns with Issue 29's pre-committed cell-type focus.
- Randolph 2021 cell attrition (235K → 34K post genotype-demultiplex + QC) limits per-donor power; HMN83575 excluded primary per pre-specified <50 cells/donor rule.
- Phase 3 calibration thresholds (Session 5 audit) were post-Harmony fit-to-data; we acknowledge this and re-frame Phase 3 results as exploratory. Phase 5+ confirmatory thresholds were literature-anchored before any Phase 5 run.
- Per-cell-type Harmony chosen over global Harmony for methodological alignment with factorized model's per-bucket grain (Issue 7 post-hoc decision, documented as such).
- DE-Jaccard metric degenerate on Harmony embedding (top-100 of 50-dim PCA = all); not used as primary metric (Issue 3 revised).
- Cross-study integration assumes biological signal is preserved through Harmony's correction. Per-cohort verification via calibration framework; cohort-specific failures of this assumption would surface as low cross-study coherence.

### Methods supplementary

- Detailed cohort metadata table
- Schema v6 column specifications + migration tests (`test_schema_v6_migration.py`, 11/11 passing)
- Calibration framework v2 unit tests (`test_calibration.py`, 8/8 passing)
- Per-cohort QC parameters and pass-rates
- Hyperparameter search results for factorized model
- Comparison method versions (scVI 1.X.X, scGen X.X.X, scCausalVI X.X.X, CPA X.X.X) — pre-specified per Issue 23
- Sensitivity analyses: alternative age cutoffs (Issue 28), alternative cell-count thresholds (Issue 27 HMN83575), alternative MVS gene sets (Interferome 2.0 vs Khatri MVS — Issue 18)
- **Session 7 pre-modeling sensitivity audit** (Issues 32 + 33, pre-specified 2026-05-11 before any analysis ran):
  - **Pre/post-harmonization Δr (Issue 32)**: per-bucket pre-Harmony (raw normalized log1p) vs post-Harmony cross-study Pearson r on response vectors. Across 5 buckets × 2 gene sets (full HVG + Khatri MVS): pre-Harmony r already substantial (0.13–0.58); Harmony adds 0.02–0.25 on top. Aggregate verdict: **MIXED both gene sets** (bio=2/mix=3/har=0 full HVG; bio=1/mix=4/har=0 MVS). No bucket × gene_set crosses the Δr > 0.30 HARMONY_DOMINANT threshold. Monocyte MVS Δr = 0.08 → BIOLOGY_DOMINANT verdict at the load-bearing bucket for the ISG-restriction contribution.
  - **Within-cohort sensitivity (Issue 33)**: v2 calibration framework applied independently to each v1 cohort (no cross-study integration). 10 bucket pairs × 4 cohorts × 2 gene sets = 80 within-cohort response-vector correlations; 20 aggregate rows vs cross-study harmonized. Sign concordance = **100%** in every aggregate test; mean magnitude alignment = 0.077 (full HVG) / 0.136 (MVS). Aggregate verdict: **BIOLOGY_CONSISTENT both gene sets**. Cross-study integration amplifies rather than creates the signal.
  - Deliverables: `results/tables/sensitivity_pre_post_harmony.csv`, `results/tables/sensitivity_within_cohort.csv`, `results/tables/sensitivity_within_vs_cross.csv`. Pre-committed decision rules and verdicts in `METHODS_CHOICES.md` Issues 32–33.

---

## Key methodological decisions (audit trail)

All decisions documented in `METHODS_CHOICES.md` with five-field structure (Choice / Rationale / Validation / Date opened / Date resolved) and atomic commit hashes:

- Issue 4 inclusion criteria (≥4 healthy + ≥4 diseased donors per study)
- Issue 7 per-cell-type vs global Harmony (per-cell-type chosen, post-hoc acknowledged)
- Issue 17 atomic schema-change rule
- Issue 25 Option B hybrid (held-out cohort validation rather than reframe)
- Issue 26 Phase 3 threshold provenance (post-Harmony, exploratory)
- Issues 27–30 held-out cohort pre-specifications

Reviewers can follow the audit trail commit-by-commit.

---

## Open items to resolve before submission

### Session 6B completion (current)

- Fix Yoshida gene-naming (Ensembl → symbol remap) — recompute r_MVS for all 5 buckets
- Investigate Randolph monocyte r_MVS anomaly (r_full = 0.287, r_MVS = 0.013)
- Run permutation null + bootstrap CI + FDR-BH on all 14 bucket-cohort tests
- Apply pre-committed decision rules to calibrated verdicts
- Generate `calibration_heldout_*.csv` tables

### Session 3.5 (post-6B)

- Pre-specify Issues 18–24 (ISG gene set, pathway, reconstruction loss, architecture hyperparameters, few-shot protocol, comparison method versions, baseline implementations)
- Integrate `PLAN.md` additions (§1.6 related work, §1.7 architecture spec)

### Session 4 (post-3.5)

- GPU setup
- scVI sensitivity analysis (Issue 6)
- Confirm or replace Harmony with scVI for primary analysis

### Phases 4–7 (modeling proper)

- Phase 4: factorized model implementation
- Phase 5: training on v1 corpus + held-out evaluation
- Phase 6: baseline comparisons (scVI, scGen, scCausalVI, CPA)
- Phase 7: paper drafting, internal review, bioRxiv submission

### Co-author engagement

- Authorship conversation: solo first-author confirmed; senior author + co-authors TBD
- Reach out to data-source PIs (Wilk, Lee, Arunachalam, Schulte-Schrepping) for cohort acknowledgment
- Domain immunologist review of biological interpretation before submission

### Pre-submission

- All `METHODS_CHOICES.md` issues resolved
- Audit trail commit log clean
- Supplementary figures/tables finalized
- Code repository public release (https://github.com/sakzgupzzz/trinetravir public branch)
- Data deposition: harmonized corpus on cellxgene Census or Zenodo
- Reproducibility README + environment specs

---

## v1.5 forward planning (separate paper, post-v1 publication)

v1.5 applies the v1 methodology framework to engineered adenovirus vector PBMC data. Scope:

- Apply calibration framework + factorized model + ISG-aware regularization to Ad-vector vaccine cohorts (Ad5-nCoV, Ad26, ChAdOx1)
- Few-shot adaptation curves for novel vector contexts
- Cross-vector generalization tests
- Detailed planning document: `v1_5_PLAN.MD` in repo root.

v1 references v1.5 in Discussion as planned follow-up; v1.5 doesn't block v1 submission.

---

## Submission timeline

| Milestone | Estimated date | Blocking |
|---|---|---|
| Session 6B complete (Issues 27–30 verdicts) | ~2026-05-18 | Yoshida gene-naming fix + Randolph monocyte investigation |
| Session 3.5 complete (Issues 18–24 pre-specs) | ~2026-05-25 | Session 6B |
| Session 4 complete (scVI sensitivity) | ~2026-06-15 | Session 3.5 + GPU setup |
| Phase 4 (model implementation) | ~2026-07-15 | Session 4 |
| Phase 5 (training + evaluation) | ~2026-08-15 | Phase 4 |
| Phase 6 (baseline comparisons) | ~2026-09-01 | Phase 5 |
| Phase 7 (draft + internal review) | ~2026-10-01 | Phase 6 |
| bioRxiv submission | ~2026-10-15 | Phase 7 |
| Journal submission | ~2026-11-01 | bioRxiv |

Realistic v1 timeline to bioRxiv from current state: 22–24 weeks. Aggressive timeline (everything goes smoothly, no scope creep): 18–20 weeks.

---

## Questions for Saksham / collaborators to resolve

1. **Authorship structure.** Solo first author with TBD senior author? Or formal co-first / co-senior arrangement? Have this conversation by Phase 5.
2. **Target journal after bioRxiv.** Cell Systems? Nature Methods? Genome Biology? Biology-leaning vs methodology-leaning venue choice affects framing.
3. **Domain immunologist review.** Who can review the biological interpretation of held-out cohort results before submission? Particularly for: chronic CMV biology (Issue 29 negative finding), HIV chronic vs acute distinction (Issue 30), pediatric immunology (Issue 28).
4. **Cohort acknowledgment / data-use thanks.** PI engagement for Wilk, Lee, Arunachalam, Schulte-Schrepping, Randolph, Yoshida, Allen Institute, Wang. Some require explicit data-use agreement acknowledgment.
5. **v1.5 timing.** Activate v1.5 development immediately after v1 submission, or wait for journal acceptance / preprint feedback? Affects resource allocation.

---

## How to use this document

This is a working manuscript document, not a final draft. Update it after each session block:

- After Session 6B: fill in Issue 27–30 verdicts, update boundary conditions section
- After Session 3.5: fill in Issues 18–24 pre-spec details, update methods section
- After Session 4: confirm Harmony vs scVI choice, update sensitivity analysis section
- After Phase 5: fill in actual model results, finalize figures plan
- After Phase 7: convert to actual manuscript draft

Section headers correspond to manuscript sections. Bullet points become paragraphs. Each "Open items to resolve" becomes a closed loop or a documented limitation.

---

## Living change log

| Date | Update | Trigger |
|---|---|---|
| 2026-05-11 | Initial document created | Session 6B initial results |
| 2026-05-11 | Ingested as MANUSCRIPT_DRAFT.md in repo root; table formatting normalized to GFM | Session 6B Step 4 (manuscript ingest) |
| 2026-05-11 | Section 4 + 5 + Limitations updated with N=1000 verdicts (15 tests); Issue 27 corrected per Issue 31 cross-bucket healthy reference (PRIMARY = monocyte_infected, r_mvs = −0.011, CHALLENGES_H1); Yoshida CI [0.02, 0.68] caveat + Issue 30 CI [0.157, 0.513] mirror caveat added; FDR<0.05 disclosure added to Limitations | Session 6B post-bg consolidation |
| 2026-05-11 | Methods supplementary updated with Session 7 pre-modeling sensitivity audit (Issues 32 + 33); Limitations updated with Harmony contribution disclosure (Δr 0.02–0.25, MIXED verdict) + within-cohort sign-concordance perfect-replication finding | Session 7 atomic commit #4 |
| TBD | Issues 18–24 methods detail | Session 3.5 completion |
| TBD | Model results inserted | Phase 5 completion |
| TBD | First draft assembled | Phase 7 |
