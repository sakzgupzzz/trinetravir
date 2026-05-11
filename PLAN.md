# Cross-Virus Generalization for Single-Cell Host Response Prediction

Project plan v1.3. Designed as the working specification for a solo computational research project producing a bioRxiv preprint plus open-source code release. Timeline target: 12–14 weeks of focused work. Compute strategy: CPU-first, escalate to GPU only when CPU becomes prohibitive.

**v1.3 changes (2026-05-11):**
- §1.1 Scientific motivation expanded from single paragraph to four subsections (§1.1.1 immunology, §1.1.2 methodological gap, §1.1.3 evaluation gap, §1.1.4 PBMC-only scope rationale). Cites Session 5+6B+7 empirical evidence anchoring the ISG-conservation finding (METHODS_CHOICES Issues 18-33) + Ahlmann-Eltze 2025 + PertEval-scFM as field-default-evaluation precedents motivating Issue 24 baseline expansion.
- §1.5 v1.5 forward planning pointer added (one-paragraph reference to `v1_5_PLAN.MD`).
- §1.6 Related work + competitive positioning added (5 buckets: existing perturbation prediction methods, foundation model baselines, cross-virus PBMC transfer learning prior work, single-cell calibration precedents, cross-cohort integration with external healthy controls).
- §1.7 Factorized model architecture (motivation-tier specification) added. References Issues 18-24 pre-specs. §4 Phase 8 retains implementation-tier pseudocode; §1.7 is motivation-tier. Frames model as secondary contribution per post-Session 7 manuscript framing.
- §1.8 Exploratory vs confirmatory evidence framing (originally added 2026-05-11 in Session 5 audit-response) retained as-is.
- Does NOT change v1 scope, deliverables, non-goals, or hypotheses. References Session 5+6B+7 empirical findings where applicable.

**v1.2 changes (2026-05-10):**
- Repo and Python package renamed `trinetravir` (was `cross-virus-scrna` / `crossvirus`).
- License: MIT.
- v1.1 scope locked to **PBMCs only** across all viruses. Multi-compartment (airway epithelium, organoids, intestinal) deferred to v2.
- §9 mini-gate executed and **PASSED**: Lee et al. 2020 PBMC (cellxgene `de2c780c`), Pearson r = 0.46 between SARS-CoV-2 and IAV response vectors, top-100 up-regulated gene Jaccard = 0.10. Signal present.
- Dataset acquisition (Phase 2) updated with cellxgene Census availability check: SARS-CoV-2 PBMC well-represented, IAV PBMC limited to Lee et al. (Census), RSV and additional IAV require GEO direct download.
- **Schema rename: `infection_status` → `donor_disease_status` (v1.1, PBMC-only).** Values changed from `infected`/`mock` → `diseased`/`healthy_control`. Reason: PBMCs from a virally-infected donor are mostly NOT directly infected at the cellular level — the signal is systemic cytokine / IFN response to viral disease, not cell-autonomous infection. The factorized cross-virus model works on this systemic signal; the framing of the eventual paper must be honest about this distinction. The PLAN §2 per-cell `infection_status` (infected/bystander/mock) semantics belong to v2 airway studies where viral reads can be assigned per cell. New obs column `label_source` (e.g. `disease_proxy`) records label origin so downstream code can branch on label semantics. `mock_control` is reserved as a future allowed value for in-vitro mock-infected studies (not used in v1). **See METHODS_CHOICES.md Issue 1 for the resolved schema decision and rationale.**

- **Methodological discipline (added 2026-05-10).** Every non-trivial methodological choice is logged in `METHODS_CHOICES.md` with required scientific rationale and validation strategy. Pre-specified rules now block scope creep at phase boundaries: study inclusion criteria (Issue 4), hyperparameter tuning policy (Issue 14), cross-virus evaluation protocol (Issue 15). Read `METHODS_CHOICES.md` at the start of every phase.
- **GEO time-budget addendum (Phase 2).** Census fetches take ~10 min per study because cellxgene already harmonized them. GEO direct downloads (RSV, second IAV, DNA virus control) need raw count matrix processing, manual metadata parsing, and ad-hoc QC because the studies were never standardized — budget **2-4 hours per GEO study**, not 10 minutes. Do not conflate the two in time estimates.

---

## 1. Context

### 1.1 Scientific motivation

#### 1.1.1 Immunology — shared antiviral programs across respiratory viruses

Type-I interferon-stimulated gene (ISG) induction is the canonical conserved component of the antiviral response across DNA and RNA viruses (Schoggins 2014 *Curr Opin Virol*; Schneider 2014 *Annu Rev Immunol*). In peripheral blood, monocytes dominate the early systemic IFN response; lymphoid compartments (B / NK / CD4T / CD8T) carry a downstream paracrine ISG signal whose magnitude and timing differ from monocytes. The biological premise of cross-virus generalization at the response-vector level is that this conserved component should transfer across viruses while virus-specific entry- and replication-associated signatures should not.

The v1 corpus empirically supports this premise. Per Session 5 audit + Session 6B held-out validation, ISG-restricted cross-study Pearson r exceeds full-HVG r by +0.06 to +0.23 across 4 of 5 v1 buckets (Khatri MVS subset; see `METHODS_CHOICES.md` Issue 18 + `references/khatri_mvs_gene_list.csv`). Pre-Harmony cross-study coherence is substantial (0.13–0.58, Session 7 Issue 32), and within-cohort effects are 100% sign-concordant with cross-study harmonized findings (Session 7 Issue 33). The ISG-conservation finding is biology with Harmony amplification, not an integration artifact.

#### 1.1.2 Methodological gap — single-virus evaluation as field default

Single-cell perturbation prediction methods (scGen, scCausalVI, CPA, CoupleVAE) and gene-expression foundation models (scGPT, Geneformer, scPRINT-2, STATE) are typically evaluated on within-virus or within-perturbation generalization: held-out cell types, held-out doses, held-out donors. Cross-virus transfer — training on virus A and predicting host response to virus B without any virus B training data — has not been benchmarked as a standardized task.

Recent work has documented that simple baselines often match or beat more complex methods on standard perturbation benchmarks when evaluation rigor is high. Ahlmann-Eltze et al. 2025 (*Nature Methods*) analyzed scPerturb-style benchmarks and showed that linear delta + per-cell-type-mean baselines are competitive with neural perturbation prediction methods. The PertEval-scFM benchmark (Tu et al. 2024) reported that gene-expression foundation models do not consistently outperform simpler baselines on cross-cell-type perturbation tasks. The implication for v1: the factorized model must clear a higher bar than "neural beats predict-mean" — it must beat simpler factorization baselines (`METHODS_CHOICES.md` Issue 24 Category 2: Sparse PCA, NMF, ISG-score regression) to justify the architectural complexity.

#### 1.1.3 Evaluation gap — calibration framework, sensitivity audit, pre-registered held-out validation

Three field gaps that v1 directly addresses, each demonstrated empirically across Sessions 5+6B+7:

- **Calibration framework.** Cross-study response-vector coherence reported in single-cell PBMC viral integration papers is typically a single Pearson r with no calibration against a permutation null, no bootstrap CI, and no FDR correction across multiple bucket-cohort tests. v1's `src/trinetravir/eval/calibration.py` (v2) implements donor-level permutation null, bootstrap CI on observed r, split-half ceiling, and FDR-BH correction; unit tests at `src/tests/test_calibration.py` (8/8 passing). All Session 5+6B+7 verdicts use this framework. Cross-study coherence claims become "passes FDR-corrected confirmatory threshold at observed effect size" rather than "r=0.X looks high".

- **Sensitivity audit.** Session 7 (Issues 32+33) quantified per-cell-type Harmony's contribution to cross-study coherence (Δr 0.02–0.25 on top of pre-Harmony baseline; aggregate MIXED both gene sets; monocyte MVS Δr=0.08 BIOLOGY_DOMINANT at the load-bearing grain) and verified within-cohort sign concordance with cross-study findings (100% across 20 aggregate bucket-pair × gene_set tests). The audit lives in the git audit trail; the manuscript discloses the Δr range honestly in Limitations.

- **Pre-registered held-out validation.** v1's cross-context, cross-age, chronic-CMV, and chronic-HIV decision rules (`METHODS_CHOICES.md` Issues 27-30) were committed to the audit trail before observing the held-out cohort data. Mechanical application of those rules produces verdicts (SUPPORTS_H1 / CHALLENGES_H1 / scope-limitation / BORDERLINE) that the manuscript reports as-is with no interpretation latitude. Issue 31 (cross-bucket healthy reference for cluster-defined subsets) was pre-specified before the corrected Randolph re-run.

These three practices are sparse in the field. v1 demonstrates them as a publishable methodology contribution; the factorized model is a secondary contribution (Phase 5/6 empirical question, Issue 24).

#### 1.1.4 PBMC-only scope rationale

v1.1 locks scope to PBMC compartments across all viruses for two reasons. First, PBMC IFN response is monocyte-dominated, fast (acute paracrine cascade within 24-72h), and clean — direct viral PAMP sensing is rare; the dominant signal is systemic IFN-α/β cascade. Airway epithelial response is slower, contaminated by bystander tissue-damage signals, and confounded by per-tissue ACE2 / TMPRSS2 expression heterogeneity for SARS-CoV-2. Second, available SARS-CoV-2 + IAV PBMC scRNA-seq cohorts pass Issue 4 inclusion criteria (≥4 healthy + ≥4 diseased per study); airway cohorts at v1 scale do not.

v1.5 extends scope to engineered adenovirus vector PBMC data (see §1.5 + `v1_5_PLAN.MD`). v2 extends to airway epithelium + organoid + intestinal compartments.

### 1.2 Project deliverables

1. A harmonized cross-virus single-cell benchmark dataset across SARS-CoV-2, influenza A, RSV, and one DNA virus control (HSV-1 or CMV), released on Zenodo with a Python loader package.
2. A benchmark evaluation of existing methods (simple baselines, scGen, scCausalVI, foundation-model-embedding heads) on within-virus vs cross-virus tasks, documenting the performance gap.
3. A factorized model architecture that explicitly decomposes predicted response into shared antiviral and virus-specific components with biological regularization, evaluated against the above baselines.
4. Few-shot adaptation experiments showing how N cells from a novel virus close the cross-virus gap.
5. A bioRxiv preprint and a public GitHub repository with reproducible code.

### 1.3 Hypotheses

- H1: Foundation model embeddings + linear heads fail to transfer across viruses zero-shot, despite within-virus performance.
- H2: Transfer is gene-specific in a biologically predictable way (ISGs transfer; entry-specific genes do not).
- H3: Transfer is asymmetric across virus pairs, reflecting biological similarity.
- H4: A factorized shared/virus-specific architecture outperforms monolithic models on cross-virus transfer.
- H5: Few-shot adaptation (N ≤ 1000 cells from target virus) closes most of the cross-virus gap.

### 1.4 Non-goals

- Wet lab validation. This is computational-only.
- Predicting clinical outcomes (severity, mortality). Out of scope for v1.
- Engineered viral vectors (AAV, oncolytic, gene therapy). Discussion section only.
- Superinfection / coinfection prediction. Different problem.
- **Multi-compartment scope (v1.1).** PBMCs only for v1.1. Airway epithelium, organoid, and intestinal datasets deferred to v2. Rationale: PBMC IFN response is monocyte-dominated, fast, clean; airway epithelial response is slower and contaminated by bystander tissue damage signals. Mixing compartments risks cross-compartment signal swamping cross-virus signal in v1 modeling.

### 1.5 v1.5 forward planning

v1.5 applies the v1 methodology framework (calibration framework + factorized model + ISG-aware regularization) to engineered adenovirus vector PBMC scRNA-seq cohorts (Ad5-nCoV, Ad26.COV2.S, ChAdOx1, heterologous prime-boost). Full specification at `v1_5_PLAN.MD` in the repo root. v1.5 is a follow-on bioRxiv preprint, not a v1 dependency; v1 references v1.5 in its Discussion as planned future work, and v1.5 does not block v1 submission.

### 1.6 Related work and competitive positioning

#### 1.6.1 Existing perturbation prediction methods

scVI (Lopez 2018 *Nat Methods*), scGen (Lotfollahi 2019 *Nat Methods*), CPA (Lotfollahi 2023 *Mol Syst Biol*), and scCausalVI (Wang 2024) are the established single-cell perturbation prediction frameworks. Each predicts per-cell response to a stimulation/perturbation, typically benchmarked on within-stimulation tasks (held-out cells from the same condition) or limited cross-stimulation tasks (e.g., transferring across drug doses). Cross-virus transfer is outside their default evaluation scope. v1 evaluates these as Issue 23 baselines against the factorized model under identical calibration framework + held-out validation.

#### 1.6.2 Foundation models (Geneformer, scGPT)

Geneformer (Theodoris 2023 *Nature*) and scGPT (Cui 2024 *Nat Methods*) are gene-expression foundation models trained on millions of single cells. v1 includes them as Issue 23 baselines per critique-document concern 2 (deep learning necessity). Foundation model fine-tuning + inference is GPU-dependent; the v1 evaluation of these baselines is gated on Session 4 GPU setup + compute envelope. PertEval-scFM (Tu 2024) found that foundation models do not consistently outperform simpler baselines on cross-cell-type perturbation tasks; v1 will report whether the cross-virus task admits a similar pattern. Foundation model checkpoints pinned by HuggingFace revision hash per Issue 23.

#### 1.6.3 Cross-virus PBMC transfer learning prior work

A direct prior literature on cross-virus PBMC scRNA-seq transfer learning is largely absent. Cross-condition perturbation prediction in PBMCs (e.g., LPS vs PolyI:C stimulation) has been explored at the bulk RNA-seq + monocyte level (Khatri MVS, Andres-Terre 2015 *Immunity*; Mostafavi 2016 *Cell*) but not at single-cell resolution with held-out viral context evaluation. This absence is the novelty hook for v1: a standardized cross-virus PBMC scRNA-seq benchmark with pre-registered held-out validation across four biological axes (cross-context IAV, cross-age SARS-CoV-2, chronic-latent CMV, chronic HIV).

#### 1.6.4 Single-cell calibration and pre-registered evaluation precedents

Pre-registered evaluation protocols and calibration frameworks (permutation null + bootstrap CI + FDR correction) are standard in GWAS and clinical-trial statistics but are sparse in single-cell methods papers. The scIB benchmark (Luecken 2022 *Nat Methods*) reports a panel of integration metrics across methods but does not pre-register thresholds for "successful integration"; the Open Problems in Single-Cell Analysis benchmark (Lance 2024) reports method rankings but does not apply formal multiple-testing correction across rank changes. v1's calibration framework v2 — donor-level permutation null + bootstrap CI + FDR-BH + literature-anchored thresholds set BEFORE Phase 5 — sits in this gap and is documented as a methodology contribution.

#### 1.6.5 Cross-cohort integration with external healthy controls

The GSE157829 chronic HIV held-out cohort (Issue 30 primary contrast) has n=1 within-cohort healthy donor, below the Issue 4 inclusion criterion. v1's design pairs GSE157829 HIV donors against the v1 corpus aggregated healthy donors as cross-cohort baseline. This design follows established field precedent for low-control chronic viral cohorts:

- *eBioMedicine* 2025 chronic-viral PBMC scRNA-seq integration paper used external SARS-CoV-2 corpus healthy aggregate as cross-cohort baseline.
- PMC10040851 (HIV-on-ART vs external healthy PBMC corpus single-cell comparison).
- PMC9434837 (inflammatory chronic disease cross-cohort GSE healthy aggregate baseline).

The cross-cohort design is documented in `METHODS_CHOICES.md` Issue 30 resolution with explicit citations to these precedents. v1's choice is methodologically defensible, not a compromise.

### 1.7 Factorized model architecture (motivation-tier specification)

The factorized model is v1's **secondary contribution** (Phase 5/6 empirical question; primary contribution is the calibration framework + methodology audit trail). Whether the model architecture adds value over simpler baselines (Issue 24 Category 2: Sparse PCA, NMF, ISG-score regression) is a question for the Phase 5+ empirical results, not a foregone conclusion of the design.

**Decomposition.** Predicted per-bucket response vector = `f_shared(baseline)` + `f_specific(baseline, virus_id)`. The shared component is a single neural encoder/decoder trained on all-virus training data; it captures the conserved antiviral component (canonical ISG response). The virus-specific component is a virus-conditioned encoder/decoder; it captures virus-specific signatures (entry-receptor effects, viral hijack pathways, type-of-virus-specific transcriptional programs).

**Virus embedding.** Each virus is represented by a learned embedding vector concatenated into the `f_specific` input. Embedding dimensionality per Issue 21 pre-specification: ∈ {8, 16, 32}. Few-shot adaptation (Issue 22) freezes `f_shared` + `f_specific` weights and trains only a new virus embedding for an unseen virus, isolating the embedding's role in cross-virus transfer.

**Regularization.**
- *ISG-aware regularization* (Issue 18): the `f_shared` component is regularized to align with the Khatri MVS gene set (`references/khatri_mvs_gene_list.csv`, 86 genes from Andres-Terre 2015 Table S2 high-confidence core subset); the `f_specific` component is penalized for nonzero output on MVS genes. The constraint encodes the empirical observation that canonical-ISG response is conserved (Session 7 Issue 32 monocyte MVS Δr=0.08 BIOLOGY_DOMINANT) and should reside in the shared component.
- *Pathway-aware regularization* (Issue 19): factor loadings on REACTOME R-HSA-913531 (interferon signaling) co-members are penalized for divergence from each other. Undirected adjacency only; no transitive expansion. Encodes pathway-level functional structure beyond identity-set membership. If pathway-aware weight tunes to ~0 at Phase 5 under held-out donor validation, the term is dropped (Issue 19 validation strategy).

**Reconstruction loss.** MSE on response vectors as primary (Issue 20; aligns the training unit with the calibration framework + cross-study coherence metric + held-out validation tests, all of which operate on per-study response vectors). NB-GLM on counts as Phase 5 sensitivity; switch headline if NB exceeds MSE by Δr ≥ 0.10 cross-study AND flips a held-out verdict (Issue 27 CHALLENGES→SUPPORTS or Issue 29 scope-limitation→appropriate-discrimination).

**Training and evaluation.** Per-bucket training per Issue 16 (lymphoid stratification): separate `f_shared` + `f_specific` per cell-type bucket (monocyte, B, NK, CD4T, CD8T). Hyperparameter search per Issue 21 (20-config budget per Issue 14; held-out donor 80/20 split). Cross-virus evaluation per Issue 15 (leave-one-virus-out, both directions reported). Comparison method versions pinned per Issue 23 at Phase 7 launch.

**Bar to beat.** Per Issue 24 Category 2: factorized model must exceed sparse PCA + NMF + ISG-score-regression by Δr ≥ 0.05 cross-study Pearson averaged across buckets. If it does not, the paper acknowledges that the model architecture does not add value beyond gene-set restriction; the methodology contribution carries the paper.

**Implementation specification** (pseudocode, hyperparameter tuning workflow, debugging notes, GPU compute envelope): see §4 Phase 8 in this document. §1.7 is motivation-tier; §4 Phase 8 is implementation-tier.

### 1.8 Exploratory vs confirmatory evidence (added v1.3, 2026-05-11)

The v1 paper distinguishes two categories of evidence in its methods + results sections. The distinction was added in Session 5 after a hostile-reviewer audit identified that Phase 3 gate thresholds were set post-Harmony (annotated as "above the pre-Harmony r"), which is fit-to-data, not pre-specification.

**Phases 1-3 produce exploratory/discovery evidence.** These phases identified which buckets have cross-study signal worth pursuing. The thresholds applied in Phase 3 (monocyte 0.60, B 0.40, NK 0.35, CD4T 0.30, CD8T 0.25) were chosen after observing Harmony output and are therefore not pre-specified. Calibrated Session 3 verdicts on the Phase 3 + Phase 3.5 data confirm which buckets have signal under permutation null + split-half ceiling criteria, but those verdicts are still computed against thresholds and on data that informed the gate design. They are exploratory and reported as such.

**Phases 4 onward produce confirmatory evidence at pre-registered thresholds.** Phase 5 thresholds will be set from external literature (Khatri MVS r≈0.45 monocyte module preservation; Pan et al. 2023 cross-virus monocyte r≈0.55-0.65) BEFORE running Phase 5. The Phase 5 pre-registration commits the v1 paper's primary claims to those thresholds in advance. Sessions 3.5 (pre-specifications for Phases 5/7/9) and Session 4 (scVI sensitivity) are scoped to make this distinction enforceable.

The methods section will state: "Phase 3 + Phase 3.5 calibrated verdicts are exploratory and were used to scope the v1 corpus. Phase 5 onward applies pre-registered thresholds set from external literature; Phase 5+ verdicts are confirmatory."

This distinction is the v1 paper's response to the post-hoc-threshold concern. See METHODS_CHOICES.md Issue 26 (Phase 3 threshold provenance) for the full audit-response acknowledgment.

---

## 2. Repository structure

```
trinetravir/
├── README.md
├── LICENSE                         # MIT
├── SETUP.md                        # Free-tier signup checklist (HUMAN TODO)
├── pyproject.toml                  # uv-managed, Python 3.11 pinned
├── .python-version
├── .gitignore
├── .pre-commit-config.yaml         # ruff, black, mypy, nbstripout
├── uv.lock
├── configs/
│   ├── datasets.yaml               # Dataset registry + per-study infection_status definitions
│   ├── models.yaml                 # Model hyperparameters
│   └── evaluation.yaml             # Eval metrics and splits
├── data/
│   ├── raw/                        # Downloaded scRNA-seq files (gitignored)
│   ├── processed/                  # Harmonized AnnData files (gitignored)
│   └── reference/                  # ISG lists, viral receptor gene sets (TRACKED)
├── src/
│   ├── trinetravir/
│   │   ├── __init__.py
│   │   ├── data/
│   │   │   ├── download.py         # Dataset download orchestration (cellxgene + GEO)
│   │   │   ├── harmonize.py        # Cross-study harmonization
│   │   │   ├── qc.py               # QC pipeline
│   │   │   └── splits.py           # Within/cross-virus split generation
│   │   ├── baselines/
│   │   │   ├── predict_mean.py
│   │   │   ├── linear_delta.py
│   │   │   └── knn.py
│   │   ├── methods/
│   │   │   ├── scgen_wrapper.py
│   │   │   ├── sccausalvi_wrapper.py
│   │   │   ├── foundation_embed.py # scGPT/Geneformer/scPRINT-2 wrappers
│   │   │   └── factorized.py       # Novel method
│   │   ├── eval/
│   │   │   ├── metrics.py          # Use Feb 2026 robust metrics
│   │   │   ├── gene_level.py       # ISG vs non-ISG decomposition
│   │   │   └── benchmark.py        # Main evaluation orchestrator
│   │   └── utils/
│   └── tests/
├── notebooks/
│   ├── 00_env_check.ipynb            # Phase 1 deliverable
│   ├── 01_data_exploration.ipynb     # §9 mini-gate (DONE — passed 2026-05-10)
│   ├── 02_sanity_check_signal.ipynb  # GATE 1 notebook (Phase 4)
│   ├── 03_baseline_results.ipynb
│   ├── 04_existing_methods.ipynb
│   ├── 05_factorized_method.ipynb
│   ├── 06_fewshot_analysis.ipynb
│   └── 07_paper_figures.ipynb
├── scripts/
│   ├── run_benchmark.py
│   └── reproduce_paper.sh
└── results/
    ├── tables/
    └── figures/
```

---

## 3. Compute infrastructure strategy

### 3.1 Principle: CPU-first, GPU last

The project is structured so that everything that can run on a personal laptop does run on a personal laptop. GPU is escalated to only when a step is either intractable on CPU (foundation model fine-tuning, dense factorized model training with many hyperparameter sweeps) or so slow on CPU that wall-clock dominates the iteration loop.

This is a deliberate ordering choice. The two most informative gates in the plan (Phase 4 signal check, Phase 7 benchmark v1) both come before any GPU work is required. If the project fails at either gate, no compute money has been spent. If it passes, you have a defensible artifact and clear justification for spinning up paid GPU.

### 3.2 Compute tier definitions

- **Tier 0 — local laptop**: 16GB+ RAM machine running Linux or macOS. Sufficient for everything through Phase 5 and parts of Phase 6 if you have patience for overnight runs.
- **Tier 1 — free CPU cloud**: for jobs you don't want to tie up your laptop with. Oracle Always Free (4 OCPU/24GB ARM Ampere) is the workhorse here, genuinely free indefinitely. AWS t3.large free tier (12 months, 2 vCPU/8GB) and Hetzner cheap CPU boxes (~€5-26/month for 4 vCPU and 8-32GB RAM) are backups.
- **Tier 2 — free GPU via providers' starter credits**: Google Cloud $300, Microsoft Azure $200, Oracle always-free tier, RunPod starter credits ($5-10). Stack these for ~$500-600 of GPU runway before any spend.
- **Tier 3 — free GPU via weekly resets**: Kaggle (30 hours/week on P100 16GB) and Google Colab free tier (T4, intermittent). Useful for ongoing small jobs throughout the project.
- **Tier 4 — paid spot/community GPU**: Vast.ai marketplace or RunPod community cloud. For the bulk of paid training time. ~$0.30-0.60/hr for RTX 4090 or A100 40GB.
- **Tier 5 — paid on-demand or research grants**: Lambda Labs research grants (apply once you have a preprint), or managed GPU clouds. Only if all of the above are exhausted.

### 3.3 What runs where

| Task | Minimum tier | Notes |
|------|------|-------|
| Data download, AnnData I/O | Tier 0 | Network-bound, not compute-bound |
| QC pipeline, scanpy preprocessing | Tier 0 | Use sparse matrix paths |
| Batch correction (Harmony) | Tier 0 | Use Harmony, not scVI, on CPU |
| Batch correction (scVI) | Tier 2/4 | scvi-tools is GPU-accelerated; CPU possible but very slow |
| Phase 4 sanity check (means, correlations) | Tier 0 | Pure numpy/scipy |
| Predict-mean, linear-delta, KNN baselines | Tier 0 | All CPU operations |
| scGen training on subset (<50k cells) | Tier 0 | Slow but viable; CPU run overnight |
| scGen / scCausalVI on full data | Tier 1 or Tier 2/3 | CPU possible; GPU 10-50x faster |
| Foundation model embedding extraction | Tier 1 or Tier 2/3 | One-time inference pass; CPU works on overnight runs for <500k cells |
| Foundation model fine-tuning (STATE, Stack) | Tier 2/4 | Not CPU-tractable |
| Factorized model training (iteration heavy) | Tier 2/4 | GPU required; many sweeps |
| Few-shot adaptation experiments | Tier 2/4 | GPU required |
| Manuscript writing, figures | Tier 0 | CPU only |

### 3.4 The honest CPU limits

Where CPU stops being practical:

- *Foundation model fine-tuning*. Extracting embeddings on CPU is fine because it's a one-time inference pass. Updating those models with gradient steps over a perturbation prediction head is not — STATE and Stack have hundreds of millions of parameters, and even tiny learning rates require many epochs.
- *Hyperparameter sweeps over the factorized model*. A single training run might fit on CPU. Twenty runs with different regularization weights, latent dimensions, and ISG mask choices is where you stop being able to iterate in human time. This is the first place you should pay for GPU.
- *Distribution-matching methods at scale*. Diffusion-based heads and Schrödinger bridge methods are GPU-native. If you decide these belong in the benchmark, plan GPU time.

### 3.5 Escalation rule

Default to laptop. Only escalate when a specific job is blocking iteration: "this needs to run overnight on CPU and I want to iterate twice today" is the trigger to move that job up a tier, not the trigger to migrate the whole project. Most of Phase 1 through Phase 6 should stay on the laptop with selected jobs offloaded to free tiers when they get long.

---

## 4. Phase-by-phase plan

Each phase ends with a defined deliverable. Phases marked GATE require a go/no-go decision before proceeding. Compute tier is annotated on each phase.

### Phase 1: Environment setup (Week 1) [Tier 0]

Tasks:
- Initialize repo with uv-managed pyproject.toml. Pin Python 3.11.
- Set up base dependencies: scanpy, anndata, scvi-tools (CPU build initially), torch (CPU-only wheel for now), pytorch-lightning, scikit-learn, numpy, pandas, pyyaml, hydra-core, wandb.
- Configure pre-commit hooks: ruff, black, mypy.
- Set up wandb project for experiment tracking.
- Sign up for the free credit accounts even if not using them yet: Google Cloud ($300, 90 days from activation), Azure ($200), Oracle Always Free, RunPod ($5-10 starter), Kaggle. Do this now so the clock isn't ticking when you actually need them.

Deliverable: working uv environment on laptop, registered accounts on relevant free tiers, "hello world" notebook importing scanpy and logging to wandb.

### Phase 2: Data acquisition (Weeks 2–3) [Tier 0]

**Scope (v1.1): PBMCs only.** Non-PBMC studies from prior plan (Chua nasopharyngeal, Ziegler lung, Triana intestinal organoids, Cao airway organoids, Steuerman mouse in vivo) are deferred to v2.

Target PBMC datasets — cellxgene Census availability verified 2026-05-10 (Census `2025-11-08` stable):

**SARS-CoV-2 PBMC (Census, plenty of options):**
| Study | Census dataset_id | n_cells | Notes |
|---|---|---|---|
| Lee et al. 2020, *Sci Immunology* | `de2c780c-1747-40bd-9ccf-9588ec186cee` | 59,572 | Already loaded; contains COVID + IAV + healthy in one design — anchor study |
| Wilk et al. 2020, *Nat Med* | `456e8b9b-f872-488b-871d-94534090a865` | 44,721 | Original publication-level dataset |
| Arunachalam et al. (Pulendran lab) | `59b69042-47c2-47fd-ad03-d21beb99818f` | 49,139 | Severity stratification |
| Schulte-Schrepping et al. 2020, *Cell* | `5e717147-0f75-4de1-8bd2-6fda01b8d75f` | 90,957 | Largest single PBMC cohort |
| Guo et al. (subset of Atlas) | `ae5341b8-60fb-4fac-86db-86e49ee66287` | 14,783 | Smaller, useful for fast iteration |
| MGH acute COVID cohort | `fa8605cf-f27e-44af-ac2a-476bee4410d3` | 59,506 | WHO severity scale, useful for stratification |

**Influenza A PBMC:**
- Census: only Lee et al. (already covered above). No standalone IAV PBMC datasets in Census.
- GEO TODO: identify additional human IAV PBMC scRNA-seq studies (e.g., Lefebvre et al. or similar live-virus in vitro PBMC stimulation cohorts). Required for cross-virus learning to be more than 2-virus.

**RSV PBMC:**
- Census: NONE (search returned empty for "RSV"|"respiratory syncytial"|"bronchiolitis").
- GEO REQUIRED. Candidate searches: pediatric bronchiolitis cohorts, Liao et al., infant PBMC studies. This is the riskiest data acquisition item; if no clean RSV PBMC dataset exists, fall back to a 3-virus benchmark (SARS-CoV-2 + IAV + DNA control) and document in limitations.

**DNA virus control (PBMC, pick one):**
- HSV-1 PBMC: Drayman et al. or follow-ups (likely GEO).
- CMV PBMC: latency/reactivation single-cell studies (GEO).
- Defer choice until SARS + IAV + RSV pipeline works.

Tasks:
- Write `src/trinetravir/data/download.py` orchestrating cellxgene Census + GEO downloads. **Central fix:** rebind `adata.var_names = adata.var["feature_name"].astype(str)` immediately after Census load (Census defaults to soma_joinid integers); assert no NaN feature names slip through.
- For each dataset, capture metadata in `configs/datasets.yaml`: virus, cell type, time post-infection, MOI or viral load proxy, donor, study, and **explicit label rule for that study** (`disease_proxy` for PBMC v1.1; `viral_read_threshold` reserved for v2 airway studies).
- Catalog viral RNA detection method per study (reads-to-viral-genome vs cell sorting vs bystander markers vs disease-status proxy).
- Use cellxgene's Python API where possible — preprocessed AnnData files save you from raw FASTQ processing.

**Output schema (v1.1 PBMC):** every persisted h5ad has obs columns `virus`, `donor_disease_status` (values: `diseased`/`healthy`), `label_source` (e.g. `disease_proxy`), `cell_type`, `donor_id`, `study_id`, plus all native cellxgene Census columns. The `donor_disease_status` label is *donor-level*, not cell-level — see v1.2 changelog at top of doc for the framing implications. Per-cell `infection_status` joins the schema in v2 when airway-epithelium datasets are added.

Deliverable: 6–10 raw PBMC datasets downloaded to `data/raw/`, with `configs/datasets.yaml` populated with full study-level provenance and per-study label rules. Storage: typically 5-30 GB total, easily fits on a laptop.

### Phase 3: Data harmonization (Week 3–4) [Tier 0, optionally Tier 1 if laptop RAM is tight]

This is the load-bearing engineering phase. Quality of harmonization determines whether the project works.

Tasks:
- Convert all datasets to AnnData with a shared schema: `obs` columns must include `virus`, `infection_status` (mock/bystander/infected), `cell_type`, `time_post_infection`, `donor_id`, `study_id`, `moi_or_viral_load`.
- Harmonize cell-type labels using a single ontology (Cell Ontology terms via OnClass or scAR). Anything left ambiguous gets dropped or relabeled.
- Define `infection_status` consistently. Use a uniform rule: a cell is "infected" if it has ≥ N viral reads (N study-specific, document choice). "Bystander" cells are uninfected cells from infected donors/wells. "Mock" cells are from mock-infected controls.
- Run a uniform QC pipeline: minimum genes/cell, mitochondrial percentage threshold, doublet detection (Scrublet).
- Use Harmony for batch correction at the integration step. Harmony runs cleanly on CPU and is faster than scVI on a single-machine setup. Use study and donor as batch keys.

Tier-1 use case: if your laptop has less than 32GB RAM and the combined dataset exceeds it, offload the harmonization to an Oracle Always Free ARM instance (4 OCPU/24GB) or a Hetzner CCX23 box (~€26/month for 4 vCPU and 16GB). Persist the harmonized AnnData files back to your laptop or to free object storage.

**Cross-study batch-correction watchpoint (added v1.2):** After Harmony correction, compute response-vector Pearson correlation *within the same virus* across studies. This must stay below ~0.7 to confirm Harmony has not over-corrected the biology. Reference floor: Lee et al. 2020 within-study, between-virus r = 0.46 (§9 mini-gate result). If post-Harmony cross-study within-virus correlations rise far above 0.7, back off Harmony `theta` (lambda parameter) or remove donor as a batch key. If they collapse below 0.46, Harmony is erasing real disease signal — also a problem.

**Phase 3 entry-gate (v1.2, settled 2026-05-10).** Pre-Harmony bulk cross-study r on 5 SARS-CoV-2 PBMC studies = mean off-diag 0.054, range [-0.58, +0.59]. mgh_acute_covid + guo_2020 excluded (1 healthy donor / 0 healthy donors; see configs/datasets.yaml). 4 clean studies remain.

Stratified diagnostic (`notebooks/03_celltype_stratified_consistency.ipynb`) confirms cell-type-composition drift is the dominant signal. Per-bucket mean off-diag r:
  - monocyte: 0.434  (8x lift from bulk)
  - B:        0.252
  - NK:       0.210
  - CD4T:     0.196
  - CD8T:     0.156

Per-study mean per-bucket r vs others: arunachalam 0.358, wilk 0.291, schulte_schrepping 0.198 (salvaged from ~0 bulk — keep), lee 0.150 (bulk inflated by lee's heavy monocyte fraction).

Severity scan: cellxgene Census strips all severity-like obs columns from these 4 studies. development_stage is age only. Severity stratification is v2 work (pull from publication metadata).

**Phase 3 harmonization strategy (v1.2, settled 2026-05-10).**
- **Per-cell-type Harmony, NOT joint Harmony with cell_type as a batch key.** Harmony's `vars_use` parameter specifies what to mix across, not what to preserve — passing `cell_type` would erase cell-type structure, the opposite of intent. Run Harmony separately per coarse bucket (monocyte / CD4T / CD8T / B / NK), each run filtered to one cell type and using `study_id` as the sole batch key. Response vectors are already a per-cell-type quantity, so a joint cross-cell-type embedding is not required for the gate.
- **donor_id default OFF.** Donors are fully nested within studies; adding donor_id without strong evidence over-corrects and can entangle with disease status when donors are unique to a condition. Add only if specific donors show as outliers in post-Harmony QC.
- For visualization (UMAP coloured by all cell types together), a separate global Harmony pass with `study_id` only is acceptable — but it is not the load-bearing harmonization for the gate.
- Implementation in `src/trinetravir/data/harmonize.py`. Pipeline per bucket: HVG selection (n_top_genes=4000, batch_key=study_id, flavor=seurat_v3) → scale → PCA (n_comps=50) → harmonypy on PCA embedding → project corrected PCA back to scaled-HVG gene space via PCA loadings → recompute mean(diseased) − mean(healthy) per study in this projected space → pairwise Pearson r across studies.

**Phase 3 exit-gate (v1.2, refined 2026-05-10).** Per-cell-type post-Harmony cross-study Pearson r thresholds — uniform r ≥ 0.5 would be a false-negative gate because T-cell repertoires are inherently more donor-specific:
  - monocyte: r ≥ 0.60   (pre 0.434; IFN-dominated, should clean up well)
  - B:        r ≥ 0.40   (pre 0.252)
  - NK:       r ≥ 0.35   (pre 0.210)
  - CD4T:     r ≥ 0.30   (pre 0.196)
  - CD8T:     r ≥ 0.25   (pre 0.156)

Plus relative-ordering check: monocytes highest, T cells lowest. If T cells suddenly correlate better than monocytes after harmonization, the correction has gone wrong.

**Phase 3 gate result (v1.2, executed 2026-05-10).** Notebook `04_phase3_harmonization.ipynb` ran per-bucket Harmony on the 4 clean studies (244,389 cells). Post-Harmony per-bucket cross-study mean off-diag r:

| bucket   | post-Harmony r | threshold | result | lift from pre-Harmony |
|----------|----------------|-----------|--------|------------------------|
| monocyte | **0.701**      | 0.60      | PASS   | +0.267                 |
| NK       | **0.384**      | 0.35      | PASS   | +0.174                 |
| CD4T     | **0.321**      | 0.30      | PASS   | +0.125                 |
| B        | 0.297          | 0.40      | FAIL   | +0.045                 |
| CD8T     | 0.169          | 0.25      | FAIL   | +0.013                 |

**3 of 5 buckets pass.** Relative ordering (monocyte > NK > CD4T > B > CD8T) preserved as predicted by the per-bucket-threshold rationale.

**Diagnostic on the two failing buckets (2026-05-10):**
1. `scripts/phase3_donor_id_retry.py` — added `donor_id` as second Harmony batch key on B + CD8T only. **Discarded.** B 0.297 → 0.135 (−0.162; schulte ↔ wilk crashed to −0.417). CD8T 0.169 → 0.110 (−0.059). Confirms PLAN rule that donor_id over-corrects when donors are fully nested in studies; entangles with disease status.
2. `scripts/phase3_lee_diagnostic.py` — Lee-specific 3-question probe. Findings:
   - Healthy donor counts per bucket: Lee 4 / arunachalam 5 / wilk 6 / schulte 21. Lee similar to arunachalam, not the bottleneck.
   - **Cell-type annotation divergence is the smoking gun.** Lee's lymphoid `cell_type` labels lack memory/naive subdivisions present in the other 3 studies. Lee's "B" bucket = exclusively memory B (IgG-neg class switched memory + IgG memory); arunachalam/schulte's "B" bucket is dominated by naive B. Lee's CD8T = `CD8-pos alpha-beta T cell` + `effector CD8`, no memory subtypes. Same coarse bucket name, different mixtures of cell populations.
   - Lee IAV-vs-SARS coarse composition Pearson r = 0.986 (ordering preserved), but absolute fractions differ substantially: IAV is 58% monocyte vs SARS 38% monocyte; B 4.7% IAV vs 15.4% SARS. **Phase 4 cross-virus eval must be cell-type stratified, not bulk** — Gate-1 r=0.46 bulk result is partly composition-driven.
3. `scripts/phase3_lee_out_retry.py` — Lee-out gate retry. CD8T 0.169 → **0.260** (PASSES 0.25 threshold). B 0.297 → 0.352 (still fails 0.40; schulte ↔ wilk = 0.021 — wilk's generic `B cell` label fails to align with schulte's fine-grained naive/memory/transitional/plasmablast). **Annotation divergence is not Lee-localized — wilk B annotation is also coarse-and-different.**

**Phase 3 commit shape (v1.2, settled 2026-05-10): monocyte-primary, lymphoid-secondary.**

- **Primary Phase 4 cross-virus benchmark = monocyte response transfer (SARS-CoV-2 ↔ IAV).** Lee is the within-study cross-virus anchor (8 SARS donors + 5 IAV donors + 4 mock donors); arunachalam/wilk/schulte provide within-virus cross-study harmonization validation. Lee's monocyte pairwise r with other studies is 0.58–0.87 — annotation divergence does not affect monocyte (label vocabulary is consistent: classical / non-classical / CD14+ / macrophage).
- **Secondary cross-virus benchmarks = CD4T and NK** (both gate-passing). Explicit caveat on Lee annotation. Sensitivity analyses, not headline.
- **Excluded from cross-virus claims = B and CD8T.** Retained in the within-virus cross-study analysis. Annotation-divergence cause documented; not a Harmony failure, not a Phase 3 strategy failure.
- **Phase 4 evaluation protocol: stratify by cell type, never bulk.** Lee's IAV-vs-SARS composition difference (58% vs 38% monocyte) means bulk response-vector correlations confound transcription with composition. Per-cell-type evaluation is mandatory. This supersedes the earlier "Cross-virus eval TODO" framing — it is a hard constraint, not a check.

**GEO acquisition still DEFERRED.** v1.5 — see `BACKLOG.md`. v1.5 must precede GEO ingest because adding RSV / 2nd IAV / DNA-control studies without unified lymphoid annotation just expands the annotation-divergence problem.

If Phase 3 gate had failed on monocyte specifically, this would be Phase 3.5 / blocker. It did not. Monocyte gate r = 0.701 with strong relative-ordering signal. Phase 4 monocyte-primary is well-founded.

Deliverable: a single harmonized AnnData per virus, stored in `data/processed/{virus}.h5ad`, plus a combined `data/processed/all_viruses.h5ad`.

### Phase 4: Sanity check (End of Week 4) — GATE 1 [Tier 0]

> **Status (2026-05-10):** §9 mini-gate executed on Lee et al. 2020 PBMC (cellxgene `de2c780c`). Pearson r = 0.46 (SARS-CoV-2 vs IAV response vectors), Spearman r = 0.31, top-100 up-regulated gene Jaccard = 0.10. **PASSED** with margin (well below 0.7 threshold). The full Phase-4 analysis below still needs to run on the harmonized multi-study dataset; the mini-gate only confirms within-study signal in one anchor cohort.

**Cross-virus evaluation protocol (METHODS_CHOICES Issue 15, pre-specified 2026-05-10).** Phase 4 uses leave-one-virus-out cross-validation per `configs/evaluation.yaml > cross_virus_protocol`. For v1 (SARS-CoV-2 + IAV), this is train-SARS / test-IAV and train-IAV / test-SARS; both directions plus the mean are reported in headline figures.

**Cell-type stratification is mandatory in Phase 4 (Phase 3 finding, 2026-05-10).** Lee's IAV-vs-SARS coarse-cell-type composition differs substantially (58% monocyte for IAV vs 38% for SARS). Bulk cross-virus response correlation conflates composition with transcription. All Phase 4 cross-virus correlations are computed per cell-type bucket (monocyte primary; CD4T / NK secondary; B / CD8T excluded from cross-virus claims per Phase 3.5 annotation-divergence findings).

Critical decision point. Before any modeling, validate that the data has the signal you assumed.

Tasks (one notebook, `02_sanity_check_signal.ipynb`):
- For each virus, compute the mean expression profile of infected cells minus mean of mock cells (the "response vector").
- Compute pairwise Pearson and Spearman correlation between response vectors across viruses.
- Identify the top differentially expressed genes per virus. Compute overlap (Jaccard) with curated ISG list (Interferome DB or Mostafavi lab ISG curation).
- Project all response vectors into 2D (UMAP or PCA) to visualize whether viruses cluster apart in response space.

Gate criteria:
- If pairwise response-vector correlations are > 0.9 across viruses, the project is in trouble — there is not enough virus-specific signal to learn, and any transfer learning is trivial. Stop and reframe.
- If correlations are < 0.5, the viruses are genuinely distinct and the project has clear signal.
- ISG overlap should be high (> 60%) across all viruses; this confirms the shared-program hypothesis.

If the gate passes, proceed. If not, the failure itself is publishable as a brief commentary preprint and saves you months.

This entire gate runs on the laptop in minutes once data is harmonized. No GPU spend before this point.

### Phase 5: Simple baselines (Week 5) [Tier 0]

Implement and evaluate three baselines on the within-virus and cross-virus tasks.

Tasks:
- `predict_mean.py`: output the mean post-infection expression vector from training, ignoring input cell.
- `linear_delta.py`: compute (post − pre) shift vector in PCA space from training; add to test cell baseline. PCA on combined mock + infected pre-projection.
- `knn.py`: for each test baseline cell, find K nearest training baseline cells (K = 25, 50, 100), return their post-infection mean.

All three baselines run on CPU in seconds-to-minutes. No GPU needed.

Deliverable: a results table showing within-virus vs cross-virus performance for each baseline on each evaluation metric.

### Phase 6: Existing methods comparison (Weeks 6–7) [Mixed: Tier 0/1 for some methods, Tier 2/3 for foundation models]

Wrap and evaluate published methods on the same task setup.

**Hyperparameter policy (METHODS_CHOICES Issue 14, pre-specified 2026-05-10).** Every benchmark method tunes via held-out donor-level validation split per `configs/evaluation.yaml > hyperparameter_policy`. Compute budget: max 20 hyperparameter configurations per method on the within-virus validation set; the configuration that maximizes the primary cross-study coherence metric is evaluated on the cross-virus test split. No published-defaults policy — that would bias against methods whose original-paper data differs most from PBMC cross-virus.

Tasks (CPU-tractable, do these first):
- `scgen_wrapper.py`: use theislab/scgen. Train on combined virus data, evaluate within-virus and cross-virus. scvi-tools runs on CPU; expect 30-90 minutes per training run on a laptop for ~100k cells. Multiple runs across virus combinations can stack overnight.
- `sccausalvi_wrapper.py`: similar story. CPU-trainable on subset; full data may want offloading to a Tier-1 cloud CPU box if laptop is also needed for other work.

Tasks (GPU-helpful but CPU-possible):
- `foundation_embed.py`: pull scGPT, Geneformer, scPRINT-2 pretrained weights. Compute zero-shot embeddings of baseline cells. This is forward-pass-only inference; CPU works but slowly. For ~100k cells through scGPT, expect 4-12 hours on CPU vs ~30 minutes on a T4 GPU. Two options:
  - *CPU path*: run overnight jobs, cache embeddings to disk, never recompute.
  - *Free-GPU path*: use Kaggle (P100 16GB, 30 hr/week free) for the embedding extraction. One Kaggle session is more than enough. This is probably the first place you actually want a GPU.
- Train a small MLP head on (embedding, virus_id) → post-infection profile. MLP head is small; trains on CPU in minutes once embeddings are cached.

Tasks (GPU-required, defer if possible):
- Fine-tuning STATE, Stack, or scPRINT-2 with a perturbation prediction head. These models are large enough that gradient updates on CPU are not practical. Use GCP credits for an A100 spot instance, or a single A100 evening on Vast.ai (~$5 of credit). Defer until after the zero-shot embedding baselines are done; the fine-tuning may not be needed if zero-shot results are already informative.

Deliverable: expanded results table including existing methods. By end of this phase, you have the headline finding of the benchmark paper.

### Phase 7: Cross-virus benchmark v1 (Week 8) — GATE 2 [Tier 0]

At this point you have enough for a publishable benchmark paper without the novel method. Decide:
- If within-virus vs cross-virus performance gaps are dramatic and consistent across methods, the benchmark paper is strong on its own. Consider submitting a short preprint now and continuing to the method as a follow-up.
- If gaps are modest or the picture is messy, the benchmark needs the novel method to make a clear contribution.

Evaluation metrics (per the Feb 2026 metrics-failure literature, avoid naive Wasserstein or Energy distance alone):
- Per-gene Pearson r between predicted and true post-infection mean expression (R^2 across all genes; separate R^2 for ISG subset vs non-ISG subset)
- Differential expression overlap: top 100 predicted DE genes vs top 100 true DE genes, Jaccard score
- Direction-of-change accuracy: for each gene, predicted vs true sign of (post − pre) shift, accuracy
- Distribution-aware metric: MMD with multiple kernels, validated against simulations
- Held-out cell type within target virus (combined cross-virus + cross-cell-type test)

This phase is all analysis and writeup. CPU only.

Deliverable: a final benchmark table that is the headline figure of the eventual paper.

### Phase 8: Factorized method (Weeks 9–10) [Tier 2/4 — GPU recommended]

Implement the novel architecture.

Design:

```
# Pseudocode for the factorized model
# response = shared_component(cell_baseline) + virus_specific_component(cell_baseline, virus_id)

class FactorizedResponseModel(nn.Module):
    """
    Decomposes virus-induced expression change into shared antiviral program
    and virus-specific program. Regularization encourages the shared component
    to dominate on a curated ISG gene set.
    """
    def __init__(self, n_genes, n_viruses, hidden_dim, isg_mask):
        # shared_encoder: cell_baseline -> shared_response_logits
        # virus_specific_encoder: (cell_baseline, virus_embedding) -> specific_response_logits
        # isg_mask: boolean mask over genes indicating ISG membership
        pass

    def forward(self, baseline, virus_id):
        shared = self.shared_encoder(baseline)
        specific = self.virus_specific_encoder(baseline, self.virus_embed(virus_id))
        return shared + specific

    def loss(self, pred, target, baseline, virus_id):
        # 1. Reconstruction: MSE between pred and target
        # 2. ISG regularization: penalize virus_specific magnitude on ISG genes
        # 3. Sparsity on virus_specific outside ISGs (encourage shared to explain
        #    as much as possible)
        # 4. Optional adversarial loss: virus_id classifier on shared output
        #    should perform poorly (shared is virus-agnostic)
        pass
```

Tasks:
- Implement in PyTorch Lightning with hydra configs.
- Initial debugging and overfitting tests: small subset on CPU is fine.
- Tune: hidden dim, regularization weights, learning rate. Use wandb sweeps. This is where GPU becomes mandatory — you'll want 10-30 sweep runs and CPU iteration is too slow.
- Recommended setup: spin up an RTX 4090 or A100 spot instance on Vast.ai (~$0.34-0.52/hr), run the sweep over a weekend, shut it down. Total cost: $10-25.
- Ablate each regularization component to show contribution.
- Compare against monolithic baseline (same architecture, no shared/specific split).

Deliverable: factorized model results in the same benchmark table.

### Phase 9: Few-shot adaptation (Week 11) [Tier 2/4 — GPU recommended]

Tasks:
- For each target virus, sample N ∈ {0, 10, 100, 1000} cells from the target virus's training data.
- Fine-tune only the virus_specific head (freeze shared backbone) on these N cells.
- Plot performance vs N. Compare against:
  - Fine-tuning a monolithic model on N cells
  - Training from scratch on N cells

Many small fine-tuning runs benefit from GPU. Same Vast.ai instance from Phase 8 can be reused.

Deliverable: few-shot curve figures for the paper.

### Phase 10: Manuscript and code release (Weeks 12–13) [Tier 0]

Tasks:
- Write the preprint. Suggested structure: intro motivating cross-virus transfer; related work positioning; methods covering benchmark construction and factorized model; results section with within-vs-cross gap, gene-level decomposition, factorized improvement, few-shot adaptation; discussion of limitations and engineered-virus extension.
- Polish the GitHub repo. README with one-command reproduction. Pin dependencies.
- Upload harmonized data to Zenodo with DOI.
- Submit to bioRxiv. Tag with relevant subject categories.

CPU only.

Deliverable: bioRxiv preprint + DOI + public repo.

---

## 5. Free credits and paid GPU options

### 5.1 Free credit stack (claim all in Week 1)

- **Google Cloud**: $300 in credits for new users, valid 90 days. Most flexible. Use this for the GPU-heavy Phase 8/9 work. Activates a 90-day countdown when you sign up, so plan around when you'll need it.
- **Microsoft Azure**: $200 for new accounts. Smaller pool but real.
- **Oracle Cloud**: Always-Free tier including 4 OCPU / 24GB ARM Ampere CPU instances. Truly free indefinitely; useful as your persistent secondary CPU box for long-running harmonization or embedding extraction jobs.
- **RunPod**: starter credits of $5-10 plus a referral bonus between $5 and $500 on first $10 spend. With community spot starting at $0.20/hr, the starter credit alone is 25-50 hours of GPU.
- **Kaggle**: 30 hours/week of P100 16GB free, indefinitely. Resets weekly. Run embedding extractions and small foundation model inference here.
- **Google Colab**: free T4 access with usage limits. Useful for prototyping. Colab Pro is $10/month if you want better availability and longer sessions.
- **Lambda Labs Research Grants**: free GPU time for published academic research. Apply once you have a draft or preprint; approval takes weeks.

Total runway from free tiers before any paid work: roughly $500-600 of effective GPU spend, plus persistent Oracle CPU compute, plus Kaggle's weekly resets.

### 5.2 Cheap paid GPU options (when free runs out)

For the bulk of training time after free credits, use spot/community pricing:

- **Vast.ai marketplace** (cheapest, less stable): A100 PCIe 40GB around $0.52/hr, A100 SXM 80GB around $0.67/hr, RTX 4090 at $0.31-0.34/hr, H100 80GB around $1.55/hr. Marketplace model — individual providers, variable reliability. Best when your training loop has good checkpoint/resume logic.
- **RunPod Community Cloud** (moderate price, better reliability): RTX 4090 from $0.34/hr, A100 80GB roughly $0.79/hr, H100 PCIe from $1.99/hr. Pre-configured templates make startup faster.
- **TensorDock** (managed marketplace): RTX 4090 from $0.35/hr, A100 from $0.75/hr, H100 SXM5 from $2.25/hr on-demand or $1.91/hr spot. Full VM control.
- **JarvisLabs**: RTX 4090 at $0.59/hr, L4 at $0.44/hr. Simple per-minute billing, persistent workspaces. Good when you want managed simplicity at near-marketplace prices.

### 5.3 What not to use

- **AWS** for GPU compute. Free tier doesn't cover GPU; on-demand pricing is ~3x neo-cloud rates; egress fees punish you for moving data out. Useful only if you need compliance, IAM integration, or AWS Activate startup credits ($1k-$100k requiring accelerator affiliation).
- **AWS, Azure, GCP egress** in general. Push your harmonized AnnData files to Hugging Face Datasets (free public hosting) or Cloudflare R2 (zero egress fees) instead of paying repeated S3 egress every time you spin up a new GPU instance.

### 5.4 Realistic total project budget

If you stay on the plan:
- Phases 1-5: $0 (laptop + free CPU tiers).
- Phase 6 foundation model embedding extraction: $0 (Kaggle weekly free, run once and cache).
- Phase 6 foundation model fine-tuning (optional): ~$10-30 if you decide it matters, on Vast.ai or GCP credits.
- Phase 8 factorized model hyperparameter sweep: ~$15-30 on Vast.ai over a weekend.
- Phase 9 few-shot adaptation: ~$10-20 on Vast.ai reusing the same instance.
- Buffer for restarts, ablations, reviewer requests: ~$30-50.

Total realistic spend: $50-150 of actual money, plus $300-500 of free credits used strategically. Possibly $0 if you're patient and disciplined about Kaggle weekly time.

### 5.5 Operational note on cheap GPU

Cheap marketplace GPUs are only cheap if your training loop survives interruption. Set up:
- Checkpoint saves every 5-10 minutes during training.
- Resume-from-checkpoint logic that's tested before you launch a long run.
- All datasets pulled to local SSD storage on the instance at session start, not streamed from object storage during training.
- wandb logging configured so if the instance dies, your experiment history persists.

Without this, the cheap GPU costs you more in restart and re-debugging time than the savings versus a managed cloud.

---

## 6. Reference resources

### 6.1 Gene sets

- ISG curation: Interferome 2.0 database (`interferome.its.monash.edu.au`)
- Alternative: Mostafavi lab ISG list (published in Mostafavi et al. 2016 Cell)
- Viral entry receptors per virus: ACE2/TMPRSS2 (SARS-CoV-2); sialic acid synthesis pathway and TMPRSS family (influenza); CX3CR1/nucleolin (RSV); HVEM/nectins (HSV-1)
- Type I IFN signaling: REACTOME R-HSA-913531; KEGG hsa04060

### 6.2 Comparison method repos

- theislab/scgen
- ArcInstitute/state
- ArcInstitute/stack
- jkobject/scPRINT
- bowang-lab/scGPT
- ctheodoris/Geneformer

### 6.3 Foundation model checkpoints

- Hugging Face: search for `state`, `scgpt`, `geneformer`, `scprint`
- Arc Virtual Cell Atlas page for STATE and Stack
- Verify license before use — most are noncommercial or CC0; some require attribution.

---

## 7. Code conventions

- Type hints required on all public functions.
- Docstrings in NumPy style.
- No emojis in code or comments.
- Comments explain *why*, not *what*. Reserve them for non-obvious decisions (e.g., why a specific QC threshold was chosen, why a regularization weight was set as it is).
- Configuration via hydra/yaml, no hardcoded paths.
- All randomness seeded; seed exposed in config.
- All experiments logged to wandb with full hyperparameter dump.
- Tests for data-loading and metric functions at minimum.
- Code should be CPU-runnable wherever possible. Use `torch.device("cuda" if torch.cuda.is_available() else "cpu")` patterns so the same scripts work on laptop and cloud.

---

## 8. Early-warning failure modes

Watch for these and stop if any occur:

1. Phase 4 gate fails (response vectors too correlated across viruses). Stop, write up the negative result, save months.
2. Datasets cannot be reasonably harmonized because of incompatible infection-status definitions across studies. Reduce scope to fewer studies per virus or fewer viruses.
3. Foundation model checkpoints turn out to be commercial-license-only and you cannot use them for this work. Drop those comparisons; the benchmark still stands with simpler models.
4. Within-virus performance of all methods is already at ceiling (R^2 > 0.95). The cross-virus problem is then less interesting because there is no gap to close. Reframe toward distribution-matching or specific gene-set recovery as the metric.
5. Laptop RAM proves insufficient for combined dataset and Oracle Always Free ARM box doesn't fit your workflow (e.g., you need x86 binaries). At this point a Hetzner CCX23 (~€26/month for 4 vCPU and 16GB) or similar is the cheapest fix.
6. **Sensitivity analyses for METHODS_CHOICES Issues 2, 3, 6, 7 produce qualitatively different cross-virus results across alternative choices.** Document the inconsistency. If headline finding depends on which bucket granularity / metric / harmonization method / per-cell-type vs global protocol is used, the conclusion is method-fragile, not biology-driven. Consider narrowing v1 scope to the subset of methodological choices under which the result is robust (e.g., monocyte primary even if B/CD8T are unstable; Pearson primary even if MMD diverges) — and report the fragility honestly rather than hiding it.

---

## 9. First week, concrete tasks

> **Status (2026-05-10):** Steps 1, 2, 4, 5, 6 DONE. Step 3 (free-tier signups) tracked in `SETUP.md` — human-only action, deferred per `SETUP.md` activation notes.

If you sit down to start tomorrow morning, do these in order:

1. ✅ `git init`, set up Python env with uv on your laptop. Install scanpy and dependencies. CPU-only torch wheel for now. — *Done; uv-managed env w/ Python 3.11 pinned, scanpy 1.11.5, anndata 0.12.13, scvi-tools 1.4.2, torch 2.11.0 (MPS available).*
2. ✅ Write a one-screen `README.md` stating the project goal in your own words. — *Done.*
3. ⬜ Sign up for free credit accounts: Google Cloud, Azure, Oracle Always Free, RunPod, Kaggle. Don't activate the GCP $300 timer yet — wait until you're closer to Phase 6 or Phase 8 when you'll actually use it. (For GCP, activating starts the 90-day clock.) — *Tracked in `SETUP.md`, awaiting human action.*
4. ✅ Open `notebooks/01_data_exploration.ipynb`. Identify one SARS-CoV-2 dataset and one influenza dataset on cellxgene that you can download in under an hour each. Read them into AnnData. — *Done; used Lee et al. 2020 (`de2c780c`) which contains both viruses in one harmonized study.*
5. ✅ Compute the mean expression of infected vs mock cells for both. Compute the Pearson correlation between the two viruses' response vectors. This is a 30-minute version of the Phase 4 gate check on a small subset of data. — *Done. Pearson r = 0.4607, Spearman r = 0.3091.*
6. ✅ If correlation < 0.7, you have signal and the project is real. Continue with the full data acquisition phase. If > 0.9, stop and reconsider. — *PROCEED. r = 0.46 << 0.7. Top-100 up-regulated gene Jaccard = 0.10 confirms dominant induced programs are virus-specific (i.e., transfer is non-trivial).*

That single sanity check is the highest-leverage hour of work in the entire project. Do it before writing any model code, before signing up for paid GPU, before harmonizing the full dataset. On your laptop, with no GPU, on freely available data. The first non-trivial compute decision of the project comes weeks later.
