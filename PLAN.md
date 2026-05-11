# Cross-Virus Generalization for Single-Cell Host Response Prediction

Project plan v1.2. Designed as the working specification for a solo computational research project producing a bioRxiv preprint plus open-source code release. Timeline target: 12–14 weeks of focused work. Compute strategy: CPU-first, escalate to GPU only when CPU becomes prohibitive.

**v1.2 changes (2026-05-10):**
- Repo and Python package renamed `trinetravir` (was `cross-virus-scrna` / `crossvirus`).
- License: MIT.
- v1.1 scope locked to **PBMCs only** across all viruses. Multi-compartment (airway epithelium, organoids, intestinal) deferred to v2.
- §9 mini-gate executed and **PASSED**: Lee et al. 2020 PBMC (cellxgene `de2c780c`), Pearson r = 0.46 between SARS-CoV-2 and IAV response vectors, top-100 up-regulated gene Jaccard = 0.10. Signal present.
- Dataset acquisition (Phase 2) updated with cellxgene Census availability check: SARS-CoV-2 PBMC well-represented, IAV PBMC limited to Lee et al. (Census), RSV and additional IAV require GEO direct download.
- **Schema rename: `infection_status` → `donor_disease_status` (v1.1, PBMC-only).** Values changed from `infected`/`mock` → `diseased`/`healthy`. Reason: PBMCs from a virally-infected donor are mostly NOT directly infected at the cellular level — the signal is systemic cytokine / IFN response, not cell-autonomous infection. The PLAN §2 per-cell `infection_status` (infected/bystander/mock) semantics belong to v2 airway studies where viral reads can be assigned per cell. New obs column `label_source` (e.g. `disease_proxy`) records label origin so downstream code can branch on label semantics.
- **GEO time-budget addendum (Phase 2).** Census fetches take ~10 min per study because cellxgene already harmonized them. GEO direct downloads (RSV, second IAV, DNA virus control) need raw count matrix processing, manual metadata parsing, and ad-hoc QC because the studies were never standardized — budget **2-4 hours per GEO study**, not 10 minutes. Do not conflate the two in time estimates.

---

## 1. Context

### 1.1 Scientific motivation

Single-cell foundation models (scGPT, Geneformer, scPRINT-2, STATE, Stack) and perturbation prediction methods (scGen, scCausalVI, CoupleVAE) are typically evaluated on within-virus or within-perturbation generalization: held-out cell types, held-out doses, held-out donors. Cross-virus transfer — training on virus A and predicting host response to virus B without any virus B data during training — has not been benchmarked as a standardized task despite the obvious biological premise that conserved antiviral programs (interferon-stimulated genes, type I IFN signaling, NF-kB inflammation) should transfer while virus-specific programs (entry receptor activation, viral hijack pathways) should not.

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

**Phase 3 entry-gate update (v1.2, added 2026-05-10).** Pre-Harmony cross-study within-virus r matrix on the 5 SARS-CoV-2 PBMC studies came out at mean off-diag r = 0.054, range [-0.58, +0.59]. mgh_acute_covid + guo_2020 are now excluded (1 healthy donor and 0 healthy donors respectively — see configs/datasets.yaml exclusion_reason). The remaining 4 clean studies show: lee↔wilk↔arunachalam form a coherent core (r ≈ 0.4–0.6), but schulte_schrepping_2020 sits near-zero with everything (r ≈ -0.06 to +0.21). Before launching Harmony, run `notebooks/03_celltype_stratified_consistency.ipynb` to recompute the cross-study r matrix per major cell type (mono / CD4T / CD8T / B / NK). Two diagnostic outcomes:
  - If per-cell-type r is meaningfully higher than bulk r → cell-type-composition drift is the problem; harmonize with `cell_type` as a covariate (or per-cell-type harmonization).
  - If per-cell-type r remains near-zero against schulte_schrepping → deeper protocol/cohort issue (10x version, capture chemistry, severity stratification, time-from-symptom-onset). Drop or down-weight schulte_schrepping rather than try to correct it away.

Also pull and report severity distribution per study during the same notebook. Severe COVID has substantially different PBMC signatures (dysregulated type I IFN, emergency myelopoiesis); if studies have systematically different severity distributions, the "study effect" in the r matrix is partly real biology and stratifying by severity in eventual benchmark splits is preferable to correcting it away.

**Phase 3 exit-gate (added v1.2).** Post-Harmony within-virus cross-study Pearson r should be **≥ 0.5 across all pairs** (5 SARS PBMC pairs after exclusions, or 6 if schulte_schrepping is retained). If yes, proceed to Phase 4 GATE 1 on the harmonized dataset. If no, diagnose before adding more data. **Do not attack GEO datasets (RSV / 2nd IAV / DNA control) until this gate passes** — more data into a broken harmonization pipeline just produces a noisier broken pipeline. GEO acquisition becomes Phase 3.5 / v1.5.

Deliverable: a single harmonized AnnData per virus, stored in `data/processed/{virus}.h5ad`, plus a combined `data/processed/all_viruses.h5ad`.

### Phase 4: Sanity check (End of Week 4) — GATE 1 [Tier 0]

> **Status (2026-05-10):** §9 mini-gate executed on Lee et al. 2020 PBMC (cellxgene `de2c780c`). Pearson r = 0.46 (SARS-CoV-2 vs IAV response vectors), Spearman r = 0.31, top-100 up-regulated gene Jaccard = 0.10. **PASSED** with margin (well below 0.7 threshold). The full Phase-4 analysis below still needs to run on the harmonized multi-study dataset; the mini-gate only confirms within-study signal in one anchor cohort.

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
