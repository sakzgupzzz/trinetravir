# Session 4 prompt — Issue 6 scVI sensitivity + GPU setup + Issue 23 compute-envelope

**Status:** Drafted 2026-05-11 after Session 3.5 closure. Reflects post-Session 3.5 state. Session 4 is the next pipeline block (block #8 of 8 per pipeline rev 3).

This session closes Issue 6 (Harmony as batch correction method) — the only methodological open issue at Session 3.5 close. Session 3.5 deferred Issue 6 to Session 4 because scVI sensitivity requires GPU compute. Session 4 also gates the Issue 23 foundation model baselines (Geneformer + scGPT) on a compute-envelope assessment performed after scVI sensitivity completes.

**Integration choice:** standalone file pattern per Session 7 + Session 3.5 precedent. Commit as `references/session_4_prompt.md`.

---

## Entry conditions

**Issues at Session 4 entry:**
- 25 fully resolved (Issues 1–5, 7–17, 25–33)
- 7 at pre-spec level (Issues 18–24, validation pending at Phase 5/7/9)
- 1 open: Issue 6 (THIS SESSION 4)

**Pipeline state:** 7/8 done. Session 4 is the final block before Phase 4 modeling launches.

**Empirical anchors from Sessions 5–7 + 3.5 that inform Session 4 decisions:**

- Per-cell-type Harmony chosen as primary per Issue 7 (validated heuristic; calibrated resolution deferred to Session 3, then absorbed into Session 5 audit framework). Global Harmony sensitivity exists on disk.
- 5 per-cell-type Harmony embeddings persisted in `data/processed/harmony_per_celltype_<bucket>.h5ad` (commit 2026-05-11). Global Harmony in `data/processed/harmony_global_embedding.h5ad`.
- **Per-bucket + global Harmony h5ads do NOT preserve raw counts** (empirically verified 2026-05-11: all per-bucket files have `layers=[]`, `raw=False`; global has only `layers['X_harmony_scaled_hvg']`). Raw counts live in per-study reannotated h5ads at `data/processed/<study_id>_reannotated.h5ad` (wilk_2020, lee_2020, arunachalam_2020, schulte_schrepping_2020). scVI input extraction script needed (Part A.1).
- Khatri MVS canonical 86-gene Table S2 high-confidence core list at `references/khatri_mvs_gene_list.csv` (Session 3.5 commit `c44d1ff`).
- Session 7 BIOLOGY_DOMINANT threshold Δr=0.08 monocyte MVS (pre-vs-post Harmony) anchors the verdict thresholds in this session.
- scvi-tools to be pinned in `configs/methods_versions.yaml` at Session 4 entry per Issue 23.

**Literature anchors informing Session 4 pre-specs:**

- scIB benchmark (Luecken 2022 *Nat Methods*): on human immune cell integration, Scanorama, FastMNN, scANVI, and Harmony perform best. scVI competitive but not dominant on immune data.
- Briefings in Bioinformatics 2022: "The performance of scVI is robust to cell populational imbalance but sensitive to magnitude of batch effects... The performance of Harmony is similar to scVI." Direct evidence scVI ≈ Harmony on immune integration.
- Hyperparameter benchmark (Feb 2026 bioRxiv): n_latent is the most sensitive scVI parameter on immune datasets; n_hidden marginal; n_layers shows BC-vs-bio trade-off.

---

## SESSION 4 — Issue 6 scVI sensitivity + GPU setup + Issue 23 compute-envelope

GPU-required. Decision-pre-specification + compute work. Roughly 4-6
hours focused work + ~7-14h GPU wall-time, distributed across 2-3
chat sessions for spec drafting + audit gate.

Prerequisite: Session 3.5 closed cleanly (33/33 issues at pre-spec
or resolved). GPU environment provisioned.

CONSTRAINT: Issue 6 scVI sensitivity + Issue 23 compute-envelope
check only. No Phase 4 modeling work. No new methodological
infrastructure. Stop and report after audit gate.

## PART A — scVI per-bucket primary sensitivity

Mirror Session 5/7 calibration discipline. Pre-spec all decisions
before compute. Apply Issue 14 hyperparameter policy (held-out
donor validation, 20-config budget per method).

---

### A.1 Input data preparation — `scripts/extract_per_bucket_counts.py` (NEW script)

Per-bucket + global Harmony h5ads do NOT preserve raw counts (verified 2026-05-11). Raw counts must be extracted from per-study reannotated h5ads and subset to bucket cells via the cell-id positional mapping. Pattern reusable from `scripts/session7_part_a_pre_post_harmony.py:compute_pre_harmony_rv_per_study` (lines ~80-145).

**Script: `scripts/extract_per_bucket_counts.py` (NEW; not a stub — write this before Part A.2 compute begins).**

Recipe:

```
INPUTS per bucket:
  data/processed/harmony_per_celltype_<bucket>.h5ad
    → obs.index format: '<int>-<study_id>' (positional row idx in study h5ad)
    → uns['hvg_genes']: 4000-gene HVG list

  data/processed/<study_id>_reannotated.h5ad (one per study)
    → wilk_2020_reannotated.h5ad
    → lee_2020_reannotated.h5ad
    → arunachalam_2020_reannotated.h5ad
    → schulte_schrepping_2020_reannotated.h5ad
    → X = raw counts (max > 20 confirms; not log-normalized)
    → var has 'feature_name' or 'gene_symbol' column for symbol resolution

ALGORITHM per bucket:
  1. Load harmony_per_celltype_<bucket>.h5ad to get hvg_genes + obs.
  2. For each study in obs['study_id']:
       a. Load data/processed/<study_id>_reannotated.h5ad.
       b. Resolve var symbols (priority: gene_symbol, feature_name, name, gene_symbols, symbol).
       c. Strip suffix '-<study_id>' from harmony obs.index → positional row idx.
       d. Slice study h5ad at those positional indices.
       e. Subset to bucket's HVG genes (intersection with study's var_names).
       f. Append obs columns: study_id, donor_id, donor_disease_status (from harmony obs).
       g. Append to per-bucket accumulator.
  3. Concat across studies → single AnnData with:
       X = raw counts (np.int32 or float32 sparse)
       obs = study_id, donor_id, donor_disease_status, original cell-id
       var = bucket's HVG (4000 genes)
  4. Write data/processed/scvi_input_<bucket>.h5ad with compression='gzip'.

OUTPUT per bucket:
  data/processed/scvi_input_<bucket>.h5ad
    → X = raw counts on 4000 HVG
    → obs columns: study_id, donor_id, donor_disease_status, cell_id
    → var index: HVG gene symbols
```

**Verification gate (post-extraction):**
- Each `scvi_input_<bucket>.h5ad` has cell count matching the bucket harmony h5ad: monocyte 68,672; CD4T 42,705; CD8T 29,855; B 26,115; NK 29,488.
- `X.max() > 20` confirms counts (not log-normalized).
- `obs['study_id'].nunique() == 4` (arunachalam, lee, schulte_schrepping, wilk).
- HVG list matches `uns['hvg_genes']` schema from source harmony h5ad.

Wall-time: ~10-15 min total for all 5 buckets. Idempotent. Atomic commit before A.2.

---

### A.2 Hyperparameter search space (pre-spec)

Open Issue 34 (scVI comparison design pre-spec) in
METHODS_CHOICES.md before compute begins. Five-field structure
per Issue 17. Decision committed atomically.

Search space (16 configurations, 100% covered by Issue 14's
20-config budget with 4 to spare):

```
n_latent     ∈ {10, 20, 30, 50}   # most sensitive on immune
                                   # (Feb 2026 bioRxiv benchmark)
n_hidden     ∈ {128, 256}         # 256 best for PBMC per
                                   # Lopez 2019 hyperopt
n_layers     ∈ {1, 2}             # captures BC-vs-bio trade-off
                                   # (Feb 2026 bioRxiv)
```

Fixed at scvi-tools defaults:

```
dropout_rate       = 0.1
learning_rate      = 1e-3
weight_decay       = 1e-6
optimizer          = Adam
gene_likelihood    = 'zinb'        # appropriate for PBMC
dispersion         = 'gene'        # per-batch is harder to
                                   # train reliably
latent_distribution = 'normal'
n_epochs_kl_warmup = 400
batch_key          = 'study_id'
```

Reproducibility:
- `random_state = 42` (numpy + torch + scvi seeds)
- `max_epochs = 400` with early stopping
- `patience = 50`
- `monitor = 'reconstruction_loss_validation'`

Per-bucket: train 16 configurations, select best on within-
bucket held-out donor validation (80/20 donor split, stratified
by donor_disease_status to preserve diseased/healthy proportions).

---

### A.3 Output space + normalization protocol (pre-spec)

Committed in Issue 34. Recipe:

```python
scvi_normalized = model.get_normalized_expression(
    library_size=1e4,
    return_numpy=True
)
scvi_log = np.log1p(scvi_normalized)
scvi_scaled = sc.pp.scale(scvi_log, zero_center=True,
                          max_value=10, copy=True)
```

Matches Harmony's `X_harmony_scaled_hvg` layer normalization protocol (log1p + scaled HVG). Without matching normalization, Δr partly measures scale difference rather than integration difference.

Rationale per scvi-tools community standard: `library_size=1e4` puts scVI output on TPM-like scale comparable to scanpy's `normalize_total(target_sum=1e4)` followed by `log1p`. Default `library_size=1` would give raw frequencies (unsuitable for direct comparison against Harmony's HVG-scaled output).

---

### A.4 Response vector computation + calibration framework

Per bucket, per selected hyperparameter configuration:

1. Compute per-(study, donor_disease_status) response vector:
   ```
   response_vec[study][status] = mean(scvi_scaled, axis=0)
                                 over cells in (study, status)
   ```

2. Mean off-diagonal Pearson r across 4-study × 2-status pairs
   (Issue 3 primary metric).

3. Apply calibration framework (Issues 8-11 + 26 v2):
     - permutation null N=1000 (Issue 8)
     - bootstrap CI N=1000 on observed r (Issue 9 v2)
     - within-study split-half ceiling N=50 (Issue 8)
     - FDR-BH correction across buckets (Issue 26)

4. Per-metric supplementary verdicts (Spearman, DE-Jaccard
   top-100, MMD-RBF median heuristic).

5. MVS-subset Pearson r using
   `references/khatri_mvs_gene_list.csv` (Issue 18 amendment) —
   this is the calibrated anchor matching Sessions 5/6B/7
   empirical defenses.

---

### A.5 Verdict computation (pre-committed thresholds)

Pre-committed in Issue 34 BEFORE any scVI training runs.

Per bucket: signed Δr (scVI minus per-bucket Harmony) on primary
metric (Pearson MVS-subset). Aggregate verdict from per-bucket
signed Δr per the four-tier rule below.

**Tier I — HARMONY_ADEQUATE:**
- `max(Δr) ≤ 0.05` across all 5 buckets
- → Harmony stays as v1 primary. No re-run.
- → Anchored against literature consensus that scVI ≈ Harmony on immune integration (scIB Luecken 2022; Briefings 2022).

**Tier II — MIXED:**
- At least one bucket Δr ∈ (0.05, 0.10], no bucket Δr > 0.10.
- → Harmony adequate for v1; scVI flagged for v2 with explicit documentation in manuscript discussion section.
- → 0.05 anchored against Session 7 BIOLOGY_CONSISTENT band; 0.10 anchored against ~1.25× Session 7 BIOLOGY_DOMINANT threshold (Δr=0.08).

**Tier III — SCVI_PREFERRED:**
- ≥3 of 5 buckets show Δr > 0.10 OR any single bucket Δr > 0.20.
- → scVI replaces Harmony as v1 primary. Sessions 5/6B/7 calibration re-run on scVI output before Phase 4. ~1-2 weeks additional work.
- → 0.20 anchored against ~2.5× Session 7 BIOLOGY_DOMINANT threshold. Would be surprising relative to scIB consensus and Briefings 2022 finding of scVI ≈ Harmony; warrants investigation.

**Tier IV — HARMONY_PREFERRED:**
- Δr < -0.10 on ≥3 of 5 buckets.
- → No action; supports current methodology choice. Document in supplementary.
- → -0.10 symmetric with Tier III's +0.10 to avoid asymmetric-threshold cognitive trap. Asymmetric Tier I (+0.05) vs Tier IV (-0.10) reflects different decision-action symmetry: Tier I avoids costly re-runs at small lifts; Tier IV needs a meaningful Harmony advantage to assert support.

**Verdict tie-break rule (rare boundary cases):**
If verdict matrix produces ambiguous classification (e.g., exactly 2 buckets with Δr > 0.10), default to the more conservative tier (Tier II MIXED over Tier III SCVI_PREFERRED). Document the boundary case in audit.

---

### A.5b Wilk sequencing-depth watchpoint (informational)

Pre-Part-A extraction (commit `f70b0ed`, 2026-05-11) revealed Wilk has X.max range 158-293 across buckets vs Arunachalam/Lee/Schulte X.max 1500-3600. Wilk's per-cell sequencing depth is ~5-10× lower than the other three studies.

**Not a blocker.** scVI's NB/ZINB likelihood handles per-cell library size internally via `size_factor`; `get_normalized_expression(library_size=1e4)` normalizes output for downstream comparison. At the response-vector level (mean across many cells per (study, status)), depth differences should normalize out.

**Verdict-review checklist** at Part A close:
- Inspect per-study response vectors (Wilk vs each of Arunachalam/Lee/Schulte) for the selected scVI configuration per bucket.
- If Wilk's response vector shows outlier pattern (e.g., dominated by housekeeping genes; orthogonal to other studies) NOT present in Harmony's per-study response vectors, flag as a Wilk-specific scVI artifact in Issue 6 resolution.
- If Wilk's response vector is qualitatively consistent with the other 3 studies (similar gene-level direction), depth difference is properly handled and Δr verdict applies as-is.

This is a watch item, not a pre-compute decision. Document outcome in Part A verdict table as an additional column or in Issue 6 resolution caveat.

### A.6 Wall-time estimate

Per-bucket (16 configs each, single A100):
- monocyte (68K cells)    ~2h
- CD4T (42K cells)        ~1.5h
- CD8T (29K cells)        ~1h
- B (26K cells)           ~1h
- NK (29K cells)          ~1h

Total Part A: ~6.5-7h single A100. Costs ~$15-25 on Lambda/RunPod pay-as-you-go. Defensible single-overnight budget.

---

## PART B — scVI global supplementary

Single scVI training run on full v1 corpus (244,389 cells).

Input: re-use raw counts extracted via `scripts/extract_per_bucket_counts.py` recipe, but concatenated across all 5 buckets rather than per-bucket. Output: `data/processed/scvi_input_global.h5ad` (244K cells × HVG union of all bucket HVG sets, or `harmony_global_embedding.h5ad` HVG list if a single global HVG list exists). Verify global HVG list source before extraction.

Hyperparameter sweep: same 16-config grid as Part A.

Output space + normalization: same recipe as A.3.

Response vector computation: per-bucket using `obs['coarse']` (NOT `cell_type_bucket` — per Issue 7 status note).

Compare to cached global Harmony per-bucket response vectors.

Verdict structure: same four-tier rule as Part A but applied to scVI_global vs Harmony_global Δr.

Reporting: supplementary table. Headline verdict from Part A; Part B addresses reviewer concern that scVI was under-resourced in per-bucket mode (since scVI's standard use is global).

Wall-time: ~2-4h additional GPU.

---

## PART C — Issue 23 compute-envelope assessment (INFERENCE-only scope)

Triggered after Part A + B complete. Pre-commit decision rule in
Issue 35 (foundation model compute-envelope decision) BEFORE
measuring foundation model wall-times.

**Scope clarification (load-bearing):** v1 Issue 23 foundation model baselines are **INFERENCE-only**. Pre-trained Geneformer + scGPT checkpoints used as-is to produce frozen embeddings; a linear head is trained per-bucket per-virus on those frozen embeddings. End-to-end fine-tuning of foundation models is **out of v1 scope**; deferred to v1.5 or v2 if reviewer feedback warrants.

Measure:
1. Actual scVI total wall-time (Parts A + B).
2. Test-run Geneformer **inference** on monocyte bucket (68K cells)
   with pretrained base checkpoint (frozen weights; embeddings only).
   Extrapolate to per-pass corpus size.
3. Test-run scGPT **inference** on monocyte bucket. Extrapolate.

**Cell-count math (per-pass):**
- v1 training corpus: 244K cells per inference pass.
- Held-out cohort eval passes (separate, one per cohort):
  - Yoshida 2022: ~168K
  - Allen Atlas monocyte: ~301K
  - GSE157829: ~36K
  - Randolph 2021: ~39K (post-Issue-31 merge with infected_monocytes)
- Per-evaluation-pass: max 301K cells (Allen Atlas) at the upper end; typical 36-244K.
- Cumulative across all eval passes ≠ what matters for single-run budgeting. A single eval cycle runs ~5-10 passes (training corpus + 4 held-out cohorts + leave-one-virus-out splits).

**Decision rule (Issue 35 pre-spec, INFERENCE-only):**

INCLUDE in v1 Phase 7 baselines if:
- Geneformer + scGPT pre-trained inference (frozen embeddings + linear head per bucket per virus) total wall-time ≤ 2 × scVI sensitivity wall-time, AND
- Total GPU cost ≤ $100.

DEFER to v1.5 if:
- Inference budget exceeded (unlikely at A100 throughput of 1-5K cells/sec; 244K cells → 1-4 min per pass), OR
- Foundation model fine-tuning is requested by reviewer (definitely exceeds v1 envelope; ~$200-1000 budget needed for full fine-tuning).

**Explicitly OUT of v1 scope:**
- End-to-end fine-tuning of Geneformer or scGPT on v1 corpus.
- Custom foundation model training from scratch.
- Foundation model architecture modifications.

v1 evaluates pre-trained models as-is. v1.5 or v2 may revisit fine-tuning if v1 reviewer feedback warrants.

**Rationale:** foundation model baselines (Geneformer, scGPT) per critique concern 2 were committed at pre-spec in Session 3.5 with the contingency that final inclusion is gated on Session 4 compute envelope. This Part operationalizes that gate.

---

## PART D — STOP AND REPORT (audit gate)

Report:

1. **Issue 34** (scVI comparison design pre-spec) opened and committed BEFORE compute. Hyperparameter grid, output space recipe, verdict thresholds, normalization protocol all in five-field structure.
2. **Part A** per-bucket scVI sensitivity:
   - 16-config sweep per bucket → 80 total runs
   - Selected configuration per bucket documented
   - Δr per bucket on primary metric (Pearson MVS-subset)
   - Calibration framework applied identically per Issues 8-11 + 26 v2
   - Verdict per Tier I/II/III/IV.
3. **Part B** scVI global supplementary:
   - Single global training run
   - Per-bucket Δr after projection to `obs['coarse']` strata
   - Supplementary verdict per same four-tier rule
4. **Part C** Issue 23 compute-envelope assessment:
   - scVI total wall-time + cost reported
   - Geneformer + scGPT extrapolated inference wall-time + cost
   - Issue 35 resolution: INCLUDE or DEFER decision
5. **Issue 6 resolution** per the verdict:
   - HARMONY_ADEQUATE → Harmony stays. No re-run.
   - MIXED → Harmony stays for v1. scVI flagged for v2.
   - SCVI_PREFERRED → scVI replaces. Re-run Sessions 5/6B/7.
   - HARMONY_PREFERRED → Harmony stays. Documented support.
6. **State after Session 4**:
   - 33 issues resolved (1-5, 7-17, 25-35) [if Issues 34, 35 opened + resolved in same session]
   - 0 open methodological issues
   - Pipeline 8/8 complete
   - Ready for Phase 4 modeling implementation

---

## CONSTRAINTS

- No Phase 4 work. No new methodological infrastructure beyond what Session 4 deliverables specify.
- Atomic schema-change rule (Issue 17) applies. Each Issue gets its own atomic commit. Spec document, Issue 34 pre-spec, scvi-tools version pin, `scripts/extract_per_bucket_counts.py` + extraction outputs, Part A results, Part B results, Issue 6 resolution, Part C compute envelope, Issue 35 resolution, audit gate — each its own commit.
- Pre-registration discipline: Issue 34 + Issue 35 pre-specs committed BEFORE any compute begins. Decision rules locked.
- No post-hoc threshold adjustment. If a verdict lands at a boundary case, apply the conservative tie-break rule. Don't re-tune thresholds.
- scvi-tools version pinned in `configs/methods_versions.yaml` per Issue 23 before Part A compute. Pin format: exact pip-installable version (e.g., `scvi-tools==1.2.0`) + matching torch + scanpy + anndata versions.
- GPU environment reproducibility: image, driver, CUDA version, Python version all documented in `configs/gpu_environment.yaml`. Pinned for Session 4 reproducibility.

---

## Issues opened in Session 4

**Issue 34 (NEW) — scVI comparison design pre-spec:**
- Decision: hyperparameter search space (16 configs), output space + normalization recipe, four-tier verdict thresholds.
- Rationale: literature-anchored per scIB Luecken 2022 + Briefings 2022 + Feb 2026 hyperparameter benchmark + Lopez 2019 hyperopt. Specific thresholds anchored against Session 7 BIOLOGY_DOMINANT Δr=0.08 internal benchmark.
- Validation: Part A + B verdict at Session 4 close.
- Status: open at pre-spec; resolution at Session 4 close.

**Issue 35 (NEW) — Foundation model compute-envelope decision (INFERENCE-only):**
- Decision rule: INCLUDE if Geneformer + scGPT pre-trained inference (frozen embeddings + linear head) total wall-time ≤ 2× scVI sensitivity wall-time AND cost ≤ $100. DEFER to v1.5 otherwise.
- Out of v1 scope: end-to-end fine-tuning of foundation models.
- Rationale: Session 3.5 Issue 23 contingency operationalized. Avoids committing v1 to compute it can't fit.
- Validation: Part C measurement at Session 4 close.
- Status: open at pre-spec; resolution at Session 4 close.

---

## Atomic commit sequence (10 logical commits; physical commits may pair tightly coupled gates)

1. `references/session_4_prompt.md` (this document)
2. Issue 34 opening in METHODS_CHOICES.md
3. Issue 35 opening in METHODS_CHOICES.md
4. `configs/methods_versions.yaml` scvi-tools pin + `configs/gpu_environment.yaml`
5. `scripts/extract_per_bucket_counts.py` + 5 per-bucket scvi_input h5ads (Part A.1)
6. Part A scVI per-bucket sensitivity results + verdict table
7. Part B scVI global supplementary results + verdict table
8. Issue 6 resolution (one of four verdicts)
9. Part C Issue 23 compute-envelope assessment + Issue 35 resolution
10. Session 4 pipeline closure (mirrors commit `ec3dfe9` from Session 3.5)

**Physical-vs-logical commit mapping (executed 2026-05-11):**
- Logical #1 = `d3b2a6b` (this spec)
- Logical #2 + #3 = `1120716` (paired physical commit; Issues 34+35 are tightly coupled pre-compute gates with same trigger — must both commit before Parts A/B/C compute; splitting would just split a cohesive section header per Issue 17 logical-unit principle).
- Logical #4 = `e8efca2` (configs)
- Logical #5 = `f70b0ed` (extract_per_bucket_counts.py)
- Logical #6 → TBD (Part A GPU)
- Logical #7 → TBD (Part B GPU)
- Logical #8 → TBD (Issue 6 resolution)
- Logical #9 → TBD (Part C GPU)
- Logical #10 → TBD (pipeline closure)

Fresh-chat Claude should reference logical commit numbers (1-10) when planning next steps; task list metadata preserves this numbering. Physical commit hashes are tracked separately and may diverge from logical numbering when tightly coupled commits are paired.

---

## Timeline estimate

- Spec drafting + Issue 34/35 opening: 1 chat session (~2h)
- GPU environment setup + scvi-tools pin + `extract_per_bucket_counts.py` script + 5 extractions: 1 chat session (~2-3h)
- Parts A + B compute: ~9-11h GPU wall-time (mostly unattended)
- Part C compute-envelope measurement: ~2-3h GPU
- Audit gate + verdict commit: 1 chat session (~1h)
- **Total: 3-4 chat sessions over 1-2 calendar days, depending on GPU availability**

---

## Pattern relation to Sessions 5/7

Session 4 follows the same pattern as Sessions 5 and 7:
- Pre-registered decision rules committed before compute begins
- Calibration framework applied identically per Issues 8–11 + 26 v2
- Per-bucket evaluation with MVS-subset as calibrated anchor
- Audit gate before pipeline progresses
- Conditional outcomes documented with action paths

---

## How to start Session 4

In fresh Claude Code session:

1. Read `SHORT_TERM_PLAN.md` (auto-loaded via memory pointer); confirm Block #8 NEXT.
2. Read `references/session_4_prompt.md` (this document).
3. Open Issue 34 + Issue 35 in `METHODS_CHOICES.md` (atomic commits per Issue 17).
4. Write `scripts/extract_per_bucket_counts.py` and run extraction (Part A.1).
5. Pin `configs/methods_versions.yaml` + `configs/gpu_environment.yaml`.
6. Launch Part A 16-config sweep per bucket.

Launch trigger: "Launch Session 4 per `references/session_4_prompt.md`. Open Issues 34 and 35 with the pre-committed decision rules verbatim before running any compute. Then write `scripts/extract_per_bucket_counts.py` and run extraction. Then begin Part A sweep."
