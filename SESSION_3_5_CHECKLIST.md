# Session 3.5 — execution checklist (DO NOT START until Session 3 closeout + Session 5 BOTH complete)

**Ordering**: Session 3 closeout → **Session 5 (audit response)** → Session 3.5 → Session 4. Session 5 spec at `SESSION_5_SPEC.md`. Session 5 OPENS Issue 25 (v1 paper framing decision) and STOPS for human decision. Session 3.5 BLOCKED until Issue 25 is resolved.

**Pre-flight gate** — before launching Session 3.5, confirm ALL of:

- [ ] `phase35_subbucket` calibration bg job (`bfs5zrjht`) completed; `results/tables/calibration_phase35_subbucket.csv` populated.
- [ ] Issue 2 resolution written to METHODS_CHOICES.md citing subbucket calibration evidence.
- [ ] Subbucket CSV + Issue 2 resolution committed in a single atomic commit (per Issue 17 rule).
- [ ] **Session 5 executed and completed** (calibration framework v2, Khatri MVS validation, Issues 3+7 revised, Issues 25+26 opened/resolved).
- [ ] **Issue 25 (v1 paper framing) resolved by human decision** (Option A reframe vs Option B acquire more data).
- [ ] If Option B was chosen, additional viral data acquired and harmonized before Session 3.5 starts.
- [ ] METHODS_CHOICES.md state: 17 issues resolved (1, 2, 3-revised, 4, 5, 7-revised, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 26). Open: Issue 6 (Session 4) + Issue 25 (human-decided). Session 3.5 may re-scope based on Issue 25 decision.
- [ ] `git status` clean (no uncommitted Session 3 / 5 work).

If any unchecked, do not start. Resume Session 3 closeout or Session 5 first.

---

## Part A — Pre-specify Issues 18-24 in METHODS_CHOICES.md (one commit)

Insert into "Open issues requiring immediate resolution" section AFTER Issue 17 stub (if any) BEFORE "Resolved Issues" section header. Each entry uses the five-field structure (status / decision / why / alternatives / validation / dates).

- [ ] **Issue 18 — ISG gene set source** (Phase 5 prereq, LOAD-BEARING for ISG-aware regularization)
  - Decision: Interferome 2.0 canonical type-I-IFN-induced genes, filtered to high-confidence subset (≥2-fold induction in PBMC studies, type I IFN as inducer).
  - Alternative: Mostafavi 2016 Cell list as supplementary sensitivity at Phase 5 gate.
  - Validation: sensitivity at Phase 5; supplementary figure showing cross-virus results robust to ISG list choice.
  - Status: open at pre-specification; resolution at Phase 5.

- [ ] **Issue 19 — Pathway gene set source** (Phase 5 prereq)
  - Decision: REACTOME R-HSA-913531 (interferon signaling) as primary source.
  - Graph construction: undirected adjacency, immediate co-members only, drop genes not in HVG space, no transitive expansion.
  - Validation: sensitivity at Phase 5; if pathway-aware weight tunes to ~0 under Issue 14 held-out validation, document and consider dropping the term.
  - Status: open at pre-specification; resolution at Phase 5.

- [ ] **Issue 20 — Reconstruction loss (LOAD-BEARING)** (Phase 5 prereq)
  - Decision: MSE on response vectors as primary; NB on counts as sensitivity.
  - Rationale: response-vector aggregation consistent with v1 pipeline; per-cell perturbed/baseline pairing for NB across studies exceeds v1 scope.
  - Validation: train both at Phase 5; supplementary comparison; switch headline if NB shows substantially better cross-virus transfer.
  - Status: open at pre-specification; resolution at Phase 5.

- [ ] **Issue 21 — Factorized model architecture hyperparameters** (Phase 5 prereq)
  - Decision: pre-specified search space within Issue 14's 20-config budget:
    - Shared latent: {16, 32, 64}
    - Virus embedding: {8, 16, 32}
    - Encoder/decoder depth: {2, 3}
    - Dropout: {0.1, 0.2, 0.3}
    - Activation: GELU (fixed)
    - Optimizer: Adam, lr ∈ {1e-3, 5e-4}, weight_decay=1e-5
  - Selection: held-out donor validation per Issue 14, donor-level split, Pearson (Issue 3) as tuning target.
  - Validation: report final hyperparameters with selection criterion in supplementary.
  - Status: open at pre-specification; resolution at Phase 5.

- [ ] **Issue 22 — Few-shot adaptation protocol (LOAD-BEARING for H5)** (Phase 9 prereq)
  - Decision:
    - Sample sizes per virus per run: 50, 100, 200, 500, 1000.
    - Seeds: 5 per (sample_size, virus).
    - Frozen: f_shared + f_specific + existing virus embeddings.
    - Trained: new virus embedding only via Adam lr=1e-3, early stop on 20% held-out fraction.
    - Selection: random sampling without replacement from target virus cells.
    - Held-out eval: remaining cells, stratified by cell-type bucket.
  - Validation: data-efficiency curves with mean±SD across seeds in Phase 9; per-bucket data-efficiency reported.
  - Status: open at pre-specification; resolution at Phase 9.

- [ ] **Issue 23 — Comparison method versions and reproducibility** (Phase 7 prereq)
  - Decision: pin exact versions in `configs/methods_versions.yaml` before Phase 7 launch. Use each method's most recent stable release. Foundation model checkpoints pinned by HuggingFace revision hash.
  - Implementation: published defaults as starting point, tuned per Issue 14 (held-out validation, 20-config budget per method).
  - Wrapper code: modifications to original training loops documented in `src/trinetravir/methods/<method>_wrapper.py` with rationale.
  - Validation: reproducibility from pinned versions + released code + released corpus.
  - Status: open at pre-specification; resolution at Phase 7.

- [ ] **Issue 24 — Baseline implementations** (Phase 5 prereq)
  - Decision:
    - Predict-mean: per-cell-type per-virus mean across training cells in HVG space.
    - Linear-delta: ridge on baseline expression + cell-type one-hot + virus one-hot; alpha cross-validated within 20-config budget.
    - KNN: cosine distance on log-normalized HVG, k=10 distance-weighted, within-virus training neighbors only for cross-virus eval per Issue 15. k sensitivity at k=5 and k=20 in supplementary.
  - Stub files (atomic with this issue's resolution per Issue 17):
    - `src/trinetravir/baselines/__init__.py` (new package)
    - `src/trinetravir/baselines/predict_mean.py` (stub + docstring)
    - `src/trinetravir/baselines/linear_delta.py` (stub + docstring)
    - `src/trinetravir/baselines/knn.py` (stub + docstring)
  - Validation: pre-specified baselines reported at Phase 5 as the bar to beat.
  - Status: open at pre-specification; resolution at Phase 5.

**Commit A**: `Pre-specify Issues 18-24 (Phase 5/7/9 prerequisites) in METHODS_CHOICES.md + stub baseline files per Issue 24.`

---

## Part B — METHODS_CHOICES.md hygiene cleanup (same commit as A or separate)

Replace `<fill in>` in `Date resolved:` for the original-open entries of these issues (each has a corresponding "Resolved Issue X" entry later with full content + 2026-05-11 date):

- [ ] Issue 8 (line ~193): `**Date resolved**: <fill in>` → `**Date resolved**: see Resolved Issue 8 entry dated 2026-05-11`
- [ ] Issue 9 (similar line near end of original entry)
- [ ] Issue 10 (similar)
- [ ] Issue 11 (similar)

Verify Issues 3, 5, 7, 12, 16 do NOT have `<fill in>` placeholders remaining (those were already updated by Session 3 with "pending Session 3" or were never opened with `<fill in>`).

**No content changes** — pointer-only fix.

---

## Part C — PLAN.md v1.3 integration (one commit, separate from Part A)

Prerequisite: read existing PLAN.md to determine if §4 Phase 5 already contains full factorized model architectural spec. If yes, fold Addition 3 content into §4 instead of creating §1.7.

- [ ] **Addition 1**: expand §1.1 Scientific motivation with §1.1.1-1.1.4 subsections (immunology, methodological gap, evaluation gap, PBMC-only scope rationale). Replace existing single-paragraph §1.1 contents.
  - Source: `plan_additions_v13.md` (if exists) — verify provided
  - If source file missing, request from user before proceeding

- [ ] **Addition 2**: new §1.6 Related work and competitive positioning. Insert after §1.5 v1.5 forward planning block, before `---` separator starting §2.

- [ ] **Addition 3 (conditional)**:
  - Check `§4 Phase 5` in PLAN.md for existing architectural spec
  - **If §4 already has full architecture detail**: fold Addition 3 content INTO §4 Phase 5; do NOT add §1.7 as standalone subsection. Note this folding in commit message.
  - **Else**: add as new §1.7 Factorized model architecture specification, after §1.6, before `---` separator starting §2.

- [ ] Update v1.2 changes block at top of PLAN.md to v1.3:
  - `v1.3 changes (2026-05-11): Expanded §1.1 scientific motivation with immunology and calibration rationale. Added §1.6 related work and competitive positioning. Added §1.7 factorized model architecture specification (or integrated with §4 if already present).`

**Commit C**: `Expand PLAN.md to v1.3. Add detailed scientific motivation, related work positioning, and factorized model architecture specification. Does not change v1 scope, deliverables, non-goals, or hypotheses.`

---

## Part D — Final report

Report items to confirm in the closing message:

1. [ ] Issues 18-24 resolutions documented in METHODS_CHOICES.md with five-field structure. Date resolved: 2026-05-11. Status: "open at pre-specification level; final validation at Phase 5 / 7 / 9".
2. [ ] METHODS_CHOICES.md hygiene cleanup applied to Issues 8-11 `<fill in>` placeholders.
3. [ ] PLAN.md v1.3 integrated. Confirmation whether §1.7 was added as new subsection or folded into §4.
4. [ ] Stub baseline files created at `src/trinetravir/baselines/` per Issue 24 decision.
5. [ ] State after Session 3.5:
   - 16 issues resolved (1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)
   - 7 issues resolved at pre-specification level (18, 19, 20, 21, 22, 23, 24)
   - 1 open: Issue 6 (Session 4 GPU work)
   - Total: 23 issues at pre-specification level or resolved; 1 remaining for Session 4
6. [ ] Recommended next step: GPU setup + Session 4 launch.

---

## Constraints (reminder)

- No calibration runs. Session 3 completed all calibration work needed for v1.
- No Phase 4 work. No Phase 5 implementation work. Stubs only for baselines.
- Atomic schema-change rule (Issue 17) applies to PLAN.md and baseline stub creation; each is its own commit.
- Pre-specifications in Part A are documentation only. Phase 5 implementation happens after Session 4 closes Issue 6.
- METHODS_CHOICES.md additions and PLAN.md additions land in separate commits.
- If §4 Phase 5 in PLAN.md already specifies factorized architecture in detail, fold Addition 3 into §4 rather than adding §1.7 standalone.

## Time budget

2-3 hours focused. Decision documentation + file integration; no compute beyond stub file creation.

## Open question for user before launching

`plan_additions_v13.md` is referenced as source for the three PLAN.md additions. Verify the file exists at repo root or is provided in a follow-up before starting Part C. If missing, request user to either paste the content into the launch prompt or add the file to the repo.
