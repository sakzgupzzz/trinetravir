# Setup

Local dev + free-tier signup checklist for Trinetravir. Tracks human-only steps the agent cannot do for you.

## 1. Local environment

Already scaffolded. To bring up a fresh checkout:

```bash
# Python 3.11 must be on PATH (pinned via .python-version)
uv sync --extra dev      # installs runtime + dev deps into .venv/
uv run pre-commit install
```

Sanity check:

```bash
uv run python -c "import scanpy, anndata, torch, scvi; print('OK', scanpy.__version__)"
```

## 2. Free-tier accounts (Phase 1, PLAN §9 — sign up but DO NOT activate paid timers yet)

| Service | Purpose | Status | Activation note |
|---|---|---|---|
| **Google Cloud** ($300, 90d) | Phase 8/9 GPU sweeps | [ ] not signed up | **DELAY activation** until Week 6+ — 90-day clock starts on signup |
| **Microsoft Azure** ($200) | Backup GPU pool | [ ] not signed up | Less time-pressure; can sign up early |
| **Oracle Cloud** (Always Free, 4 OCPU/24GB ARM) | Persistent CPU box for harmonization / overnight embedding jobs | [ ] not signed up | Sign up now; truly free indefinitely |
| **RunPod** ($5-10 starter + referral) | Cheap community GPU | [ ] not signed up | Sign up + verify; redeem credits when needed |
| **Kaggle** (30 hr/wk P100 16GB) | Foundation model embedding extraction | [ ] not signed up | First place GPU is genuinely needed (Phase 6) |
| **Google Colab** (free T4) | Prototyping | [ ] not signed up | Optional |
| **Weights & Biases** | Experiment tracking | [ ] not signed up | Required before Phase 5+. Free for personal/academic |
| **Hugging Face** | Foundation model checkpoints + dataset hosting | [ ] not signed up | Free; needed for scGPT/Geneformer/scPRINT-2 |
| **Zenodo** | Final dataset DOI release | [ ] not signed up | Phase 10 only |
| **bioRxiv** | Preprint submission | [ ] not signed up | Phase 10 only |

After creating the wandb account:

```bash
uv run wandb login
```

## 3. Compute escalation triggers (PLAN §3.5)

Stay on laptop until a specific job blocks iteration. Escalation order:
1. Laptop (Tier 0)
2. Oracle Always Free CPU (Tier 1) — for long jobs you don't want tying up the laptop
3. Kaggle weekly P100 (Tier 3) — for foundation model embedding extraction
4. Vast.ai community spot RTX 4090 / A100 (Tier 4) — for Phase 8/9 sweeps
5. GCP $300 credits (Tier 2) — fallback when stability matters more than price

## 4. Curated reference data to download (Phase 2 prep, NOT yet)

These go into `data/reference/` (git-tracked, small):
- Interferome 2.0 ISG list — `interferome.its.monash.edu.au`
- Mostafavi et al. 2016 ISG curation
- Viral entry receptor gene sets per virus (ACE2/TMPRSS2, sialic acid + TMPRSS family, CX3CR1/nucleolin, HVEM/nectins)
- REACTOME R-HSA-913531 (Type I IFN signaling), KEGG hsa04060
