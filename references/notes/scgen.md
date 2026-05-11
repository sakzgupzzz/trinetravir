# scGen

**Paper**: scGen predicts single-cell perturbation responses, Lotfollahi M., Wolf F.A., Theis F.J., 2019, *Nature Methods*
**Repo**: https://github.com/theislab/scgen
**License**: GPL-3.0 per `LICENSE` file at repo root. Note: `pyproject.toml` lists `license = "MIT"` — internal inconsistency; treat as GPL-3.0 (the more restrictive of the two) unless maintainers clarify.

## One-paragraph summary
scGen is a generative model that predicts how cells respond to a perturbation in conditions where you have not directly observed that response. It learns a low-dimensional latent space from control + perturbed cells across multiple cell types using a vanilla VAE, then performs *vector arithmetic* in latent space — `latent(predicted_perturbed) = latent(observed_control) + (mean_latent(perturbed) − mean_latent(control))` — and decodes back to gene-expression space. The promise is that the perturbation-effect vector transfers across cell types you only saw in the control condition, enabling predictions for held-out cell types and species. It is implemented on top of `scvi-tools` and is the canonical "simple but strong" perturbation-prediction baseline that Trinetravir must beat.

## Architecture
A single VAE (no condition-specific branches, no adversarial heads). Inputs are normalized + log1p-transformed gene-expression vectors; default config is HVG-filtered to ~7000 genes per the README.

- Encoder `z_encoder`: scvi-tools `Encoder` with `n_layers=2`, `n_hidden=800`, `LeakyReLU` activation, batch-norm enabled, dropout `0.2`. Outputs `(qz_m, qz_v)` for a Gaussian variational posterior.
- Latent: `n_latent=100` by default (PLAN.md models.yaml uses this value).
- Decoder `DecoderSCGEN`: mirror of encoder, `n_hidden=800`, `n_layers=2`, `LeakyReLU`. Reconstructs gene-expression vector.
- Input shape: `(n_cells, n_genes)`. Latent shape: `(n_cells, 100)`. Output shape: `(n_cells, n_genes)`.
- Prediction (post-training): `predict()` does latent-space arithmetic across cell types, with optional `balancer()` that subsamples cells to equalize cell-type counts before computing the perturbation-direction vector.

## Loss function
Standard ELBO with a heavily weighted reconstruction term. From `scgen/_scgenvae.py:loss`:

`loss = 0.5 * reconstruction_loss + 0.5 * (kl_divergence * kl_weight)`

where `kl_weight = 5e-5` by default — i.e., the KL term is downweighted by ~10000 relative to reconstruction. `reconstruction_loss` is MSE on log1p-normalized data (not raw counts; not negative binomial). KL is computed against `Normal(0, 1)` prior. No adversarial loss, no auxiliary classifier, no perturbation-conditional branch.

## Hyperparameters (from paper or defaults in repo)
- learning rate: scvi-tools default (1e-3, AdamW)
- batch size: scvi-tools default (128)
- epochs: 100 in the paper experiments; `train()` accepts `max_epochs` argument
- latent dimensions: 100 (configurable via `n_latent`)
- hidden dim: 800 (configurable via `n_hidden`)
- n_layers: 2
- dropout: 0.2 (encoder/decoder)
- KL weight: 5e-5 (heavily downweighted; reconstruction-dominated)
- Recommended preprocessing: `sc.pp.normalize_total` + `sc.pp.log1p` + top ~7000 HVGs
- other: balancer subsamples cell-type-specific cells before computing the perturbation direction (rebalances cell-type composition between control and stimulated populations to avoid composition bias in the direction vector)

## Evaluation protocol used in paper
- Datasets: PBMC IFN-β stimulation (Kang et al. 2018), HPoly + Salmonella intestinal infection (Haber et al. 2017), influenza species-transfer (Plasschaert et al. 2018, mouse → human airway).
- Setup: train on (control + stimulated) for some cell types and (control only) for held-out cell types; predict the held-out cell type's stimulated profile.
- Metrics: per-gene Pearson r (predicted vs true post-perturbation mean), per-DE-gene Pearson r (top 100), R² on means + variances, qualitative UMAPs.
- Comparisons: CVAE, MMD-CVAE, scVI baselines, plus a no-prediction "control mean" floor.

## Relevance to Trinetravir
scGen is one of the three Phase 6 "existing methods" we will benchmark per PLAN.md §4.6 (alongside scCausalVI and the foundation-model embedding heads). It is the strongest pure-generative baseline that does not rely on causal disentanglement or large pretraining; if our factorized model cannot beat scGen on cross-virus transfer, the novel architecture has not earned its keep. **Known limitation for cross-virus generalization:** scGen's latent-arithmetic assumes a single global perturbation direction transfers across cell types and species. For cross-virus PBMC, this implies the SARS-CoV-2 → IAV mapping is a single latent shift, ignoring our hypothesis that responses decompose into a shared antiviral component plus a virus-specific component (PLAN H4). scGen will likely interpolate between the two virus directions rather than additively decompose them — that is exactly the failure mode we expect, and the comparison should make that visible.

## Implementation notes
- Built on `scvi-tools` (`scvi.model.base.BaseModelClass` plus `VAEMixin` + `UnsupervisedTrainingMixin`). Inherits scvi's training loop, AnnData manager, and lightning-based fit. Familiar API if you already use scvi.
- `setup_anndata()` requires categorical `condition_key` (e.g. control / stimulated) and `cell_type_key` registered before `train()`. Easy to overlook — `train()` will fail with an opaque error otherwise.
- The `predict()` API is awkward: takes `ctrl_key` and `stim_key` strings plus either `celltype_to_predict` OR `adata_to_predict` (XOR — both raises). The choice between modes is not obvious from the docstring.
- Pinned dep: `scvi-tools >= 0.20.0`. Our env has `scvi-tools 1.4.2` which is far ahead; the older API has churned (e.g. `BaseModelClass`, `REGISTRY_KEYS`, `LossOutput` signatures). May need a wrapper that pins a compatible scvi-tools version, or vendor a thin re-implementation if upstream is too out of date for our environment.
- License inconsistency: `LICENSE` file at root is GPL-3.0; `pyproject.toml` claims MIT. We will treat as GPL-3.0 for redistribution purposes — flag for our own LICENSE-compliance check before any code from this repo is vendored into Trinetravir.
- Repo last meaningfully updated 2023; the Travis CI badge in README is dead. Maintenance status is "stable but quiet."
- The `balancer()` utility (`scgen/_utils.py`) is the unsung hero — composition-bias correction in the perturbation direction is genuinely important for our cross-study PBMC setting where cell-type proportions vary substantially (notebook 03 confirmed this).
