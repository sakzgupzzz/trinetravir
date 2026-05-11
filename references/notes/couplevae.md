# CoupleVAE

**Paper**: CoupleVAE: coupled variational autoencoders for predicting perturbational single-cell RNA sequencing data, Wu Y., Liu J., Xiao Y., Zhang S., Li L., 2025, *Briefings in Bioinformatics* 26(2):bbaf126
**Repo**: https://github.com/LiminLi-xjtu/CoupleVAE
**License**: NO LICENSE file in the repo as of clone date. Default = "all rights reserved" under copyright; **cannot be redistributed or vendored without explicit permission from the authors**. We may run it as a baseline locally; we cannot ship its code as part of Trinetravir.

## One-paragraph summary
CoupleVAE is a perturbation-prediction model that pairs two parallel VAEs — one for the control (unperturbed) condition and one for the perturbed condition — and ties their latent spaces together with two learned "coupler" networks that map control-latent ↔ perturbed-latent. Inference path: encode a control cell → couple to the perturbed-latent space → decode through the perturbed decoder → synthetic post-perturbation cell. The system is trained simultaneously on direct reconstruction (control→control, perturbed→perturbed) and translation reconstruction (control→perturbed→control roundtrip and vice versa) plus an alignment loss that pulls the coupler output toward the encoded distribution of the opposite condition. Useful as a third generative-baseline class for Trinetravir, complementary to scGen (latent arithmetic) and scCausalVI (structural causal).

## Architecture
A pair of mirrored VAEs plus two couplers, all simple MLPs (no scvi-tools, no normalization-flow, no negative-binomial decoder).

- `Encoderc` (control) and `Encoderp` (perturbed): identical architecture. Two hidden layers of 800 units each, BatchNorm + LayerNorm + LeakyReLU + Dropout. Output `(mean, log_var)` of latent dim `z_dim=16`.
- `Decoderc`, `Decoderp`, `Decodercp` (control→perturbed translation), `Decoderpc` (perturbed→control translation): four parallel decoders, each two hidden layers of 800 units with LeakyReLU + Dropout, final ReLU on output (so the model implicitly assumes nonnegative inputs — they recommend log1p+scale preprocessing).
- `Couplerc`, `Couplerp`: small MLPs that map a latent vector through one hidden layer (`z_dim → z_dim` with BatchNorm + ReLU) to produce `(mean, log_var)` of the *other* condition's latent. The coupler IS the cross-condition mapping.
- Forward pass given a paired batch `(x_0, x_1)` (control, perturbed):
  1. `mu_0, log_var_0 = encoder_c(x_0)`, sample `z_c`
  2. `mu_1, log_var_1 = encoder_p(x_1)`, sample `z_p`
  3. `mu_p_via_c, lv_p_via_c = coupler_c(z_c)` → sample `z_1` (control's prediction of where it would land in perturbed-latent)
  4. `mu_c_via_p, lv_c_via_p = coupler_p(z_p)` → sample `z_0` (perturbed's prediction of where it would land in control-latent)
  5. Reconstructions: `x_hat_0 = decoder_c(z_c)`, `x_hat_1 = decoder_p(z_p)`, `x_hat_cp = decoder_cp(z_1)` (control through coupler then through perturbed-decoder), `x_hat_pc = decoder_pc(z_0)`.
- Input shape: `(batch, n_genes)` for both conditions, paired (requires control and perturbed to be loaded in equal-size batches via the `Trainer` wrapper).
- Latent shape: `(batch, 16)`.
- Output shape: 4 reconstructions, each `(batch, n_genes)`.

## Loss function
All four loss terms operate per-cell, weighted by `0.25` constants in the source so they sum to ~1 across the four conditions. Combined loss per cell:

```
kl_loss     = 0.25 * KL(N(mu_0, log_var_0) || N(0,1)) + 0.25 * KL(N(mu_1, log_var_1) || N(0,1))
recon_loss  = 0.25 * ||x_0 - x_hat_0||^2 + 0.25 * ||x_1 - x_hat_1||^2          # direct reconstruction
trans_loss  = 0.25 * ||x_0 - x_hat_pc||^2 + 0.25 * ||x_1 - x_hat_cp||^2        # cross-translation reconstruction
coupl_loss  = 0.25 * ||mu_c - mu_0||^2 + 0.25 * ||mu_p - mu_1||^2              # coupler outputs match the other encoder's mu

vae_loss = mean( recon_loss + trans_loss + alpha * kl_loss + beta * coupl_loss )
```

Where `alpha = 0.01` (KL weight), `beta = 1.0` (coupler-alignment weight) by default. All reconstruction losses are MSE (Gaussian likelihood on log1p+scaled inputs); no negative binomial.

## Hyperparameters (from paper or defaults in repo)
- learning rate: 1e-3 (Adam)
- batch size: 32 (Trainer default)
- epochs: 25 (Trainer default; paper experiments use up to 100s)
- latent dim `z_dim`: 16
- hidden dim: 800 (hard-coded, not configurable via constructor)
- dropout: 0.2
- alpha (KL weight): 0.01
- beta (coupler weight): 1.0
- Recommended preprocessing: filter, normalize_per_cell, log1p, scale (full standardization)
- Early stopping: patience=20 in Trainer

## Evaluation protocol used in paper
- Datasets in the manuscript: PBMC IFN-β stimulation (Kang 2018), Salmonella + HPoly intestinal infection (Haber 2017), and a cross-species PBMC dataset (per the README, also COVID-19).
- Setup: paired control / perturbed cells from the same cell type used for training; evaluation predicts the perturbed profile of held-out cell types from their control profiles.
- Metrics: Pearson r on per-gene means, Pearson r on top DE genes, R² on gene means, qualitative UMAPs.
- Comparisons: scGen, CPA, scPreGAN, basic CVAE.

## Relevance to Trinetravir
CoupleVAE rounds out the generative-baseline trio for Phase 6 (PLAN.md §4.6) alongside scGen and scCausalVI. Its couplers learn a *bidirectional* control↔perturbed mapping, which is a different inductive bias than scGen's single linear shift or scCausalVI's per-treatment causal encoders. **Known limitation for cross-virus generalization:** the architecture is fundamentally pairwise (one control, one perturbed condition). For our 3-virus PBMC benchmark we would either need to train a separate CoupleVAE per (mock, virus_X) pair (no parameter sharing across viruses, no chance of zero-shot transfer to a held-out virus) or extend the architecture to multi-condition (which is paper-worthy but means re-implementing rather than wrapping). Worth running as a strict pairwise baseline (within-virus and cross-virus pairs separately) so reviewers cannot say we ignored it.

## Implementation notes
- Pure PyTorch, no scvi-tools dependency. The simplest of the three baselines to wrap.
- BUT: requires *paired* input batches (`x_0`, `x_1` both passed at the same time and assumed cell-aligned). PBMC data is unpaired between donors; the dataloader (`load_h5ad_to_dataloader` in `couplevae/model/util.py`) appears to randomly pair control and perturbed cells per batch. This is a known weakness of pairwise-VAE designs and partially defeats the "translation" semantics — there is no meaningful per-cell pairing in cross-donor scRNA-seq. Worth flagging in our writeup.
- Hard-coded `device='cuda'` in `to_latent()` and `predict()` — will fail on Mac MPS or CPU-only environments. Patch needed.
- Hidden dim (800), number of hidden layers (2 enc, 2 dec), and the architecture in general are not parameterizable through the constructor; would require source edits to ablate.
- `requirements.txt` pins very specific versions (anndata 0.10.7, scanpy 1.9.2, torch 2.3.1). Our env (anndata 0.12, scanpy 1.11, torch 2.11) is ahead; minor API changes likely. Defer dependency-resolution work until we actually wrap it.
- `Decoder.fc3` ends in ReLU, forcing nonnegative outputs. Together with the recommended `sc.pp.scale` preprocessing (which produces zero-mean unit-variance inputs that ARE often negative), this looks like a bug or at least a dataset-specific assumption; the model will systematically under-reconstruct negative values. Worth a sanity check before trusting any results.
- License gap: no `LICENSE` file in the repo. We must obtain explicit permission from the authors before vendoring or redistributing any of their code in Trinetravir. Running it locally as a baseline is fine.
- Maintenance: last meaningful commits early-2024 around publication; not actively maintained.
