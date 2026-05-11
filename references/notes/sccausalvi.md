# scCausalVI

**Paper**: scCausalVI disentangles single-cell perturbation responses with causality-aware generative model, Shaokun An et al., 2025, *Cell Systems* (preprint *bioRxiv* 10.1101/2025.02.02.636136)
**Repo**: https://github.com/ShaokunAn/scCausalVI
**License**: MIT (Copyright 2023 Shaokun An)

## One-paragraph summary
scCausalVI is a structural-causal generative model for case-control single-cell RNA-seq that explicitly factorizes each cell's representation into (a) a *background* latent capturing the cell's intrinsic state (cell type, donor, batch effects) and (b) a *treatment-effect* latent that is forced to be zero for control cells and nonzero only for cells from a treatment condition. A negative-binomial decoder reconstructs counts from the sum. Training pulls the background distributions of control and treated cells together via MMD so the only systematic difference is the treatment-effect latent. This gives an in-silico counterfactual: take a treated cell, swap its treatment-effect latent to zero, decode → "what this cell would have looked like without the treatment." The architecture is the closest existing analogue to Trinetravir's planned factorized model (shared antiviral component + virus-specific component) and is therefore the most important head-to-head baseline.

## Architecture
Multi-encoder, multi-decoder VAE on top of `scvi-tools`. Shapes for `n_genes` input genes:

- `control_background_encoder`: `(x, batch_index) → (qbg_m, qbg_v, z_bg)`. Produces background latent of dim `n_background_latent=10` for control cells. Uses `scvi.nn.Encoder` with `n_layers=1`, `n_hidden=128`, dropout `0.1`, batch covariate injected.
- One `treatment_background_encoder` PER treatment condition: same architecture as the control encoder; each produces a `z_bg` for cells of its specific treatment label.
- One `treatment_effect_encoder` PER treatment condition: produces treatment-effect latent of dim `n_te_latent=10`. Forced to zero for control cells in the inference step (`z_t = zeros`).
- Optional `library_encoder` if `use_observed_lib_size=False`; otherwise library size is computed directly from input counts.
- Decoder `scvi.nn.DecoderSCVI`: input is the concatenation `[z_bg, z_t]`, conditioned on batch. Produces `(px_rate, px_dropout)` for the ZINB likelihood. `px_r` (per-gene dispersion) is a global learnable parameter.
- Input shape `(n_cells, n_genes)`. Background latent `(n_cells, 10)`. Treatment latent `(n_cells, 10)`. Output: ZINB parameters over `n_genes`.

The structural-causal part is the *invariance constraint*: at training time, the model is told which cells are control vs treated, and the generative path requires `z_t = 0` for control cells. Treatment effects must therefore live entirely in the `z_t` channel, and `z_bg` becomes the "what this cell would have been without treatment" representation.

## Loss function
Per-cell contributions, all summed and meaned:

- `recon_loss` = `−ZINB(mu=px_rate, theta=px_r, zi_logits=px_dropout).log_prob(x).sum(-1)` — ZINB negative log-likelihood on raw counts (the model takes counts as input and applies `log(x+1)` internally before encoding).
- `kl_bg` = `KL( N(qbg_m, qbg_v) || N(0, 1) )` summed over latent dims, computed for both the control and the treatment paths.
- `kl_library` = library-size KL against per-batch log-normal prior, only when `use_observed_lib_size=False`; otherwise zero.
- `loss_mmd` = sum over treatment labels of `MMD(z_bg_control, z_bg_treatment_t)` using a multi-bandwidth Gaussian kernel — pulls the background distributions of control and each treatment condition together so that systematic background differences cannot leak into `z_t`.
- `loss_norm` = `norm_weight * ||z_t||_2^2` summed per cell — L2 sparsity on the treatment-effect latent for treated cells, pushing the model to use `z_t` only when it must.

`total_loss = recon_loss.mean() + kl_bg.mean() + kl_library.mean() + mmd_weight * loss_mmd + norm_weight * loss_norm.mean()`

## Hyperparameters (from paper or defaults in repo)
- learning rate: scvi-tools default (~1e-3, AdamW)
- batch size: scvi-tools default
- epochs: not pinned in module; tutorial uses ~400 epochs
- background latent dim: 10
- treatment-effect latent dim: 10
- hidden dim: 128
- n_layers: 1
- dropout: 0.1
- mmd_weight: 1.0 (default)
- norm_weight: 1.0 (default)
- gammas: multi-bandwidth Gaussian kernel (must be supplied when `use_mmd=True`)
- use_observed_lib_size: True (default, recommended)
- Recommended preprocessing: raw counts in `.X`; the model internally applies `log(x+1)` and uses observed library size.

## Evaluation protocol used in paper
- Datasets in the manuscript benchmarks: PBMC IFN-β stimulation (Kang 2018), HPoly + Salmonella intestinal (Haber 2017), LPS + IFNγ macrophage stimulation, an unpublished cross-condition cancer dataset.
- Tasks evaluated: counterfactual prediction (predict treated profile from control cell's `z_bg`), held-out cell-type generalization, identification of treatment-responsive subpopulations via clustering on `z_t`, and gene-level treatment-effect inference.
- Metrics: per-gene Pearson r between predicted and true post-treatment mean, R² on means + variances, AUROC for responsive-cell classification.
- Comparisons: scGen, CPA, CINEMA-OT, scVI baseline, plus an unsupervised differential-expression baseline.

## Relevance to Trinetravir
scCausalVI is the closest published analogue to our planned factorized model (PLAN.md §4.8). Both decompose the cell representation into a perturbation-invariant background and a perturbation-specific component, and both use a regularizer (MMD here, ISG/sparsity in our design) to keep the decomposition honest. Running scCausalVI on the cross-virus task is the head-to-head we most need to win — if scCausalVI's per-treatment encoders generalize to unseen treatments out of the box, our explicit ISG-mask regularization may not buy much. **Known limitation for cross-virus generalization:** scCausalVI's treatment-effect encoders are *per-treatment* — there is one set of weights per training condition. Predicting on a held-out virus requires either (a) treating the held-out virus as one of the training treatments at training time, which defeats the cross-virus setup, or (b) assuming the trained treatment-effect encoders' average applies, which has no causal justification. We need a wrapper that handles "novel treatment at inference" properly; the upstream code does not.

## Implementation notes
- Built on `scvi-tools` (`scvi.module.base.BaseModuleClass`, `scvi.distributions.ZeroInflatedNegativeBinomial`, `scvi.nn.Encoder`, `scvi.nn.DecoderSCVI`). Same scvi-tools version-skew risk as scGen — repo's `requirements.txt` should be checked against our `1.4.2`.
- The custom `SCCAUSALVI_REGISTRY_KEYS` (in `scCausalVI/model/base/`) extends scvi's `REGISTRY_KEYS` with a `CONDITION_KEY` for the treatment label. `setup_anndata()` in this repo registers that field and the standard batch field.
- `condition2int` is a hard-coded dict that the user must supply at `__init__`, mapping each condition string to an integer index. The control condition's name must also be passed explicitly. Easy footgun: if `condition2int` has the wrong control key, the entire causal invariance breaks silently (control cells get treatment-effect encoders).
- `use_mmd=True` (the paper's headline setting) requires `gammas` (kernel bandwidths) to be supplied — the constructor raises if they are missing. The repo does not document a sensible default; tutorials use a list like `[1, 2, 4, 8, 16]`.
- The `_generic_inference` loop over unique treatment labels per batch is `O(n_treatments)` and does Python-level masking — likely slow for our 4-virus + 1-mock setup but tractable.
- Last commit recent (2025); maintained by the original author.
- Re-implementation risk: about half a day to wrap for our pipeline, mostly because of the `condition2int` bookkeeping and the per-treatment encoder construction loop. Lower risk than scGen because the API uses current scvi-tools idioms.
