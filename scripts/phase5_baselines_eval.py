"""Phase 5 baselines evaluation harness (PLAN.md Phase 5 deliverable).

Per Issue 24 Category 1 + Issue 15 cross-virus protocol + Issue 14 hyperparameter
policy. Evaluates 3 trivial baselines (predict_mean, linear_delta, knn) on Lee
2020 within-study cross-virus task (only cross-virus data source in v1 corpus
per Issue 16).

Lee 2020 has 3 virus states: sars_cov_2 (31,463 cells), iav (10,519 cells),
mock (17,590 cells). Mock = baseline; sars + iav = response targets.

EVALUATION PROTOCOL (per Issue 15 + 14):
  Within-virus (leave-one-donor-out per virus):
    For each virus v in {sars_cov_2, iav}:
      For each Lee donor d that has both mock + v-infected cells in bucket b:
        Train on other donors (mock + v-infected cells, same bucket)
        Predict held-out donor d's v-infected response vector
        Compute Pearson r between predicted and observed

  Cross-virus (leave-one-virus-out):
    Hold out IAV donors' iav cells:
      Train on (all mock cells + sars_cov_2 cells from all donors)
      For each Lee donor with iav cells in bucket b:
        Predict iav response vector from donor's mock cells + virus_id='iav'
        Compute Pearson r
    Hold out SARS donors' sars cells:
      Symmetric.

OUTPUT:
  results/tables/phase5_baselines_eval.csv
    rows: (baseline, bucket, virus, split_type, n_donors, mean_r_pearson, sd_r_pearson)

Baselines operate at per-(donor, bucket, virus) response-vector level
(consistent with v1 calibration framework grain).
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from trinetravir.baselines.knn import KNNBaseline
from trinetravir.baselines.linear_delta import LinearDeltaBaseline
from trinetravir.baselines.predict_mean import PredictMeanBaseline

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"

LEE_FILE = PROC / "lee_2020_reannotated.h5ad"
HARMONY_GLOBAL = PROC / "harmony_global_embedding.h5ad"

BUCKETS = ("monocyte", "B", "NK", "CD4T", "CD8T")
VIRUSES = ("sars_cov_2", "iav")
MIN_CELLS_PER_DONOR_VIRUS = 50  # Issue 4 + Issue 24 floor


def load_lee_with_buckets() -> ad.AnnData:
    """Load Lee 2020 reannotated h5ad; map cells to coarse buckets via harmony_global.

    Returns: AnnData with cells in any of the 5 buckets + valid virus + donor.
    """
    logger.info("loading harmony_global for bucket assignment...")
    hg = ad.read_h5ad(HARMONY_GLOBAL, backed="r")
    hg_obs = hg.obs[["study_id", "coarse"]].copy()
    hg_lee = hg_obs[hg_obs["study_id"] == "lee_2020"].copy()
    # cell IDs in harmony_global are '<positional>-lee_2020'; strip suffix
    hg_lee["positional_idx"] = hg_lee.index.str.replace("-lee_2020", "", regex=False).astype(int)
    pos_to_coarse = dict(zip(hg_lee["positional_idx"], hg_lee["coarse"], strict=False))
    logger.info("harmony_global Lee cells: %d", len(hg_lee))

    logger.info("loading Lee reannotated h5ad...")
    lee = ad.read_h5ad(LEE_FILE)
    logger.info("Lee shape: %s", lee.shape)

    # Map cells to coarse bucket
    lee_coarse = pd.Series(
        [pos_to_coarse.get(i) for i in range(lee.shape[0])],
        index=lee.obs.index,
        dtype="object",
    )
    lee.obs["coarse"] = lee_coarse
    keep_mask = (lee.obs["coarse"].isin(BUCKETS) & lee.obs["virus"].isin(["mock", *VIRUSES])).values
    lee_kept = lee[keep_mask].copy()
    logger.info("Lee post-filter (in-bucket + valid virus): %d cells", lee_kept.shape[0])
    logger.info("  coarse counts: %s", dict(lee_kept.obs["coarse"].value_counts()))
    logger.info("  virus counts: %s", dict(lee_kept.obs["virus"].value_counts()))
    logger.info("  donor counts: %d unique", lee_kept.obs["donor_id"].nunique())
    return lee_kept


def compute_response_vectors(
    lee: ad.AnnData, hvg_genes: list[str]
) -> dict[tuple[str, str, str], np.ndarray]:
    """Compute per-(donor, bucket, virus_state) mean log-normalized expression.

    Returns dict keyed by (donor_id, bucket, virus) → mean expression on 4000 HVG.
    Virus state 'mock' = baseline; 'sars_cov_2' + 'iav' = post-infection.
    """
    logger.info("subsetting Lee to harmony_global HVG (4000 genes)...")
    var_sym = lee.var.index.astype(str)
    common = pd.Index(hvg_genes).intersection(var_sym)
    logger.info("  HVG intersect: %d / %d", len(common), len(hvg_genes))
    lee_hvg = lee[:, common].copy()

    logger.info("normalizing (target_sum=1e4, log1p)...")
    sc.pp.normalize_total(lee_hvg, target_sum=1e4)
    sc.pp.log1p(lee_hvg)

    X = lee_hvg.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)

    rvs = {}
    obs = lee_hvg.obs.copy()
    for (donor, bucket, virus), idx in obs.groupby(
        ["donor_id", "coarse", "virus"], observed=True
    ).groups.items():
        if len(idx) < MIN_CELLS_PER_DONOR_VIRUS:
            continue
        positional = obs.index.get_indexer(idx)
        rvs[(str(donor), str(bucket), str(virus))] = X[positional].mean(axis=0)
    logger.info("computed %d (donor, bucket, virus) response vectors", len(rvs))

    # Reindex to harmony_global HVG order (fill missing with 0)
    n_genes = len(hvg_genes)
    sym_to_idx = {g: i for i, g in enumerate(common)}
    aligned = {}
    for key, vec in rvs.items():
        full = np.zeros(n_genes, dtype=np.float32)
        for j, g in enumerate(hvg_genes):
            if g in sym_to_idx:
                full[j] = vec[sym_to_idx[g]]
        aligned[key] = full
    return aligned


def pearson_per_gene(pred: np.ndarray, true: np.ndarray) -> float:
    """Cell-wise Pearson r between predicted and observed expression vector."""
    if pred.std() < 1e-12 or true.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(pred, true)[0, 1])


def cohort_mock_baseline(
    rvs: dict[tuple[str, str, str], np.ndarray], bucket: str
) -> np.ndarray | None:
    """Mean of all mock donors' response vectors for given bucket.

    Lee donors are unpaired (each donor is either mock OR infected, not both).
    Cohort-level baseline = mean of all mock donor rvs for that bucket.
    """
    mock_rvs = [v for (d, b, virus), v in rvs.items() if b == bucket and virus == "mock"]
    if not mock_rvs:
        return None
    return np.stack(mock_rvs).mean(axis=0)


def within_virus_eval(
    rvs: dict[tuple[str, str, str], np.ndarray],
    baseline_cls,
    baseline_name: str,
    **bl_kwargs,
) -> list[dict]:
    """Leave-one-infected-donor-out within-virus evaluation.

    Donors are unpaired (separate mock vs infected donor groups). For each
    (bucket, virus): hold out one infected donor's response vector, train on
    other infected donors (+ cohort mock baseline), predict held-out donor's rv.
    """
    rows = []
    for bucket in BUCKETS:
        mock_baseline = cohort_mock_baseline(rvs, bucket)
        if mock_baseline is None:
            continue
        for virus in VIRUSES:
            infected_donors = sorted({d for (d, b, v) in rvs if b == bucket and v == virus})
            if len(infected_donors) < 3:
                continue
            r_values = []
            for held_donor in infected_donors:
                train_donors = [d for d in infected_donors if d != held_donor]
                X_train_baseline = np.stack([mock_baseline] * len(train_donors))
                X_train_post = np.stack([rvs[(d, bucket, virus)] for d in train_donors])
                train_cell_type = [bucket] * len(train_donors)
                train_virus = [virus] * len(train_donors)

                model = baseline_cls(**bl_kwargs)
                if baseline_name == "predict_mean":
                    model.fit(X_train_post, train_cell_type, train_virus)
                    pred = model.predict(mock_baseline.reshape(1, -1), [bucket], [virus])[0]
                elif baseline_name == "linear_delta":
                    X_train_response = X_train_post - X_train_baseline
                    model.fit(X_train_baseline, X_train_response, train_cell_type, train_virus)
                    pred = model.predict(mock_baseline.reshape(1, -1), [bucket], [virus])[0]
                elif baseline_name == "knn":
                    model.fit(X_train_baseline, X_train_post, train_virus)
                    pred = model.predict(mock_baseline.reshape(1, -1), [virus])[0]
                else:
                    continue

                true = rvs[(held_donor, bucket, virus)]
                r = pearson_per_gene(pred, true)
                if not np.isnan(r):
                    r_values.append(r)

            if not r_values:
                continue
            rows.append(
                {
                    "baseline": baseline_name,
                    "split_type": "within_virus",
                    "bucket": bucket,
                    "virus": virus,
                    "n_donors": len(infected_donors),
                    "n_evals": len(r_values),
                    "mean_r_pearson": round(float(np.mean(r_values)), 4),
                    "sd_r_pearson": round(float(np.std(r_values, ddof=1)), 4)
                    if len(r_values) > 1
                    else 0.0,
                    "min_r_pearson": round(float(np.min(r_values)), 4),
                    "max_r_pearson": round(float(np.max(r_values)), 4),
                }
            )
            logger.info(
                "  within %s %s: n=%d r_mean=%.4f (sd=%.4f)",
                bucket,
                virus,
                len(r_values),
                np.mean(r_values),
                np.std(r_values, ddof=1) if len(r_values) > 1 else 0.0,
            )
    return rows


def cross_virus_eval(
    rvs: dict[tuple[str, str, str], np.ndarray],
    baseline_cls,
    baseline_name: str,
    **bl_kwargs,
) -> list[dict]:
    """Leave-one-virus-out: train on N-1 viruses' donors, predict held-out virus donors.

    Per Issue 15 cross-virus protocol with cohort mock baseline.
    """
    rows = []
    for bucket in BUCKETS:
        mock_baseline = cohort_mock_baseline(rvs, bucket)
        if mock_baseline is None:
            continue
        for held_virus in VIRUSES:
            train_viruses = [v for v in VIRUSES if v != held_virus]
            train_donor_virus = [(d, v) for (d, b, v) in rvs if b == bucket and v in train_viruses]
            if len(train_donor_virus) < 2:
                continue

            X_train_baseline = np.stack([mock_baseline] * len(train_donor_virus))
            X_train_post = np.stack([rvs[(d, bucket, v)] for d, v in train_donor_virus])
            train_cell_type = [bucket] * len(train_donor_virus)
            train_virus = [v for _d, v in train_donor_virus]

            model = baseline_cls(**bl_kwargs)
            if baseline_name == "predict_mean":
                model.fit(X_train_post, train_cell_type, train_virus)
            elif baseline_name == "linear_delta":
                X_train_response = X_train_post - X_train_baseline
                model.fit(X_train_baseline, X_train_response, train_cell_type, train_virus)
            elif baseline_name == "knn":
                model.fit(X_train_baseline, X_train_post, train_virus)
            else:
                continue

            held_donors = sorted({d for (d, b, v) in rvs if b == bucket and v == held_virus})
            r_values = []
            for d in held_donors:
                if baseline_name == "knn":
                    pred = model.predict(mock_baseline.reshape(1, -1), [held_virus])[0]
                else:
                    pred = model.predict(mock_baseline.reshape(1, -1), [bucket], [held_virus])[0]
                true = rvs[(d, bucket, held_virus)]
                r = pearson_per_gene(pred, true)
                if not np.isnan(r):
                    r_values.append(r)
            if not r_values:
                continue
            rows.append(
                {
                    "baseline": baseline_name,
                    "split_type": "cross_virus",
                    "bucket": bucket,
                    "virus": held_virus,
                    "n_donors": len(held_donors),
                    "n_evals": len(r_values),
                    "mean_r_pearson": round(float(np.mean(r_values)), 4),
                    "sd_r_pearson": round(float(np.std(r_values, ddof=1)), 4)
                    if len(r_values) > 1
                    else 0.0,
                    "min_r_pearson": round(float(np.min(r_values)), 4),
                    "max_r_pearson": round(float(np.max(r_values)), 4),
                }
            )
            logger.info(
                "  cross %s held-out %s: n=%d r_mean=%.4f (sd=%.4f)",
                bucket,
                held_virus,
                len(r_values),
                np.mean(r_values),
                np.std(r_values, ddof=1) if len(r_values) > 1 else 0.0,
            )
    return rows


def main() -> int:
    # Load harmony_global HVG list (4000 genes; canonical input space per Issue 38 + Part B)
    logger.info("loading harmony_global HVG list...")
    hg = ad.read_h5ad(HARMONY_GLOBAL, backed="r")
    hvg_genes = list(hg.var.index)
    logger.info("HVG: %d genes", len(hvg_genes))

    lee = load_lee_with_buckets()
    rvs = compute_response_vectors(lee, hvg_genes)

    all_rows = []
    for cls, name, kwargs in [
        (PredictMeanBaseline, "predict_mean", {}),
        (LinearDeltaBaseline, "linear_delta", {"alpha": 1.0}),
        (KNNBaseline, "knn", {"k": 5, "within_virus_only": True}),
    ]:
        logger.info("=== %s within-virus ===", name)
        all_rows.extend(within_virus_eval(rvs, cls, name, **kwargs))
        logger.info("=== %s cross-virus ===", name)
        all_rows.extend(cross_virus_eval(rvs, cls, name, **kwargs))

    df = pd.DataFrame(all_rows)
    out = TABLES / "phase5_baselines_eval.csv"
    df.to_csv(out, index=False)
    logger.info("wrote %s (%d rows)", out.name, len(df))

    # Headline summary: mean r per (baseline, split_type) averaged across buckets/viruses
    logger.info("\n=== Phase 5 baseline headline ===")
    summary = (
        df.groupby(["baseline", "split_type"], observed=True)["mean_r_pearson"]
        .agg(["mean", "min", "max", "count"])
        .round(4)
    )
    logger.info("%s", summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
