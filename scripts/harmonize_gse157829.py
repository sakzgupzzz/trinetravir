"""Session 6A Part C harmonization for GSE157829 (Wang 2020 HIV exhaustion atlas).

Inputs:
  data/raw/gse157829/GSM*_*barcodes.tsv.gz
  data/raw/gse157829/GSM*_*genes.tsv.gz
  data/raw/gse157829/GSM*_*matrix.mtx.gz

Outputs:
  data/processed/gse157829_processed_v6.h5ad

Pipeline (per SESSION_6A_CHECKLIST C-pre.6):
  1. Load 7 samples (C1 healthy + Q1, Q2, Q3, Q4, Q5, Q7 HIV) into per-sample
     AnnData; concatenate.
  2. Apply CellTypist Immune_All_Low (per Issue 12) on combined matrix.
  3. Map to v1 5-bucket scheme.
  4. donor_disease_status: C1 -> healthy_control; Q* -> diseased.
  5. Apply schema_v6_migration:
     - donor_response_design = cross_sectional
     - exposure_type = retroviral_infection
     - infection_state = chronic_latent (HIV+) | naive (C1 healthy)
     - donor_serostatus = positive (HIV+) | negative (C1)
     - age_group_category = adult (per Issue 30; chronic HIV adult cohort)
  6. Save data/processed/gse157829_processed_v6.h5ad.

Note (Issue 30 amendment 9f6b79e): C1 is a single healthy donor in GEO
deposit (paper's '4 healthy' includes 3 EXTERNAL public 10X datasets not
in GEO). Cross-cohort integration design used: v1 corpus's 41 aggregated
healthy donors serve as primary baseline for Session 6B calibration.
Within-GSE157829 C1 = supplementary sanity check.
"""

from __future__ import annotations

import gzip
import logging
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio

from trinetravir.data.annotate import annotate_unified
from trinetravir.data.schema_v6_migration import V6_CATEGORIES

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw" / "gse157829"
OUT = REPO / "data" / "processed" / "gse157829_processed_v6.h5ad"

# Sample -> donor mapping. C1 = healthy. Q* = HIV donors.
# Per GSE157829 GEO + PMC7646563.
SAMPLE_TO_DONOR_DISEASE = {
    "C1": ("C1", "healthy_control"),
    "Q1": ("Q1", "diseased"),
    "Q2": ("Q2", "diseased"),
    "Q3": ("Q3", "diseased"),
    "Q4": ("Q4", "diseased"),
    "Q5": ("Q5", "diseased"),
    "Q7": ("Q7", "diseased"),
}


def load_sample(barcodes_p: Path, features_p: Path, matrix_p: Path) -> ad.AnnData:
    """Load one 10x sample (barcodes + features + matrix gz). Use ENSG IDs as var index."""
    with gzip.open(matrix_p, "rb") as f:
        X = sio.mmread(f).tocsr().T.tocsr()  # genes-by-cells -> cells-by-genes
    with gzip.open(barcodes_p, "rt") as f:
        barcodes = [line.strip() for line in f]
    with gzip.open(features_p, "rt") as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    n_features = X.shape[1]
    if rows and len(rows[0]) >= 2:
        gene_id_raw = [r[0] for r in rows[:n_features]]
        gene_sym_raw = [r[1] for r in rows[:n_features]]
    else:
        gene_id_raw = [r[0] for r in rows[:n_features]]
        gene_sym_raw = gene_id_raw

    # Some samples (e.g. GSM4775594_Q7) prefix gene IDs with "hg19_" because
    # cellranger was run against a custom hg19 + viral reference. Strip the
    # hg19_ prefix and drop any non-human entries (those don't start with hg19_
    # but also aren't bare ENSG — drop). After stripping, the gene IDs align
    # with the other samples' bare ENSG IDs.
    keep_idx = []
    gene_id = []
    gene_sym = []
    for i, (gid, gs) in enumerate(zip(gene_id_raw, gene_sym_raw, strict=False)):
        if gid.startswith("hg19_"):
            keep_idx.append(i)
            gene_id.append(gid[len("hg19_") :])
            gene_sym.append(gs[len("hg19_") :] if gs.startswith("hg19_") else gs)
        elif gid.startswith("ENSG"):
            keep_idx.append(i)
            gene_id.append(gid)
            gene_sym.append(gs)
        # else: viral or non-human contig (e.g. "hiv_*"); drop

    if len(keep_idx) != X.shape[1]:
        X = X[:, np.asarray(keep_idx)]
    var = pd.DataFrame({"gene_symbol": gene_sym}, index=gene_id)

    # Dedup index
    dup_mask = var.index.duplicated(keep="first")
    if dup_mask.any():
        keep = ~dup_mask
        var = var.loc[keep]
        X = X[:, keep]
    obs = pd.DataFrame(index=barcodes)
    return ad.AnnData(X=X, obs=obs, var=var)


def main() -> int:
    # Group files by sample
    file_groups: dict[str, dict[str, Path]] = {}
    for p in RAW.iterdir():
        m = re.match(r"^GSM\d+_(C1|Q\d)(barcodes|genes|features|matrix)\.(tsv|mtx)\.gz$", p.name)
        if not m:
            continue
        sample, kind = m.group(1), m.group(2)
        # normalize kind: gene/features -> features
        kind_norm = "features" if kind in ("genes", "features") else kind
        file_groups.setdefault(sample, {})[kind_norm] = p

    parts: list[ad.AnnData] = []
    for sample, files in sorted(file_groups.items()):
        b = files.get("barcodes")
        f = files.get("features")
        m = files.get("matrix")
        if not (b and f and m):
            logger.warning("sample %s missing files: have %s", sample, list(files.keys()))
            continue
        logger.info("loading sample %s", sample)
        a = load_sample(b, f, m)
        a.obs["sample_id"] = sample
        donor, disease = SAMPLE_TO_DONOR_DISEASE[sample]
        a.obs["donor_id"] = donor
        a.obs["donor_disease_status"] = disease
        # Use first cellxgene-style label as cell_type placeholder for annotate_unified
        a.obs["cell_type"] = "unknown"
        a.obs_names = [f"{sample}_{x}" for x in a.obs_names]
        parts.append(a)

    if not parts:
        raise RuntimeError("no GSE157829 samples loaded")

    logger.info("concatenating %d samples", len(parts))
    combined = ad.concat(parts, join="inner", merge="first", uns_merge="first")
    logger.info("combined shape: %s", combined.shape)

    # Make obs_names unique + restore proper categorical dtypes
    combined.obs["sample_id"] = combined.obs["sample_id"].astype("category")
    combined.obs["donor_id"] = combined.obs["donor_id"].astype("category")
    combined.obs["donor_disease_status"] = pd.Categorical(
        combined.obs["donor_disease_status"].astype(str).values,
        categories=["diseased", "healthy_control"],
    )

    # ---- CellTypist Immune_All_Low (per Issue 12) ----
    # CellTypist requires gene symbols as var_names. Swap from ENSG -> gene_symbol.
    if "gene_symbol" in combined.var.columns:
        gs = combined.var["gene_symbol"].astype(str).values
        # Dedup gene_symbol (some symbols map to multiple ENSG)
        new_idx = pd.Index(gs)
        dup = new_idx.duplicated(keep="first")
        if dup.any():
            keep = ~dup
            combined = combined[:, keep].copy()
            gs = gs[keep]
        combined.var_names = pd.Index(gs)
        combined.var.index.name = "gene_symbol"
    logger.info("running CellTypist Immune_All_Low")
    combined = annotate_unified(combined, model_name="Immune_All_Low.pkl", majority_voting=True)
    bucket_counts = combined.obs["cell_type_bucket_unified"].value_counts()
    other_count = int(bucket_counts.get("other", 0))
    total = int(combined.n_obs)
    other_pct = 100.0 * other_count / total
    logger.info("bucket counts: %s | 'other' pct: %.2f%%", dict(bucket_counts), other_pct)

    # ---- Schema v6 obs ----
    n = combined.n_obs
    combined.obs["study_id"] = "gse157829"
    combined.obs["study_id"] = combined.obs["study_id"].astype("category")

    combined.obs["donor_response_design"] = pd.Categorical(
        ["cross_sectional"] * n, categories=V6_CATEGORIES["donor_response_design"]
    )
    combined.obs["exposure_pair_id"] = pd.array([""] * n, dtype="string")
    combined.obs["exposure_type"] = pd.Categorical(
        ["retroviral_infection"] * n, categories=V6_CATEGORIES["exposure_type"]
    )
    combined.obs["exposure_duration_hours"] = np.full(n, np.nan, dtype=np.float64)
    combined.obs["age_years"] = np.full(n, np.nan, dtype=np.float64)
    combined.obs["age_group_category"] = pd.Categorical(
        ["adult"] * n, categories=V6_CATEGORIES["age_group_category"]
    )

    disease_str = combined.obs["donor_disease_status"].astype(str)
    combined.obs["infection_state"] = pd.Categorical(
        np.where(disease_str == "diseased", "chronic_latent", "naive"),
        categories=V6_CATEGORIES["infection_state"],
    )
    combined.obs["donor_serostatus"] = pd.Categorical(
        np.where(disease_str == "diseased", "positive", "negative"),
        categories=V6_CATEGORIES["donor_serostatus"],
    )

    combined.uns["annotation_source"] = "celltypist_immune_all_low"
    combined.uns["celltypist_model"] = "Immune_All_Low.pkl"
    combined.uns["celltypist_majority_voting"] = True

    OUT.parent.mkdir(parents=True, exist_ok=True)
    logger.info("writing %s", OUT)
    combined.write_h5ad(OUT)
    logger.info(
        "done. processed shape %s. n_donors: %d HIV + %d healthy",
        combined.shape,
        sum(1 for v in SAMPLE_TO_DONOR_DISEASE.values() if v[1] == "diseased"),
        sum(1 for v in SAMPLE_TO_DONOR_DISEASE.values() if v[1] == "healthy_control"),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
