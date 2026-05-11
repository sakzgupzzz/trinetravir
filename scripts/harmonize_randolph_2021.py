"""Session 6A Part C harmonization for Randolph 2021 (GSE162632 IAV ex vivo).

Multiplexed cohort: 30 10x pools, each containing multiple donors with hashtag
oligo (HTO) demultiplexing. Per-cell donor identity is in the published Seurat
.rds files (Zenodo 4273999 inputs.tar.gz, downloaded + selectively extracted).

Inputs:
  data/raw/randolph_2021/GSM*_*_barcodes.tsv.gz, *_features.tsv.gz, *_matrix.mtx.gz  (GEO raw 10x)
  data/raw/randolph_2021/zenodo_inputs/inputs/1_calculate_pseudobulk/*_cluster_singlets.rds  (per-cell-type Seurat singlets with demuxed donor IDs)

Output:
  data/processed/randolph_2021_processed_v6.h5ad

Pipeline:
  1. Read per-cell-type Seurat singlets .rds via rdata library; extract per-cell
     metadata (barcode, donor_id from MS_indiv_ID, infection_status, celltype).
  2. Consolidate barcode -> donor_id + condition + celltype mapping.
  3. Load GEO raw 10x mtx files per pool (B1_c1 etc), build per-pool AnnData.
  4. Filter to demuxed singlet barcodes only.
  5. Map cells to v1's 5-bucket scheme via Randolph celltype labels.
  6. Apply schema v6 (paired_within_donor design, ex_vivo_challenge, 6h duration).
  7. Save processed h5ad.

Note: paired_within_donor design — each donor has both mock (NI) + IAV-infected
cells in same Randolph experiment. exposure_pair_id = donor_id so paired test
can pair mock vs IAV per donor.
"""

from __future__ import annotations

import gzip
import logging
import re
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import rdata
import scipy.io as sio

from trinetravir.data.schema_v6_migration import V6_CATEGORIES

warnings.filterwarnings("ignore", category=UserWarning, module="rdata")

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
GEO_RAW = REPO / "data" / "raw" / "randolph_2021"
ZENODO_INPUTS = (
    REPO / "data" / "raw" / "randolph_2021" / "zenodo_inputs" / "inputs" / "1_calculate_pseudobulk"
)
OUT = REPO / "data" / "processed" / "randolph_2021_processed_v6.h5ad"

# Randolph cell types -> v1 5-bucket
RANDOLPH_CELLTYPE_TO_BUCKET = {
    "monocytes": "monocyte",
    "infected_monocytes": "monocyte",
    "CD4_T": "CD4T",
    "CD8_T": "CD8T",
    "B": "B",
    "NK": "NK",
    "NK_high_response": "NK",
    "DC": "other",
    "neutrophils": "other",
    "NKT": "other",
}


def load_rds_meta(rds_path: Path) -> pd.DataFrame:
    """Extract Seurat object's meta.data via rdata. Returns pandas DataFrame."""
    logger.info("reading rds: %s (%.1f GB)", rds_path.name, rds_path.stat().st_size / 1e9)
    parsed = rdata.parser.parse_file(str(rds_path))
    conv = rdata.conversion.convert(parsed)
    md = getattr(conv, "meta.data")
    # Coerce column names from np.str_ to str
    md.columns = [str(c) for c in md.columns]
    return md


def build_barcode_to_donor_map() -> pd.DataFrame:
    """Concatenate per-cell-type demux CSVs (pre-extracted via extract_randolph_demux.py).

    CD4T + CD8T cluster_singlets .rds files (10GB + 2.3GB compressed) cause OOM
    during rdata parsing on a 16GB laptop. Randolph harmonization in Session 6A
    uses only monocyte + B + NK buckets, sufficient for the Issue 27 primary
    biological test (monocyte cross-context conserved-component). CD4T + CD8T
    are deferred to v1.5 when streaming-Seurat-parse is available.

    Returns DataFrame indexed by cell barcode with columns:
    donor_id, infection_status, celltype, bucket, pool_capture.
    """
    csv_files = {
        "monocyte": GEO_RAW / "demux_meta_monocyte.csv",
        "B": GEO_RAW / "demux_meta_B.csv",
        "NK": GEO_RAW / "demux_meta_NK.csv",
    }
    parts = []
    for bucket, path in csv_files.items():
        if not path.exists():
            logger.warning("missing %s", path)
            continue
        df = pd.read_csv(path)
        donor_id = df["MS_indiv_ID"].astype(str)
        cond = df["MS_infection_status"].astype(str)
        pool = df["orig.ident"].astype(str)
        keep_mask = (donor_id != "nan") & (donor_id != "") & (cond != "nan")
        sub = pd.DataFrame(
            {
                "donor_id": donor_id.values,
                "infection_status": cond.values,
                "celltype": df.get("celltype", pd.Series([bucket] * len(df))).astype(str).values,
                "bucket": bucket,
                "pool_capture": pool.values,
            },
            index=df["barcode"].astype(str).values,
        )
        sub = sub.loc[keep_mask.values]
        parts.append(sub)
        logger.info("  %s: %d cells, %d unique donors", bucket, len(sub), sub["donor_id"].nunique())
    if not parts:
        raise RuntimeError("no demux CSVs found; run extract_randolph_demux.py first")
    return pd.concat(parts, axis=0)


def load_geo_pool(pool_name: str) -> ad.AnnData | None:
    """Load one GEO 10x pool from data/raw/randolph_2021/GSM*_<pool>_*.

    pool_name from demux CSV is e.g. "B1_c1"; GEO filenames use dashes:
    "GSM4955739_B1-c1-10X_barcodes.tsv.gz". Convert underscore->dash for match.
    """
    geo_pool = pool_name.replace("_", "-")
    pat = re.compile(rf"^GSM\d+_{re.escape(geo_pool)}-10X_(barcodes|features|matrix)\.")
    files = {}
    for p in GEO_RAW.iterdir():
        m = pat.match(p.name)
        if m:
            files[m.group(1)] = p
    if not all(k in files for k in ("barcodes", "features", "matrix")):
        return None
    with gzip.open(files["matrix"], "rb") as f:
        X = sio.mmread(f).tocsr().T.tocsr()  # genes x cells -> cells x genes
    with gzip.open(files["barcodes"], "rt") as f:
        barcodes = [line.strip() for line in f]
    with gzip.open(files["features"], "rt") as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    n_features = X.shape[1]
    if rows and len(rows[0]) >= 2:
        gene_id = [r[0] for r in rows[:n_features]]
        gene_sym = [r[1] for r in rows[:n_features]]
    else:
        gene_id = [r[0] for r in rows[:n_features]]
        gene_sym = gene_id
    var = pd.DataFrame({"gene_symbol": gene_sym}, index=gene_id)
    dup = var.index.duplicated(keep="first")
    if dup.any():
        var = var.loc[~dup]
        X = X[:, ~dup]
    obs = pd.DataFrame(index=[f"{pool_name}_{b.split('-')[0]}" for b in barcodes])
    obs["pool"] = pool_name
    return ad.AnnData(X=X, obs=obs, var=var)


def main() -> int:
    # ---- Step 1: build barcode -> donor mapping from Seurat metadata ----
    logger.info("extracting demux metadata from Seurat .rds files")
    demux = build_barcode_to_donor_map()
    logger.info("total demuxed cells across cell types: %d", len(demux))
    logger.info("unique donors: %d", demux["donor_id"].nunique())
    logger.info("condition counts: %s", dict(demux["infection_status"].value_counts()))

    # ---- Step 2: load each GEO pool + match to demux ----
    pools = sorted(demux["pool_capture"].astype(str).unique().tolist())
    logger.info("loading %d GEO pools", len(pools))
    parts = []
    for pool in pools:
        a = load_geo_pool(pool)
        if a is None:
            logger.warning("missing GEO data for pool %s", pool)
            continue
        # match demux subset for this pool
        sub_demux = demux[demux["pool_capture"] == pool]
        common_idx = a.obs.index.intersection(sub_demux.index)
        if len(common_idx) == 0:
            logger.warning("pool %s: no cells matched between GEO + demux", pool)
            continue
        a = a[common_idx].copy()
        a.obs["donor_id"] = sub_demux.loc[common_idx, "donor_id"].values
        a.obs["infection_status"] = sub_demux.loc[common_idx, "infection_status"].values
        a.obs["randolph_celltype"] = sub_demux.loc[common_idx, "celltype"].values
        a.obs["bucket"] = sub_demux.loc[common_idx, "bucket"].values
        parts.append(a)
        logger.info("pool %s: %d demuxed cells", pool, a.n_obs)
    if not parts:
        raise RuntimeError("no pools loaded")
    logger.info("concatenating %d pools", len(parts))
    combined = ad.concat(parts, join="inner", merge="first", uns_merge="first")
    logger.info("combined: %d cells x %d genes", combined.n_obs, combined.n_vars)

    # ---- Step 3: map to 5-bucket scheme via Randolph celltype ----
    bk = combined.obs["randolph_celltype"].astype(str).map(RANDOLPH_CELLTYPE_TO_BUCKET)
    bk = bk.fillna("other")
    combined.obs["cell_type_bucket_unified"] = pd.Categorical(
        bk.values, categories=["monocyte", "CD4T", "CD8T", "B", "NK", "other"]
    )
    other_pct = 100.0 * int((bk == "other").sum()) / len(bk) if len(bk) else 0.0
    logger.info(
        "bucket counts: %s | 'other' pct: %.2f%%",
        dict(combined.obs["cell_type_bucket_unified"].value_counts()),
        other_pct,
    )

    # ---- Step 4: schema v6 ----
    n = combined.n_obs
    # infection_status: NI = mock (healthy_control), flu = diseased (IAV-infected)
    cond_str = combined.obs["infection_status"].astype(str)
    combined.obs["donor_disease_status"] = pd.Categorical(
        np.where(cond_str == "flu", "diseased", "healthy_control"),
        categories=["diseased", "healthy_control"],
    )

    combined.obs["study_id"] = "randolph_2021"
    combined.obs["study_id"] = combined.obs["study_id"].astype("category")
    combined.obs["donor_id"] = combined.obs["donor_id"].astype("category")

    combined.obs["donor_response_design"] = pd.Categorical(
        ["paired_within_donor"] * n, categories=V6_CATEGORIES["donor_response_design"]
    )
    # paired ID = donor_id (each donor contributes both NI + flu samples)
    combined.obs["exposure_pair_id"] = pd.array(
        combined.obs["donor_id"].astype(str).values, dtype="string"
    )
    combined.obs["exposure_type"] = pd.Categorical(
        ["ex_vivo_challenge"] * n, categories=V6_CATEGORIES["exposure_type"]
    )
    combined.obs["exposure_duration_hours"] = np.full(n, 6.0, dtype=np.float64)
    combined.obs["age_years"] = np.full(n, np.nan, dtype=np.float64)
    combined.obs["age_group_category"] = pd.Categorical(
        ["adult"] * n, categories=V6_CATEGORIES["age_group_category"]
    )
    combined.obs["infection_state"] = pd.Categorical(
        np.where(cond_str == "flu", "acute", "naive"),
        categories=V6_CATEGORIES["infection_state"],
    )
    combined.obs["donor_serostatus"] = pd.Categorical(
        ["unknown"] * n, categories=V6_CATEGORIES["donor_serostatus"]
    )

    combined.uns["annotation_source"] = "randolph_2021_seurat_singlets_celltype"
    combined.uns["bucket_column"] = "cell_type_bucket_unified"
    combined.uns["demux_source"] = "Zenodo_4273999_inputs.tar.gz_MS_indiv_ID"

    # Issue 4 verification per condition
    donor_per_cond = combined.obs.groupby("donor_disease_status", observed=True)[
        "donor_id"
    ].nunique()
    logger.info("donor counts per condition: %s", donor_per_cond.to_dict())

    OUT.parent.mkdir(parents=True, exist_ok=True)
    logger.info("writing %s", OUT)
    combined.write_h5ad(OUT)
    logger.info("done. shape: %s", combined.shape)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
