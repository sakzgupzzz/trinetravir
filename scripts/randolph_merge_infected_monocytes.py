"""Merge Randolph infected_monocytes_cluster_singlets.rds into randolph_2021_processed_v6.h5ad.

The bystander 'monocyte' bucket already in v6 (15,531 cells from GEO mtx) is preserved.
Adds 4,964 cluster-8 infected-monocyte cells as new bucket 'monocyte_infected' so that
Issue 27 primary test can be re-run on the responsive subpopulation.

Pipeline:
  1. Parse rds via rdata (Seurat S4 object).
  2. Build sparse CSC counts matrix from dgCMatrix (i, p, x).
  3. Map rds gene symbols → v6 Ensembl IDs via v6.var.gene_symbol; pad zeros for symbols
     not present in v6's Ensembl space.
  4. Build per-cell obs aligned to v6 schema v6 columns; donor_id ← SOC_indiv_ID;
     donor_disease_status from SOC_infection_status (flu→diseased, NI→healthy_control).
     Copy donor-level v6 schema cols (exposure_*, age_*, infection_state, donor_serostatus)
     from existing v6 per donor_id; defaults for donors not in existing v6.
  5. Set bucket='monocyte_infected' for all new cells (NOT 'monocyte').
  6. Concat to existing v6 → randolph_2021_processed_v6.h5ad (in-place augment).

Output: data/processed/randolph_2021_processed_v6.h5ad rewritten with 39,240 cells
        (34,276 existing + 4,964 new).
Backup: data/processed/randolph_2021_processed_v6.pre_infected_backup.h5ad written first.
"""

from __future__ import annotations

import logging
import shutil
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore", category=UserWarning, module="rdata")

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
RDS_PATH = (
    REPO
    / "data/raw/randolph_2021/zenodo_inputs/inputs/1_calculate_pseudobulk/infected_monocytes_cluster_singlets.rds"
)
V6_PATH = REPO / "data/processed/randolph_2021_processed_v6.h5ad"
BACKUP_PATH = REPO / "data/processed/randolph_2021_processed_v6.pre_infected_backup.h5ad"


def parse_rds() -> tuple[sp.csc_matrix, list[str], list[str], pd.DataFrame]:
    """Parse Seurat rds. Returns (counts_csc 19248×4964, gene_symbols, cell_barcodes, meta_df)."""
    import rdata

    logger.info("parsing %s", RDS_PATH.name)
    parsed = rdata.parser.parse_file(str(RDS_PATH))
    conv = rdata.conversion.convert(parsed)
    rna = conv.assays[np.str_("RNA")]
    raw = rna.counts
    n_features, n_cells = int(raw.Dim[0]), int(raw.Dim[1])
    # dgCMatrix: column-compressed sparse. i = row indices, p = column pointers, x = values.
    i = np.asarray(raw.i, dtype=np.int64)
    p = np.asarray(raw.p, dtype=np.int64)
    x = np.asarray(raw.x, dtype=np.float32)
    counts = sp.csc_matrix((x, i, p), shape=(n_features, n_cells))
    gene_symbols = [str(g) for g in raw.Dimnames[0]]
    cell_barcodes = [str(c) for c in raw.Dimnames[1]]
    md = getattr(conv, "meta.data")
    # rdata returns columns as np.str_ keys; normalize to plain str
    md.columns = [str(c) for c in md.columns]
    logger.info("parsed: counts %dx%d, meta.data %s", n_features, n_cells, md.shape)
    return counts, gene_symbols, cell_barcodes, md


def align_to_v6_gene_space(
    counts: sp.csc_matrix,
    gene_symbols: list[str],
    v6_var: pd.DataFrame,
) -> sp.csr_matrix:
    """Map rds counts (gene-symbol rows) into v6 Ensembl-indexed gene space; pad zeros for missing."""
    # v6.var.index = Ensembl, v6.var['gene_symbol'] = symbol
    sym_to_ens: dict[str, str] = {}
    for ens, sym in v6_var["gene_symbol"].astype(str).items():
        if sym not in sym_to_ens:
            sym_to_ens[sym] = ens
    n_v6_genes = len(v6_var)
    ens_to_pos = {ens: i for i, ens in enumerate(v6_var.index)}
    n_cells = counts.shape[1]
    mapped = 0
    missing = 0
    counts_csr = counts.tocsr()
    new_counts = sp.lil_matrix((n_v6_genes, n_cells), dtype=np.float32)
    for row_idx, sym in enumerate(gene_symbols):
        ens = sym_to_ens.get(sym)
        if ens is None:
            missing += 1
            continue
        pos = ens_to_pos.get(ens)
        if pos is None:
            missing += 1
            continue
        row = counts_csr.getrow(row_idx)
        if row.nnz == 0:
            continue
        new_counts[pos, :] = row.toarray().flatten()
        mapped += 1
    logger.info(
        "gene mapping: %d/%d symbols mapped to v6 Ensembl; %d missing",
        mapped,
        len(gene_symbols),
        missing,
    )
    return new_counts.tocsr()


def build_obs(
    meta_df: pd.DataFrame, v6_obs: pd.DataFrame, cell_barcodes: list[str]
) -> pd.DataFrame:
    """Build obs DataFrame for new cells with full v6 schema."""
    obs = pd.DataFrame(index=cell_barcodes)
    # Map per-cell columns
    obs["pool"] = meta_df["orig.ident"].astype(str).values
    obs["donor_id"] = meta_df["SOC_indiv_ID"].astype(str).values
    obs["infection_status"] = meta_df["SOC_infection_status"].astype(str).values
    obs["randolph_celltype"] = "infected_monocytes"
    obs["bucket"] = "monocyte_infected"
    obs["cell_type_bucket_unified"] = "monocyte_infected"
    obs["donor_disease_status"] = obs["infection_status"].map(
        {"flu": "diseased", "NI": "healthy_control"}
    )
    obs["study_id"] = "randolph_2021"

    # Per-donor schema v6 columns: copy from existing v6 by donor_id
    donor_v6 = (
        v6_obs[
            [
                "donor_id",
                "donor_response_design",
                "exposure_pair_id",
                "exposure_type",
                "exposure_duration_hours",
                "age_years",
                "age_group_category",
                "infection_state",
                "donor_serostatus",
            ]
        ]
        .drop_duplicates("donor_id")
        .set_index("donor_id")
    )

    for col in [
        "donor_response_design",
        "exposure_pair_id",
        "exposure_type",
        "exposure_duration_hours",
        "age_years",
        "age_group_category",
        "infection_state",
        "donor_serostatus",
    ]:
        obs[col] = obs["donor_id"].map(donor_v6[col]).astype(donor_v6[col].dtype)

    # infection_state at cell-level: NI cells are 'naive', flu cells are 'acute' (per v6 convention)
    obs.loc[obs["infection_status"] == "NI", "infection_state"] = "naive"
    obs.loc[obs["infection_status"] == "flu", "infection_state"] = "acute"

    n_unmapped = obs["donor_response_design"].isna().sum()
    if n_unmapped:
        logger.warning(
            "%d cells have donor_id not in existing v6 obs; setting defaults", n_unmapped
        )
        obs["donor_response_design"] = obs["donor_response_design"].fillna("paired_within_donor")
        obs["exposure_type"] = obs["exposure_type"].fillna("ex_vivo_iav_6h_moi_0_5")
        obs["donor_serostatus"] = obs["donor_serostatus"].fillna("unknown")
    return obs


def main() -> int:
    if not RDS_PATH.exists():
        logger.error("rds not found: %s", RDS_PATH)
        return 1
    if not V6_PATH.exists():
        logger.error("v6 h5ad not found: %s", V6_PATH)
        return 1

    if not BACKUP_PATH.exists():
        logger.info("creating backup at %s", BACKUP_PATH.name)
        shutil.copy(V6_PATH, BACKUP_PATH)
    else:
        logger.info("backup already exists at %s; skipping", BACKUP_PATH.name)

    counts, gene_symbols, cell_barcodes, meta_df = parse_rds()
    v6 = ad.read_h5ad(V6_PATH)
    logger.info("v6 existing: %s", v6.shape)

    new_counts = align_to_v6_gene_space(counts, gene_symbols, v6.var)
    new_obs = build_obs(meta_df, v6.obs, cell_barcodes)
    logger.info(
        "new infected cells: %d cells × %d genes; %d flu (diseased) / %d NI (healthy)",
        new_counts.shape[1],
        new_counts.shape[0],
        (new_obs["donor_disease_status"] == "diseased").sum(),
        (new_obs["donor_disease_status"] == "healthy_control").sum(),
    )

    new_ann = ad.AnnData(
        X=new_counts.T.tocsr(),  # AnnData wants cells × genes
        obs=new_obs,
        var=v6.var.copy(),
    )

    # Concat (v6 + new); use unique cell barcode index — prefix new cells to avoid collision
    new_ann.obs_names = [f"infected_{n}" for n in new_ann.obs_names]
    merged = ad.concat([v6, new_ann], axis=0, join="outer", merge="same", index_unique=None)
    logger.info("merged shape: %s", merged.shape)
    logger.info(
        "merged bucket counts:\n%s",
        merged.obs["cell_type_bucket_unified"].value_counts().to_string(),
    )

    merged.write_h5ad(V6_PATH, compression="gzip")
    logger.info("wrote %s", V6_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
