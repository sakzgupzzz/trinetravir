"""Phase 3.5 validation: produces all tables Session 3 needs to decide model
choice and granularity. Reads <study>_reannotated_{low,high}.h5ad and emits:

  results/tables/phase35_bucket_sizes_<model>.csv
  results/tables/phase35_confusion_<study>_<model>.csv
  results/tables/phase35_label_vocab_overlap_<model>.csv
  results/tables/phase35_lee_lymphoid_breakdown.csv
  results/tables/phase35_bucket_conflict_per_cell.csv
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import anndata as ad
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
PROCESSED = REPO / "data" / "processed"
TABLES = REPO / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

STUDIES = ("wilk_2020", "arunachalam_2020", "lee_2020", "schulte_schrepping_2020")
MODELS = ("low", "high")


def load(study: str, model: str) -> ad.AnnData:
    p = PROCESSED / f"{study}_reannotated_{model}.h5ad"
    return ad.read_h5ad(p)


def per_study_bucket_sizes() -> None:
    rows = []
    for model in MODELS:
        for study in STUDIES:
            a = load(study, model)
            obs = a.obs
            for bucket in sorted(obs["cell_type_bucket_unified"].astype(str).unique()):
                mask = obs["cell_type_bucket_unified"].astype(str) == bucket
                d = obs[mask & (obs["donor_disease_status"].astype(str) == "diseased")]
                h = obs[mask & (obs["donor_disease_status"].astype(str) == "healthy_control")]
                rows.append(
                    {
                        "model": model,
                        "study": study,
                        "bucket": bucket,
                        "n_cells_diseased": int(len(d)),
                        "n_cells_healthy": int(len(h)),
                        "n_donors_diseased": int(d["donor_id"].astype(str).nunique()),
                        "n_donors_healthy": int(h["donor_id"].astype(str).nunique()),
                    }
                )
    df = pd.DataFrame(rows)
    out = TABLES / "phase35_bucket_sizes_combined.csv"
    df.to_csv(out, index=False)
    print(f"wrote {out}")
    for model in MODELS:
        sub = df[df["model"] == model]
        out = TABLES / f"phase35_bucket_sizes_{model}.csv"
        sub.to_csv(out, index=False)
        print(f"wrote {out}")


def confusion_matrices() -> None:
    for model in MODELS:
        for study in STUDIES:
            a = load(study, model)
            ct_orig = a.obs["cell_type_original"].astype(str)
            ct_bucket = a.obs["cell_type_bucket_unified"].astype(str)
            xtab = pd.crosstab(ct_orig, ct_bucket)
            out = TABLES / f"celltypist_confusion_{study}_{model}.csv"
            xtab.to_csv(out)
            print(f"wrote {out}: {xtab.shape}")


def label_vocab_overlap() -> None:
    """Per model, list which study uses each unified label."""
    for model in MODELS:
        vocab = {
            study: set(load(study, model).obs["cell_type_unified"].astype(str).unique())
            for study in STUDIES
        }
        all_labels = sorted({lbl for s in vocab.values() for lbl in s})
        rows = []
        for lbl in all_labels:
            present = [s for s in STUDIES if lbl in vocab[s]]
            rows.append(
                {
                    "label": lbl,
                    "n_studies": len(present),
                    "studies": ",".join(present),
                }
            )
        df = pd.DataFrame(rows).sort_values(["n_studies", "label"], ascending=[False, True])
        out = TABLES / f"phase35_label_vocab_overlap_{model}.csv"
        df.to_csv(out, index=False)
        print(f"wrote {out}: {len(df)} labels")


def lee_lymphoid_breakdown() -> None:
    """Lee-specific: does Immune_All_Low surface naive B + Tcm/Tem CD8 sublineages?"""
    rows = []
    for model in MODELS:
        a = load("lee_2020", model)
        for bucket in ("B", "CD8T"):
            sub = a.obs[a.obs["cell_type_bucket_unified"].astype(str) == bucket]
            counts = Counter(sub["cell_type_unified"].astype(str))
            for lbl, n in counts.most_common():
                rows.append(
                    {
                        "model": model,
                        "bucket": bucket,
                        "unified_label": lbl,
                        "n_cells": int(n),
                    }
                )
        # also under High, B and T (no CD8T):
        if model == "high":
            for bucket in ("T",):
                sub = a.obs[a.obs["cell_type_bucket_unified"].astype(str) == bucket]
                counts = Counter(sub["cell_type_unified"].astype(str))
                for lbl, n in counts.most_common():
                    rows.append(
                        {
                            "model": model,
                            "bucket": bucket,
                            "unified_label": lbl,
                            "n_cells": int(n),
                        }
                    )
    out = TABLES / "phase35_lee_lymphoid_breakdown.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"wrote {out}")


def cross_model_bucket_conflict() -> None:
    """Per study, count cells where Low and High disagree on coarse bucket
    (after collapsing Low CD4T+CD8T -> T for parity)."""
    rows = []
    for study in STUDIES:
        a_low = load(study, "low")
        a_high = load(study, "high")
        # Align cells by index. Low has more buckets; collapse to Low_collapsed
        # comparable to High: CD4T, CD8T -> T; NK -> nk_only_in_low; others as is.
        if not a_low.obs.index.equals(a_high.obs.index):
            print(f"WARNING {study}: index mismatch between Low and High")
            # Reindex to intersection.
            common = a_low.obs.index.intersection(a_high.obs.index)
            low_bk = a_low.obs.loc[common, "cell_type_bucket_unified"].astype(str)
            high_bk = a_high.obs.loc[common, "cell_type_bucket_unified"].astype(str)
        else:
            low_bk = a_low.obs["cell_type_bucket_unified"].astype(str)
            high_bk = a_high.obs["cell_type_bucket_unified"].astype(str)
        low_collapsed = low_bk.replace({"CD4T": "T", "CD8T": "T", "NK": "NK_only_low"})
        xtab = pd.crosstab(low_collapsed, high_bk)
        out = TABLES / f"phase35_bucket_conflict_{study}.csv"
        xtab.to_csv(out)
        print(f"wrote {out}: {xtab.shape}")
        # Headline rows: total cells per (low_collapsed, high) pair.
        for low_label, high_row in xtab.iterrows():
            for high_label, n in high_row.items():
                if n > 0:
                    rows.append(
                        {
                            "study": study,
                            "low_bucket_collapsed": low_label,
                            "high_bucket": high_label,
                            "n_cells": int(n),
                        }
                    )
    out = TABLES / "phase35_bucket_conflict_per_cell.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"wrote {out}")


def main() -> None:
    per_study_bucket_sizes()
    confusion_matrices()
    label_vocab_overlap()
    lee_lymphoid_breakdown()
    cross_model_bucket_conflict()


if __name__ == "__main__":
    main()
