"""Map CellTypist Immune_All_Low predicted labels to our 5 coarse buckets.

Phase 3.5 mapping (2026-05-10). The Immune_All_Low model exposes 98
predicted labels. We collapse them to monocyte / CD4T / CD8T / B / NK.

Anything in the model's vocabulary that we cannot cleanly assign — bone-
marrow progenitors not expected in PBMC, dendritic cells, granulocytes,
erythrocytes, mast cells, gamma-delta T cells, NKT cells, cycling T
cells of unknown lineage, tissue-resident macrophages — maps to
``other`` and is filtered out before the per-bucket Harmony / gate
computation. ``other`` is also the fallback for any label CellTypist
returns that is not in our explicit mapping.

Distinct from ``trinetravir.data.harmonize.coarse_cell_type``, which
maps the original cellxgene ``cell_type`` strings using substring rules.
This module is mapping-table-driven because CellTypist labels are a
controlled vocabulary of known short strings.
"""

from __future__ import annotations

from typing import Final

# Explicit per-label mapping. Source: CellTypist Immune_All_Low.pkl
# vocabulary as dumped from `Model.cell_types` (n=98, listed 2026-05-10).
LABEL_TO_BUCKET: Final[dict[str, str]] = {
    # ---- monocyte ------------------------------------------------------
    "Classical monocytes": "monocyte",
    "Non-classical monocytes": "monocyte",
    "Monocytes": "monocyte",
    "Cycling monocytes": "monocyte",
    "Mono-mac": "monocyte",
    "Monocyte precursor": "monocyte",
    "MNP": "monocyte",  # Mononuclear phagocyte progenitor
    "Macrophages": "monocyte",
    "Alveolar macrophages": "monocyte",
    "Erythrophagocytic macrophages": "monocyte",
    "Hofbauer cells": "monocyte",
    "Intermediate macrophages": "monocyte",
    "Intestinal macrophages": "monocyte",
    "Kidney-resident macrophages": "monocyte",
    "Kupffer cells": "monocyte",
    # ---- B -------------------------------------------------------------
    "Naive B cells": "B",
    "Memory B cells": "B",
    "Age-associated B cells": "B",
    "Cycling B cells": "B",
    "Follicular B cells": "B",
    "Germinal center B cells": "B",
    "Proliferative germinal center B cells": "B",
    "Plasma cells": "B",
    "Plasmablasts": "B",
    "Transitional B cells": "B",
    "B cells": "B",
    # Note: Pre-B / Pro-B / Pre-pro-B are bone-marrow precursors; not
    # expected in PBMC. They map to "other" via the unknown fallback
    # below so we can detect mis-annotation rather than absorb it.
    # ---- NK ------------------------------------------------------------
    "NK cells": "NK",
    "CD16+ NK cells": "NK",
    "CD16- NK cells": "NK",
    "Cycling NK cells": "NK",
    "Transitional NK": "NK",
    # NKT cells excluded (hybrid; ambiguous lineage). Mapped to "other".
    # ---- CD4T ----------------------------------------------------------
    "Tcm/Naive helper T cells": "CD4T",
    "Tem/Effector helper T cells": "CD4T",
    "Tem/Effector helper T cells PD1+": "CD4T",
    "Follicular helper T cells": "CD4T",
    "Type 1 helper T cells": "CD4T",
    "Type 17 helper T cells": "CD4T",
    "Regulatory T cells": "CD4T",
    "Treg(diff)": "CD4T",
    "Memory CD4+ cytotoxic T cells": "CD4T",
    # ---- CD8T ----------------------------------------------------------
    "Tcm/Naive cytotoxic T cells": "CD8T",
    "Tem/Temra cytotoxic T cells": "CD8T",
    "Tem/Trm cytotoxic T cells": "CD8T",
    "Trm cytotoxic T cells": "CD8T",
    "CD8a/a": "CD8T",
    "CD8a/b(entry)": "CD8T",
    "MAIT cells": "CD8T",
}

# All Immune_All_Low labels we explicitly route to "other". Stored as a
# frozenset for fast contains. Used only by the validation tooling to
# distinguish "deliberately excluded" from "unmapped / unexpected".
EXPLICIT_OTHER: Final[frozenset[str]] = frozenset(
    {
        "NKT cells",
        "Cycling T cells",  # CD4 vs CD8 ambiguous
        "T(agonist)",
        "CRTAM+ gamma-delta T cells",
        "gamma-delta T cells",
        "Cycling gamma-delta T cells",
        "DC",
        "DC1",
        "DC2",
        "DC3",
        "DC precursor",
        "Migratory DCs",
        "Cycling DCs",
        "Transitional DC",
        "pDC",
        "pDC precursor",
        "ILC",
        "ILC1",
        "ILC2",
        "ILC3",
        "ILC precursor",
        "Mast cells",
        "Granulocytes",
        "Neutrophils",
        "Myelocytes",
        "Promyelocytes",
        "Neutrophil-myeloid progenitor",
        "Megakaryocytes/platelets",
        "Megakaryocyte precursor",
        "Megakaryocyte-erythroid-mast cell progenitor",
        "Early MK",
        "HSC/MPP",
        "MEMP",
        "GMP",
        "CMP",
        "ELP",
        "ETP",
        "Early lymphoid/T lymphoid",
        "Early erythroid",
        "Mid erythroid",
        "Late erythroid",
        "Erythrocytes",
        "Endothelial cells",
        "Epithelial cells",
        "Fibroblasts",
        "Double-negative thymocytes",
        "Double-positive thymocytes",
        "Large pre-B cells",
        "Small pre-B cells",
        "Pre-pro-B cells",
        "Pro-B cells",
    }
)


def map_label_to_bucket(label: str) -> str:
    """Map one CellTypist Immune_All_Low label to one of the 5 buckets, or 'other'.

    Parameters
    ----------
    label
        Predicted label string from CellTypist (cell_type_unified).

    Returns
    -------
    bucket
        ``monocyte`` / ``CD4T`` / ``CD8T`` / ``B`` / ``NK`` / ``other``.
        Labels not in either ``LABEL_TO_BUCKET`` or ``EXPLICIT_OTHER``
        also return ``other`` and should be logged as unexpected by the
        caller.
    """
    return LABEL_TO_BUCKET.get(label, "other")


def is_unexpected_label(label: str) -> bool:
    """Return True for labels that are neither mapped nor in EXPLICIT_OTHER.

    These are CellTypist outputs we did not anticipate. The caller should
    log them so the bucket-map can be reviewed.
    """
    return label not in LABEL_TO_BUCKET and label not in EXPLICIT_OTHER


# ---------------------------------------------------------------------------
# Sub-bucket map for Immune_All_Low (finer granularity, Issue 2 sensitivity).
# Splits monocyte/B/CD4T/CD8T/NK into biologically meaningful sub-types so
# the cross-study + cross-virus gate can be run at finer granularity to
# show the headline result is robust to bucket choice.
# Labels not present in LABEL_TO_SUBBUCKET_LOW fall through to map_label_to_bucket
# (i.e., they keep the coarse bucket name as their sub-bucket).
# ---------------------------------------------------------------------------
LABEL_TO_SUBBUCKET_LOW: Final[dict[str, str]] = {
    # monocyte sub-buckets
    "Classical monocytes": "mono_classical",
    "Non-classical monocytes": "mono_nonclassical",
    "Cycling monocytes": "mono_cycling",
    "Monocytes": "mono_unspecified",
    "Mono-mac": "mono_unspecified",
    "Monocyte precursor": "mono_unspecified",
    "MNP": "mono_unspecified",
    "Macrophages": "mono_macrophage",
    "Alveolar macrophages": "mono_macrophage",
    "Erythrophagocytic macrophages": "mono_macrophage",
    "Hofbauer cells": "mono_macrophage",
    "Intermediate macrophages": "mono_macrophage",
    "Intestinal macrophages": "mono_macrophage",
    "Kidney-resident macrophages": "mono_macrophage",
    "Kupffer cells": "mono_macrophage",
    # B sub-buckets
    "Naive B cells": "B_naive",
    "Memory B cells": "B_memory",
    "Age-associated B cells": "B_memory",
    "Plasma cells": "B_plasma",
    "Plasmablasts": "B_plasma",
    "Transitional B cells": "B_transitional",
    "Follicular B cells": "B_memory",
    "Germinal center B cells": "B_gc",
    "Proliferative germinal center B cells": "B_gc",
    "Cycling B cells": "B_cycling",
    "B cells": "B_unspecified",
    # NK sub-buckets
    "CD16+ NK cells": "NK_cd16pos",
    "CD16- NK cells": "NK_cd16neg",
    "NK cells": "NK_unspecified",
    "Cycling NK cells": "NK_cycling",
    "Transitional NK": "NK_transitional",
    # CD4T sub-buckets
    "Tcm/Naive helper T cells": "CD4T_naive_cm",
    "Tem/Effector helper T cells": "CD4T_em",
    "Tem/Effector helper T cells PD1+": "CD4T_em",
    "Follicular helper T cells": "CD4T_tfh",
    "Type 1 helper T cells": "CD4T_th1",
    "Type 17 helper T cells": "CD4T_th17",
    "Regulatory T cells": "CD4T_treg",
    "Treg(diff)": "CD4T_treg",
    "Memory CD4+ cytotoxic T cells": "CD4T_cyt",
    # CD8T sub-buckets
    "Tcm/Naive cytotoxic T cells": "CD8T_naive_cm",
    "Tem/Temra cytotoxic T cells": "CD8T_em_temra",
    "Tem/Trm cytotoxic T cells": "CD8T_em_trm",
    "Trm cytotoxic T cells": "CD8T_trm",
    "CD8a/a": "CD8T_cd8aa",
    "CD8a/b(entry)": "CD8T_cd8ab",
    "MAIT cells": "CD8T_mait",
}


def map_label_to_subbucket_low(label: str) -> str:
    """Map a CellTypist Immune_All_Low label to a finer sub-bucket.

    Sub-buckets prefix the coarse bucket (e.g., ``mono_classical``,
    ``B_naive``) so a downstream consumer can recover the 5-bucket level
    by string-prefix. Labels not in ``LABEL_TO_SUBBUCKET_LOW`` but in
    ``LABEL_TO_BUCKET`` keep the coarse bucket name as their sub-bucket
    (no further subdivision available for that label). Labels not in
    either map return ``other``.
    """
    if label in LABEL_TO_SUBBUCKET_LOW:
        return LABEL_TO_SUBBUCKET_LOW[label]
    if label in LABEL_TO_BUCKET:
        return LABEL_TO_BUCKET[label]
    return "other"


# ---------------------------------------------------------------------------
# Immune_All_High bucket map (Issue 12 model-choice sensitivity).
# Immune_All_High has only 32 labels and is *coarser* than Immune_All_Low.
# It does NOT split CD4 vs CD8 (only "T cells") and has no NK label. So the
# 5-bucket scheme collapses to:
#   monocyte (Monocytes / Mono-mac / Monocyte precursor / Macrophages / MNP)
#   B        (B cells / B-cell lineage / Plasma cells)
#   T        (T cells)  <- CD4 and CD8 NOT distinguishable
# The "NK" bucket is unrecoverable from Immune_All_High; cells assigned
# "ILC" are mapped to "other" because ILC at the High level conflates
# NK with other innate lymphoid lineages. This is the documented
# methodological caveat for the Issue 12 comparison.
# ---------------------------------------------------------------------------
LABEL_TO_BUCKET_HIGH: Final[dict[str, str]] = {
    # monocyte
    "Monocytes": "monocyte",
    "Mono-mac": "monocyte",
    "Monocyte precursor": "monocyte",
    "Macrophages": "monocyte",
    "MNP": "monocyte",
    # B
    "B cells": "B",
    "B-cell lineage": "B",
    "Plasma cells": "B",
    # T (no CD4/CD8 distinction in High)
    "T cells": "T",
    "Cycling cells": "T",  # rough; mostly T in PBMC. Caveat in methods.
}

EXPLICIT_OTHER_HIGH: Final[frozenset[str]] = frozenset(
    {
        "DC",
        "DC precursor",
        "pDC",
        "pDC precursor",
        "Granulocytes",
        "Myelocytes",
        "Promyelocytes",
        "Mast cells",
        "ILC",
        "ILC precursor",
        "HSC/MPP",
        "ETP",
        "Erythrocytes",
        "Erythroid",
        "Early MK",
        "Megakaryocytes/platelets",
        "Megakaryocyte precursor",
        "Double-negative thymocytes",
        "Double-positive thymocytes",
        "Endothelial cells",
        "Epithelial cells",
        "Fibroblasts",
    }
)


def map_label_to_bucket_high(label: str) -> str:
    """Map one CellTypist Immune_All_High label to monocyte / B / T / other.

    The High model cannot resolve CD4 vs CD8 vs NK; this map collapses
    the 5-bucket scheme to 3 buckets plus ``other``. Used only for the
    Issue 12 model-choice sensitivity comparison; downstream consumers
    must NOT mix High and Low bucket assignments without explicit
    documentation of the asymmetry.
    """
    return LABEL_TO_BUCKET_HIGH.get(label, "other")


def is_unexpected_label_high(label: str) -> bool:
    """Return True for High labels not in LABEL_TO_BUCKET_HIGH or EXPLICIT_OTHER_HIGH."""
    return label not in LABEL_TO_BUCKET_HIGH and label not in EXPLICIT_OTHER_HIGH
