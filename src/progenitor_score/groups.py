"""Proposed grouping of the 55 BoneMarrowMap states. A proposal: see docs/bm_state_groups_v1.csv (per-state held-out recall and precision)."""
GROUPS = {
    "myeloid_progenitor": ["Early GMP", "GMP-Cycle", "GMP-Mono"],
    "monocyte_precursor": ["Early ProMono", "Late ProMono"],
    "dc_precursor": ["Pre-cDC", "Pre-pDC", "Pre-pDC Cycling"],
    "granulocyte_primed": ["GMP-Neut"],
    "multipotent": ["HSC", "MPP-MyLy", "MPP-MkEry", "LMPP", "Cycling Progenitor"],
    "lymphoid_shared": ["MLP", "MLP-II"],
    "mature_mnp": ["CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "ASDC", "pDC"],
}
DECISION = {"myeloid_progenitor": "in", "dc_precursor": "in", "monocyte_precursor": "in (monocyte continuum)", "mature_mnp": "in",
            "granulocyte_primed": "review", "multipotent": "review", "lymphoid_shared": "review", "other": "out"}
IN_PROGENITOR = ("myeloid_progenitor", "dc_precursor")
REVIEW = ("granulocyte_primed", "multipotent", "lymphoid_shared")
_STATE_TO_GROUP = {s: g for g, ss in GROUPS.items() for s in ss}

def group_of(state: str) -> str:
    return _STATE_TO_GROUP.get(state, "other")
