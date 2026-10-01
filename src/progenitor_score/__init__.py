"""progenitor-score: call hematopoietic progenitor types by scoring cells against the BoneMarrowMap atlas (Zeng et al. 2025)."""
from .groups import GROUPS, DECISION, IN_PROGENITOR, REVIEW, group_of
from .model import Model, load_model, DEFAULT_MODEL
from .score import normalise_counts, score_matrix, score_counts, score_anndata, annotate, decide
from .cite import citation

__version__ = "0.1.0"
