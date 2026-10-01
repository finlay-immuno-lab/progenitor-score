from __future__ import annotations
import numpy as np, pandas as pd, scipy.sparse as sp
from .groups import group_of, IN_PROGENITOR, REVIEW
from .model import Model, load_model

CALLS = ["progenitor_in", "progenitor_review", "monocyte_precursor", "mature_mnp", "uncertain", "not_progenitor", "not_scored"]

def normalise_counts(X, target_sum: float = 1e4):
    """log1p(CP10k) from raw counts (cells x genes). The library size is the sum over the columns given, so pass only the model
    genes present in the query: that is how the model was trained."""
    X = sp.csr_matrix(X, dtype=np.float64) if sp.issparse(X) else np.asarray(X, dtype=np.float64)
    tot = np.asarray(X.sum(axis=1)).ravel(); f = np.divide(target_sum, tot, out=np.zeros_like(tot), where=tot > 0)
    X = X.multiply(f[:, None]).tocsr() if sp.issparse(X) else X * f[:, None]
    return X.log1p() if sp.issparse(X) else np.log1p(X)

def score_matrix(decision: np.ndarray, states) -> pd.DataFrame:
    """From decision values to: bm_state, runner_up, top1, margin, bm_group. Use margin/top1, not probabilities (saturated)."""
    states = np.asarray(states); order = np.argsort(decision, axis=1); r = np.arange(len(decision))
    i1, i2 = order[:, -1], order[:, -2]
    return pd.DataFrame({"bm_state": states[i1], "runner_up": states[i2], "top1": decision[r, i1], "margin": decision[r, i1] - decision[r, i2],
                         "bm_group": [group_of(s) for s in states[i1]]})

def decide(scores: pd.DataFrame, min_margin: float = 100.0, min_top1: float = 0.0) -> pd.Series:
    """progenitor_in / progenitor_review need a clear margin and a good match; a progenitor best state without them is 'uncertain'.
    monocyte_precursor (Early/Late ProMono) is part of the monocyte continuum, not a progenitor call."""
    g = scores["bm_group"]; clear = (scores["margin"] >= min_margin) & (scores["top1"] >= min_top1)
    prog = g.isin([*IN_PROGENITOR, *REVIEW])
    d = np.select([g.isin(IN_PROGENITOR) & clear, g.isin(REVIEW) & clear, prog & ~clear, g == "monocyte_precursor", g == "mature_mnp"],
                  ["progenitor_in", "progenitor_review", "uncertain", "monocyte_precursor", "mature_mnp"], "not_progenitor")
    return pd.Series(d, index=scores.index, name="progenitor_call")

def score_counts(X, genes, cell_names=None, model: Model | None = None, qc=None, min_margin: float = 100.0, min_top1: float = 0.0) -> pd.DataFrame:
    """Score any RAW-count matrix (cells x genes; dense, scipy sparse or DataFrame) whose columns are gene SYMBOLS `genes`.
    `qc`: optional boolean array; False cells get 'not_scored'. Raises if fewer than half of the model genes are found."""
    model = model or load_model()
    if isinstance(X, pd.DataFrame): cell_names = X.index if cell_names is None else cell_names; genes = X.columns if genes is None else genes; X = X.to_numpy()
    genes = np.asarray(genes, dtype=str); cell_names = pd.RangeIndex(X.shape[0]) if cell_names is None else pd.Index(cell_names)
    s = pd.Series(np.arange(len(genes)), index=genes); s = s[~s.index.duplicated(keep=False)]
    feat = [g for g in model.genes if g in s.index]
    if len(feat) < 0.5 * len(model.genes):
        raise ValueError(f"only {len(feat)} of {len(model.genes)} model genes present: are gene symbols used?")
    keep = np.ones(X.shape[0], bool) if qc is None else np.asarray(qc, bool)
    out = pd.DataFrame(index=cell_names, columns=["bm_state", "runner_up", "top1", "margin", "bm_group", "progenitor_call"], dtype=object)
    out["progenitor_call"] = "not_scored"
    if keep.any():
        sub = X[np.where(keep)[0]][:, s[feat].to_numpy()]
        sc = score_matrix(model.decision(normalise_counts(sub), feat), model.states); sc.index = cell_names[keep]
        sc["progenitor_call"] = decide(sc, min_margin, min_top1)
        out.loc[sc.index, sc.columns] = sc.values
    for c in ("top1", "margin"): out[c] = pd.to_numeric(out[c])
    return out

def score_anndata(adata, model: Model | None = None, symbol_col: str | None = None, qc_col: str | None = None,
                  min_margin: float = 100.0, min_top1: float = 0.0, layer: str | None = None) -> pd.DataFrame:
    """Raw counts in adata.X (or adata.layers[layer]). Gene symbols from var[symbol_col] (or var_names). Cells with qc_col False are 'not_scored'."""
    sym = adata.var[symbol_col].astype(str).to_numpy() if symbol_col else np.asarray(adata.var_names, dtype=str)
    qc = adata.obs[qc_col].astype(bool).to_numpy() if qc_col else None
    return score_counts(adata.X if layer is None else adata.layers[layer], sym, adata.obs_names, model, qc, min_margin, min_top1)

def annotate(adata, model: Model | None = None, symbol_col: str | None = None, qc_col: str | None = None, prefix: str = "bmps_",
             min_margin: float = 100.0, min_top1: float = 0.0, layer: str | None = None):
    """Score and write the result into adata.obs (columns prefix+bm_state, runner_up, top1, margin, bm_group, progenitor_call) in place; provenance and
    the citation go to adata.uns['progenitor_score']. Returns adata, so it chains."""
    from .cite import SHORT
    model = model or load_model()
    sc = score_anndata(adata, model, symbol_col, qc_col, min_margin, min_top1, layer)
    for c in sc.columns:
        adata.obs[prefix + c] = pd.Categorical(sc[c]) if sc[c].dtype == object else sc[c].to_numpy()
    adata.uns["progenitor_score"] = {"model_id": str(model.meta.get("model_id")), "rules_version": str(model.meta.get("rules_version")),
        "min_margin": float(min_margin), "min_top1": float(min_top1), "prefix": prefix, "citation": SHORT}
    return adata
