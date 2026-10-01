"""The scoring model as plain tables: a linear classifier (coef, intercept) with the training scaler (mean, scale).
No pickle is read. Tables are verified against the sha256 values in model.json before use."""
from __future__ import annotations
import hashlib, json
from dataclasses import dataclass
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp

DEFAULT_MODEL = Path(__file__).resolve().parent / "models" / "bone_marrow_zeng2025_v1"

@dataclass
class Model:
    genes: np.ndarray          # (n_genes,) symbols
    states: list               # (n_states,)
    coef: np.ndarray           # (n_states, n_genes)
    intercept: np.ndarray      # (n_states,)
    mean: np.ndarray           # (n_genes,)
    scale: np.ndarray          # (n_genes,)
    meta: dict

    def decision(self, X, genes=None, chunk: int = 20000) -> np.ndarray:
        """Decision values (cells x states). X is log1p(CP10k), cells x genes (dense or scipy sparse); `genes` names its columns
        (default: the model's genes). Query genes missing from the model are ignored; model genes missing from the query contribute 0,
        exactly as CellTypist does."""
        genes = self.genes if genes is None else np.asarray(genes, dtype=str)
        pos = {g: i for i, g in enumerate(self.genes)}
        cols = [j for j, g in enumerate(genes) if g in pos]
        if not cols:
            raise ValueError("no query genes overlap the model genes; gene symbols are required")
        idx = np.array([pos[genes[j]] for j in cols]); W = self.coef[:, idx]; mu, sd = self.mean[idx], self.scale[idx]
        out = np.empty((X.shape[0], len(self.states)))
        for a in range(0, X.shape[0], chunk):
            blk = X[a:a + chunk][:, cols]; blk = blk.toarray() if sp.issparse(blk) else np.asarray(blk)
            blk = (blk - mu) / sd; np.minimum(blk, 10, out=blk)       # CellTypist clips above 10 only
            out[a:a + chunk] = blk @ W.T + self.intercept
        return out

    def to_celltypist(self):
        """Rebuild a celltypist.models.Model from the tables (needs the 'celltypist' extra); no pickle involved."""
        import celltypist
        from sklearn.linear_model import SGDClassifier
        from sklearn.preprocessing import StandardScaler
        clf = SGDClassifier(loss="log_loss"); clf.classes_ = np.array(self.states); clf.coef_ = self.coef.copy(); clf.intercept_ = self.intercept.copy()
        clf.n_features_in_ = len(self.genes); clf.features = np.array(self.genes)
        sc = StandardScaler(); sc.mean_ = self.mean.copy(); sc.scale_ = self.scale.copy(); sc.var_ = self.scale ** 2
        sc.n_features_in_ = len(self.genes); sc.n_samples_seen_ = 1
        desc = {k: str(self.meta.get(k, "")) for k in ("model_id", "source", "training")}; desc["number_celltypes"] = str(len(self.states))
        return celltypist.models.Model(clf, sc, desc)

def _sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_model(path=DEFAULT_MODEL, verify: bool = True) -> Model:
    path = Path(path); meta = json.loads((path / "model.json").read_text())
    if verify:
        for f, want in meta["files"].items():
            if _sha256(path / f) != want:
                raise ValueError(f"checksum mismatch for {f}: the model tables were changed")
    coef = pd.read_csv(path / "coef.csv.gz", index_col=0); sc = pd.read_csv(path / "scaler.csv"); st = pd.read_csv(path / "states.csv")
    assert list(coef.index) == list(sc["gene"]) and list(coef.columns) == list(st["state"]), "table order mismatch"
    return Model(genes=coef.index.to_numpy(dtype=str), states=list(coef.columns), coef=coef.to_numpy().T.copy(), intercept=st["intercept"].to_numpy(),
                 mean=sc["mean"].to_numpy(), scale=sc["scale"].to_numpy(), meta=meta)
