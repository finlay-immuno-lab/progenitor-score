"""Pseudobulk state profiles (training-time only; needs the BoneMarrowMap h5ad) used by examples/simulate_counts.py.
usage: python training/make_sim_profiles.py /path/to/BoneMarrowMap_Annotated_Dataset_expandedFeatures.h5ad
Output: examples/data/sim_state_profiles.csv.gz, genes x states, parts-per-million of counts over the 2,997 model genes (3 significant digits).
Aggregate profiles only: no single-cell data is stored."""
import sys, numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, "src"); import progenitor_score as ps
STATES = ["HSC", "MPP-MyLy", "Early GMP", "GMP-Mono", "Pre-cDC", "cDC2", "Pre-pDC", "pDC", "Early ProMono", "Late ProMono", "CD14 Mono", "CD4 Naive", "Orthochromatic Erythroblast"]
m = ps.load_model(); bb = ad.read_h5ad(sys.argv[1], backed="r"); rng = np.random.default_rng(0)
st = bb.obs.CellType.astype(str).str.replace(r"\s+", " ", regex=True).str.strip().to_numpy()
sym = pd.Index(bb.var_names.astype(str)); uniq = ~np.asarray(sym.duplicated(keep=False))
cols = np.array([sym.get_loc(g) for g in m.genes if g in sym and uniq[sym.get_loc(g)]]); genes = [g for g in m.genes if g in sym and uniq[sym.get_loc(g)]]
out = {}
for s in STATES:
    idx = np.sort(rng.choice(np.where(st == s)[0], min(400, (st == s).sum()), replace=False))
    X = bb[idx].X[:, cols]
    tot = np.asarray(X.sum(0)).ravel().astype(float); out[s] = tot / tot.sum() * 1e6; print(s, len(idx), flush=True)
pd.DataFrame(out, index=genes).apply(lambda c: c.map(lambda v: float(f"{v:.3g}"))).to_csv("examples/data/sim_state_profiles.csv.gz")
