import os, sys, time, numpy as np, pandas as pd, h5py, anndata as ad, scipy.sparse as sp
os.environ.setdefault("CELLTYPIST_FOLDER", os.path.abspath("celltypist_home"))
import scanpy as sc
from scipy.special import softmax
from sklearn.metrics import roc_auc_score
sys.path.insert(0, "."); import bm_progenitor_scoring as bps
R = os.environ.get("BMM_DIR", ".") + "/"          # folder holding the BoneMarrowMap files (set BMM_DIR)
BMP = R + "BoneMarrowMap_Annotated_Dataset_expandedFeatures.h5ad"
T0 = time.time(); log = lambda *a: print(f"[{time.time()-T0:5.0f}s]", *a, flush=True)
D = pd.read_csv("out/holdout_decision.csv.gz", index_col=0); P = pd.read_csv("out/holdout_prob.csv.gz", index_col=0); O = pd.read_csv("out/holdout_obs.csv.gz", index_col=0)
y = O.state.values; ytrue_g = pd.Series(y).map(bps.group_of).values
S = pd.DataFrame(softmax(D.values, axis=1), index=D.index, columns=D.columns)          # mutually exclusive weights across the 55 states
pred = D.columns[D.values.argmax(1)]; pg = pd.Series(pred).map(bps.group_of).values
# ---- 1. group-level confusion by argmax
rows = []
for g in [*bps.GROUPS, "other"]:
    m = ytrue_g == g; c = pg == g
    rows.append(dict(group=g, decision=bps.DECISION[g], n_true=int(m.sum()), recall=(c & m).sum() / max(m.sum(), 1), precision=(c & m).sum() / max(c.sum(), 1)))
pd.DataFrame(rows).to_csv("out/validation_group_argmax.csv", index=False)
cg = pd.crosstab(pd.Series(ytrue_g, name="true"), pd.Series(pg, name="called"), normalize="index").round(3); cg.to_csv("out/validation_group_confusion.csv")
# ---- 2. in-progenitor call vs softmax confidence
pin_state = [c for c in D.columns if bps.group_of(c) in ("myeloid_committed", "dc_precursor")]
p_in = S[pin_state].sum(axis=1).values; truth = np.isin(ytrue_g, ["myeloid_committed", "dc_precursor"])
cal = []
for t in (0.3, 0.5, 0.7, 0.8, 0.9):
    call = p_in >= t; cal.append(dict(softmax_threshold=t, n_called=int(call.sum()), recall=call[truth].mean(), precision=truth[call].mean() if call.any() else np.nan,
                                      false_calls_from_mature_mnp=int((call & (ytrue_g == "mature_mnp")).sum()), false_calls_from_multipotent=int((call & (ytrue_g == "multipotent")).sum())))
pd.DataFrame(cal).to_csv("out/validation_in_progenitor_softmax.csv", index=False); print(pd.DataFrame(cal).round(3).to_string())
# ---- 3. leakage-free marker route: markers re-derived from the TRAINING cells only
tr_ids = pd.read_csv("out/train_cells.csv.gz", index_col=0).iloc[:, 0].values
bb = ad.read_h5ad(BMP, backed="r"); sym = pd.Index(bb.var_names.astype(str)); uniq = ~np.asarray(sym.duplicated(keep=False))
model_genes = pd.read_csv("out/holdout_decision.csv.gz", nrows=1).shape  # placeholder to keep shapes explicit
import celltypist; m = celltypist.models.Model.load("out/holdout_model.pkl"); genes = list(m.features); gset = set(genes)
gcols = np.array([i for i, s in enumerate(sym) if s in gset and uniq[i]]); gnames = sym[gcols]
def read_ids(ids):
    pos = bb.obs_names.get_indexer(ids); rows = np.sort(pos[pos >= 0])
    with h5py.File(BMP, "r") as f:
        ip = f["X/indptr"][:]; Dd, I = f["X/data"], f["X/indices"]; n = len(rows)
        data = np.empty(int((ip[rows + 1] - ip[rows]).sum()), np.float32); ind = np.empty(len(data), np.int32); ptr = np.zeros(n + 1, np.int64); pos_ = 0; r0 = 0
        for k in range(1, n + 1):
            if k == n or rows[k] != rows[k - 1] + 1:
                a, z = ip[rows[r0]], ip[rows[k - 1] + 1]; L = int(z - a); data[pos_:pos_ + L] = Dd[a:z]; ind[pos_:pos_ + L] = I[a:z]
                for r in range(r0, k): ptr[r + 1] = ptr[r] + (ip[rows[r] + 1] - ip[rows[r]])
                pos_ += L; r0 = k
    X = sp.csr_matrix((data, ind, ptr), shape=(n, len(sym)))[:, gcols].tocsr()
    return rows, X
rows, X = read_ids(tr_ids); a = ad.AnnData(X, obs=bb.obs.iloc[rows][["CellType"]].astype(str).copy(), var=pd.DataFrame(index=gnames))
a.obs["state"] = a.obs.CellType.str.replace(r"\s+", " ", regex=True).str.strip(); a = bps.normalise(a); log("train cells", a.shape)
sc.tl.rank_genes_groups(a, "state", method="t-test", n_genes=50, use_raw=False); rg = a.uns["rank_genes_groups"]["names"]
mk_train = {s: [g for g in rg[s]] for s in rg.dtype.names}
# held-out expression
te_ids = O.index.values; rows, X = read_ids(te_ids); t = ad.AnnData(X, obs=O.loc[bb.obs_names[rows]].copy() if False else pd.DataFrame(index=bb.obs_names[rows]), var=pd.DataFrame(index=gnames)); t = bps.normalise(t)
order = pd.Index(t.obs_names).get_indexer(O.index); assert (order >= 0).all(); t = t[order].copy()
ms = {}
for s, gl in mk_train.items():
    sc.tl.score_genes(t, gl, score_name="_m", use_raw=False); ms[s] = t.obs["_m"].values.copy()
MS = pd.DataFrame(ms, index=O.index); MSz = (MS - MS.mean()) / MS.std(); out = []
for s in sorted(set(y)):
    mm = y == s
    if mm.sum() >= 5 and s in MS: out.append(dict(state=s, auc_marker_train_only=roc_auc_score(mm, MS[s].values), auc_classifier_softmax=roc_auc_score(mm, S[s].values)))
pd.DataFrame(out).to_csv("out/validation_auc_by_state_fair.csv", index=False)
g_out = []
for g, ss in bps.GROUPS.items():
    mm = np.isin(y, ss); cls = S[[c for c in S.columns if c in ss]].sum(axis=1).values; mk = MSz[[c for c in ss if c in MSz]].max(axis=1).values
    g_out.append(dict(group=g, decision=bps.DECISION[g], n_test=int(mm.sum()), auc_classifier_softmax=roc_auc_score(mm, cls), auc_marker_train_only_max=roc_auc_score(mm, mk)))
pd.DataFrame(g_out).to_csv("out/validation_auc_by_group_fair.csv", index=False)
pd.to_pickle(mk_train, "out/markers_train_only.pkl")
print(pd.DataFrame(g_out).round(3).to_string()); print(pd.DataFrame(rows if False else out).round(3).to_string()); log("DONE")
