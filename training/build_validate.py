import os, sys, json, time, numpy as np, pandas as pd, h5py, anndata as ad, scipy.sparse as sp
os.environ.setdefault("CELLTYPIST_FOLDER", os.path.abspath("celltypist_home")); os.makedirs(os.environ["CELLTYPIST_FOLDER"], exist_ok=True)
import scanpy as sc, celltypist
from sklearn.metrics import roc_auc_score
sys.path.insert(0, "."); import bm_progenitor_scoring as bps
R = os.environ.get("BMM_DIR", ".") + "/"          # folder holding the BoneMarrowMap files (set BMM_DIR)
BMP = R + "BoneMarrowMap_Annotated_Dataset_expandedFeatures.h5ad"
XLS = R + "bcd-24-0342_supplementary_tables_s1-s21_suppst1.xlsx"
CAP, TEST_CAP, SEED = 300, 300, 0; T0 = time.time(); os.makedirs("out", exist_ok=True)
log = lambda *a: print(f"[{time.time()-T0:5.0f}s]", *a, flush=True)
bb = ad.read_h5ad(BMP, backed="r")
obs = bb.obs[["CellType", "CellType_Broad", "Study", "Donor", "Sorting"]].copy()
obs["state"] = obs.CellType.astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
sym = pd.Index(bb.var_names.astype(str)); uniq = ~np.asarray(sym.duplicated(keep=False))
# gene universe for the model: BM HVGs + top-50 markers of every state (Supp Table S2)
s2 = pd.read_excel(XLS, sheet_name="S2 - Cell State Markers"); log("S2 shape", s2.shape, list(s2.columns[:4]))
markers = {str(c).strip(): [g for g in s2[c].dropna().astype(str).tolist()][:50] for c in s2.columns}
hvg = set(sym[bb.var["HVG_intersect3000"].values.astype(bool)]) if "HVG_intersect3000" in bb.var else set()
mk_all = set(g for v in markers.values() for g in v)
universe = sorted((hvg | mk_all) & set(sym[uniq])); log("model gene universe", len(universe), "| hvg", len(hvg), "| markers", len(mk_all))
_u = set(universe); gcols = np.array([i for i, s in enumerate(sym) if s in _u and uniq[i]]); gnames = sym[gcols]
def read_rows(rows):
    rows = np.sort(rows)
    with h5py.File(BMP, "r") as f:
        ip = f["X/indptr"][:]; D, I = f["X/data"], f["X/indices"]; n = len(rows)
        data = np.empty(int((ip[rows + 1] - ip[rows]).sum()), np.float32); ind = np.empty(len(data), np.int32); ptr = np.zeros(n + 1, np.int64); pos = 0; r0 = 0
        for k in range(1, n + 1):
            if k == n or rows[k] != rows[k - 1] + 1:
                a, z = ip[rows[r0]], ip[rows[k - 1] + 1]; L = int(z - a); data[pos:pos + L] = D[a:z]; ind[pos:pos + L] = I[a:z]
                for r in range(r0, k): ptr[r + 1] = ptr[r] + (ip[rows[r] + 1] - ip[rows[r]])
                pos += L; r0 = k
    X = sp.csr_matrix((data, ind, ptr), shape=(n, len(sym)))[:, gcols].tocsr()
    return rows, X
def cap_sample(mask, by, cap, rng):
    idx = []
    for _, g in obs[mask].groupby(by, observed=True):
        ii = np.where(obs.index.isin(g.index))[0] if False else g.index.map(obs.index.get_loc).values
        idx.append(ii if len(ii) <= cap else rng.choice(ii, cap, replace=False))
    return np.concatenate(idx)
def make_adata(rows):
    rows, X = read_rows(rows); a = ad.AnnData(X, obs=obs.iloc[rows][["state", "CellType_Broad", "Study", "Donor", "Sorting"]].astype(str).copy(), var=pd.DataFrame(index=gnames))
    a.obs["group"] = a.obs.state.map(bps.group_of); return a
rng = np.random.default_rng(SEED); donors = np.array(sorted(obs.Donor.astype(str).unique())); hold = rng.choice(donors, 9, replace=False)
is_hold = obs.Donor.astype(str).isin(hold).values; log("held-out donors", list(hold))
tr_rows = cap_sample(~is_hold, ["state", "Study"], CAP, rng); te_rows = cap_sample(is_hold, ["state"], TEST_CAP, rng)
log("train", len(tr_rows), "test", len(te_rows))
tr = make_adata(tr_rows); te = make_adata(te_rows); log("loaded", tr.shape, te.shape)
tr_n = bps.normalise(tr); te_n = bps.normalise(te)
t = time.time(); model = celltypist.train(tr_n, labels="state", n_jobs=8, check_expression=True, feature_selection=False, use_SGD=True, alpha=1e-4, max_iter=200,
                                         date="2026-09-30", details="BoneMarrowMap states (Zeng 2025) holdout-donor validation model", source="Zeng et al. Blood Cancer Discov 2025")
log("trained holdout model in", round(time.time() - t), "s; classes", len(model.cell_types))
res = celltypist.annotate(te_n, model=model, majority_voting=False)
model.write('out/holdout_model.pkl'); P = res.probability_matrix; res.decision_matrix.to_csv('out/holdout_decision.csv.gz'); P.to_csv('out/holdout_prob.csv.gz'); te.obs.to_csv('out/holdout_obs.csv.gz'); pd.Series(list(tr.obs.index)).to_csv('out/train_cells.csv.gz'); pred = res.predicted_labels["predicted_labels"].values; y = te.obs.state.values
# ---- per-state metrics
rows = []
for s in sorted(set(y)):
    m = y == s; tp = (pred[m] == s).sum(); fp = ((pred == s) & ~m).sum()
    rows.append(dict(state=s, group=bps.group_of(s), n_test=int(m.sum()), recall=tp / m.sum(), precision=tp / max(tp + fp, 1)))
st = pd.DataFrame(rows); st["f1"] = 2 * st.precision * st.recall / (st.precision + st.recall).replace(0, np.nan)
# ---- crude route: mean of S2 top-50 markers (scanpy score_genes) vs classifier probability, one-vs-rest AUC
log("scoring marker sets")
auc = []
for s in sorted(set(y)):
    gl = [g for g in markers.get(s, []) if g in te_n.var_names]
    m = (y == s)
    if len(gl) < 5 or m.sum() < 5 or m.all(): continue
    sc.tl.score_genes(te_n, gl, score_name="_ms", use_raw=False)
    auc.append(dict(state=s, n_markers_used=len(gl), auc_marker_score=roc_auc_score(m, te_n.obs["_ms"].values), auc_classifier=roc_auc_score(m, P[s].values) if s in P.columns else np.nan))
st = st.merge(pd.DataFrame(auc), on="state", how="left"); st.to_csv("out/validation_by_state.csv", index=False)
# ---- group-level: classifier summed probability vs best marker score
grp_rows = []
ms = {}
for s in sorted(set(y)):
    gl = [g for g in markers.get(s, []) if g in te_n.var_names]
    if len(gl) >= 5: sc.tl.score_genes(te_n, gl, score_name="_ms", use_raw=False); ms[s] = te_n.obs["_ms"].values.copy()
MS = pd.DataFrame(ms); MSz = (MS - MS.mean()) / MS.std()
for g, ss in bps.GROUPS.items():
    m = np.isin(y, ss); cls = P[[c for c in P.columns if c in ss]].sum(axis=1).values; mk = MSz[[c for c in ss if c in MSz]].max(axis=1).values
    row = dict(group=g, decision=bps.DECISION[g], n_test=int(m.sum()), auc_classifier=roc_auc_score(m, cls), auc_marker_max=roc_auc_score(m, mk))
    for thr in (0.5, 0.7, 0.9):
        call = cls >= thr; row[f"recall_p{thr}"] = call[m].mean(); row[f"precision_p{thr}"] = call[m].sum() / max(call.sum(), 1)
    grp_rows.append(row)
pd.DataFrame(grp_rows).to_csv("out/validation_by_group.csv", index=False)
# in-progenitor call (myeloid_committed + dc_precursor)
scr = bps.score  # not used here; direct sum below
pin = P[[c for c in P.columns if bps.group_of(c) in ("myeloid_committed", "dc_precursor")]].sum(axis=1).values
truth = np.isin(pd.Series(y).map(bps.group_of), ["myeloid_committed", "dc_precursor"])
cal = [dict(threshold=t, recall=(pin >= t)[truth].mean(), precision=(truth[pin >= t]).mean() if (pin >= t).any() else np.nan, n_called=int((pin >= t).sum())) for t in (0.3, 0.5, 0.7, 0.8, 0.9, 0.95)]
pd.DataFrame(cal).to_csv("out/validation_in_progenitor_calls.csv", index=False)
# confusion within the myeloid / DC part of the tree
keep = bps.GROUPS["myeloid_progenitor"] + bps.GROUPS["dc_precursor"] + bps.GROUPS["granulocyte_primed"] + bps.GROUPS["mature_mnp"] + ["HSC", "MPP-MyLy", "LMPP", "MLP"]
conf = pd.crosstab(pd.Series(y, name="true"), pd.Series(pred, name="pred")); conf = conf.reindex(index=[k for k in keep if k in conf.index], columns=[k for k in keep if k in conf.columns], fill_value=0)
(conf.div(conf.sum(1), axis=0)).round(3).to_csv("out/validation_confusion_myeloid_fraction.csv")
log("validation done"); print(st.round(3).to_string()); print(pd.DataFrame(grp_rows).round(3).to_string()); print(pd.DataFrame(cal).round(3).to_string())
if os.path.exists('out/BoneMarrowMap_CellTypist_v1.pkl'): log('final model exists; DONE'); sys.exit(0)
# ---- final model: all donors
tr2 = make_adata(cap_sample(np.ones(len(obs), bool), ["state", "Study"], CAP, np.random.default_rng(SEED))); tr2n = bps.normalise(tr2)
t = time.time(); final = celltypist.train(tr2n, labels="state", n_jobs=8, check_expression=True, feature_selection=False, use_SGD=True, alpha=1e-4, max_iter=200, date="2026-09-30",
        details=f"BoneMarrowMap 55 states, {tr2.n_obs} cells (cap {CAP} per state x study), {len(gnames)} genes (BM HVGs + Supp S2 markers); logistic regression", source="Zeng et al. Blood Cancer Discov 2025;6:307", version="v1", url="https://github.com/andygxzeng/BoneMarrowMap")
final.write("out/BoneMarrowMap_CellTypist_v1.pkl"); log("final model", tr2.n_obs, "cells,", len(final.cell_types), "states, trained in", round(time.time() - t), "s")
json.dump(dict(cap=CAP, seed=SEED, holdout_donors=list(hold), n_train_holdout=int(len(tr_rows)), n_test=int(len(te_rows)), n_train_final=int(tr2.n_obs), n_genes=int(len(gnames))), open("out/run_info.json", "w"), indent=1)
log("DONE")
