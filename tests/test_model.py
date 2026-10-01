import json, shutil
from pathlib import Path
import numpy as np, pandas as pd, pytest
import progenitor_score as ps
from progenitor_score.model import DEFAULT_MODEL

T = Path(__file__).resolve().parent
X = pd.read_csv(T / "golden_X.csv.gz", index_col=0); FULL = pd.read_csv(T / "golden_decision_full.csv.gz", index_col=0)
MISS = pd.read_csv(T / "golden_decision_missing.csv.gz", index_col=0); DROPPED = set((T / "golden_dropped_genes.txt").read_text().split())

def test_tables_verify_and_tamper_is_refused(tmp_path):
    ps.load_model()
    shutil.copytree(DEFAULT_MODEL, tmp_path / "m"); p = tmp_path / "m" / "scaler.csv"; p.write_text(p.read_text() + " ")
    with pytest.raises(ValueError, match="checksum"): ps.load_model(tmp_path / "m")

def test_numpy_decision_equals_original_celltypist_full_genes():
    m = ps.load_model(); D = m.decision(X.to_numpy(), list(X.columns))
    assert np.allclose(D, FULL.to_numpy(), rtol=1e-4, atol=1e-2) and (D.argmax(1) == FULL.to_numpy().argmax(1)).all()

def test_numpy_decision_equals_original_with_missing_genes():
    m = ps.load_model(); keep = [g for g in X.columns if g not in DROPPED]
    D = m.decision(X[keep].to_numpy(), keep); assert len(keep) == len(X.columns) - 300
    assert np.allclose(D, MISS.to_numpy(), rtol=1e-4, atol=1e-2)

def test_the_clip_at_10_is_exercised_and_matters():
    m = ps.load_model(); S = (X.to_numpy() - m.mean) / m.scale
    assert (S > 10).sum() > 0
    unclipped = S @ m.coef.T + m.intercept; assert not np.allclose(unclipped, FULL.to_numpy(), rtol=1e-4, atol=1e-2)

def test_rebuilt_celltypist_model_matches_original():
    celltypist = pytest.importorskip("celltypist"); import anndata as ad
    cm = ps.load_model().to_celltypist()
    a = ad.AnnData(X.to_numpy().astype(np.float32), obs=pd.DataFrame(index=X.index), var=pd.DataFrame(index=X.columns))
    res = celltypist.annotate(a, model=cm, majority_voting=False)
    assert np.allclose(res.decision_matrix.to_numpy(), FULL.to_numpy(), rtol=1e-4, atol=1e-2)
    assert (res.decision_matrix.columns == FULL.columns).all()

def test_normalise_counts_gives_cp10k():
    C = np.random.default_rng(1).poisson(1.0, size=(20, 50)).astype(float); C[0] = 0
    L = ps.normalise_counts(C); s = np.expm1(L).sum(1)
    assert np.allclose(s[1:], 1e4) and s[0] == 0

def test_decide_rules_and_groups():
    sc = pd.DataFrame({"bm_group": ["myeloid_progenitor", "myeloid_progenitor", "multipotent", "monocyte_precursor", "mature_mnp", "other"],
                       "margin": [300, 20, 250, 500, 500, 500], "top1": [50, 50, 50, 50, 50, 50]})
    assert list(ps.decide(sc)) == ["progenitor_in", "uncertain", "progenitor_review", "monocyte_precursor", "mature_mnp", "not_progenitor"]
    assert ps.group_of("Late ProMono") == "monocyte_precursor" and ps.group_of("CD4 Naive") == "other"

def test_score_anndata_not_scored_and_wrong_gene_ids():
    import anndata as ad, scipy.sparse as sp
    m = ps.load_model(); rng = np.random.default_rng(0)
    C = sp.csr_matrix(rng.poisson(0.4, size=(30, len(m.genes))).astype(np.float32)); obs = pd.DataFrame({"qc": [True] * 20 + [False] * 10}, index=[f"c{i}" for i in range(30)])
    a = ad.AnnData(C, obs=obs, var=pd.DataFrame({"symbol": m.genes}, index=[f"ENSG{i}" for i in range(len(m.genes))]))
    out = ps.score_anndata(a, m, symbol_col="symbol", qc_col="qc")
    assert (out.progenitor_call.iloc[20:] == "not_scored").all() and set(out.progenitor_call) <= set(ps.score.CALLS)
    with pytest.raises(ValueError, match="model genes present"): ps.score_anndata(a, m)       # Ensembl-style var_names, no symbol_col


# ---- end-to-end entry points ----------------------------------------------------------------------------------------------------------------
E2E = pd.read_csv(T / "e2e_counts.csv.gz", index_col=0); E2E_QC = pd.read_csv(T / "e2e_qc.csv", index_col=0)["qc"]; E2E_EXP = pd.read_csv(T / "e2e_expected.csv", index_col=0)

def test_score_counts_dataframe_and_sparse_agree_with_fixture():
    import scipy.sparse as sp
    a = ps.score_counts(E2E.T.to_numpy(), E2E.index, E2E.columns, qc=E2E_QC.to_numpy())
    b = ps.score_counts(sp.csr_matrix(E2E.T.to_numpy()), E2E.index, E2E.columns, qc=E2E_QC.to_numpy())
    for r in (a, b):
        assert (r.progenitor_call.to_numpy() == E2E_EXP.progenitor_call.to_numpy()).all() and (r.bm_state.dropna().to_numpy() == E2E_EXP.bm_state.dropna().to_numpy()).all()
        assert np.allclose(r.top1.dropna(), E2E_EXP.top1.dropna(), rtol=1e-6, atol=1e-6)
    assert (a.progenitor_call[~E2E_QC.to_numpy()] == "not_scored").all()

def test_annotate_writes_obs_and_provenance():
    import anndata as ad
    a = ad.AnnData(E2E.T.to_numpy().astype(np.float32), obs=pd.DataFrame({"qc": E2E_QC.to_numpy()}, index=E2E.columns), var=pd.DataFrame(index=E2E.index))
    out = ps.annotate(a, qc_col="qc")
    assert out is a and {"bmps_bm_state", "bmps_margin", "bmps_progenitor_call"} <= set(a.obs.columns)
    assert (a.obs["bmps_progenitor_call"].astype(str).to_numpy() == E2E_EXP.progenitor_call.to_numpy()).all()
    assert "Zeng" in a.uns["progenitor_score"]["citation"] and a.uns["progenitor_score"]["model_id"]

def test_cli_cite_and_run(tmp_path, capsys):
    import anndata as ad
    from progenitor_score.cli import main
    main(["--cite"]); o = capsys.readouterr().out
    assert "Zeng" in o and "CC BY-NC 4.0" in o
    a = ad.AnnData(E2E.T.to_numpy().astype(np.float32), obs=pd.DataFrame({"qc": E2E_QC.to_numpy()}, index=E2E.columns), var=pd.DataFrame(index=E2E.index))
    a.write_h5ad(tmp_path / "x.h5ad"); main([str(tmp_path / "x.h5ad"), "-o", str(tmp_path / "s.csv"), "--qc-col", "qc"])
    s = pd.read_csv(tmp_path / "s.csv", index_col=0); assert (s.progenitor_call.to_numpy() == E2E_EXP.progenitor_call.to_numpy()).all()

def test_r_package_model_copy_is_identical():
    src = DEFAULT_MODEL; cp = T.parent / "inst" / "extdata" / "models" / src.name
    if not cp.exists(): pytest.skip("not in a repo checkout")
    for f in src.iterdir(): assert (cp / f.name).read_bytes() == f.read_bytes(), f"{f.name} differs from the R copy: re-copy inst/extdata/models"

def test_attribution_and_licence_files_present():
    root = T.parent
    if not (root / "LICENSE").exists(): pytest.skip("not in a repo checkout")
    for f in ("LICENSE", "LICENSE-MODEL.md", "ATTRIBUTION.md", "CITATION.cff"): assert (root / f).exists()
    assert "CC BY-NC 4.0" in ps.load_model().meta["licence"] and "Zeng" in ps.load_model().meta["attribution"]
