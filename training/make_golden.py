"""Build tests/golden_*.csv.gz from the ORIGINAL CellTypist pickle (training-time only; the pickle is not in the repo).
Inputs are synthetic counts, so no real cell data is redistributed. Run from the repo root with the pickle path as argument."""
import sys, numpy as np, pandas as pd, anndata as ad, celltypist
sys.path.insert(0, "src"); from progenitor_score import load_model, normalise_counts
pkl = sys.argv[1]; m = load_model(); rng = np.random.default_rng(0); n = 40
C = rng.poisson(0.6, size=(n, len(m.genes))).astype(float)
for i in range(0, n, 8): C[i, rng.choice(len(m.genes), 5, replace=False)] = rng.integers(200, 900, 5)      # spikes: exercise the clip at 10
X = normalise_counts(C); cells = [f"golden{i}" for i in range(n)]
orig = celltypist.models.Model.load(pkl)
def annotate(Xs, genes):
    a = ad.AnnData(Xs.astype(np.float32), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
    return celltypist.annotate(a, model=orig, majority_voting=False).decision_matrix
fmt = "%.17g"; genes = list(m.genes)
pd.DataFrame(X, index=cells, columns=genes).to_csv("tests/golden_X.csv.gz", float_format=fmt)
annotate(X, genes).to_csv("tests/golden_decision_full.csv.gz", float_format=fmt)
drop = sorted(rng.choice(len(genes), 300, replace=False)); keep = [j for j in range(len(genes)) if j not in set(drop)]
pd.Series([genes[j] for j in drop]).to_csv("tests/golden_dropped_genes.txt", index=False, header=False)
annotate(X[:, keep], [genes[j] for j in keep]).to_csv("tests/golden_decision_missing.csv.gz", float_format=fmt)
print("clipped entries:", int((((X - m.mean) / m.scale) > 10).sum()))
