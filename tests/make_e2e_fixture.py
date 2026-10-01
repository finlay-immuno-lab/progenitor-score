"""Synthetic RAW counts (genes x cells, shaped like a real object) + the expected Python output, for the end-to-end tests of both languages.
Includes: genes outside the model, a duplicated symbol, a zero-count cell, three cells with qc = FALSE, shuffled gene order."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "src"); import progenitor_score as ps
m = ps.load_model(); rng = np.random.default_rng(7); n = 40
genes = list(m.genes) + [f"EXTRA{i}" for i in range(200)] + [m.genes[0]]          # last = duplicate of the first model gene -> that gene is dropped
C = rng.poisson(0.03, size=(len(genes), n)).astype(int)
for j in range(0, n, 5): C[rng.choice(len(m.genes), 8, replace=False), j] = rng.integers(100, 600, 8)
# marker-driven cells: each target state gets cells carrying that state's strongest genes (so every call type occurs)
tgt = ["GMP-Mono", "Early GMP", "Pre-cDC", "HSC", "MLP", "GMP-Neut", "Early ProMono", "CD14 Mono", "cDC2", "Pre-pDC", "CD4 Naive"]
for k, st in enumerate(tgt):
    top = np.argsort(-m.coef[m.states.index(st)])[:40]
    for j in (14 + 2 * k, 15 + 2 * k): C[top, j] = rng.integers(3, 15, len(top))
C[:, 3] = 0; perm = rng.permutation(len(genes)); genes = [genes[i] for i in perm]; C = C[perm]
cells = [f"cell{i}" for i in range(n)]; qc = np.array([True] * n); qc[[10, 11, 12]] = False
pd.DataFrame(C, index=genes, columns=cells).to_csv("tests/e2e_counts.csv.gz")
pd.Series(qc, index=cells, name="qc").to_csv("tests/e2e_qc.csv")
exp = ps.score_counts(C.T, genes, cells, m, qc); exp.to_csv("tests/e2e_expected.csv", float_format="%.10g")
print(exp.progenitor_call.value_counts().to_dict())
