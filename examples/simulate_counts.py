"""Simulate a small 10x-style RAW count matrix with KNOWN labels, to demonstrate and sanity-check progenitor-score.

How it is built
  * Each population draws counts from the pseudobulk profile of one BoneMarrowMap state (examples/data/sim_state_profiles.csv.gz, parts per
    million over the model genes), blended 80:20 with a neighbouring state for the transitional populations, plus 5% ambient RNA.
  * Cell depth is log-normal (median ~1,500 UMIs, a 20-fold spread); counts are gamma-Poisson (per-gene gamma noise), so cells of one
    population differ from each other. About 25 cells fall below 500 UMIs and are flagged qc_pass = False.
  * 'Not in reference' cells use the CD14 Mono profile with the gene labels shuffled: a cell type the atlas does not contain.
What it shows and what it does not
  * It checks the plumbing (input formats, gene matching, QC flag, output columns, thresholds, R == Python) and how calls behave for clean,
    transitional, out-of-reference and low-depth cells.
  * It is NOT an accuracy estimate: the profiles are averages of the atlas the model was trained on, and real cells are noisier and more
    varied. Accuracy on 9 held-out real donors is in docs/validation/.
Run from the repo root:  python examples/simulate_counts.py
"""
import gzip, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.io import mmwrite
sys.path.insert(0, "src"); import progenitor_score as ps

OUT = Path("examples/data"); OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(2026); m = ps.load_model()
prof = pd.read_csv(OUT / "sim_state_profiles.csv.gz", index_col=0)                    # genes x states, ppm
n_extra = 500
genes = list(prof.index) + [f"SIMGENE{i:04d}" for i in range(n_extra)]; G = len(genes)
def vec(state):                                                                         # profile over all simulated genes, sums to 1
    v = np.concatenate([prof[state].to_numpy(float), np.full(n_extra, prof[state].mean() * 0.2)]); return v / v.sum()
allstates = np.mean([vec(s) for s in prof.columns], axis=0)

# population: (label, own state, blended-in neighbour or None, n cells)
POP = [("HSC", "HSC", "MPP-MyLy", 50), ("Early GMP", "Early GMP", "HSC", 70), ("GMP-Mono", "GMP-Mono", "Early ProMono", 70),
       ("Pre-cDC", "Pre-cDC", "cDC2", 50), ("Pre-pDC", "Pre-pDC", "pDC", 40), ("Early ProMono", "Early ProMono", "Late ProMono", 70),
       ("Late ProMono", "Late ProMono", "CD14 Mono", 70), ("CD14 Mono", "CD14 Mono", None, 100), ("cDC2", "cDC2", None, 50),
       ("CD4 Naive", "CD4 Naive", None, 60), ("Erythroblast", "Orthochromatic Erythroblast", None, 40), ("Not in reference", None, None, 60)]
rows, cols, vals, meta = [], [], [], []
c = 0
for label, st, nb, n in POP:
    if st is None: mu = rng.permutation(vec("CD14 Mono"))                               # same expression levels, scrambled gene identity
    else: mu = vec(st) * (0.8 if nb else 1.0) + (vec(nb) * 0.2 if nb else 0.0)
    for _ in range(n):
        depth = max(int(rng.lognormal(np.log(1500), 0.6)), 120)
        lam = depth * (0.95 * mu + 0.05 * allstates) * rng.gamma(2.0, 0.5, G)           # ambient 5%, gamma-Poisson noise
        x = rng.poisson(lam); nz = np.nonzero(x)[0]
        rows += list(nz); cols += [c] * len(nz); vals += list(x[nz])
        meta.append({"cell": f"cell{c:04d}", "true_population": label, "true_state": st or "none", "true_group": ps.group_of(st) if st else "not_in_reference", "umis": int(x.sum())})
        c += 1
M = sp.csr_matrix((vals, (rows, cols)), shape=(G, c), dtype=np.int32)
meta = pd.DataFrame(meta); meta["qc_pass"] = meta.umis >= 500
mmwrite(str(OUT / "sim_counts.mtx"), M)
with open(OUT / "sim_counts.mtx", "rb") as f, gzip.open(OUT / "sim_counts.mtx.gz", "wb", compresslevel=9) as g: g.write(f.read())
(OUT / "sim_counts.mtx").unlink()
pd.Series(genes).to_csv(OUT / "sim_genes.tsv", index=False, header=False)
meta.to_csv(OUT / "sim_cells.csv", index=False)
print(f"{G} genes x {c} cells, {M.nnz} non-zeros; qc fail {int((~meta.qc_pass).sum())}; median UMIs {int(meta.umis.median())}")
print(meta.groupby("true_population", sort=False).umis.median().astype(int).to_dict())
