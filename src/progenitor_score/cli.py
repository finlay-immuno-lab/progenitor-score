import argparse, json, sys
from pathlib import Path
import anndata as ad
from .model import load_model, DEFAULT_MODEL
from .score import score_anndata, CALLS
from .cite import citation, SHORT

def main(argv=None):
    ap = argparse.ArgumentParser(prog="progenitor-score", description="Score cells against BoneMarrowMap and call progenitor types.")
    ap.add_argument("h5ad", nargs="?", help="AnnData with RAW counts in X"); ap.add_argument("--layer", default=None, help="use adata.layers[LAYER] instead of X"); ap.add_argument("--cite", action="store_true", help="print the references to cite and exit"); ap.add_argument("-o", "--out", default="progenitor_scores.csv.gz")
    ap.add_argument("--symbol-col", default=None, help="var column holding gene symbols (default: var_names)"); ap.add_argument("--qc-col", default=None, help="obs bool column; False cells are not scored")
    ap.add_argument("--model", default=str(DEFAULT_MODEL)); ap.add_argument("--min-margin", type=float, default=None); ap.add_argument("--min-top1", type=float, default=None)
    a = ap.parse_args(argv)
    if a.cite: print(citation()); return
    if not a.h5ad: ap.error("an input .h5ad is required (or --cite)")
    m = load_model(a.model); th = m.meta["default_thresholds"]
    mm = th["min_margin"] if a.min_margin is None else a.min_margin; mt = th["min_top1"] if a.min_top1 is None else a.min_top1
    sc = score_anndata(ad.read_h5ad(a.h5ad), m, a.symbol_col, a.qc_col, mm, mt, a.layer); sc.to_csv(a.out)
    print(SHORT, file=sys.stderr)
    print(json.dumps({c: int((sc.progenitor_call == c).sum()) for c in CALLS}))

if __name__ == "__main__":
    main()
