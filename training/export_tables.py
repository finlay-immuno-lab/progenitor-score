"""Export a CellTypist pickle to the plain tables in src/progenitor_score/models/<id>/. usage: python training/export_tables.py model.pkl <model_id>"""
import json, hashlib, sys, numpy as np, pandas as pd, celltypist
from pathlib import Path
sys.path.insert(0, "src"); from progenitor_score.groups import group_of, DECISION
pkl, mid = sys.argv[1], sys.argv[2]; m = celltypist.models.Model.load(pkl)
D = Path("src/progenitor_score/models") / mid; D.mkdir(parents=True, exist_ok=True); fmt = "%.17g"
genes = np.asarray(m.features, dtype=str); states = [str(s) for s in m.classifier.classes_]
pd.DataFrame(m.classifier.coef_.T, index=genes, columns=states).rename_axis("gene").to_csv(D / "coef.csv.gz", float_format=fmt)
pd.DataFrame({"gene": genes, "mean": m.scaler.mean_, "scale": m.scaler.scale_}).to_csv(D / "scaler.csv", index=False, float_format=fmt)
pd.DataFrame({"state": states, "intercept": m.classifier.intercept_, "group": [group_of(s) for s in states], "decision": [DECISION[group_of(s)] for s in states]}).to_csv(D / "states.csv", index=False, float_format=fmt)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
meta = {"model_id": mid, "rules_version": "progenitor_rules_v1", "n_states": len(states), "n_genes": len(genes), "training": m.description,
        "celltypist_version": celltypist.__version__, "source_pickle_sha256": sha(pkl), "default_thresholds": {"min_margin": 100.0, "min_top1": 0.0},
        "attribution": "Model weights derived from the BoneMarrowMap atlas: Zeng et al. Blood Cancer Discov 2025;6:307-24, doi:10.1158/2643-3230.BCD-24-0342. Trained by the Conor Finlay Lab, Trinity College Dublin.", "licence": "CC BY-NC 4.0 (see LICENSE-MODEL.md)",
        "files": {f: sha(D / f) for f in ["coef.csv.gz", "scaler.csv", "states.csv"]}}
json.dump(meta, open(D / "model.json", "w"), indent=1)
