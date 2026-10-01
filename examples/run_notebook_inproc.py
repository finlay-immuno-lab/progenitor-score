"""Execute a notebook cell by cell in ONE process and embed the outputs. For sandboxes that cannot open the sockets a Jupyter kernel needs.
In Jupyter just use Run All. usage (from the repo root): python examples/run_notebook_inproc.py examples/progenitor_score_demo.ipynb"""
import sys, ast, io, base64, contextlib, traceback, time, os, nbformat
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
path = os.path.abspath(sys.argv[1]); nb = nbformat.read(path, as_version=4); cur = []
os.chdir(os.path.dirname(os.path.abspath(path)))                       # notebooks run with their own folder as working directory, as in Jupyter
sys.path.insert(0, "..")
def show(*a, **k):
    for n in plt.get_fignums():
        f = plt.figure(n); b = io.BytesIO(); f.savefig(b, format="png", dpi=110, bbox_inches="tight")
        cur.append(nbformat.v4.new_output("display_data", data={"image/png": base64.b64encode(b.getvalue()).decode(), "text/plain": "<Figure>"})); plt.close(f)
plt.show = show; ns = {"__name__": "__main__"}; k = 0
for c in nb.cells:
    if c.cell_type != "code": continue
    k += 1; c.execution_count = k; cur.clear(); buf = io.StringIO(); t0 = time.time(); err = None; res = None
    try:
        tree = ast.parse(c.source); last = tree.body.pop() if tree.body and isinstance(tree.body[-1], ast.Expr) else None
        with contextlib.redirect_stdout(buf):
            exec(compile(tree, f"<cell {k}>", "exec"), ns)
            if last is not None: res = eval(compile(ast.Expression(last.value), f"<cell {k}>", "eval"), ns)
    except BaseException:
        err = traceback.format_exc()
    outs = []
    if buf.getvalue(): outs.append(nbformat.v4.new_output("stream", name="stdout", text=buf.getvalue()))
    outs += list(cur)
    if res is not None:
        d = {"text/plain": repr(res)}
        if hasattr(res, "_repr_html_"): d["text/html"] = res._repr_html_()
        outs.append(nbformat.v4.new_output("execute_result", data=d, execution_count=k))
    if err: outs.append(nbformat.v4.new_output("error", ename="Error", evalue=err.strip().splitlines()[-1], traceback=err.splitlines()))
    c.outputs = outs; print(f"cell {k} {'ERROR' if err else 'ok'} {time.time() - t0:.1f}s", flush=True)
    if err: print(err); nbformat.write(nb, path); sys.exit(1)
nb.metadata["execution_note"] = "Executed by examples/run_notebook_inproc.py (single process, outputs embedded)."
nbformat.write(nb, path); print("executed ->", path)
