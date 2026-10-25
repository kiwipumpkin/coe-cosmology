"""Extrema of A/A_all(l) over widened input ranges. Model imported unchanged
via extrema-vs-inputs.py; every local max/min found at each sweep value is
plotted (no branch tracking needed, so merges/new extrema show as they are)."""
import importlib.util, sys, warnings
import numpy as np, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
warnings.filterwarnings("ignore")
NOTES = Path("/workspace/cosmology-course/notes")
spec = importlib.util.spec_from_file_location("ex", NOTES / "extrema-vs-inputs.py")
ex = importlib.util.module_from_spec(spec); sys.modules["ex"] = ex; spec.loader.exec_module(ex)
plt.rcParams.update({"mathtext.fontset": "dejavusans", "font.size": 14})
N = 400
SW = [  # key, label, lo, hi, base, log, held overrides
    ("x", r"$D/L$", 1.01, 15.0, 5.9, False, {}),
    ("alpha", r"$\alpha$", 0.01, 8.0, 1.0, False, {}),
    ("kappa", r"$\kappa$", 0.001, 2.0, ex.lv.KAPPA0, False, {}),
    ("eps", r"$d/D$", 0.005, 0.999, 0.10, False, {"alpha": 1.3}),
]
def fam(kind, l):
    if kind == "max":
        return 0 if (l < 100 or l > 300) else 1
    return 0 if l < 180 else 1
COL = {("max", 0): "#d62728", ("max", 1): "#ff7f0e", ("min", 0): "#1f77b4", ("min", 1): "#2ca02c"}
data = []
for key, label, lo, hi, base, lg, over in SW:
    vals = np.geomspace(lo, hi, N) if lg else np.linspace(lo, hi, N)
    pts = []
    for v in vals:
        kw = dict(ex.HOLD); kw.update(over); kw[key] = float(v)
        for kind, l, A in ex.extrema(kw):
            pts.append((v, kind, l, A))
    data.append(pts)
    print(key, "extrema counts at ends:", sum(p[0] == vals[0] for p in pts), sum(p[0] == vals[-1] for p in pts))

fig, axs = plt.subplots(2, 4, figsize=(16, 8.6), sharex="col",
                        gridspec_kw=dict(hspace=0.08, wspace=0.28))
for c, ((key, label, lo, hi, base, lg, over), pts) in enumerate(zip(SW, data)):
    al, aa = axs[0, c], axs[1, c]
    for (kind, f), col in COL.items():
        sel = [p for p in pts if p[1] == kind and fam(kind, p[2]) == f]
        if not sel: continue
        v = np.array([p[0] for p in sel]); l = np.array([p[2] for p in sel]); A = np.array([p[3] for p in sel])
        mk = "o" if kind == "max" else "s"
        al.scatter(v, l, s=5, color=col, marker=mk, lw=0)
        aa.scatter(v, A, s=5, color=col, marker=mk, lw=0)
    for a in (al, aa):
        a.axvline(base, color="k", ls=":", lw=1.8)
        a.grid(alpha=0.35, which="both"); a.tick_params(labelsize=13)
        if lg: a.set_xscale("log")
        a.set_xlim(lo, hi)
    aa.axhline(1.0, color="#888888", lw=0.9)
    al.set_ylim(0, 360); al.set_yticks(np.arange(0, 361, 60))
    aa.set_ylim(0, 2.2)
    al.set_title(f"vary {label}" + (r" ($\alpha$=1.3 held)" if key == "eps" else ""), fontsize=15)
    aa.set_xlabel(label + (" (log scale)" if lg else ""), fontsize=15)
axs[0, 0].set_ylabel(r"$\ell$ of extremum [deg]", fontsize=15)
axs[1, 0].set_ylabel(r"$A/A_{\mathrm{all}}$ at extremum", fontsize=15)
hs = [Line2D([], [], color=COL[("max", 0)], marker="o", ls="none", ms=8, label="max near 10°"),
      Line2D([], [], color=COL[("max", 1)], marker="o", ls="none", ms=8, label="max near 170°"),
      Line2D([], [], color=COL[("min", 0)], marker="s", ls="none", ms=8, label="min near 95°"),
      Line2D([], [], color=COL[("min", 1)], marker="s", ls="none", ms=8, label="min near 266°"),
      Line2D([], [], color="k", lw=1.8, ls=":", label="baseline")]
fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=5,
           fontsize=14, frameon=False, columnspacing=1.4)
fig.subplots_adjust(top=0.87, bottom=0.10, left=0.07, right=0.99)
fig.savefig(NOTES / "extrema-vs-inputs-wide-linear.png", dpi=200)
print("ok")
