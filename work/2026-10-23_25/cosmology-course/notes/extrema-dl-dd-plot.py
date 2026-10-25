"""Figure: the four extremum angles and d(ell)/dd vs d/D (other inputs at baseline, alpha = 1).
Exact no-boxcar curve with tilt vs first-order card formula vs arccos form (minima). Not a fit."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from nobox_dl_dd_curves import formulas
plt.rcParams.update({"font.size": 14, "mathtext.fontset": "cm", "font.family": "serif"})
z = np.load("nobox_dl_dd_curves.npz"); es, ex = z["es"], z["ex"]
lin = z["lin"]; arc = z["arc"]
dex = np.gradient(ex, es, axis=0)
F = [formulas(e) for e in es]
dlin_min = np.array([[f[2][2], f[2][3]] for f in F]); darc = np.array([f[3] for f in F])
# first-order d(ell_max)/dd from the explicit formula (sympy forms)
import sympy as sp
from nobox_dl_dd_forms import forms
from nobox_dl_dd_sym import D, L, d, M, al, G, H, f, a, ph, full
from nobox_check import K0
sub = {D: 1.0, L: 1/5.9, M: K0, al: 1, G: 1, H: 1, f: 0.94, a: -0.06, ph: 0.35}
dmaxf = [sp.lambdify(d, full(forms(nm)[d]).subs(sub), "numpy") for nm in ("lmax1", "lmax2")]
dlin_max = np.array([[fn(e) for fn in dmaxf] for e in es])
def cut_pole(y, jump):
    y = y.copy(); bad = np.where(np.abs(np.diff(y)) > jump)[0]; y[bad + 1] = np.nan; return y
fig, axs = plt.subplots(4, 2, figsize=(12, 16), sharex=True)
titles = [r"$\ell_{\min,1}$", r"$\ell_{\min,2}$", r"$\ell_{\max,1}$", r"$\ell_{\max,2}$"]
rows = [(2, 0, 0), (3, 1, 1), (0, None, 0), (1, None, 1)]   # (exact col, arccos col, d-col)
for r, (k, ka, j) in enumerate(rows):
    A, Bx = axs[r, 0], axs[r, 1]
    A.plot(es, ex[:, k], "k-", lw=2.4, label="exact (no boxcar, with tilt)")
    A.plot(es, cut_pole(lin[:, k], 20), color="tab:blue", ls="--", lw=2, label="first-order card formula")
    if ka is not None:
        A.plot(es, arc[:, ka], color="tab:red", ls="-.", lw=2, label=r"arccos form: $\cos\ell_{\min}=-3d/4D$ + tilt term")
        Bx.plot(es, dex[:, k], "k-", lw=2.4); Bx.plot(es, dlin_min[:, j], color="tab:blue", ls="--", lw=2)
        Bx.plot(es, darc[:, j], color="tab:red", ls="-.", lw=2)
    else:
        Bx.plot(es, dex[:, k], "k-", lw=2.4); Bx.plot(es, cut_pole(dlin_max[:, j], 50), color="tab:blue", ls="--", lw=2)
    lo, hi = np.nanmin(ex[:, k]), np.nanmax(ex[:, k]); pad = 0.15*(hi - lo) + 1
    A.set_ylim(lo - pad, hi + pad)
    dl, dh = np.nanmin(dex[:, k]), np.nanmax(dex[:, k]); dp = 0.15*(dh - dl) + 2
    Bx.set_ylim(dl - dp, dh + dp)
    if k == 0:
        for ax in (A, Bx): ax.axvline(0.399, color="tab:blue", lw=1, ls=":")
        A.text(0.41, lo - 0.6*pad, "formula pole, $d/D=0.399$", color="tab:blue", fontsize=13)
    A.set_ylabel(titles[r] + "  (deg)"); Bx.set_ylabel(r"$d$" + titles[r] + r"$/d(d/D)$  (deg)")
    for ax in (A, Bx): ax.grid(alpha=0.3); ax.set_xlim(0, 0.9)
axs[3, 0].set_xlabel(r"$d/D$"); axs[3, 1].set_xlabel(r"$d/D$")
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=1, fontsize=15, frameon=False, bbox_to_anchor=(0.5, 0.995))
fig.text(0.5, 0.905, r"Other inputs at baseline: $D=1$, $L=D/5.9$, $M=0.2306$, $\alpha=1$; Amp$_A=-0.06$, $\phi_A=0.35$ rad, $f_u=0.94$",
         ha="center", fontsize=14)
fig.tight_layout(rect=(0, 0, 1, 0.895))
fig.savefig("extrema-dl-dd-plot.png", dpi=200)
print("saved")
