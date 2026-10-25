"""loop_Agraph-style A(l) chart for the Newtonian Omega_s = 0.049 refit (universe_newt_omb.json, sampling on the 14.3 Gpc sphere)."""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
from bubble_core import AFfour, ML, MA, MLO, MHI, TH1, IMIG, r2
J = json.load(open("universe_newt_omb.json")); d = J["Rs14.3|sine"]; dg = J["Rs14.3|geom"]
A = np.array(d["A"]); Ag = np.array(dg["A"]); Aold = np.load("gr_ltb/refit_mass/newton_curve.npy")[1]
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
th2 = np.append(TH1, 360); cl = lambda y: np.append(y, y[0]); Af2 = AFfour(th2); Af = AFfour(TH1)
fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, 5.25), sharex=True, gridspec_kw={"height_ratios": [2.0, 1], "hspace": 0.06})
for L, col in ((29.5, "#0097b2"), (122.6, "#6fa800"), (187.5, "#e0107a"), (264.4, "#2166d4")):
    for ax_ in (a, b): ax_.axvspan(L - 9, L + 9, color=col, alpha=0.22, lw=0, zorder=0)
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
a.plot(th2, Af2, color="C0", lw=1.5, label="$A_F$")
a.plot(th2, cl(Aold), color="0.35", lw=0.9, ls="--", label="previous fit, $\\Omega_s$=1.01 ($R^2$=%.3f)" % r2(MA, Aold[IMIG]))
a.plot(th2, cl(A), color="C3", lw=1.6, label="Newtonian, $\\Omega_s$=0.049 forced: $D$=%.1f, $R_{\\rm in}$=%.1f, $R_{\\rm out}$=%.1f Gpc, $\\ell_0$=%.0f° ($R^2$=%.3f)" % (d["D"], d["R_in"], d["R_out"], d["l0"], d["R2"]))
b.plot(th2, cl(100 * (A - Af) / Af), color="C3", lw=1.4, label="model")
b.plot(th2, cl(100 * (Aold - Af) / Af), color="0.35", lw=0.9, ls="--", label="previous")
a.set_ylabel("$A/A_{\\rm all}$"); a.set_ylim(0.62, 1.22); a.legend(fontsize=7.0, loc="lower left", frameon=False)
b.axhline(0, color="k", lw=0.6); b.plot(ML, 100 * (MA - AFfour(ML)) / AFfour(ML), "o", ms=2.5, color="0.45", label="data")
b.set_ylabel("$100\\,(\\cdot-A_F)/A_F$ [%]", fontsize=8); b.set_xlabel("Galactic longitude $\\ell$ (deg)")
b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45)); b.legend(fontsize=7, frameon=False, loc="upper right", ncol=3)
fig.savefig("newt_omb_Agraph.png", dpi=200, bbox_inches="tight", pad_inches=0.04); plt.close(fig)
im = Image.open("newt_omb_Agraph.png").convert("RGB"); g = np.asarray(im.convert("L")) < 250; ys, xs = np.nonzero(g)
im = im.crop((max(xs.min() - 12, 0), max(ys.min() - 12, 0), min(xs.max() + 13, im.width), min(ys.max() + 13, im.height))); im.save("newt_omb_Agraph.png")
q = im.copy(); q.thumbnail((1700, 1700)); q.save("newt_omb_Agraph_web.jpg", quality=88)
