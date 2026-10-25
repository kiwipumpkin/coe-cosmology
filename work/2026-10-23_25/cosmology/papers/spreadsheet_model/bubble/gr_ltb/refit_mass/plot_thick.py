"""Widened best fits: (1) Figure-2 style gr_fig2_R_migkas_wide.png; (2) loop_Agraph style gr_Agraph_wide.png
(same figsize, fonts, panel ratios, bands +-9 deg at alpha 0.22, trimming and web export as loops/loop_stack.py)."""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
import rm_core as C
from bubble_core import AFfour, ML, MA, MLO, MHI, TH1, IMIG, r2
W = json.load(open("thick.json"))
from thick import valid, OS
best = {}
for c in ("median", "zmax"):
    V = [o for o in W["grid"] + W["refined"] if o.get("case") == c and valid(o)]; best[c] = max(V, key=lambda o: o["R2_s"])
curves = {}; rows = {}
for c, o in best.items():
    C.ZL = (o["z"],); e = C.evaluate(o["Mc_kg"], o["D"], o["h"], Os=OS, full=True); f = e["fits"][str(o["z"])]
    Ad = np.array(f["Ad"]); curves[c] = C.curve(Ad, f["l0_s"], f["amp"], f["phi"])
    rows[c] = {k: v for k, v in e.items() if k != "fits"}; rows[c].update({k: v for k, v in f.items() if k != "Ad"}, z=o["z"], mu=o["mu"], h=o["h"])
json.dump(rows, open("best_rows_thick.json", "w"), indent=1, default=float)
Nn = np.load("newton_curve.npy")[1]; th2 = np.append(TH1, 360); cl = lambda y: np.append(y, y[0]); Af2 = AFfour(th2)
LAB = {c: "LTB GR, thick shell $\\Omega_s$=0.049, sampled $z$=%g: $D$=%.0f, $R_{\\rm in}$=%.0f, $R_{\\rm out}$=%.0f Gpc ($R^2$=%.3f)" % (best[c]["z"], best[c]["D"], best[c]["D"] - best[c]["h"], best[c]["D"] + best[c]["h"], rows[c]["R2_s"]) for c in best}
STY = {"median": dict(color="k", lw=1.4), "zmax": dict(color="#2ca02c", lw=1.1, ls="--")}
def trim(im, pad=12):
    g = np.asarray(im.convert("L")) < 250; ys, xs = np.nonzero(g)
    return im.crop((max(xs.min() - pad, 0), max(ys.min() - pad, 0), min(xs.max() + pad + 1, im.width), min(ys.max() + pad + 1, im.height)))
def draw(fname, bands):
    plt.rcdefaults(); plt.rcParams.update({"font.size": 9, "font.family": "serif"})
    fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, 5.25), sharex=True, gridspec_kw={"height_ratios": [2.0, 1], "hspace": 0.06})
    if bands:
        for L, col in ((29.5, "#0097b2"), (122.6, "#6fa800"), (187.5, "#e0107a"), (264.4, "#2166d4")):
            for ax_ in (a, b): ax_.axvspan(L - 9, L + 9, color=col, alpha=0.22, lw=0, zorder=0)
    a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
    a.plot(th2, Af2, color="C0", lw=1.5, label="$A_F$")
    a.plot(th2, cl(Nn), color="C3", lw=1.6, label="Newtonian thin shell + core, $R$=14.3 Gpc ($R^2$=%.3f)" % r2(MA, Nn[IMIG]))
    b.plot(th2, cl(100 * (Nn - AFfour(TH1)) / AFfour(TH1)), color="C3", lw=1.4, label="Newtonian")
    for c in ("median", "zmax"):
        a.plot(th2, cl(curves[c]), label=LAB[c], **STY[c]); b.plot(th2, cl(100 * (curves[c] - AFfour(TH1)) / AFfour(TH1)), label="LTB $z$=%g" % best[c]["z"], **STY[c])
    a.set_ylabel("$A/A_{\\rm all}$"); a.set_ylim(0.62, 1.22); a.legend(fontsize=7.2 if bands else 6.8, loc="lower left", frameon=False)
    b.axhline(0, color="k", lw=0.6); b.plot(ML, 100 * (MA - AFfour(ML)) / AFfour(ML), "o", ms=2.5, color="0.45", label="data")
    b.set_ylabel("$100\\,(\\cdot-A_F)/A_F$ [%]", fontsize=8); b.set_xlabel("Galactic longitude $\\ell$ (deg)")
    b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45)); b.legend(fontsize=7, frameon=False, loc="upper right", ncol=2)
    fig.savefig(fname, dpi=200, bbox_inches="tight", pad_inches=0.04); plt.close(fig)
    im = trim(Image.open(fname).convert("RGB")); im.save(fname); q = im.copy(); q.thumbnail((1700, 1700)); q.save(fname.replace(".png", "_web.jpg"), quality=88)
draw("gr_Agraph_omb_thick.png", True)
for c in rows: print(c, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[c].items()})
