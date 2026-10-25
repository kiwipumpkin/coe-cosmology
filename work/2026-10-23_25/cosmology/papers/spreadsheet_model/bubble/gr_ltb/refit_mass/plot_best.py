"""Figure 2-style plot: Migkas points, A_F, Newtonian red curve (0.909), LTB best fit (R = median Migkas distance), residual panel."""
import json, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
import rm_core as C
from bubble_core import AFfour, ML, MA, MLO, MHI, TH1, IMIG, r2
plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})
res = json.load(open("rgrid.json")); rows = {}
for case in ("median", "zmax"):
    V = [o for o in res if o["case"] == case and "R2_s" in o and o["twoE_min"] > -1]
    rows[case] = max(V, key=lambda o: o["R2_s"])
curves = {}
for case, o in rows.items():
    C.ZL = (o["z"],); e = C.evaluate(o["Mc_kg"], o["D"], o["R"], Os=o["Os"], full=True); f = e["fits"][str(o["z"])]
    Ad = np.array(f["Ad"]); curves[case] = (C.curve(Ad, f["l0_s"], f["amp"], f["phi"]), C.curve(Ad, f["l0_ns"]), e)
Nn = np.load("newton_curve.npy")[1]
Af = AFfour(TH1); th2 = np.append(TH1, 360); cl = lambda y: np.append(y, y[0])
fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, 4.6), sharex=True, gridspec_kw={"height_ratios": [2.0, 1], "hspace": 0.06})
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
a.plot(th2, cl(Af), color="C0", lw=1.0, label="$A_F$")
a.plot(th2, cl(Nn), color="#d62728", lw=1.5, label="Newtonian best fit, $R=14.3$ Gpc (fit $R^2$=%.3f)" % r2(MA, Nn[IMIG]))
Ag = curves["median"][0]; o = rows["median"]
a.plot(th2, cl(Ag), color="k", lw=1.3, label="LTB GR, $R=0.32$ Gpc ($z=0.075$), with sine ($R^2$=%.3f)" % r2(MA, Ag[IMIG]))
Ag2 = curves["zmax"][0]
a.plot(th2, cl(Ag2), color="#2ca02c", lw=1.0, ls="--", label="LTB GR, $R=1.19$ Gpc ($z=0.3$), with sine ($R^2$=%.3f)" % r2(MA, Ag2[IMIG]))
a.set_ylabel("$A/A_{\\rm all}$"); a.set_ylim(0.62, 1.22); a.legend(fontsize=6.8, loc="lower left", frameon=False)
b.axhline(0, color="k", lw=0.6); b.plot(ML, 100 * (MA - AFfour(ML)) / AFfour(ML), "o", ms=2.5, color="0.45", label="data")
for y, c, ls, lab in ((Nn, "#d62728", "-", "Newtonian"), (Ag, "k", "-", "LTB $z=0.075$"), (Ag2, "#2ca02c", "--", "LTB $z=0.3$")):
    b.plot(th2, cl(100 * (y - Af) / Af), color=c, ls=ls, lw=1.1, label=lab)
b.set_ylabel("$100\\,(\\cdot-A_F)/A_F$ [%]", fontsize=8); b.set_xlabel("Galactic longitude $\\ell$ (deg)"); b.set_xticks(range(0, 361, 45)); b.set_xlim(0, 360)
b.legend(fontsize=6.5, ncol=4, loc="lower center", frameon=False)
fig.savefig("gr_fig2_R_migkas.png", dpi=200, bbox_inches="tight")
im = Image.open("gr_fig2_R_migkas.png").convert("RGB"); w, h = im.size; s = min(1.0, 1400 / w)
im.resize((int(w * s), int(h * s)), Image.LANCZOS).save("gr_fig2_R_migkas_web.jpg", quality=88)
out = {}
for case, (As, Ans, e) in curves.items():
    f = e["fits"][str(e["fits"] and rows[case]["z"])]
    out[case] = {k: v for k, v in e.items() if k != "fits"}; out[case].update({kk: vv for kk, vv in f.items() if kk != "Ad"}, z=rows[case]["z"], mu=rows[case]["mu"])
json.dump(out, open("best_rows.json", "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
