"""Figure-2-top-panel style plot: Migkas points, old thin-shell fit with sine (red), refit without sine."""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
import thin_nosine as N, analysis_curve as AC
from bubble_core import ML, MA, MLO, MHI, TH1, IMIG, r2
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
plt.rcParams.update({"mathtext.fontset": "dejavusans", "axes.linewidth": 0.8, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.labelsize": 9, "axes.titlesize": 12})
J = json.load(open("thin_nosine.json")); old = json.load(open("thin_b_free_hug.json")); new = J["refit"]
A_old = AC.curve("ptshell", old); A_new = N.model_nosine(new)
R_old, R_new = r2(MA, A_old[IMIG]), r2(MA, A_new[IMIG]); assert abs(R_old - old["R2"]) < 1e-6 and abs(R_new - new["R2"]) < 1e-6
th2 = np.append(TH1, 360); c2 = lambda A: np.append(A, A[0])
fig, a = plt.subplots(1, 1, figsize=(6.6, 5.25 * 2.0 / 3.06 + 0.45))
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
a.plot(th2, c2(A_old), color="C3", ls="-", lw=1.6, label="thin shell + core, with sine ($R^2$=%.3f)" % R_old)
a.plot(th2, c2(A_new), color="C0", ls="-", lw=1.6, label="thin shell + core, no sine, refit ($R^2$=%.3f)" % R_new)
a.set_ylabel("$A/A_{\\rm all}$"); a.set_ylim(0.62, 1.22); a.legend(fontsize=7.2, loc="lower left", frameon=False)
a.set_xlabel("Galactic longitude $\\ell$ (deg)"); a.set_xlim(0, 360); a.set_xticks(range(0, 361, 45))
print("new curve range", A_new.min(), A_new.max())
fig.savefig("../thin_nosine.png", dpi=200, bbox_inches="tight")
im = Image.open("../thin_nosine.png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../thin_nosine_web.jpg", quality=88); print(im.size)
