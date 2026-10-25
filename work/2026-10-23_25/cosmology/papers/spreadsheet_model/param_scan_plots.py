import sys, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "."); from param_scan_core import *
plt.rcParams.update({"font.size": 8.5, "font.family": "serif"})
P = json.load(open("param_scan_profiles.json")); G_ = json.load(open("param_scan_results.json"))
gb = G_["global"]["mig"]["params"]; gl = G_["global"]["mig_l0"]["params"]
lab = dict(DR="D/R: CoE distance / shell radius (C23)", LD="L/D: density scale / CoE distance (D23)", M="M: mass inside D (B27) [kg]",
           amp="sine amplitude (BA57)", phi="sine phase (BA59) [rad]", l0=r"optional axis offset $\ell_0$ [deg] ($\theta=\ell-\ell_0$)")
fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6)); axs = axs.ravel()
for ax, n in zip(axs, ["DR", "LD", "M", "amp", "phi", "l0"]):
    d = P[n]; g = np.array(d["grid"])
    ax.plot(g, d["profile"], color="C3", lw=1.5, label="others re-optimised")
    ax.plot(g, d["fixed"], color="C0", lw=1.2, ls="--", label="others at sheet values")
    ax.axvline(CUR[n], color="0.4", lw=0.8, ls=":")
    bv = gl["l0"] if n == "l0" else gb[n]
    ax.axvline(bv, color="C3", lw=0.8, ls="-.")
    if n in ("LD", "M"): ax.set_xscale("log")
    ax.set_ylim(0.0, 0.95); ax.set_xlabel(lab[n], fontsize=7.5); ax.grid(alpha=0.25)
    if n in ("DR", "amp"): ax.set_ylabel(r"$R^2$ vs 71 Migkas points")
axs[0].legend(fontsize=7, frameon=False, loc="lower right")
fig.text(0.5, 1.0, r"R = 14.3 Gpc.  Dotted grey: sheet value;  dash-dot red: best fit", ha="center", fontsize=8)
fig.tight_layout(); fig.savefig("param_scan_profiles.png", dpi=200, bbox_inches="tight"); plt.close(fig)

m = np.load("param_scan_map.npz")
fig, ax = plt.subplots(figsize=(4.6, 3.6))
Z = np.clip(m["Z"], 0, None)
cs = ax.pcolormesh(m["DR"], m["LD"], Z, shading="auto", cmap="viridis", vmin=0, vmax=0.8)
c = ax.contour(m["DR"], m["LD"], m["Z"], levels=[0.5, 0.6, 0.7, 0.74], colors="w", linewidths=0.7)
ax.clabel(c, fontsize=6.5, fmt="%.2f")
ax.plot(CUR["DR"], CUR["LD"], "o", mfc="none", mec="r", ms=7, label="sheet (5, 0.18)")
ax.plot(gb["DR"], gb["LD"], "*", color="r", ms=9, label="best (%.2f, %.2f)" % (gb["DR"], gb["LD"]))
ax.set_yscale("log"); ax.set_xlabel("D/R: CoE distance / shell radius (C23)"); ax.set_ylabel("L/D: density scale length / CoE distance (D23)")
fig.colorbar(cs, ax=ax, label=r"$R^2$ (M, amp, phase re-optimised)")
ax.legend(fontsize=7, loc="upper right", frameon=True)
fig.tight_layout(); fig.savefig("param_scan_DR_LD.png", dpi=200, bbox_inches="tight"); plt.close(fig)

th = np.arange(0, 361.0)
bh_c = BH(th, CUR["DR"], CUR["LD"], CUR["M"], CUR["amp"], CUR["phi"])
bh_b = BH(th, gb["DR"], gb["LD"], gb["M"], gb["amp"], gb["phi"])
bh_l = BH(th - gl["l0"], gl["DR"], gl["LD"], gl["M"], gl["amp"], gl["phi"])
Af = AFfour(th)
fig, (a, b) = plt.subplots(2, 1, figsize=(6.5, 5.0), sharex=True, gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.06})
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
a.plot(th, Af, color="C0", lw=1.5, label=r"$A_F$")
a.plot(th, bh_c, color="0.3", lw=1.1, ls="--", label=r"BH, sheet values ($R^2$=%.3f)" % G_["current"]["mig"])
a.plot(th, bh_b, color="C3", lw=1.6, label=r"BH, best fit ($R^2$=%.3f)" % G_["global"]["mig"]["R2"])
a.plot(th, bh_l, color="C2", lw=1.2, ls="-.", label=r"BH, best fit with $\ell_0$=%.1f$^\circ$ ($R^2$=%.3f)" % (gl["l0"], G_["global"]["mig_l0"]["R2"]))
a.set_ylabel(r"$A/A_{\rm all}$"); a.set_ylim(0.74, 1.22); a.legend(fontsize=7, loc="lower left", frameon=False)
b.axhline(0, color="k", lw=0.6)
b.plot(ML, MA - AFfour(ML), "o", ms=2.5, color="0.45", label=r"data$-A_F$")
b.plot(th, bh_c - Af, color="0.3", lw=1.0, ls="--"); b.plot(th, bh_b - Af, color="C3", lw=1.4); b.plot(th, bh_l - Af, color="C2", lw=1.1, ls="-.")
b.set_ylabel(r"model$-A_F$"); b.set_xlabel(r"$\ell$ (deg);  sheet $\theta=\ell$ unless $\ell_0$ shown"); b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45))
b.legend(fontsize=7, frameon=False, loc="upper left")
fig.savefig("param_scan_bestfit.png", dpi=200, bbox_inches="tight")

from PIL import Image
for f in ["param_scan_profiles", "param_scan_DR_LD", "param_scan_bestfit"]:
    im = Image.open(f + ".png").convert("RGB"); im.thumbnail((1700, 1700)); im.save(f + "_web.jpg", quality=88)
