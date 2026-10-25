import sys, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "."); from integ_scan import *
sys.path.insert(0, ".."); import param_scan_core as PS
plt.rcParams.update({"font.size": 8.5, "font.family": "serif"})
res = json.load(open("integ_results.json")); P = json.load(open("integ_profiles.json"))
gb, gl, gs = res["mig"], res["mig_l0"], res["mig_shoot"]; hisM = res["his_geom_Mfit"]["mc"]
CUR = dict(DR=5.0, LD=0.18, mc=hisM, amp=-0.06, phi=0.35, l0=0.0)
lab = dict(DR="D/R: CoE distance / shell radius", LD="L/D: density scale / CoE distance", mc=r"M/M$_{\rm crit}$ (mass inside D)",
           amp="sine amplitude", phi="sine phase [rad]", l0=r"optional axis offset $\ell_0$ [deg]")
fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6)); axs = axs.ravel()
for ax, n in zip(axs, ["DR", "LD", "mc", "amp", "phi", "l0"]):
    d = P[n]; g = np.array(d["grid"]); pr = np.array(d["profile"], float); fx = np.array(d["fixed"], float)
    pr[~np.isfinite(pr)] = np.nan; fx[~np.isfinite(fx)] = np.nan
    ax.plot(g, pr, color="C3", lw=1.5, label="others re-optimised")
    ax.plot(g, fx, color="C0", lw=1.2, ls="--", label="others at your values\n(M re-fitted)")
    ax.axvline(CUR[n], color="0.4", lw=0.8, ls=":")
    ax.axvline(gl["l0"] if n == "l0" else gb[n], color="C3", lw=0.8, ls="-.")
    if n in ("LD", "mc"): ax.set_xscale("log")
    ax.set_ylim(0.0, 0.95); ax.set_xlabel(lab[n], fontsize=7.5); ax.grid(alpha=0.25)
    if n in ("DR", "amp"): ax.set_ylabel(r"$R^2$ vs 71 Migkas points")
axs[0].legend(fontsize=6.5, frameon=False, loc="lower right")
fig.text(0.5, 1.0, r"Integrated model, fixed profile, R = 14.3 Gpc.  Dotted grey: your value;  dash-dot red: best fit", ha="center", fontsize=8)
fig.tight_layout(); fig.savefig("../integrated_profiles.png", dpi=200, bbox_inches="tight"); plt.close(fig)

m = np.load("integ_map.npz")
fig, ax = plt.subplots(figsize=(4.6, 3.6))
Z = np.where(np.isfinite(m["Z"]), m["Z"], np.nan)
cs = ax.pcolormesh(m["DR"], m["LD"], np.clip(Z, 0, None), shading="auto", cmap="viridis", vmin=0, vmax=0.8)
c = ax.contour(m["DR"], m["LD"], Z, levels=[0.5, 0.6, 0.7, 0.74], colors="w", linewidths=0.7); ax.clabel(c, fontsize=6.5, fmt="%.2f")
ax.plot(CUR["DR"], CUR["LD"], "o", mfc="none", mec="r", ms=7, label="your values (5, 0.18)")
ax.plot(gb["DR"], gb["LD"], "*", color="r", ms=9, label="best (%.2f, %.3f)" % (gb["DR"], gb["LD"]))
ax.plot(gs["DR"], gs["LD"], "^", color="orange", ms=6, label="best, shot (%.2f, %.3f)" % (gs["DR"], gs["LD"]))
ax.set_yscale("log"); ax.set_xlabel("D/R: CoE distance / shell radius"); ax.set_ylabel("L/D: density scale length / CoE distance")
fig.colorbar(cs, ax=ax, label=r"$R^2$ (M, amp, phase re-optimised)"); ax.legend(fontsize=6.5, loc="upper right", frameon=True)
fig.tight_layout(); fig.savefig("../integrated_DR_LD.png", dpi=200, bbox_inches="tight"); plt.close(fig)

th = TH1; ib = json.load(open("../param_scan_results.json"))["global"]["mig"]
ip = ib["params"]
A_imp = PS.BH(th, ip["DR"], ip["LD"], ip["M"], ip["amp"], ip["phi"])
A_int = A_from(af_grid(gb["DR"], gb["LD"], gb["mc"]/MCRIT_FAC), th, gb["amp"], gb["phi"])
afl = shifted(af_grid(gl["DR"], gl["LD"], gl["mc"]/MCRIT_FAC), gl["l0"]); A_l = A_from(afl, th - gl["l0"], gl["amp"], gl["phi"])
A_s = A_from(af_grid(gs["DR"], gs["LD"], gs["mc"]/MCRIT_FAC, dict(CFG, shoot=6)), th, gs["amp"], gs["phi"])
Af = AFfour(th)
fig, (a, b) = plt.subplots(2, 1, figsize=(6.5, 5.0), sharex=True, gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.06})
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
a.plot(th, Af, color="C0", lw=1.5, label=r"$A_F$")
a.plot(th, A_imp, color="k", lw=1.1, ls="--", label=r"impulse best fit ($R^2$=%.3f)" % ib["R2"])
a.plot(th, A_int, color="C3", lw=1.6, label=r"integrated best fit ($R^2$=%.3f)" % gb["R2"])
a.plot(th, A_s, color="orange", lw=1.1, ls=":", label=r"integrated, shot to your geometry ($R^2$=%.3f)" % gs["R2"])
a.plot(th, A_l, color="C2", lw=1.2, ls="-.", label=r"integrated with $\ell_0$=%.1f$^\circ$ ($R^2$=%.3f)" % (gl["l0"], gl["R2"]))
a.set_ylabel(r"$A/A_{\rm all}$"); a.set_ylim(0.66, 1.22); a.legend(fontsize=6.8, loc="lower center", ncol=2, frameon=False)
b.axhline(0, color="k", lw=0.6)
b.plot(ML, MA - AFfour(ML), "o", ms=2.5, color="0.45", label=r"data$-A_F$")
b.plot(th, A_imp - Af, color="k", lw=1.0, ls="--"); b.plot(th, A_int - Af, color="C3", lw=1.4)
b.plot(th, A_s - Af, color="orange", lw=1.0, ls=":"); b.plot(th, A_l - Af, color="C2", lw=1.1, ls="-.")
b.set_ylabel(r"model$-A_F$"); b.set_xlabel(r"Galactic longitude $\ell$ (deg);  $\theta=\ell$ unless $\ell_0$ shown"); b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45))
b.legend(fontsize=7, frameon=False, loc="upper left")
fig.savefig("../integrated_bestfit.png", dpi=200, bbox_inches="tight")
from PIL import Image
for f in ["integrated_profiles", "integrated_DR_LD", "integrated_bestfit"]:
    im = Image.open("../" + f + ".png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../" + f + "_web.jpg", quality=88)
