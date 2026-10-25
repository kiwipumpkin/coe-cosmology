import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
plt.rcParams.update({"font.size": 8.5, "font.family": "serif"})
res = json.load(open("bound_results.json")); P = res["prof"]
MCF = 1.8470
sty = {"mean|mig|ob=0": ("C3", "-", "vs Migkas"), "mean|mig_l0|ob=0": ("C3", "--", r"vs Migkas, with $\ell_0$"),
       "mean|four|ob=0": ("C0", "-", r"vs $A_F$"), "mean|four_l0|ob=0": ("C0", "--", r"vs $A_F$, with $\ell_0$"),
       "mean|mig|ob=0.3": ("C2", "-", r"vs Migkas, uniform background"), "mean|mig_l0|ob=0.3": ("C2", "--", r"vs Migkas, background, $\ell_0$"),
       "kR|mig|ob=0": ("k", ":", "vs Migkas, your $kR$ normalisation")}
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.9))
for k, (c, ls, lab) in sty.items():
    rows = [r for r in P[k] if r.get("R2") is not None and r["R2"] > (-1 if k.startswith("kR") else 0.5)]
    mnow = np.array([2 * r["mc"] / MCF / r["Ht"] ** 2 for r in rows]); R2 = np.array([r["R2"] for r in rows])
    v2 = np.array([r["vesc2"] for r in rows]); age = np.array([r["age_Gyr"] for r in rows]); b = v2 <= 1
    for ax, x in ((axs[0], mnow), (axs[1], v2)):
        ax.plot(x, R2, color=c, ls=ls, lw=1.2, label=lab)
        ax.plot(x[b], R2[b], "o", color=c, ms=3); ax.plot(x[~b], R2[~b], "o", mfc="none", color=c, ms=3)
    axs[2].plot(mnow, age, color=c, ls=ls, lw=1.2); axs[2].plot(mnow[b], age[b], "o", color=c, ms=3)
for ax in axs: ax.set_xscale("log"); ax.grid(alpha=0.25)
axs[0].set_xlabel(r"$M/M_{\rm crit}$ inside $D$ (today, $H_0=70$)"); axs[0].set_ylabel(r"best $R^2$")
axs[1].set_xlabel(r"Earth $(v/v_{\rm esc})^2$ today"); axs[1].axvline(1, color="0.3", lw=0.8)
axs[1].text(0.6, 0.05, "bound", ha="right", fontsize=7.5); axs[1].text(1.7, 0.05, "unbound", ha="left", fontsize=7.5)
for ax in axs[:2]: ax.set_ylim(0, 1.0)
axs[2].axhspan(0, 12, color="C3", alpha=0.08); axs[2].axhline(12, color="C3", lw=0.8)
axs[2].text(0.012, 11.6, "younger than 12 Gyr", color="C3", fontsize=7, va="top")
axs[2].set_xlabel(r"$M/M_{\rm crit}$ inside $D$"); axs[2].set_ylabel(r"age $t_0$ [Gyr] at $H_0=70$"); axs[2].set_ylim(9, 14.5)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=4, fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.1))
fig.tight_layout(); fig.savefig("../bound_fit_R2.png", dpi=200, bbox_inches="tight")
im = Image.open("../bound_fit_R2.png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../bound_fit_R2_web.jpg", quality=88)
