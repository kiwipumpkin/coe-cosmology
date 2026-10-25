import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
from PIL import Image
from bubble_core import AFfour, ML, MA, MLO, MHI, TH1, RKM, GPC, IMIG, r2
import analysis_curve as AC
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
L = lambda k: json.load(open("thin_%s.json" % k.replace("|", "_")))
d = L("b|free|hug"); Rg = RKM / GPC; Dg = d["D_Gpc"]; l0 = d["l0"]
fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax.add_patch(Wedge((0, 0), Dg + Rg, -40, 40, width=2 * Rg, fc="0.78", ec="0.35", lw=1.0, zorder=1))
ax.add_patch(Circle((Dg, 0), Rg, fc="C2", alpha=0.25, ec="C2", lw=1.6, zorder=3))
ax.plot(0, 0, "o", color="k", ms=5, zorder=5); ax.text(0, -4, "CoE\n(core $M_c$ = %.2g M$_\\odot$)" % d["Mc_sun"], ha="center", va="top", fontsize=8)
ax.plot(Dg, 0, "o", color="C0", mec="k", ms=6, zorder=6); ax.text(Dg, 2, "Earth", ha="center", va="bottom", fontsize=8.5, zorder=7)
ax.annotate("", xy=(Dg, -44), xytext=(0, -44), arrowprops=dict(arrowstyle="<->", lw=0.9)); ax.text(Dg / 2, -46, "$D$ = %.1f Gpc" % Dg, ha="center", va="top")
t = np.radians(60); ax.annotate("", xy=(Dg + Rg * np.cos(t), Rg * np.sin(t)), xytext=(Dg, 0), arrowprops=dict(arrowstyle="->", color="C2", lw=1.1))
ax.text(Dg + 8, 15, "$R$ = %.1f Gpc" % Rg, color="C2", fontsize=8.5)
ax.text(Dg - Rg - 6, 38, "shell, constant density\n$\\Omega_s$ = %.2f, from $D-R$ to $D+R$" % d["Os"], ha="right", fontsize=8, color="0.25")
ax.text(Dg * 0.42, 14, "empty\n(hole, mass $-\\frac{4\\pi}{3}\\rho_s R_{\\rm in}^3$ as seen\nfrom inside the shell)", ha="center", fontsize=8, color="0.3")
ax.text(Dg - Rg, -3, "$D-R$\n%.1f" % (Dg - Rg), ha="right", va="top", fontsize=7.5); ax.text(Dg + Rg + 1, -3, "$D+R$\n%.1f" % (Dg + Rg), ha="left", va="top", fontsize=7.5)
ax.set_xlim(-12, Dg + Rg + 14); ax.set_ylim(-52, 52); ax.set_aspect("equal")
ax.set_xlabel("x [Gpc]"); ax.set_ylabel("y [Gpc]")
ax.set_title("Thin shell hugging the visible sphere, best fit (core + shell, $\\Omega_s$ free), to scale", fontsize=9)
fig.savefig("../thin_shell_drawing.png", dpi=200, bbox_inches="tight"); plt.close(fig)
fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, 4.8), sharex=True, gridspec_kw={"height_ratios": [2.0, 1], "hspace": 0.06})
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.78", lw=0.8, label="Migkas Fig. 6")
Af = AFfour(TH1); a.plot(TH1, Af, color="C0", lw=1.5, label="$A_F$")
for k, lab, c, ls in [("b|free|hug", "core + shell, $\\Omega_s$=%.2f", "C3", "-"), ("b|0.05|hug", "core + shell, $\\Omega_s$=%.2f (baryons)", "C1", "--"),
                      ("a|free|hug", "shell only, $\\Omega_s$=%.2f", "C4", "-."), ("a|0.05|hug", "shell only, $\\Omega_s$=%.2f (baryons)", "C7", ":")]:
    dd = L(k); A = AC.curve("ptshell", dd)
    chk = r2(MA, A[IMIG]); assert abs(chk - dd["R2_mig"]) < 1e-6, (k, chk)
    lab2 = (lab % dd["Os"]) + ", $\\ell_0$=%.0f° ($R^2$=%.3f)" % (dd["l0"], dd["R2_mig"])
    a.plot(TH1, A, color=c, ls=ls, lw=1.4, label=lab2); b.plot(TH1, 100 * (A - Af) / Af, color=c, ls=ls, lw=1.2)
a.set_ylabel("$A/A_{\\rm all}$"); a.set_ylim(0.62, 1.22); a.legend(fontsize=6.8, loc="lower left", frameon=False, ncol=2)
b.axhline(0, color="k", lw=0.6); b.plot(ML, 100 * (MA - AFfour(ML)) / AFfour(ML), "o", ms=2.5, color="0.45")
b.set_ylabel("$100\\,(\\cdot-A_F)/A_F$ [%]", fontsize=8); b.set_xlabel("Galactic longitude $\\ell$ (deg)"); b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45))
fig.savefig("../thin_shell_Achart.png", dpi=200, bbox_inches="tight")
for f in ["thin_shell_drawing", "thin_shell_Achart"]:
    im = Image.open("../" + f + ".png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../" + f + "_web.jpg", quality=88)
print("ok")
