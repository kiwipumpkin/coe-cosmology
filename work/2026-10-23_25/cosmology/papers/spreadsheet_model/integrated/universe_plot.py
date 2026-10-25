import sys, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
sys.path.insert(0, "."); from integ_scan import *
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
p = json.load(open("integ_results.json"))["mig_l0"]
Rg = RKM / GPC; Dg = p["DR"] * Rg; Lg = p["LD"] * Dg; l0 = p["l0"]
Dkm, Lkm = p["DR"] * RKM, p["LD"] * p["DR"] * RKM
M = p["mc"] / MCRIT_FAC * (Dkm * 1e3) ** 3 * K ** 2 / GSI
lpk = (90 - np.degrees(p["phi"]) + l0) % 360          # max of 1 - Amp sin(theta+phi), Amp<0
E = np.array([Dg, 0.0])
def sky(l, r): t = np.radians(l - l0); return E + r * np.array([np.cos(t), np.sin(t)])
fig, ax = plt.subplots(figsize=(7.2, 6.3))
# density shading e^{-r/L}
xx, yy = np.meshgrid(np.linspace(-48, 72, 980), np.linspace(-47, 47, 770))
rho = np.exp(-np.hypot(xx, yy) / Lg)
im0 = ax.imshow(np.log10(np.maximum(rho, 1e-12)), extent=(-48, 72, -47, 47), origin="lower", cmap="Greys", vmin=-8, vmax=0, alpha=0.9, zorder=0)
ax.add_patch(Circle((0, 0), Lg, fill=False, color="C1", lw=1.3, zorder=3))
ax.add_patch(Circle((0, 0), Dg, fill=False, color="C0", lw=1.2, ls="--", zorder=2))
ax.add_patch(Circle(E, Rg, fill=True, fc="C2", alpha=0.12, ec="C2", lw=1.5, zorder=2))
ax.plot(0, 0, "k+", ms=10, mew=1.5, zorder=5); ax.text(-1.5, -3.2, "CoE", ha="right", va="top", fontsize=9)
ax.plot(*E, "o", color="C0", ms=6, zorder=6, mec="k"); ax.text(Dg - 0.8, 1.2, "Earth", ha="right", va="bottom", fontsize=9, zorder=7)
# D and L labels
ax.annotate("", xy=(Dg, -20.5), xytext=(0, -20.5), arrowprops=dict(arrowstyle="<->", lw=0.9))
ax.text(Dg / 2, -22.0, r"$D$ = %.1f Gpc" % Dg, ha="center", va="top", fontsize=9)
ax.plot([0, 0], [-19, -22], "k", lw=0.6); ax.plot([Dg, Dg], [-19, -22], "k", lw=0.6)
ax.annotate(r"$r=L$ = %.2f Gpc" % Lg, xy=(Lg * np.cos(2.2), Lg * np.sin(2.2)), xytext=(-24, 9), fontsize=9, color="C1",
            arrowprops=dict(arrowstyle="-", color="C1", lw=0.8))
ax.text(-21, 42.5, r"circle of radius $D$ (mass $M$ inside)", color="C0", fontsize=8.5, ha="center")
# R radius, drawn toward l = 270 side (down-ish)
t = np.radians(125); ax.annotate("", xy=E + Rg * np.array([np.cos(t), np.sin(t)]), xytext=E, arrowprops=dict(arrowstyle="->", color="C2", lw=1.0))
q = E + Rg * np.array([np.cos(t), np.sin(t)]); ax.text(q[0] - 0.8, q[1] + 0.8, r"$R$ = %.1f Gpc" % Rg, color="C2", fontsize=9, ha="right", va="bottom")
# Earth velocity (radially away from CoE) and direction to CoE
ax.annotate("", xy=E + [9, 0], xytext=E, arrowprops=dict(arrowstyle="-|>", color="C3", lw=1.8), zorder=6)
ax.text(Dg + 5.5, -1.2, "Earth's\nvelocity", color="C3", fontsize=7.5, ha="center", va="top")
ax.annotate("", xy=E - [9, 0], xytext=E - [1, 0], arrowprops=dict(arrowstyle="-|>", color="k", lw=1.0, ls="-"), zorder=6)
ax.text(Dg - 5.5, -1.2, "to CoE", fontsize=7.5, ha="center", va="top")
# sky directions in Galactic longitude on the shell
for l, txt, c in [(0, r"$\ell$=0°", "0.25"), (90, r"$\ell$=90°", "0.25"), (180, r"$\ell$=180°", "0.25"), (270, r"$\ell$=270°", "0.25")]:
    a, b = sky(l, Rg - 0.9), sky(l, Rg + 0.9); ax.plot([a[0], b[0]], [a[1], b[1]], color=c, lw=1.2)
    q = sky(l, Rg + 3.3); ax.text(*q, txt, ha="center", va="center", fontsize=8, color=c)
q = sky(l0, Rg + 7.4); ax.text(*q, "anti-CoE\n$\\ell$=%.1f°" % l0, ha="center", va="center", fontsize=8, color="C3")
q = sky(l0 + 180, Rg + 1.0); ax.text(q[0], q[1] - 0.8, "CoE\n$\\ell$=%.1f°" % ((l0 + 180) % 360), ha="right", va="top", fontsize=8, color="w")
a, b = sky(lpk, 0), sky(lpk, Rg + 5)
ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color="C4", lw=1.3, ls="--"))
q = sky(lpk, Rg + 5.3); ax.text(q[0] + 0.3, q[1], "sine-term peak\n$\\ell$=%.1f°" % lpk, ha="left", va="bottom", fontsize=8, color="C4")
cb = fig.colorbar(im0, ax=ax, shrink=0.6, pad=0.02); cb.set_label(r"$\log_{10}$ density / central density ($e^{-r/L}$)", fontsize=8)
ax.set_xlim(-48, 72); ax.set_ylim(-47, 47); ax.set_aspect("equal")
ax.set_xlabel("x [Gpc] (away from CoE through Earth)"); ax.set_ylabel("y [Gpc]")
fig.tight_layout(); fig.savefig("../universe_drawing.png", dpi=200, bbox_inches="tight"); plt.close(fig)
# scale check: data-coordinate distances used above
chk = dict(D=float(np.hypot(*E)), R=float(np.hypot(*(sky(123.0, Rg) - E))), L=Lg, D_over_R=float(np.hypot(*E) / Rg), L_over_D=Lg / Dg,
           Dkm=Dkm, Lkm=Lkm, M=M, lpk=lpk, antiCoE=l0, CoE=(l0 + 180) % 360, rcore_min=(Dg - Rg) / Lg)
print(chk); json.dump(chk, open("universe_check.json", "w"), indent=1)

# A vs angle chart
th = TH1; afl = shifted(af_grid(p["DR"], p["LD"], p["mc"] / MCRIT_FAC), l0); A_l = A_from(afl, th - l0, p["amp"], p["phi"])
Af = AFfour(th); th2 = np.append(th, 360); A2 = np.append(A_l, A_l[0]); Af2 = AFfour(th2)
fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, 4.3), sharex=True, gridspec_kw={"height_ratios": [2.0, 1], "hspace": 0.06})
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
a.plot(th2, Af2, color="C0", lw=1.5, label=r"$A_F$")
a.plot(th2, A2, color="C3", lw=1.6, label=r"integrated model A, $\ell_0$=%.1f° ($R^2$=%.3f)" % (l0, p["R2"]))
a.set_ylabel(r"$A/A_{\rm all}$"); a.set_ylim(0.74, 1.22); a.legend(fontsize=7.5, loc="lower left", frameon=False)
b.axhline(0, color="k", lw=0.6)
b.plot(ML, 100 * (MA - AFfour(ML)) / AFfour(ML), "o", ms=2.5, color="0.45", label=r"data")
b.plot(th2, 100 * (A2 - Af2) / Af2, color="C3", lw=1.4, label="model")
b.set_ylabel(r"$100\,(\cdot-A_F)/A_F$ [%]", fontsize=8); b.set_xlabel(r"Galactic longitude $\ell$ (deg)")
b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45)); b.legend(fontsize=7, frameon=False, loc="upper right", ncol=2)
fig.savefig("../universe_Achart.png", dpi=200, bbox_inches="tight")
from PIL import Image
for f in ["universe_drawing", "universe_Achart"]:
    im = Image.open("../" + f + ".png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../" + f + "_web.jpg", quality=88)
