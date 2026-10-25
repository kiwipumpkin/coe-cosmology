import sys, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
from bubble_core import *
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
key = sys.argv[1]; p = json.load(open("bubble_results.json"))[key]
Rg = RKM / GPC; Dg = p["D_Gpc"]; Lg = p["L_Gpc"]; l0 = p["l0"]; Rbg = 1.5 * (Dg + Rg)
rr = run(p["DR"], p["LD"], p["mc"], p["Ob"])
lpk = (90 - np.degrees(p["phi"]) + l0) % 360
E = np.array([Dg, 0.0])
def sky(l, r): t = np.radians(l - l0); return E + r * np.array([np.cos(t), np.sin(t)])
XL, YL = (-48, 72), (-47, 47)
fig, ax = plt.subplots(figsize=(7.2, 6.3))
xx, yy = np.meshgrid(np.linspace(*XL, 1000), np.linspace(*YL, 800))
rgd = np.hypot(xx, yy) / Dg
rel = rho(rgd, rr["mu"], p["LD"], rr["mub"]) / (3 * rr["mub"] / (4 * np.pi))     # total density / background density
im0 = ax.imshow(np.log10(rel), extent=XL + YL, origin="lower", cmap="magma_r", vmin=0, vmax=2, zorder=0)
cs = ax.contour(xx, yy, rel, levels=[1.01, 1.1], colors=["0.35", "0.2"], linewidths=0.7, linestyles=":")
ax.clabel(cs, fmt={1.01: "+1%", 1.1: "+10%"}, fontsize=7, manual=[(-20, -20), (-14, -14)])
ax.add_patch(Circle((0, 0), Lg, fill=False, color="C1", lw=1.3, zorder=3))
ax.add_patch(Circle((0, 0), Dg, fill=False, color="C0", lw=1.2, ls="--", zorder=2))
ax.add_patch(Circle(E, Rg, fill=False, ec="C2", lw=1.8, zorder=4))
ax.plot(0, 0, "k+", ms=10, mew=1.5, zorder=5); ax.text(-1.5, -3.2, "CoE", ha="right", va="top", fontsize=9, color="w")
ax.plot(*E, "o", color="C0", ms=6, zorder=6, mec="k"); ax.text(Dg - 0.8, 1.2, "Earth", ha="right", va="bottom", fontsize=9, zorder=7)
ax.annotate("", xy=(Dg, -21.5), xytext=(0, -21.5), arrowprops=dict(arrowstyle="<->", lw=0.9))
ax.text(Dg / 2, -23.0, r"$D$ = %.1f Gpc" % Dg, ha="center", va="top", fontsize=9)
ax.annotate(r"$r=L$ = %.2f Gpc" % Lg, xy=(Lg * np.cos(2.2), Lg * np.sin(2.2)), xytext=(-26, 10), fontsize=9, color="C1",
            arrowprops=dict(arrowstyle="-", color="C1", lw=0.8))
ax.text(-21, 42.5, r"circle of radius $D$ (CoE mass $M$ inside)", color="C0", fontsize=8.5, ha="center")
t = np.radians(125); q = E + Rg * np.array([np.cos(t), np.sin(t)])
ax.annotate("", xy=q, xytext=E, arrowprops=dict(arrowstyle="->", color="C2", lw=1.0))
ax.text(q[0] - 0.8, q[1] + 0.8, r"$R$ = %.1f Gpc" % Rg, color="C2", fontsize=9, ha="right", va="bottom")
ax.annotate("", xy=E + [9, 0], xytext=E, arrowprops=dict(arrowstyle="-|>", color="C3", lw=1.8), zorder=6)
ax.text(Dg + 5.5, -1.2, "Earth's\nvelocity", color="C3", fontsize=7.5, ha="center", va="top")
ax.annotate("", xy=E - [9, 0], xytext=E - [1, 0], arrowprops=dict(arrowstyle="-|>", color="k", lw=1.0), zorder=6)
ax.text(Dg - 5.5, -1.2, "to CoE", fontsize=7.5, ha="center", va="top")
for l in (0, 90, 180, 270):
    a, b = sky(l, Rg - 0.9), sky(l, Rg + 0.9); ax.plot([a[0], b[0]], [a[1], b[1]], color="0.25", lw=1.2)
    q = sky(l, Rg + 3.3); ax.text(*q, r"$\ell$=%d°" % l, ha="center", va="center", fontsize=8, color="0.25")
q = sky(l0, Rg + 7.4); ax.text(*q, "anti-CoE\n$\\ell$=%.1f°" % l0, ha="center", va="center", fontsize=8, color="C3")
q = sky(l0 + 180, Rg + 1.0); ax.text(q[0], q[1] - 0.8, "CoE\n$\\ell$=%.1f°" % ((l0 + 180) % 360), ha="right", va="top", fontsize=8)
a, b = sky(lpk, 0), sky(lpk, Rg + 5)
ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color="C4", lw=1.3, ls="--"))
q = sky(lpk, Rg + 5.3); ax.text(q[0] + 0.3, q[1], "sine-term peak\n$\\ell$=%.1f°" % lpk, ha="left", va="bottom", fontsize=8, color="C4")
ax.text(Dg - 2, -Rg - 11.5, "density contrast across bubble: %.2f%%\n(at $\\pm$3.4 Gpc: %.3f%%)" % (100 * p["contrast"], 100 * p["c34"]), ha="center", va="top", fontsize=8.5, color="C2",
        bbox=dict(fc="w", ec="C2", lw=0.6, alpha=0.85))
ax.text(71, -45.5, r"uniform background ($\Omega_b$=%.2f) fills a sphere of radius %.0f Gpc" % (p["Ob"], Rbg), ha="right", va="bottom", fontsize=7.5)
cb = fig.colorbar(im0, ax=ax, shrink=0.6, pad=0.02); cb.set_label(r"$\log_{10}$ (total density / background density)", fontsize=8)
ax.set_xlim(*XL); ax.set_ylim(*YL); ax.set_aspect("equal")
ax.set_xlabel("x [Gpc] (away from CoE through Earth)"); ax.set_ylabel("y [Gpc]")
fig.tight_layout(); fig.savefig("../uniform_bubble_drawing.png", dpi=200, bbox_inches="tight"); plt.close(fig)
im = Image.open("../uniform_bubble_drawing.png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../uniform_bubble_drawing_web.jpg", quality=88)
print("scale check D", np.hypot(*E), "R", np.hypot(*(sky(77, Rg) - E)), "L", Lg, "DR", np.hypot(*E) / Rg, "LD", Lg / Dg)
