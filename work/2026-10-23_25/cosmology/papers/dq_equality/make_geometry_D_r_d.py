"""Geometry for d < r < D: CoE beyond the vision shell, sources inside at depth d.
(a) schematic (not to scale); (b) true scale in Gpc with a zoom inset for d."""
import numpy as np, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc, FancyArrowPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

mpl.rcParams.update({"font.size": 10, "mathtext.fontset": "cm", "savefig.dpi": 220})
C_COE, C_OBS, C_SRC, C_SH, C_AX = "#D55E00", "#0072B2", "#009E73", "#555555", "#CC79A7"
d_true = 2.99792458e5 * 0.075 / 70.0 / 1e3   # 0.321 Gpc
r_true = 14.3                                 # Gpc
D_ill = 20.0                                  # Gpc, illustrative only (D > r; not fixed by the fit)

def darrow(ax, p, q, color, lw=1.4, style="<->"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=11, color=color, lw=lw,
                                 shrinkA=0, shrinkB=0, zorder=4))

fig = plt.figure(figsize=(13.0, 6.0))
gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.08)

# ---------------- (a) schematic
ax = fig.add_subplot(gs[0])
r, D, d, th = 1.0, 1.75, 0.42, np.radians(60)
O, coe = np.array([0.0, 0.0]), np.array([-D, 0.0])
S = d * np.array([np.cos(th), np.sin(th)])
RS = np.linalg.norm(S - coe)
ax.add_patch(Circle(O, r, fill=True, fc="#f4f4f4", ec=C_SH, lw=1.8, zorder=0))
ax.text(0.0, -1.04, r"vision shell (limit of vision), radius $r\approx14.3$ Gpc, centred on O",
        color=C_SH, ha="center", va="top", fontsize=9.5)
# r
ang = np.radians(-50)
ax.plot([0, r * np.cos(ang)], [0, r * np.sin(ang)], color=C_SH, lw=1.1)
ax.text(0.5 * np.cos(ang) + 0.07, 0.5 * np.sin(ang) + 0.02, r"$r$", color=C_SH, fontsize=14)
# anti-CoE axis
ax.plot([-D, 2.2], [0, 0], color=C_AX, lw=1.0, ls="--", zorder=1)
ax.annotate("", xy=(2.22, 0), xytext=(2.05, 0), arrowprops=dict(arrowstyle="-|>", color=C_AX, lw=1.3))
ax.text(1.05, 0.05, "anti-CoE axis\n" r"($\mu=+1$, $\ell_0\approx15^\circ$)", color=C_AX, ha="left",
        va="bottom", fontsize=9)
# CoE direction
ax.annotate("toward CoE\n" r"($\ell\approx195^\circ$)", xy=(-1.05, -0.62), xytext=(-0.35, -0.62),
            color=C_COE, fontsize=9, ha="left", va="center",
            arrowprops=dict(arrowstyle="-|>", color=C_COE, lw=1.2))
# D
yD = -0.22
darrow(ax, (coe[0], yD), (0, yD), C_COE)
for xx in (coe[0], 0.0):
    ax.plot([xx, xx], [yD - 0.04, -0.03], color=C_COE, lw=0.7)
ax.text(-D / 2, yD - 0.05, r"$D$", color=C_COE, fontsize=14, ha="center", va="top")
# d
darrow(ax, O, S, C_SRC, lw=1.8)
nperp = np.array([-np.sin(th), np.cos(th)])
ax.text(*(S / 2 + 0.07 * nperp), r"$d$", color=C_SRC, fontsize=14, ha="center", va="center")
# R_S
ax.plot([coe[0], S[0]], [coe[1], S[1]], color="k", lw=1.0, ls=":", zorder=2)
mid = (coe + S) / 2
ax.text(mid[0] - 0.02, mid[1] + 0.06, r"$R_S$", fontsize=14, ha="center")
# theta
ax.add_patch(Arc(O, 0.36, 0.36, theta1=0, theta2=np.degrees(th), color="k", lw=1.0))
ax.text(0.24 * np.cos(th / 2), 0.24 * np.sin(th / 2), r"$\theta$", fontsize=13, ha="center", va="center")
# points
ax.plot(*coe, marker="*", ms=20, color=C_COE, zorder=5)
ax.text(coe[0], 0.1, "CoE", color=C_COE, ha="center", va="bottom", fontsize=11, weight="bold")
ax.plot(*O, marker="o", ms=8, color=C_OBS, zorder=5)
ax.text(0.04, -0.04, "O (us)", color=C_OBS, ha="left", va="top", fontsize=10.5, weight="bold")
ax.plot(*S, marker="s", ms=8, color=C_SRC, zorder=5)
ax.text(S[0] + 0.05, S[1], "S (X-ray clusters)", color=C_SRC, fontsize=10, ha="left", va="center")
ax.text(0.995, 0.99, r"$d<r<D$" "\n" "CoE outside shell,\nsources inside",
        transform=ax.transAxes, fontsize=11.5, va="top", ha="right",
        bbox=dict(boxstyle="round", fc="white", ec=C_COE))
ax.text(0.01, 0.01, r"$R_S^2=D^2+d^2+2Dd\cos\theta$,   $\mu=\cos\theta$" "\n"
        r"NOT TO SCALE ($d$ drawn $\approx19\times$ too large relative to $r$)",
        transform=ax.transAxes, fontsize=9.5, va="bottom", ha="left")
ax.set_xlim(-2.05, 2.3); ax.set_ylim(-1.55, 1.12); ax.set_aspect("equal"); ax.axis("off")
ax.set_title(r"(a) geometry, $d<r<D$ (schematic)", fontsize=11.5)

# ---------------- (b) true scale
bx = fig.add_subplot(gs[1])
bx.add_patch(Circle((0, 0), r_true, fill=True, fc="#f4f4f4", ec=C_SH, lw=1.5, zorder=0))
bx.plot([-D_ill, 16], [0, 0], color=C_AX, lw=0.9, ls="--", zorder=1)
bx.plot(-D_ill, 0, "*", color=C_COE, ms=14, zorder=5)
bx.text(-D_ill, 0.9, "CoE", color=C_COE, ha="center", va="bottom", fontsize=10, weight="bold")
bx.plot(0, 0, "o", color=C_OBS, ms=3, zorder=5)
bx.add_patch(Circle((0, 0), d_true, fill=True, fc=C_SRC, ec=C_SRC, lw=0.5, zorder=4))
darrow(bx, (-D_ill, -1.6), (0, -1.6), C_COE, lw=1.1)
bx.text(-D_ill / 2, -2.0, r"$D>r$ (drawn at 20 Gpc;" "\n" "value not fixed by the fit)", color=C_COE,
        fontsize=8.5, ha="center", va="top")
ang = np.radians(55)
bx.plot([0, r_true * np.cos(ang)], [0, r_true * np.sin(ang)], color=C_SH, lw=1.0)
bx.text(0.5 * r_true * np.cos(ang) - 1.8, 0.5 * r_true * np.sin(ang) + 0.6, r"$r=14.3$",
        color=C_SH, fontsize=11, rotation=55, ha="center", va="center")
bx.set_xlim(-22, 16); bx.set_ylim(-16, 16); bx.set_aspect("equal")
bx.set_xlabel("Gpc"); bx.tick_params(labelsize=8)
bx.set_title("(b) true scale (Gpc)", fontsize=11.5)
ins = inset_axes(bx, width="36%", height="36%", loc="lower right", borderpad=1.0)
ins.add_patch(Circle((0, 0), d_true, fill=False, ec=C_SRC, lw=1.2, ls="-"))
ins.plot(0, 0, "o", color=C_OBS, ms=4)
Sx = d_true * np.array([np.cos(th), np.sin(th)])
ins.annotate("", xy=Sx, xytext=(0, 0), arrowprops=dict(arrowstyle="<->", color=C_SRC, lw=1.2, shrinkA=0, shrinkB=0))
ins.plot(*Sx, "s", color=C_SRC, ms=4)
ins.text(-0.02, 0.2, r"$d$", color=C_SRC, fontsize=11)
ins.plot([-0.6, 0.6], [0, 0], color=C_AX, lw=0.8, ls="--")
ins.set_xlim(-0.6, 0.6); ins.set_ylim(-0.6, 0.6); ins.set_aspect("equal")
ins.tick_params(labelsize=6)
ins.set_title(r"zoom: $d=0.32$ Gpc ($d/r=0.022$)", fontsize=7.5)
mark_inset(bx, ins, loc1=2, loc2=1, fc="none", ec="0.5", lw=0.6)
bx.text(0.02, 0.98, r"$d:r=0.32:14.3$ Gpc" "\n" r"$d/D<d/r=0.022$", transform=bx.transAxes, fontsize=9,
        va="top", ha="left", bbox=dict(boxstyle="round", fc="white", ec="0.7"))

fig.savefig("geometry_D_r_d.png", bbox_inches="tight")
fig.savefig("geometry_D_r_d.pdf", bbox_inches="tight")
print(f"d={d_true:.4f} Gpc, d/r={d_true/r_true:.4f}")
