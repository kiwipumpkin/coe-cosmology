import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.optimize import brentq
from PIL import Image
from bubble_core import AFfour, ML, MA, MLO, MHI, TH1, H0, GSI, RKM, GPC
import analysis_curve as AC
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
d = json.load(open("ptshell_results.json"))["inside|Os=0.05|mig_l0|c=0.02"]
c = 2.99792458e8; Msun = 1.989e30
Rg = RKM / GPC; Dg = d["DR"] * Rg; l0 = d["l0"]
Dm = Dg * GPC * 1e3
rhoc = 3 * H0**2 / (8 * np.pi * GSI); rhos = d["Os"] * rhoc
Mcrit = 4 * np.pi / 3 * rhoc * Dm**3
Mc = d["mcc"] * Mcrit; rin, rout = d["rin"] * Dm, d["rout"] * Dm
def Mshell(r): return 4 * np.pi / 3 * rhos * (np.clip(r, rin, rout)**3 - rin**3)
def Menc(r): return Mc + Mshell(r)
Mtot = Menc(Dm); Msh = Mshell(Dm)
rs = lambda M: 2 * GSI * M / c**2
rs_c, rs_D = rs(Mc) / (GPC * 1e3), rs(Mtot) / (GPC * 1e3)
def Phi(r):
    if r < rin: ps = -GSI * 2 * np.pi * rhos * (rout**2 - rin**2)
    elif r < rout: ps = -GSI * (4 * np.pi / 3 * rhos * (r**3 - rin**3) / r + 2 * np.pi * rhos * (rout**2 - r**2))
    else: ps = -GSI * Mshell(rout) / r
    return -GSI * Mc / r + ps
vesc2 = 2 * abs(Phi(Dm)); v2 = d["vesc2"] * vesc2; E = v2 / 2 + Phi(Dm)
rt = brentq(lambda r: Phi(r) - E, Dm, 1e4 * Dm)
Mt = Menc(rt); rs_t = rs(Mt)
vE = np.sqrt(v2)
chk = dict(D_Gpc=Dg, R_Gpc=Rg, Lc_Gpc=d["Lc"] * Dg, Rin_Gpc=rin / GPC / 1e3, Rout_Gpc=rout / GPC / 1e3, Mc_kg=Mc, Mc_sun=Mc / Msun, Msh_kg=Msh, Mtot_kg=Mtot,
           Mtot_json=d["M_kg"], rhoc=rhoc, rhos=rhos, rs_c=rs_c, rs_D=rs_D, rs_c_D=rs_c / Dg, rs_D_D=rs_D / Dg, vE_over_c=vE / c, rt_Gpc=rt / GPC / 1e3,
           rt_D=rt / Dm, rs_t_Gpc=rs_t / GPC / 1e3, Mt_kg=Mt, R2=d["R2_mig"], R2F=d["R2_four"], age=d["age"], vesc2=d["vesc2"], amp=d["amp"], phi=d["phi"], l0=l0,
           DR=d["DR"], Dkm=d["DR"] * RKM, Lckm=d["Lc"] * d["DR"] * RKM, Rinkm=d["rin"] * d["DR"] * RKM, Routkm=d["rout"] * d["DR"] * RKM, t_turn=d["t_turn"], Mhole=d["Mhole"])
print(json.dumps(chk, indent=1)); json.dump(chk, open("universe_real.json", "w"), indent=1)
lpk = (90 - np.degrees(d["phi"]) + l0) % 360
Ev = np.array([Dg, 0.0])
def sky(l, r): t = np.radians(l - l0); return Ev + r * np.array([np.cos(t), np.sin(t)])
Rin, Rout = chk["Rin_Gpc"], chk["Rout_Gpc"]
def shade(ax, ext, n):
    xx, yy = np.meshgrid(np.linspace(ext[0], ext[1], n), np.linspace(ext[2], ext[3], n))
    rr = np.hypot(xx, yy); z = np.where((rr > Rin) & (rr < Rout), 0.28, 0.0)
    return ax.imshow(z, extent=ext, origin="lower", cmap="Greys", vmin=0, vmax=1, zorder=0, interpolation="nearest")
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.4, 3.95), gridspec_kw=dict(width_ratios=[1.35, 1], wspace=0.22))
ext = (-22, 80, -38, 38); shade(ax, ext, 1200)
ax.add_patch(Circle((0, 0), Rin, fill=False, color="0.3", lw=0.8, ls=":", zorder=2))
ax.add_patch(Circle((0, 0), Dg, fill=False, color="C0", lw=1.2, ls="--", zorder=2))
ax.add_patch(Circle(Ev, Rg, fill=True, fc="C2", alpha=0.15, ec="C2", lw=1.5, zorder=2))
ax.add_patch(Circle((0, 0), rs_c, fill=False, color="C3", lw=1.4, zorder=3))
ax.plot(0, 0, "o", color="k", ms=3.5, zorder=5)
ax.annotate("point core", xy=(0.3, -0.3), xytext=(6, -9), fontsize=7, zorder=6, arrowprops=dict(arrowstyle="-", lw=0.7))
ax.annotate("hole, $R_{\\rm in}$ = %.1f Gpc" % Rin, xy=(Rin * np.cos(4.0), Rin * np.sin(4.0)), xytext=(-20, -22), fontsize=7, color="0.2",
            arrowprops=dict(arrowstyle="-", color="0.3", lw=0.7))
ax.annotate("$r_s$ of core = %.1f Gpc (%.2f D)" % (rs_c, rs_c / Dg), xy=(rs_c * np.cos(2.0), rs_c * np.sin(2.0)), xytext=(-20, 20), fontsize=7, color="C3",
            arrowprops=dict(arrowstyle="-", color="C3", lw=0.7))
ax.plot(*Ev, "o", color="C0", ms=5, zorder=6, mec="k"); ax.text(Dg - 0.8, 1.2, "Earth", ha="right", va="bottom", fontsize=8, zorder=7)
ax.annotate("", xy=(Dg, -27), xytext=(0, -27), arrowprops=dict(arrowstyle="<->", lw=0.9))
ax.text(Dg / 2, -28.5, r"$D$ = %.1f Gpc" % Dg, ha="center", va="top", fontsize=8)
ax.annotate("", xy=Ev + [8, 0], xytext=Ev, arrowprops=dict(arrowstyle="-|>", color="C3", lw=1.6), zorder=6)
t = np.radians(125); q = Ev + Rg * np.array([np.cos(t), np.sin(t)])
ax.annotate("", xy=q, xytext=Ev, arrowprops=dict(arrowstyle="->", color="C2", lw=1.0)); ax.text(q[0] - 0.5, q[1] + 0.6, r"$R$ = %.1f Gpc" % Rg, color="C2", fontsize=8, ha="right", va="bottom")
for l in (0, 90, 180, 270):
    a, b = sky(l, Rg - 0.9), sky(l, Rg + 0.9); ax.plot([a[0], b[0]], [a[1], b[1]], color="0.2", lw=1.1)
    q = sky(l, Rg + 3.6); ax.text(*q, r"$\ell$=%d°" % l, ha="center", va="center", fontsize=7, color="0.2")
q = sky(l0, Rg + 5.0); ax.text(q[0] + 1.0, q[1] + 6.5, "anti-CoE\n$\\ell$=%.1f°" % l0, ha="center", va="center", fontsize=7, color="C3")
a, b = sky(lpk, 0), sky(lpk, Rg + 4.5); ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color="C4", lw=1.2, ls="--"))
q = sky(lpk, Rg + 4.8); ax.text(q[0] + 0.3, q[1], "sine peak\n$\\ell$=%.1f°" % lpk, ha="left", va="bottom", fontsize=7, color="C4")
ax.text(-20, 35, "circle of radius $D$", color="C0", fontsize=7.5, va="top")
ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")
ax.set_xlabel("x [Gpc] (CoE to Earth)"); ax.set_ylabel("y [Gpc]"); ax.set_title("(a) near view", fontsize=9)
E2 = 580; shade(bx, (-E2, E2, -E2, E2), 900)
bx.add_patch(Circle((0, 0), Rout, fill=False, color="k", lw=1.0, zorder=2))
bx.add_patch(Circle((0, 0), Dg, fill=False, color="C0", lw=1.0, ls="--", zorder=2))
bx.add_patch(Circle(Ev, Rg, fill=True, fc="C2", alpha=0.3, ec="C2", lw=1.0, zorder=2))
bx.add_patch(Circle((0, 0), rs_D, fill=False, color="C3", lw=1.4, ls="-", zorder=3))
bx.add_patch(Circle((0, 0), rt / GPC / 1e3, fill=False, color="C4", lw=1.1, ls="-.", zorder=3))
bx.plot(0, 0, "o", color="k", ms=2.5)
bx.text(0, Rout + 12, "front $R_{\\rm out}$\n%.0f Gpc (%.2f D)" % (Rout, Rout / Dg), ha="center", va="bottom", fontsize=6.3, linespacing=1.0)
bx.text(0, -rs_D - 14, "$r_s(<D)$ = %.0f Gpc (%.1f D)" % (rs_D, rs_D / Dg), ha="center", va="top", fontsize=7, color="C3")
bx.text(0, rt / GPC / 1e3 - 45, "Earth's turnaround\n%.0f Gpc (%.1f D)" % (rt / GPC / 1e3, rt / Dm), ha="center", va="top", fontsize=6.5, color="C4")
bx.set_xlim(-E2, E2); bx.set_ylim(-E2, E2); bx.set_aspect("equal"); bx.set_xlabel("x [Gpc]"); bx.set_title("(b) zoomed out", fontsize=9)
fig.savefig("../universe_real_drawing.png", dpi=200, bbox_inches="tight"); plt.close(fig)
A = AC.curve("ptshell", d); Af = AFfour(TH1); th2 = np.append(TH1, 360); A2 = np.append(A, A[0]); Af2 = AFfour(th2)
fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, float(__import__("os").environ.get("AH", "6.15"))), sharex=True, gridspec_kw={"height_ratios": [2.0, 1], "hspace": 0.06})
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
a.plot(th2, Af2, color="C0", lw=1.5, label=r"$A_F$")
a.plot(th2, A2, color="C3", lw=1.6, label=r"point + shell (baryons only), $\ell_0$=%.1f° ($R^2$=%.3f; vs $A_F$ %.3f)" % (l0, d["R2_mig"], d["R2_four"]))
a.set_ylabel(r"$A/A_{\rm all}$"); a.set_ylim(0.66, 1.22); a.legend(fontsize=7.5, loc="lower left", frameon=False)
b.axhline(0, color="k", lw=0.6)
b.plot(ML, 100 * (MA - AFfour(ML)) / AFfour(ML), "o", ms=2.5, color="0.45", label="data")
b.plot(th2, 100 * (A2 - Af2) / Af2, color="C3", lw=1.4, label="model")
b.set_ylabel(r"$100\,(\cdot-A_F)/A_F$ [%]", fontsize=8); b.set_xlabel(r"Galactic longitude $\ell$ (deg)")
b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45)); b.legend(fontsize=7, frameon=False, loc="upper right", ncol=2)
fig.savefig("../universe_real_Achart.png", dpi=200, bbox_inches="tight")
from bubble_core import r2
from bubble_core import IMIG
print("R2 check", r2(MA, A[IMIG]), r2(Af, A))
for f in ["universe_real_drawing", "universe_real_Achart"]:
    im = Image.open("../" + f + ".png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../" + f + "_web.jpg", quality=88)
