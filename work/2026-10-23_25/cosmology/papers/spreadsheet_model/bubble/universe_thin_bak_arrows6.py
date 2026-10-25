import json, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
from scipy.optimize import brentq
from PIL import Image
from bubble_core import AFfour, ML, MA, MLO, MHI, TH1, H0, GSI, RKM, GPC, IMIG, r2
import analysis_curve as AC
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
c = 2.99792458e8; Msun = 1.989e30; Gm = GPC * 1e3
def phys(k):
    d = json.load(open("thin_%s.json" % k.replace("|", "_")))
    Dm = d["D_Gpc"] * Gm; rhoc = 3 * H0**2 / (8 * np.pi * GSI); rhos = d["Os"] * rhoc; Mcrit = 4 * np.pi / 3 * rhoc * Dm**3
    Mc = d["mcc"] * Mcrit; rin = d["rin"] * Dm; rout = d["rout"] * Dm
    Msh = lambda r: 4 * np.pi / 3 * rhos * (np.clip(r, rin, rout)**3 - rin**3)
    def Phi(r):
        if r < rin: ps = -GSI * 2 * np.pi * rhos * (rout**2 - rin**2)
        elif r < rout: ps = -GSI * (4 * np.pi / 3 * rhos * (r**3 - rin**3) / r + 2 * np.pi * rhos * (rout**2 - r**2))
        else: ps = -GSI * Msh(rout) / r
        return -GSI * Mc / r + ps
    v2 = d["vesc2"] * 2 * abs(Phi(Dm)); E = v2 / 2 + Phi(Dm)
    rt = brentq(lambda r: Phi(r) - E, Dm, 1e6 * Dm) if E < 0 else np.inf
    rs = lambda M: 2 * GSI * M / c**2 / Gm
    Mt = Mc + Msh(min(rt, 1e300))
    return dict(d=d, D=d["D_Gpc"], R=RKM / GPC, DR=d["DR"], Os=d["Os"], rhos=rhos, Mc=Mc, Mcrit=Mcrit, mcc=d["mcc"], Msh=Msh(Dm), Mtot=Mc + Msh(Dm),
                Mhole=4 * np.pi / 3 * rhos * rin**3, Rin=rin / Gm, Rout=rout / Gm, rs_c=rs(Mc), rs_D=rs(Mc + Msh(Dm)), rt=rt / Gm, rs_t=rs(Mt),
                vE=np.sqrt(v2) / c, vesc2=d["vesc2"], age=d["age"], t_turn=d["t_turn"], amp=d["amp"], phi=d["phi"], l0=d["l0"], R2=d["R2_mig"], R2F=d["R2_four"],
                Lc=d["Lc"] * d["D_Gpc"], Dkm=d["DR"] * RKM, Lckm=d["Lc"] * d["DR"] * RKM)
P1, P2 = phys("b|free|hug"), phys("b|0.05|hug")
json.dump({k: {kk: vv for kk, vv in p.items() if kk != "d"} for k, p in (("free", P1), ("baryons", P2))}, open("universe_thin.json", "w"), indent=1, default=float)
for p in (P1, P2): print({k: (round(v, 4) if isinstance(v, float) and abs(v) < 1e6 else v) for k, v in p.items() if k != "d"})
p = P1; Dg, Rg, l0 = p["D"], p["R"], p["l0"]
lpk = (90 - np.degrees(p["phi"]) + l0) % 360
Ev = np.array([Dg, 0.0])
def sky(l, r): t = np.radians(l - l0); return Ev + r * np.array([np.cos(t), np.sin(t)])
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.4, 3.95), gridspec_kw=dict(width_ratios=[1.6, 1], wspace=0.2))
ext = (-10, 127, -47, 47)
ax.add_patch(Circle((0, 0), p["Rout"], fc="0.82", ec="0.35", lw=0.9, zorder=0))
ax.add_patch(Circle((0, 0), p["Rin"], fc="white", ec="0.35", lw=0.9, zorder=0.5))
ax.add_patch(Circle((0, 0), Dg, fill=False, color="C0", lw=1.0, ls="--", zorder=2))
ax.add_patch(Circle(Ev, Rg, fill=True, fc="C2", alpha=0.18, ec="C2", lw=1.4, zorder=2))
ax.add_patch(Circle((0, 0), p["rs_c"], fill=False, color="C3", lw=1.4, zorder=3))
ax.plot(0, 0, "o", color="k", ms=4, zorder=5)
ax.text(0, -3, "CoE core", ha="center", va="top", fontsize=7.5)
ax.text(36, 18, "empty interior", ha="center", fontsize=7.5, color="0.35")
ax.text(p["Rin"] - 4, 29, "shell, $D-R$ to $D+R$", ha="right", fontsize=7, color="0.3")
ax.text(-8, 44, "red arc: core $r_s$ = %.0f Gpc (%.2f D)" % (p["rs_c"], p["rs_c"] / Dg), ha="left", va="top", fontsize=6.8, color="C3")
ax.plot(*Ev, "o", color="C0", ms=4.5, zorder=6, mec="k"); ax.text(Dg - 1.8, -2.2, "Earth", ha="right", va="top", fontsize=7, zorder=7)
ax.annotate("", xy=(Dg, -32), xytext=(0, -32), arrowprops=dict(arrowstyle="<->", lw=0.9))
ax.text(Dg / 2, -33.5, "$D$ = %.1f Gpc" % Dg, ha="center", va="top", fontsize=8)
ax.annotate("", xy=(p["rs_c"] - 2.5, 0), xytext=Ev + [0.8, 0], arrowprops=dict(arrowstyle="simple,head_width=0.9,head_length=0.9,tail_width=0.38", fc="#0055d4", ec="#0055d4", lw=0.5), zorder=8)
t = np.radians(125); q = Ev + Rg * np.array([np.cos(t), np.sin(t)])
ax.annotate("", xy=q, xytext=Ev, arrowprops=dict(arrowstyle="->", color="C2", lw=1.0)); ax.text(q[0] - 0.8, q[1] + 0.5, "$R$ = %.1f Gpc" % Rg, color="C2", fontsize=7.5, ha="right", va="bottom")
for l in (0, 90, 180, 270):
    a, b = sky(l, Rg - 1.2), sky(l, Rg + 1.2); ax.plot([a[0], b[0]], [a[1], b[1]], color="0.15", lw=1.0)
    q = sky(l, Rg + 4.6); ax.text(*q, "$\\ell$=%d°" % l, ha="center", va="center", fontsize=6.3, color="0.15", zorder=7)
ax.text(Dg + 11, -31, "velocity $v\\approx$%.1f$c$ outward\n(anti-CoE, $\\ell$=%.1f°)" % (p["vE"], l0), ha="center", va="bottom", fontsize=6.6, color="#0055d4", fontweight="bold", bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.8))
ax.annotate("", xy=(Dg + 12, -1.8), xytext=(Dg + 11, -23.5), arrowprops=dict(arrowstyle="-", color="#0055d4", lw=0.7))
a, b = sky(lpk, 0), sky(lpk, Rg + 3); ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color="C4", lw=1.2, ls="--"))
ax.text(Dg + 13, 26, "sine peak\n$\\ell$=%.1f°" % lpk, ha="center", va="bottom", fontsize=6.5, color="C4")
ax.annotate("", xy=b + [0, 1], xytext=(Dg + 13, 25.5), arrowprops=dict(arrowstyle="-", color="C4", lw=0.6))
ax.text(-8, -45, "circle of radius $D$ (dashed)", color="C0", fontsize=6.8, va="bottom")
ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")
ax.set_xlabel("x [Gpc] (CoE to Earth)"); ax.set_ylabel("y [Gpc]"); ax.set_title("(a) near view", fontsize=9)
E2 = 1150
bx.add_patch(Circle((0, 0), p["Rout"], fc="0.6", ec="0.3", lw=0.6, zorder=1)); bx.add_patch(Circle((0, 0), p["Rin"], fc="white", ec="0.3", lw=0.4, zorder=1.5))
bx.add_patch(Circle(Ev, Rg, fc="C2", alpha=0.5, ec="C2", lw=0.6, zorder=2))
bx.add_patch(Circle((0, 0), p["rs_c"], fill=False, color="C3", lw=1.2, zorder=3))
bx.add_patch(Circle((0, 0), p["rt"], fill=False, color="C4", lw=1.1, ls="-.", zorder=3))
bx.text(0, p["rs_c"] + 25, "core $r_s$", ha="center", va="bottom", fontsize=6.5, color="C3")
bx.text(0, p["rt"] - 45, "Earth's turnaround\n%.0f Gpc (%.1f D)" % (p["rt"], p["rt"] / Dg), ha="center", va="top", fontsize=6.5, color="C4")
bx.text(0, -250, "shell %.0f–%.0f Gpc" % (p["Rin"], p["Rout"]), ha="center", va="top", fontsize=6.5, color="0.2")
bx.text(0, -p["rt"] + 60, "$r_s(<D)$ = %s Gpc (%.0f D),\noff scale" % ("{:,.0f}".format(p["rs_D"]), p["rs_D"] / Dg), ha="center", va="bottom", fontsize=6.5, color="C3")
for ang in (0, 30, 150, 180, 210, 330):
    u = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
    bx.annotate("", xy=u * 800, xytext=u * 200, arrowprops=dict(arrowstyle="simple,head_width=0.7,head_length=0.8,tail_width=0.25", fc="#0055d4", ec="#0055d4", lw=0.4), zorder=4)
bx.text(0, 600, "outward velocity\n(decelerating)", ha="center", va="center", fontsize=6.8, color="#0055d4", fontweight="bold")
bx.set_xlim(-E2, E2); bx.set_ylim(-E2, E2); bx.set_aspect("equal"); bx.set_xlabel("x [Gpc]"); bx.set_title("(b) zoomed out", fontsize=9)
fig.savefig("../universe_thin_drawing.png", dpi=200, bbox_inches="tight"); plt.close(fig)
AH = float(os.environ.get("AH", "5.25"))
Af = AFfour(TH1); th2 = np.append(TH1, 360); Af2 = AFfour(th2)
fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, AH), sharex=True, gridspec_kw={"height_ratios": [2.0, 1], "hspace": 0.06})
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
a.plot(th2, Af2, color="C0", lw=1.5, label="$A_F$")
for pp, col, ls, lab in [(P1, "C3", "-", "thin shell + core, $\\Omega_s$=%.2f"), (P2, "C1", "--", "thin shell + core, $\\Omega_s$=%.2f (baryons)")]:
    A = AC.curve("ptshell", pp["d"]); assert abs(r2(MA, A[IMIG]) - pp["R2"]) < 1e-6
    A2 = np.append(A, A[0])
    a.plot(th2, A2, color=col, ls=ls, lw=1.6 if ls == "-" else 1.3, label=(lab % pp["Os"]) + ", $\\ell_0$=%.1f° ($R^2$=%.3f; vs $A_F$ %.3f)" % (pp["l0"], pp["R2"], pp["R2F"]))
    b.plot(th2, 100 * (A2 - Af2) / Af2, color=col, ls=ls, lw=1.4 if ls == "-" else 1.2, label="model" if ls == "-" else "baryons")
a.set_ylabel("$A/A_{\\rm all}$"); a.set_ylim(0.62, 1.22); a.legend(fontsize=7.2, loc="lower left", frameon=False)
b.axhline(0, color="k", lw=0.6); b.plot(ML, 100 * (MA - AFfour(ML)) / AFfour(ML), "o", ms=2.5, color="0.45", label="data")
b.set_ylabel("$100\\,(\\cdot-A_F)/A_F$ [%]", fontsize=8); b.set_xlabel("Galactic longitude $\\ell$ (deg)")
b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45)); b.legend(fontsize=7, frameon=False, loc="upper right", ncol=3)
fig.savefig("../universe_thin_Achart.png", dpi=200, bbox_inches="tight")
for f in ["universe_thin_drawing", "universe_thin_Achart"]:
    im = Image.open("../" + f + ".png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../" + f + "_web.jpg", quality=88)
