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
# ---- relativistic (exact LTB) thick ordinary-matter shell, z = 0.3 case (gr_ltb/refit_mass) ----
def rel():
    o = json.load(open("universe_rel_mcfree.json"))   # all free: D, Mc, R_in, R_out, l0, Amp, phi (old Mc-constrained fit: best_rows_thick.json)
    D = o["D"]; Dm = D * Gm; rhoc = 3 * H0**2 / (8 * np.pi * GSI); rhos = o["Os"] * rhoc
    rs = lambda M: 2 * GSI * M / c**2 / Gm; Rv = RKM / GPC; Rin, Rout = o["R_in"], o["R_out"]
    cur = np.array(o["A"])
    return dict(D=D, R=Rv, DR=D / Rv, Os=o["Os"], rhos=rhos, Mc=o["Mc_kg"], Mcrit=4 * np.pi / 3 * rhoc * Dm**3, mcc=o["Mc_kg"] / (4 * np.pi / 3 * rhoc * Dm**3),
                Msh=o["M_D_kg"] - o["Mc_kg"], Mtot=o["M_D_kg"], Mshtot=4 * np.pi / 3 * rhos * ((Rout * Gm)**3 - (Rin * Gm)**3), Mhole=4 * np.pi / 3 * rhos * (Rin * Gm)**3,
                Mvis=4 * np.pi / 3 * rhos * (Rv * Gm)**3, Rin=Rin, Rout=Rout, rs_c=rs(o["Mc_kg"]), rs_D=rs(o["M_D_kg"]), vE=H0 * Dm / c, twoE=o["twoE"],
                age=o["age_Gyr"], amp=o["amp"], phi=o["phi"], l0=o["l0_s"], R2=o["R2_s"], R2ns=o["R2_ns"], R2F=r2(AFfour(TH1), cur), H_r=o["H_r"], H_perp=o["H_perp"],
                z=o["z"], Lc=0.005 * D, Dkm=D * GPC, Lckm=0.005 * D * GPC, cur=cur)
P3 = rel()
def newt_omb():   # Newtonian, Omega_s = 0.049 forced (newt_omb.py)
    o = json.load(open("universe_newt_omb.json"))["Rs14.3|sine"]; D = o["D"]; Dm = D * Gm; rhoc = 3 * H0**2 / (8 * np.pi * GSI); rhos = o["Os"] * rhoc
    rs = lambda M: 2 * GSI * M / c**2 / Gm; Rv = RKM / GPC; cur = np.array(o["A"])
    return dict(D=D, R=Rv, DR=D / Rv, Os=o["Os"], rhos=rhos, Mc=o["Mc_kg"], Mcrit=o["Mc_kg"] / o["mcc"], mcc=o["mcc"], Msh=o["shell_inside_D_kg"], Mtot=o["M_inside_D_kg"],
                Mshtot=o["shell_total_kg"], Mhole=4 * np.pi / 3 * rhos * (o["R_in"] * Gm)**3, Mvis=4 * np.pi / 3 * rhos * (Rv * Gm)**3, Rin=o["R_in"], Rout=o["R_out"],
                rs_c=rs(o["Mc_kg"]), rs_D=rs(o["M_inside_D_kg"]), vE=o["VO_over_c"], vesc2=o["vesc2_full_potential"], age=o["age_Gyr"], amp=o["amp"], phi=o["phi"], l0=o["l0"],
                R2=o["R2"], R2geo=json.load(open("universe_newt_omb.json"))["Rs14.3|geom"]["R2"], R2F=r2(AFfour(TH1), cur), Lc=0.005 * D, Dkm=D * GPC, Lckm=0.005 * D * GPC, cur=cur)
P4 = newt_omb()
for q in (P1, P2): q["Mshtot"] = 4 * np.pi / 3 * q["rhos"] * ((q["Rout"] * Gm)**3 - (q["Rin"] * Gm)**3); q["Mvis"] = 4 * np.pi / 3 * q["rhos"] * (q["R"] * Gm)**3
json.dump({k: {kk: vv for kk, vv in q.items() if kk not in ("d", "cur")} for k, q in (("free", P1), ("baryons", P2), ("rel", P3), ("newt_omb", P4))}, open("universe_thin.json", "w"), indent=1, default=float)
print("REL", {k: (round(v, 4) if isinstance(v, float) and abs(v) < 1e6 else v) for k, v in P3.items() if k != "cur"})
p = P3; Dg, Rg, l0 = p["D"], p["R"], p["l0"]
lpk = (90 - np.degrees(p["phi"]) + l0) % 360
Ev = np.array([Dg, 0.0])
def sky(l, r): t = np.radians(l - l0); return Ev + r * np.array([np.cos(t), np.sin(t)])
# ---------------- Figure 1: one figure, three axes, unified style ----------------
from matplotlib.gridspec import GridSpec
from itertools import combinations
FS_A, FS_L, FS_T = 7.5, 8.5, 9            # annotations, axis labels/ticks, titles
plt.rcParams.update({"mathtext.fontset": "cm", "axes.linewidth": 0.8, "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.labelsize": FS_L, "axes.titlesize": FS_T,
                     "xtick.major.width": 0.8, "ytick.major.width": 0.8})
CRED, CBLU, CDC, CGRN, CPUR, CSH = "#d62728", "#0055d4", "#1f77b4", "#2ca02c", "#9467bd", "0.82"
VEL = dict(arrowstyle="-|>", color=CBLU, lw=1.3, mutation_scale=9, shrinkA=0, shrinkB=0)
POS = lambda ls="-": dict(arrowstyle="-|>", color="0.2", lw=0.9, ls=ls, mutation_scale=8, shrinkA=0, shrinkB=0)
exA = (-0.6, 5.6, -2.6, 3.4); exB = (-7, 82.7, -34.5, 34.5); E2 = 230
ka = (exA[1] - exA[0]) / (exA[3] - exA[2]); kb = (exB[1] - exB[0]) / (exB[3] - exB[2])
FW = 7.4; L0, G1, G2, R0 = 0.04, 0.58, 0.54, 0.10
H = (FW - L0 - G1 - G2 - R0) / (ka + kb + 1); BOT, TOP = 0.42, 0.24; FH = H + BOT + TOP
fig = plt.figure(figsize=(FW, FH))
gs = GridSpec(1, 5, figure=fig, width_ratios=[ka, G1 / H, kb, G2 / H, 1], wspace=0, left=L0 / FW, right=1 - R0 / FW, bottom=BOT / FH, top=1 - TOP / FH)
vx, ax, bx = fig.add_subplot(gs[0]), fig.add_subplot(gs[2]), fig.add_subplot(gs[4])
def T(axx, *a, **k): k.setdefault("fontsize", FS_A); return axx.text(*a, **k)
# (a) Hubble velocity vectors (schematic)
C0_, O_, P_ = np.array([0, 0.9]), np.array([2.6, 0.9]), np.array([3.6, 2.3]); h_ = 0.4
VO, VP = h_ * (O_ - C0_), h_ * (P_ - C0_); VH = VP - VO
vx.annotate("", xy=O_, xytext=C0_, arrowprops=POS()); vx.annotate("", xy=P_, xytext=C0_, arrowprops=POS("--")); vx.annotate("", xy=P_, xytext=O_, arrowprops=POS(":"))
for base, vec in ((O_, VO), (P_, VP), (P_, VH)): vx.annotate("", xy=base + vec, xytext=base, arrowprops=VEL)
vx.plot(*C0_, "o", color=CRED, ms=5, zorder=6); vx.plot(*O_, "o", color=CBLU, mec="k", mew=0.6, ms=5, zorder=6); vx.plot(*P_, "o", color="0.3", ms=4.5, zorder=6)
T(vx, *(C0_ + [0, -0.18]), "CoE", color=CRED, ha="center", va="top")
T(vx, *(O_ + [0.02, -0.18]), "us ($O$)", ha="left", va="top")
T(vx, *(P_ + [0.1, -0.1]), "galaxy", ha="left", va="top")
T(vx, *((C0_ + O_) / 2 + [0, -0.12]), "$\\vec D$", ha="center", va="top")
T(vx, *((C0_ + P_) / 2 + [-0.12, 0.12]), "$\\vec D+\\vec r$", ha="right", va="bottom")
T(vx, *((O_ + P_) / 2 + [0.1, -0.05]), "$\\vec r$", ha="left", va="center")
T(vx, *(O_ + VO + [0.08, 0]), "$\\vec V_O$", color=CBLU, ha="left", va="bottom")
T(vx, *(P_ + VP + [0.04, -0.06]), "$\\vec V_P$", color=CBLU, ha="left", va="top")
T(vx, *(P_ + VH + [0.07, -0.02]), "$\\vec V_H$", color=CBLU, ha="left", va="bottom")
A_ = np.array([1.1, -1.45]); s_ = 1.6
for base, vec in ((A_, s_ * VO), (A_, s_ * VP), (A_ + s_ * VO, s_ * VH)): vx.annotate("", xy=base + vec, xytext=base, arrowprops=VEL)
T(vx, *(A_ + s_ * VO / 2 + [0, -0.1]), "$\\vec V_O$", color=CBLU, ha="center", va="top")
T(vx, *(A_ + s_ * VP / 2 + [-0.1, 0.1]), "$\\vec V_P$", color=CBLU, ha="right", va="bottom")
T(vx, *(A_ + s_ * VO + s_ * VH / 2 + [0.1, 0]), "$\\vec V_H$", color=CBLU, ha="left", va="center")
T(vx, 2.5, -2.55, "$\\vec V_H=\\vec V_P-\\vec V_O=H\\vec r$", ha="center", va="bottom")
vx.set_xlim(exA[0], exA[1]); vx.set_ylim(exA[2], exA[3]); vx.set_aspect("equal"); vx.set_xticks([]); vx.set_yticks([])
vx.set_title("(a) Hubble velocity vectors")
# (b) Hubble bubble scale
ax.add_patch(Circle((0, 0), p["Rout"], fc=CSH, ec="0.35", lw=0.8, zorder=0)); ax.add_patch(Circle((0, 0), p["Rin"], fc="white", ec="0.35", lw=0.8, zorder=0.5))
ax.add_patch(Circle((0, 0), Dg, fill=False, color=CDC, lw=1.0, ls="--", zorder=2))
ax.add_patch(Circle(Ev, Rg, fc=CGRN, alpha=0.18, ec=CGRN, lw=1.2, zorder=2))
ax.add_patch(Circle((0, 0), p["rs_c"], fill=False, color=CRED, lw=1.4, zorder=3))
ax.plot(0, 0, "o", color=CRED, ms=5, zorder=6); T(ax, -5.5, 2.4, "CoE core", color=CRED, ha="left", va="bottom")
ax.plot(*Ev, "o", color=CBLU, mec="k", mew=0.6, ms=5, zorder=7); T(ax, p["Rin"] - 1.4, -1.0, "Earth", ha="right", va="top", zorder=7)
T(ax, 5.5, 11.7, "empty\ninterior", color="0.4", ha="center", va="center", linespacing=1.1)
T(ax, 0.5 * (Dg + p["Rout"]) * np.cos(np.radians(42)), 0.5 * (Dg + p["Rout"]) * np.sin(np.radians(42)), "shell", color="0.25", ha="center", va="center", rotation=-48, rotation_mode="anchor")
tq = np.radians(125); q = Ev + Rg * np.array([np.cos(tq), np.sin(tq)])
ax.annotate("", xy=q, xytext=Ev, arrowprops=dict(arrowstyle="-|>", color=CGRN, lw=1.0, mutation_scale=8, shrinkA=0, shrinkB=0)); T(ax, (Ev[0] + q[0]) / 2 - 1.0, (Ev[1] + q[1]) / 2 - 1.0, "$R$", color=CGRN, ha="right", va="center")
td = np.radians(-38); qd = Dg * np.array([np.cos(td), np.sin(td)])
ax.annotate("", xy=qd, xytext=(0, 0), arrowprops=dict(arrowstyle="<|-|>", color=CDC, lw=0.9, mutation_scale=7, shrinkA=0, shrinkB=0)); T(ax, *(qd * 0.58 + [1.4, -1.4]), "$D$", color=CDC, ha="left", va="top")
T(ax, 43.5, 33, "core $r_s$ = %.0f Gpc\n(%.1f D), off scale" % (p["rs_c"], p["rs_c"] / Dg), color=CRED, ha="left", va="top", linespacing=1.1)
VX1 = 81; ax.annotate("", xy=(VX1, 0), xytext=Ev, arrowprops=VEL, zorder=8)
vlab = T(ax, VX1, -2.8, "$V_O\\!\\approx\\!%.1fc$" % p["vE"], color=CBLU, ha="right", va="top")
ax.set_xlim(exB[0], exB[1]); ax.set_ylim(exB[2], exB[3]); ax.set_aspect("equal")
ax.set_xlabel("$x$ [Gpc]"); ax.set_ylabel("$y$ [Gpc]"); ax.set_title("(b) Hubble bubble scale")
# (c) Cosmic shell scale
bx.add_patch(Circle((0, 0), p["Rout"], fc=CSH, ec="0.35", lw=0.6, zorder=1)); bx.add_patch(Circle((0, 0), p["Rin"], fc="white", ec="0.35", lw=0.5, zorder=1.5))
bx.add_patch(Circle((0, 0), p["rs_c"], fill=False, color=CRED, lw=1.4, zorder=3))
bx.add_patch(Circle(Ev, Rg, fc=CGRN, alpha=0.18, ec=CGRN, lw=1.2, zorder=2))   # vision sphere R, to scale (as in b)
T(bx, 0, 0.87 * p["rs_c"], "core $r_s$\n= %.0f Gpc" % p["rs_c"], color=CRED, ha="center", va="top", linespacing=1.1); T(bx, p["Rout"] + 6, 25, "shell", color="0.25", ha="left", va="center")
T(bx, 0, -0.996 * E2, "$r_s(<D)$ = %s Gpc\n(%.1f D), off scale" % ("{:,.0f}".format(p["rs_D"]), p["rs_D"] / Dg), color=CRED, ha="center", va="bottom", linespacing=1.0)
for ang in (45, 135, 225, 315):
    u = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))]); bx.annotate("", xy=u * 0.87 * p["rs_c"] * 0.75, xytext=u * p["Rout"] * 1.05, arrowprops=VEL, zorder=4)
bx.set_xlim(-E2, E2); bx.set_ylim(-E2, E2); bx.set_aspect("equal"); bx.set_xticks([-200, 0, 200]); bx.set_yticks([-200, 0, 200])
bx.set_xlabel("$x$ [Gpc]"); bx.set_title("(c) Cosmic shell scale")
# checks
fig.canvas.draw(); rr = fig.canvas.get_renderer()
for A_x in (vx, ax, bx):
    bbs = [(t.get_text()[:22], t.get_window_extent(rr)) for t in A_x.texts if t.get_text()]
    for (n1, b1), (n2, b2) in combinations(bbs, 2):
        if b1.overlaps(b2): print("OVERLAP", repr(n1), "|", repr(n2))
    W_ = A_x.get_window_extent(rr)
    for n1, b1 in bbs:
        if b1.x0 < W_.x0 or b1.x1 > W_.x1 or b1.y0 < W_.y0 or b1.y1 > W_.y1: print("OUTSIDE", repr(n1))
for A_x, circ, pad in ((ax, (p["Rin"], p["Rout"], p["rs_c"], Dg), 1.2), (bx, (p["Rin"], p["Rout"], p["rs_c"]), 17)):
    inv_ = A_x.transData.inverted()
    for t in A_x.texts:
        (x0, y0_), (x1, y1) = inv_.transform(t.get_window_extent(rr)); xs, ys = np.meshgrid(np.linspace(x0, x1, 60), np.linspace(y0_, y1, 20)); rr_ = np.hypot(xs, ys)
        for R_ in circ:
            if rr_.min() - pad < R_ < rr_.max() + pad: print("EDGE", repr(t.get_text()[:20]), round(R_, 1), np.round([x0, y0_, x1, y1], 1))
print("axes heights (in)", [round(a_.get_window_extent(rr).height / fig.dpi, 3) for a_ in (vx, ax, bx)], "H", round(H, 3), "tight W", round(fig.get_tightbbox(rr).width, 3))
fig.savefig("../universe_thin_drawing.png", dpi=200, bbox_inches="tight", pad_inches=0.02); plt.close(fig)
plt.rcParams.update({"mathtext.fontset": "dejavusans", "axes.linewidth": 0.8, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.labelsize": 9, "axes.titlesize": 12})
AH = float(os.environ.get("AH", "7.58"))
Af = AFfour(TH1); th2 = np.append(TH1, 360); Af2 = AFfour(th2)
fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, AH), sharex=True, gridspec_kw={"height_ratios": [2.0, 1], "hspace": 0.06})
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
a.plot(th2, Af2, color="C0", lw=1.5, label="$A_F$")
for pp, col, ls, lab in [(P1, "C3", "-", "Newtonian thin shell + core, $\\Omega_s$=%.2f"), (P4, "C2", ":", "Newtonian shell + core, $\\Omega_s$=%.3f forced"), (P3, "C1", "--", "relativistic (LTB) shell + core, $\\Omega_s$=%.3f, $z$=0.3")]:
    A = AC.curve("ptshell", pp["d"]) if "d" in pp else pp["cur"]; assert abs(r2(MA, A[IMIG]) - pp["R2"]) < 2e-3, (r2(MA, A[IMIG]), pp["R2"])
    A2 = np.append(A, A[0])
    a.plot(th2, A2, color=col, ls=ls, lw=1.6 if ls != "--" else 1.3, label=(lab % pp["Os"]) + ", $\\ell_0$=%.1f° (fit $R^2$=%.3f; vs $A_F$ %.3f)" % (pp["l0"], pp["R2"], pp["R2F"]))
    b.plot(th2, 100 * (A2 - Af2) / Af2, color=col, ls=ls, lw=1.4 if ls == "-" else 1.2, label={"-": "Newtonian", ":": "Newt. $\\Omega_s$=0.049", "--": "relativistic"}[ls])
a.set_ylabel("$A/A_{\\rm all}$"); a.set_ylim(0.62, 1.22); a.legend(fontsize=7.2, loc="lower left", frameon=False)
b.axhline(0, color="k", lw=0.6); b.plot(ML, 100 * (MA - AFfour(ML)) / AFfour(ML), "o", ms=2.5, color="0.45", label="data")
b.set_ylabel("$100\\,(\\cdot-A_F)/A_F$ [%]", fontsize=8); b.set_xlabel("Galactic longitude $\\ell$ (deg)")
b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45)); b.legend(fontsize=7, frameon=False, loc="upper right", ncol=4)
fig.savefig("../universe_thin_Achart.png", dpi=200, bbox_inches="tight")
for f in ["universe_thin_drawing", "universe_thin_Achart"]:
    im = Image.open("../" + f + ".png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../" + f + "_web.jpg", quality=88)
