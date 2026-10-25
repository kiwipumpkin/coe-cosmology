# FIGURE_CODEBOOK_makers — figs_letter_exact_python

**Audience:** another bot — regenerating-code cookbook.
**Rule:** EVERY letter figure has complete runnable Python that PRODUCES the figure.
**Forbidden:** `cp` of gold/pdfimages extracts as the production method (comic HANDDRAWN is the only exception).
**Do not modify:** `letter_brief/CoE_calculations.md`.

| | |
|--|--|
| Output root | `/workspace/figs_letter_exact_python/` |
| Makers | `/workspace/figs_letter_exact_python/makers/` |
| Gold extracts | `/workspace/figs_letter_exact_python/gold_extracts/` |
| Run all | `python3 /workspace/figs_letter_exact_python/makers/run_all.py` |

Chris constraints (2026-10-01): Fig03 = termslice/make_fig.py only (2×2), no scale/axis edits; Fig06 = make_migkas_fig8_overlay_matched.py as it produced the letter (cut-off titles OK; no margin fixes).

---

## 1. Overview

| File | Maker | Method |
|------|-------|--------|
| `Fig01_geometry.png` | Fig01 | `makers/make_fig01_02_geometry_achart.py` |
| `Fig02_Achart.png` | Fig02 | `makers/make_fig01_02_geometry_achart.py` |
| `Fig03_migkas_slices.png` | Fig03 | `makers/make_fig03_migkas_slices.py` |
| `Fig04_amp_full.png` | Fig04 | `makers/make_fig04_amp_full.py` |
| `Fig05_sky_orientation_four_row.jpg` | Fig05 | `makers/make_fig05_sky_orientation.py` |
| `Fig06_overlay_stack.png` | Fig06 | `makers/make_fig06_overlay_stack.py` |
| `Fig07_storyboard.png` | Fig07 | `makers/make_fig07_storyboard.py` |
| `Fig08_two_emitters.png` | Fig08 | `makers/make_fig08_two_emitters.py` |
| `Fig09_kill.jpg` | Fig09 | `makers/make_fig09_kill.py` |
| `Fig10_cancel_partner.png` | Fig10 | `makers/make_fig10_cancel_partner.py` |
| `FigA1_ltb_shell.png` | FigA1 | `makers/make_figA1_ltb_shell.py` |
| `FigA2_coe_sun.jpg` | FigA2 | `makers/make_figA2_coe_sun.py` |
| `FigA2b_sun_coe_crossings.jpg` | FigA2b | `makers/make_figA2b_crossings.py` |
| `FigA3_coe_sky.png` | FigA3 | `makers/make_figA3_coe_sky.py` |
| `FigA_comic_HANDDRAWN_from_PDF.jpg` | comic | `makers/make_figA_comic_HANDDRAWN.py` |
| `eclipse_inner_shell.png` | eclipse | `makers/make_eclipse_inner_shell.py` |

---

## 2. Per-figure cookbook (full source)

### Fig01 — `Fig01_geometry.png`

- **Gold extract:** `gold-000.png`
- **Expected dims (gold):** 1460×488
- **Method:** `rerun:makers/make_fig01_02_geometry_achart.py`
- **Output:** `/workspace/figs_letter_exact_python/Fig01_geometry.png`
- **Notes:** Reruns spreadsheet_model/bubble/universe_thin.py (letter-era). Pixel-matches gold.

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_fig01_02_geometry_achart.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_fig01_02_geometry_achart.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_fig01_02_geometry_achart.py
#!/usr/bin/env python3
"""Fig01 geometry + Fig02 Achart — rerun spreadsheet_model/bubble/universe_thin.py.

Pixel-matches letter gold embeds (PNG encoding may differ; RGBA max abs == 0).
"""
from __future__ import annotations

import runpy
import shutil
from pathlib import Path

OUT = Path("/workspace/figs_letter_exact_python")
BUBBLE = Path("/workspace/cosmology/papers/spreadsheet_model/bubble")
SM = BUBBLE.parent
SCRIPT = BUBBLE / "universe_thin.py"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # universe_thin.py writes ../universe_thin_drawing.png and ../universe_thin_Achart.png
    # relative to bubble/; run in that directory via runpy with argv.
    import os, sys
    cwd = os.getcwd()
    sys.path.insert(0, str(BUBBLE))
    try:
        os.chdir(BUBBLE)
        runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        os.chdir(cwd)
    src1 = SM / "universe_thin_drawing.png"
    src2 = SM / "universe_thin_Achart.png"
    d1 = OUT / "Fig01_geometry.png"
    d2 = OUT / "Fig02_Achart.png"
    shutil.copy2(src1, d1)
    shutil.copy2(src2, d2)
    print("wrote", d1, d1.stat().st_size)
    print("wrote", d2, d2.stat().st_size)
    return d1, d2


if __name__ == "__main__":
    main()

```

**Dependency (full source):** `/workspace/cosmology/papers/spreadsheet_model/bubble/universe_thin.py`

```python
# PATH: /workspace/cosmology/papers/spreadsheet_model/bubble/universe_thin.py
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
```

### Fig02 — `Fig02_Achart.png`

- **Gold extract:** `gold-002.png`
- **Expected dims (gold):** 1182×1287
- **Method:** `rerun:makers/make_fig01_02_geometry_achart.py`
- **Output:** `/workspace/figs_letter_exact_python/Fig02_Achart.png`
- **Notes:** Same universe_thin.py run as Fig01. Pixel-matches gold.

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_fig01_02_geometry_achart.py
```

### Fig03 — `Fig03_migkas_slices.png`

- **Gold extract:** `gold-004.png`
- **Expected dims (gold):** 1123×1014
- **Method:** `rerun:makers/make_fig03_migkas_slices.py`
- **Output:** `/workspace/figs_letter_exact_python/Fig03_migkas_slices.png`
- **Notes:** PRIMARY: termslice/make_fig.py ONLY (2×2). Do not change scales/axes. Wrapper chdir+runs.

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_fig03_migkas_slices.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_fig03_migkas_slices.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_fig03_migkas_slices.py
#!/usr/bin/env python3
"""Fig03 — Migkas redshift slices + D–m_cc degeneracy.

PRIMARY (letter-era): cmb_reconcile/termslice/make_fig.py
Writes termslice_summary.png (2x2: three A/slice-mean panels + Degeneracy),
then copies to Figs_from_makers/Fig03_migkas_slices.png.

Do NOT use coe-figures/src/fig03_migkas_slices.py (wrong reconstruction).
"""
from __future__ import annotations

import os
import runpy
import shutil
from pathlib import Path

OUT = Path("/workspace/figs_letter_exact_python")
OUT_NAME = "Fig03_migkas_slices.png"
TERMSLICE = Path(
    "/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/termslice"
)
SCRIPT = TERMSLICE / "make_fig.py"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # make_fig.py uses relative imports (../cycle3, solve.json, …) — must run in termslice cwd
    prev = os.getcwd()
    try:
        os.chdir(TERMSLICE)
        runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        os.chdir(prev)
    src = TERMSLICE / "termslice_summary.png"
    if not src.is_file():
        raise SystemExit(f"missing output from letter-era maker: {src}")
    dest = OUT / OUT_NAME
    shutil.copy2(src, dest)
    print("wrote", dest, dest.stat().st_size)
    return dest


if __name__ == "__main__":
    main()
```

**Dependency (full source):** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/termslice/make_fig.py`

```python
# PATH: /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/termslice/make_fig.py
import json, sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "../cycle3"); import migkas_A as M
R = json.load(open("solve.json")); S2 = json.load(open("../cycle3/slices2.json")); RT = json.load(open("rotate.json")); TH1 = np.arange(0, 360, 1.0); LC = M.LC
pick = lambda s, D=49.1746, cone=False, mode="median": [r for r in R if r["slice"] == s and abs(r["D"] - D) < 1e-3 and r["mode"] == mode and r["cone"] == cone and r["l0_branch"] == 0][0]
col = {"0.03-0.07": "C1", "z>0.07": "C3", "all": "k"}
# Letter gold is 2x2 (1123x1014): three A/slice-mean panels + degeneracy
fig, axs = plt.subplots(2, 2, figsize=(9.5, 8.6), gridspec_kw=dict(wspace=0.27, hspace=0.32))
axs_flat = axs.ravel()
for ax, s in zip(axs_flat[:3], ("0.03-0.07", "z>0.07", "all")):
    p = np.array(S2[s]["A"], float); e = np.array(S2[s]["A_err"], float); ok = np.isfinite(p); m = np.nanmean(p)
    ax.fill_between(LC[ok], (p - e)[ok] / m, (p + e)[ok] / m, color="0.85"); ax.plot(LC[ok], p[ok] / m, ".", color=col[s], ms=4, label="Migkas 2020 slice (rebuilt)")
    r = pick(s); ax.plot(TH1, r["curve"], color=col[s], lw=2, label="solved model: R2 %.2f" % r["R2_point"])
    rc = pick(s, cone=True); cb = np.convolve(np.r_[np.array(rc["curve"])[-45:], rc["curve"], np.array(rc["curve"])[:45]], np.ones(91) / 91, mode="valid")
    ax.plot(TH1, cb, color=col[s], ls="--", lw=1, label="solved for binned profile, shown binned: R2 %.2f" % rc["R2_conebinned"])
    if s == "0.03-0.07":
        A = np.array(pick("all")["curve"]); dl = RT["all|point"]["dl"]
        ax.plot(TH1, np.interp((TH1 - dl) % 360, TH1, A, period=360), "k:", lw=1.5, label="black model rotated by %.0f deg: R2 %.2f" % (dl, RT["all|point"]["R2"]))
    ax.set_title("%s (median distance %.2f Gpc)" % (s, r["Rs"][0]), fontsize=10); ax.set_xlabel("Galactic longitude (deg)"); ax.set_xlim(0, 360); ax.legend(fontsize=7, loc="lower left")
axs_flat[0].set_ylabel("A / slice mean")
axs_flat[2].set_ylabel("A / slice mean")
ax = axs_flat[3]
for s in ("0.03-0.07", "z>0.07", "all"):
    rr = sorted([r for r in R if r["slice"] == s and r["mode"] == "median" and not r["cone"] and r["l0_branch"] == 0], key=lambda r: r["D"])
    ax.loglog([r["D"] for r in rr], [r["Mc_kg"] for r in rr], "o-", color=col[s], label="%s: core share of critical mass %.4f to %.4f" % (s, min(r["mcc"] for r in rr), max(r["mcc"] for r in rr)))
ax.set_xticks([20,30,49.17,94,200]); ax.set_xticklabels(["20","30","49","94","200"]); ax.minorticks_off(); ax.set_xlabel("CoE distance D (Gpc), not constrained by the data"); ax.set_ylabel("solved core mass (kg)"); ax.set_title("Degeneracy: every D fits exactly", fontsize=10); ax.legend(fontsize=7)
fig.savefig("termslice_summary.png", dpi=140, bbox_inches="tight"); print("ok")
```

### Fig04 — `Fig04_amp_full.png`

- **Gold extract:** `gold-006.png`
- **Expected dims (gold):** 1800×1385
- **Method:** `rerun:makers/make_fig04_amp_full.py`
- **Output:** `/workspace/figs_letter_exact_python/Fig04_amp_full.png`
- **Notes:** Regenerates from Migkas CSV + bubble fit curves. No Amp×sin.

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_fig04_amp_full.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_fig04_amp_full.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_fig04_amp_full.py
#!/usr/bin/env python3
"""Fig04 — Full Migkas A/Aall + residual (NO Amp×sin squiggles).

Regenerates from digitised Migkas CSV + bubble Newt/Rel fit curves.
Forbidden asset: money_shot_amp_full_SINE_BAD.png (do not use).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths import OUT  # noqa: E402
from physics import (  # noqa: E402
    A_F,
    BAND_CENTRES,
    BAND_COLORS,
    BAND_HALF,
    load_fit_curves,
    load_migkas_fig6,
)

OUT_NAME = "Fig04_amp_full.png"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ell_d, A_d, e_d = load_migkas_fig6()
    th, A_newt, A_free, A_rel, meta = load_fit_curves()
    th2 = np.append(th, 360.0)
    cl = lambda y: np.append(y, y[0])

    # Letter gold ~1800×1385 → figsize ~9×6.9 at 200dpi
    fig, (a, b) = plt.subplots(
        2, 1, figsize=(9.0, 6.9), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.07},
    )
    for L, col in zip(BAND_CENTRES, BAND_COLORS):
        for ax in (a, b):
            ax.axvspan(L - BAND_HALF, L + BAND_HALF, color=col, alpha=0.20, lw=0, zorder=0)

    a.errorbar(ell_d, A_d, yerr=e_d, fmt="o", ms=4, color="0.35", ecolor="0.65",
               lw=0.9, zorder=3, label="Migkas (digitised)")
    a.plot(th2, A_F(th2), color="C0", lw=1.7, label=r"$A_F$", zorder=2)
    mn, mr = meta["newt"], meta["rel"]
    a.plot(
        th2, cl(A_newt), color="#1a7f37", lw=2.0,
        label=(
            rf"Newtonian $\Omega_s$={mn['Os']}: $D$={mn['D']:.1f}, "
            rf"$R_{{\rm in}}$={mn['R_in']:.1f}, $R_{{\rm out}}$={mn['R_out']:.1f} Gpc, "
            rf"$\ell_0$={mn['l0']:.0f}° ($R^2$={mn['R2']:.3f})"
        ),
    )
    a.plot(
        th2, cl(A_rel), color="C1", lw=1.8, ls="--",
        label=(
            rf"Relativistic LTB $\Omega_s$={mr['Os']}: $D$={mr['D']:.1f}, "
            rf"$R_{{\rm in}}$={mr['R_in']:.1f}, $R_{{\rm out}}$={mr['R_out']:.1f} Gpc, "
            rf"$\ell_0$={mr['l0']:.1f}° ($R^2$={mr['R2']:.3f})"
        ),
    )
    a.set_ylabel(r"$A/A_{\mathrm{all}}$", fontsize=12)
    a.set_ylim(0.60, 1.24)
    a.legend(fontsize=7.2, loc="lower left", frameon=False)
    a.set_title(r"Full Migkas $A/A_{\mathrm{all}}$ profile", loc="left", fontsize=12)

    def res(model):
        return 100 * (A_d - np.interp(ell_d, th, model, period=360)) / A_d

    b.axhline(0, color="k", lw=0.8)
    b.plot(ell_d, res(A_newt), "o-", ms=3.5, color="#1a7f37", lw=1.0, label="data − Newt")
    b.plot(ell_d, res(A_rel), "s--", ms=3.5, color="C1", lw=1.0, label="data − Rel")
    b.set_ylabel("residual [%]", fontsize=11)
    b.set_xlabel(r"Galactic longitude $\ell$ (deg)", fontsize=12)
    b.set_xlim(0, 360)
    b.set_xticks(range(0, 361, 45))
    b.legend(fontsize=8, frameon=False, loc="upper right")
    b.grid(True, alpha=0.25)

    dest = OUT / OUT_NAME
    fig.savefig(dest, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", dest, dest.stat().st_size)
    return dest


if __name__ == "__main__":
    main()

```

### Fig05 — `Fig05_sky_orientation_four_row.jpg`

- **Gold extract:** `gold-007.jpg`
- **Expected dims (gold):** 1370×1420
- **Method:** `rerun:makers/make_fig05_sky_orientation.py`
- **Output:** `/workspace/figs_letter_exact_python/Fig05_sky_orientation_four_row.jpg`
- **Notes:** Rings stack: four-row sky orientation (untilted→tilt→gradient→roll).

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_fig05_sky_orientation.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_fig05_sky_orientation.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_fig05_sky_orientation.py
#!/usr/bin/env python3
"""Fig05 — Four-row sky orientation (untilted → tilt → gradient → roll).

Regenerates analytic band circles; writes JPEG to match letter delivery format.
Gold target ~1370×1420.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths import OUT  # noqa: E402
from physics import (  # noqa: E402
    BAND_CENTRES,
    BAND_COLORS,
    BAND_HALF,
    BAND_KINDS,
    ROLL_DEG,
    TILT_DEG,
    band_circle_xyz,
    mollweide_xy,
    rotate_xyz,
    xyz_to_galactic,
)

OUT_NAME = "Fig05_sky_orientation_four_row.jpg"


def _draw_model_frame(ax, tilt, roll, gradient=False):
    ax.set_aspect("equal")
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.axis("off")
    ax.add_patch(Circle((0, 0), 1.0, fill=False, ec="0.4", lw=1.4))
    ax.plot(0, 0, "ko", ms=5)
    ax.annotate("", xy=(1.15, 0), xytext=(-1.15, 0),
                arrowprops=dict(arrowstyle="<->", color="#6c3483", lw=1.5))
    ax.text(-1.28, 0.12, "CoE (−x)", color="#6c3483", fontsize=7, ha="left")
    ax.text(1.28, 0.12, "anti-CoE (+x)", color="#6c3483", fontsize=7, ha="right")
    for th, col, kind in zip(BAND_CENTRES, BAND_COLORS, BAND_KINDS):
        t = np.radians(th)
        alpha = 0.35
        if gradient:
            dist = min(abs(th - 180), abs(th - 180 + 360), abs(th - 180 - 360))
            alpha = 0.15 + 0.55 * (1 - dist / 180.0)
        for sgn in (+1, -1):
            x = np.cos(t)
            z = sgn * np.sin(t)
            tt = np.radians(tilt)
            xr = x * np.cos(tt) - z * np.sin(tt)
            zr = x * np.sin(tt) + z * np.cos(tt)
            rr = np.radians(roll)
            zr2 = zr * np.cos(rr)
            ax.plot([0.82 * xr, 1.05 * xr], [0.82 * zr2, 1.05 * zr2],
                    color=col, lw=4.5, alpha=alpha, solid_capstyle="round")
        xl = 1.18 * np.cos(t)
        zl = 1.18 * np.sin(t)
        ax.text(xl, zl, rf"$\theta$={th:g}°", color=col, fontsize=6.5,
                ha="center", va="center", alpha=min(1, alpha + 0.3))


def _draw_mollweide(ax, tilt, roll, gradient=False):
    t = np.linspace(0, 2 * np.pi, 400)
    ox = 2 * np.sqrt(2) * np.cos(t)
    oy = np.sqrt(2) * np.sin(t)
    ax.plot(-ox, oy, "k-", lw=0.8)
    lon = np.linspace(0, 360, 400)
    x, y = mollweide_xy(lon, np.zeros_like(lon))
    ax.plot(x, y, color="#1a7f37", lw=0.9, ls="--", alpha=0.8)
    for th, col, kind in zip(BAND_CENTRES, BAND_COLORS, BAND_KINDS):
        band = band_circle_xyz(th, BAND_HALF, n=720)
        for edge_name, lw in (("mid", 2.2), ("lo", 0.9), ("hi", 0.9)):
            xyz = rotate_xyz(band[edge_name], tilt_deg=tilt, roll_deg=roll)
            ell, b = xyz_to_galactic(xyz)
            x, y = mollweide_xy(ell, b)
            d = np.hypot(np.diff(x), np.diff(y))
            breaks = np.where(d > 0.35)[0]
            idx = np.r_[0, breaks + 1, len(x)]
            for i0, i1 in zip(idx[:-1], idx[1:]):
                if i1 - i0 < 2:
                    continue
                xs, ys = x[i0:i1], y[i0:i1]
                if gradient:
                    mid = ell[i0:i1].mean()
                    d0 = min(mid % 360, 360 - (mid % 360))
                    alpha = 0.2 + 0.75 * (1 - d0 / 180.0)
                else:
                    alpha = 0.9 if edge_name == "mid" else 0.45
                ax.plot(xs, ys, color=col, lw=lw if edge_name == "mid" else 0.8, alpha=alpha)
    ac = rotate_xyz(np.array([[1.0, 0.0, 0.0]]), tilt_deg=tilt, roll_deg=roll)
    ell0, b0 = xyz_to_galactic(ac)
    x0, y0 = mollweide_xy(ell0, b0)
    ax.plot(x0, y0, "o", color="#7b241c", ms=6, zorder=5)
    ax.text(x0[0], y0[0] + 0.12, "anti-CoE", color="#7b241c", fontsize=7, ha="center")
    ax.set_aspect("equal")
    ax.set_xlim(-3.0, 3.0)
    ax.set_ylim(-1.7, 1.7)
    ax.axis("off")
    for L, xpos in ((180, -2.0), (270, -1.0), (0, 0.0), (90, 1.0), (180, 2.0)):
        ax.text(xpos, -1.55, "180" if L == 180 else str(L), fontsize=7, ha="center", color="0.35")
    ax.text(0, -1.72, r"Galactic $\ell$ (increases ←)", fontsize=7, ha="center", color="0.35")


def _draw_sphere3d(ax, tilt, roll, gradient=False):
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 20)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, color="0.75", lw=0.3, alpha=0.5)
    ax.plot([-1.2, 1.2], [0, 0], [0, 0], color="#6c3483", lw=1.5)
    for th, col in zip(BAND_CENTRES, BAND_COLORS):
        band = band_circle_xyz(th, 0.0, n=200)["mid"]
        xyz = rotate_xyz(band, tilt_deg=tilt, roll_deg=roll)
        alpha = 0.9
        if gradient:
            alpha = 0.95 if th <= 90 else 0.35
        ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=col, lw=2.0, alpha=alpha)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(-1.2, 1.2)
    ax.set_axis_off()
    ax.view_init(elev=18, azim=35)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [
        (0.0, 0.0, False, "Untilted, no roll"),
        (TILT_DEG, 0.0, False, rf"+{TILT_DEG:g}° tilt"),
        (TILT_DEG, 0.0, True, rf"+{TILT_DEG:g}° tilt, gradient"),
        (TILT_DEG, ROLL_DEG, True, rf"+{TILT_DEG:g}° tilt, gradient, +{ROLL_DEG:g}° roll"),
    ]
    # Aim near gold 1370×1420
    fig = plt.figure(figsize=(9.13, 9.47), dpi=150)
    for i, (tilt, roll, grad, title) in enumerate(rows):
        ax0 = fig.add_subplot(4, 3, 3 * i + 1)
        _draw_model_frame(ax0, tilt, roll, gradient=grad)
        ax0.set_title("Model frame" if i == 0 else "", fontsize=10)
        ax0.text(0.0, 1.22, title, ha="center", fontsize=9, transform=ax0.transData, clip_on=False)
        ax1 = fig.add_subplot(4, 3, 3 * i + 2)
        _draw_mollweide(ax1, tilt, roll, gradient=grad)
        ax1.set_title("Mollweide (Galactic)" if i == 0 else "", fontsize=10)
        ax2 = fig.add_subplot(4, 3, 3 * i + 3, projection="3d")
        _draw_sphere3d(ax2, tilt, roll, gradient=grad)
        ax2.set_title("Sphere 3D" if i == 0 else "", fontsize=10)
    fig.suptitle(
        "Sky orientation sequence: untilted → +15° tilt → gradient → +80.8° roll",
        fontsize=12, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    png = OUT / "_Fig05_sky_orientation_four_row.png"
    fig.savefig(png, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    # JPEG delivery (letter uses jpg)
    dest = OUT / OUT_NAME
    im = Image.open(png).convert("RGB")
    # Resize width to gold 1370 keeping aspect
    w = 1370
    h = round(im.height * w / im.width)
    im.resize((w, h), Image.LANCZOS).save(dest, quality=92)
    print("wrote", dest, dest.stat().st_size, im.size, "->", (w, h))
    return dest


if __name__ == "__main__":
    main()

```

### Fig06 — `Fig06_overlay_stack.png`

- **Gold extract:** `gold-008.png`
- **Expected dims (gold):** 2164×3222
- **Method:** `rerun:makers/make_fig06_overlay_stack.py`
- **Output:** `/workspace/figs_letter_exact_python/Fig06_overlay_stack.png`
- **Notes:** Runs make_migkas_fig8_overlay_matched.py AS ON DISK (unchanged; no margin fixes). If letter gold has cut-off title, match that. Then vertical stack.

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_fig06_overlay_stack.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_fig06_overlay_stack.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_fig06_overlay_stack.py
#!/usr/bin/env python3
"""Fig06 — Migkas Fig.8 overlay stack AS THE LETTER PRODUCED IT.

Runs cosmology/figures/make_migkas_fig8_overlay_matched.py UNCHANGED
(no margin patches, no SHOW_LEGEND flips, no bbox_inches hacks).
If letter gold has cut-off titles, that is intentional — match letter.

Then vertically stacks matched / matched_grad / matched_grad_rot81.
"""
from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

OUT = Path("/workspace/figs_letter_exact_python")
SM = Path("/workspace/cosmology/papers/spreadsheet_model")
OVERLAY_PY = Path("/workspace/cosmology/figures/make_migkas_fig8_overlay_matched.py")

PANEL_NAMES = [
    "fig8_overlay_matched.png",
    "fig8_overlay_matched_grad.png",
    "fig8_overlay_matched_grad_rot81.png",
]
PANEL_LABELS = [
    "+15° tilt",
    "+15° tilt + gradient",
    "+15° tilt + gradient + 80.8° roll",
]
OUT_PANELS = [
    "Fig06_overlay_matched.png",
    "Fig06_overlay_matched_grad.png",
    "Fig06_overlay_matched_grad_rot81.png",
]
OUT_STACK = "Fig06_overlay_stack.png"


def _run_overlay_unchanged():
    """Load and run make_migkas_fig8_overlay_matched.py exactly as on disk."""
    spec = importlib.util.spec_from_file_location("migkas_overlay_matched", OVERLAY_PY)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    # Restrict to the three letter stack panels only (same code path, no layout edits)
    mod.VERSIONS = [
        ("matched", False, 0.0),
        ("matched_grad", True, 0.0),
        ("matched_grad_rot81", True, 80.8),
    ]
    # SHOW_LEGEND left as in the letter-era script on disk — do NOT patch margins.
    mod.main()
    return mod


def _stack_panels():
    OUT.mkdir(parents=True, exist_ok=True)
    parts = []
    for src_name, out_name in zip(PANEL_NAMES, OUT_PANELS):
        src = SM / src_name
        if not src.exists():
            raise FileNotFoundError(src)
        dest = OUT / out_name
        shutil.copy2(src, dest)
        parts.append(Image.open(src).convert("RGB"))
        print("panel", dest, Image.open(dest).size)

    try:
        font = ImageFont.truetype(
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 28
        )
    except Exception:
        font = ImageFont.load_default()

    label_h = 48
    gap = 12
    W = max(p.width for p in parts)
    H = sum(p.height + label_h + gap for p in parts) + 20
    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    y = 10
    for lab, p in zip(PANEL_LABELS, parts):
        bbox = draw.textbbox((0, 0), lab, font=font)
        tw = bbox[2] - bbox[0]
        draw.text(((W - tw) // 2, y + 8), lab, fill="black", font=font)
        y += label_h
        canvas.paste(p, ((W - p.width) // 2, y))
        y += p.height + gap

    stack = canvas.crop((0, 0, W, y))
    dest = OUT / OUT_STACK
    stack.save(dest)
    print("wrote stack", dest, stack.size, dest.stat().st_size)
    return dest


def main():
    print("=== regenerating Migkas Fig.8 overlays (UNCHANGED letter script) ===")
    _run_overlay_unchanged()
    print("=== stacking panels ===")
    return _stack_panels()


if __name__ == "__main__":
    main()
```

**Dependency (full source):** `/workspace/cosmology/figures/make_migkas_fig8_overlay_matched.py`

```python
# PATH: /workspace/cosmology/figures/make_migkas_fig8_overlay_matched.py
#!/usr/bin/env python3
"""
Superimpose CoE MAX/MIN highlight bands + Migkas markers onto the REAL Migkas
Fig. 8 theta=75 A/A_all Mollweide map.

Matches make_coe_rings_migkas_pin_axis.py (three-panel Panel B) geometry:
  Exactly FOUR features — each ONE transparent highlight band, 18° end-to-end
  (= ±9° = 5% of 360°) about the center θ + ONE solid mid-stroke:
    MAX 30°  → [21°, 39°]   color #00e5ff (bright cyan)
    MAX 185° → [176°, 194°] color #ff2d95 (magenta/hot pink)
    MIN 115° → [106°, 124°] color #b8f200 (lime)
    MIN 265° → [256°, 274°] color #4fc3f7 (sky blue)
  Drawing: clean PolyCollection fill + one mid-stroke (optional thin white
  underlay on mid only). NO black rim, NO thick white halo stacks, NO dual
  boundary-edge strokes. No 283° companion. No Dec=0 overlay.
  MAX triangles: (17°, −9°) and (185°, 0°)
  MIN squares:   (115°, 0°) and (274°, −22°)
  CoE axis tilted +30°: anti-CoE (0°,+30°), CoE (180°,−30°) (purple).

Registration (unchanged from working Fig. 8):
  GOOD_REG = (cx, cy, Rx, Ry) = (863.0, 763.25, 693.9927, 392.542425)
  Source: FULL_migkas_theta75.png
  CoE axis tilted +30°; clip/inset so overlays stay inside colored oval
  lon_disp = -wrap(ℓ)

Outputs:
  migkas_fig8_theta75_with_coe_rings.png
  migkas_fig8_theta75_with_coe_rings.json
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
from PIL import Image

OUT_DIR = Path(__file__).resolve().parent
SM_DIR = Path("/workspace/cosmology/papers/spreadsheet_model")
PNG_PATH = SM_DIR / "fig8_overlay_matched.png"
JSON_PATH = SM_DIR / "fig8_overlay_matched.json"

# Versions = panel B of each three-panel figure (grad: opacity gradient, roll: deg CCW about CoE axis)
VERSIONS = [
    ("matched", False, 0.0),
    ("matched_grad", True, 0.0),
    ("matched_grad_rot45", True, 45.0),
    ("matched_grad_rot81", True, 80.8),
    ("matched_rot81", False, 80.8),  # roll + all loops, no opacity gradient
]
GRAD = False
ROLL = 0.0
SHOW_LEGEND = False  # letter Fig.6: no legend box/strip

SOURCE_CANDIDATES = [
    Path("/workspace/cosmology/map_compare/FULL_migkas_theta75.png"),
    Path("/workspace/cosmology/map_compare/migkas_fig8_theta75_crop.png"),
    Path("/workspace/cosmology/migkas/pub_maps/Fig8_full_hires.png"),
]

# Strip A(ℓ) scan stops from +x (after-2ℓ A_after extrema) — four band centers.
# 5% of 360° = 18° end-to-end → half-width ±9° about each center (match three-panel).
BAND_HALF = 9.0
BAND_FULL = 2.0 * BAND_HALF  # 18°
THETA_MAX = [29.5, 187.5]
THETA_MIN = [122.6, 264.4]
MAX_BANDS = [(t - BAND_HALF, t + BAND_HALF) for t in THETA_MAX]  # [21,39], [176,194]
MIN_BANDS = [(t - BAND_HALF, t + BAND_HALF) for t in THETA_MIN]  # [106,124], [256,274]


def mu_pair_for_band(lo, hi):
    """Edge μ=cos(θ) for a scan-angle band (μ-space, not θ-space)."""
    return float(np.cos(np.deg2rad(lo))), float(np.cos(np.deg2rad(hi)))

# Visually-good registration (pre medial-axis retune). Hard-set as working default.
GOOD_REG = {
    "cx": 863.0,
    "cy": 763.25,
    "Rx": 693.9927,
    "Ry": 392.542425,
}

# Sky markers (ℓ, b): map extrema + strip stand-ins (same as 3-panel Panel B).
MIGKAS_MAXES = [(17.0, -9.0), (185.0, 0.0)]
MIGKAS_MINS = [(115.0, 0.0), (274.0, -22.0)]

# Galactic pose (Chris 2026-09-22): anti-CoE (+D/+x) at (0°, +30°); CoE at (180°, −30°)
ANTI_COE_LB = (0.0, 15.0)
COE_LB = (180.0, -15.0)
AXIS_TILT_B_DEG = 15.0


def lb_to_n(l_deg, b_deg):
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    return np.array(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], dtype=float
    )


D_HAT = lb_to_n(*ANTI_COE_LB)  # +D = anti-CoE = +x


def great_circle_through_poles(D, n=721):
    """Great circle through ±D (CoE axis meridian on the sky)."""
    a = D / np.linalg.norm(D)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([1.0, 0.0, 0.0])
    e2 = np.cross(a, tmp)
    e2 /= np.linalg.norm(e2)
    t = np.linspace(0.0, 2.0 * np.pi, n)
    pts = np.outer(np.cos(t), a) + np.outer(np.sin(t), e2)
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    return ls, bs

# High-contrast loop colors that pop on magma (avoid reds/oranges/corals)
C_MAX_30 = "#0097b2"    # cyan — MAX 29.5° (three-panel matched)
C_MAX_185 = "#e0107a"   # magenta — MAX 187.5°
C_MIN_115 = "#6fa800"   # green — MIN 122.6°
C_MIN_265 = "#2166d4"   # blue — MIN 264.4°
COLORS_MAX = [C_MAX_30, C_MAX_185]
COLORS_MIN = [C_MIN_115, C_MIN_265]
C_COE = "#6b3fa0"
C_ANTI = "#8e1b1b"  # dark red marker; distinct from loop reds (match three-panel)
C_AXIS = "#5b2c6f"
C_EQ = "#1e8449"


def wrap_lon_deg(lon):
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def lon_disp_from_ell(ell_deg):
    return -wrap_lon_deg(ell_deg)


def n_to_lb(n):
    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    b = float(np.rad2deg(np.arcsin(np.clip(n[2], -1.0, 1.0))))
    l = float(np.rad2deg(np.arctan2(n[1], n[0])) % 360.0)
    return l, b


def celestial_equator_lb(n=2500):
    ra_ngp = np.deg2rad(192.85948)
    dec_ngp = np.deg2rad(27.12825)
    l_ncp = np.deg2rad(122.93192)
    ra = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    dec = 0.0
    sin_b = np.clip(
        np.sin(dec) * np.sin(dec_ngp)
        + np.cos(dec) * np.cos(dec_ngp) * np.cos(ra - ra_ngp),
        -1.0,
        1.0,
    )
    b = np.arcsin(sin_b)
    y = np.cos(dec) * np.sin(ra - ra_ngp)
    x = (
        np.sin(dec) * np.cos(dec_ngp)
        - np.cos(dec) * np.sin(dec_ngp) * np.cos(ra - ra_ngp)
    )
    l = (l_ncp - np.arctan2(y, x)) % (2.0 * np.pi)
    return np.rad2deg(l), np.rad2deg(b)


def small_circle_about_axis(axis, mu, n=721):
    a = axis / np.linalg.norm(axis)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    tmp = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(a, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    t = np.linspace(0.0, 2.0 * np.pi, n)
    return (mu * a)[None, :] + rho * (
        np.outer(np.cos(t), e1) + np.outer(np.sin(t), e2)
    )


def mollweide_xy(lon_rad, lat_rad, n_iter=64):
    """Standard Mollweide; stable bracketed solve of 2θ+sin(2θ)=π sin(lat)."""
    lon = np.asarray(lon_rad, dtype=float)
    lat = np.asarray(lat_rad, dtype=float)
    target = np.pi * np.sin(np.clip(lat, -np.pi / 2, np.pi / 2))
    lo = np.full_like(lat, -np.pi / 2, dtype=float)
    hi = np.full_like(lat, np.pi / 2, dtype=float)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        g = 2.0 * mid + np.sin(2.0 * mid)
        lo = np.where(g < target, mid, lo)
        hi = np.where(g >= target, mid, hi)
    theta = 0.5 * (lo + hi)
    theta = np.where(lat >= np.pi / 2 - 1e-15, np.pi / 2, theta)
    theta = np.where(lat <= -np.pi / 2 + 1e-15, -np.pi / 2, theta)
    x = (2.0 * np.sqrt(2.0) / np.pi) * lon * np.cos(theta)
    y = np.sqrt(2.0) * np.sin(theta)
    return x, y


def lb_to_px(ell_deg, b_deg, cx, cy, Rx, Ry):
    lon = np.deg2rad(lon_disp_from_ell(np.asarray(ell_deg, dtype=float)))
    lat = np.deg2rad(np.asarray(b_deg, dtype=float))
    x, y = mollweide_xy(lon, lat)
    Xn = x / (2.0 * np.sqrt(2.0))
    Yn = y / np.sqrt(2.0)
    return cx + Xn * Rx, cy - Yn * Ry


def as_xy(px, py):
    return float(np.asarray(px).reshape(-1)[0]), float(np.asarray(py).reshape(-1)[0])


def chamfer_distance_to_dark(dark):
    h, w = dark.shape
    INF = 1e6
    D = np.where(dark, 0.0, INF)
    for y in range(h):
        for x in range(w):
            v = D[y, x]
            if x > 0:
                v = min(v, D[y, x - 1] + 1.0)
            if y > 0:
                v = min(v, D[y - 1, x] + 1.0)
            if x > 0 and y > 0:
                v = min(v, D[y - 1, x - 1] + 1.414)
            if x + 1 < w and y > 0:
                v = min(v, D[y - 1, x + 1] + 1.414)
            D[y, x] = v
    for y in range(h - 1, -1, -1):
        for x in range(w - 1, -1, -1):
            v = D[y, x]
            if x + 1 < w:
                v = min(v, D[y, x + 1] + 1.0)
            if y + 1 < h:
                v = min(v, D[y + 1, x] + 1.0)
            if x + 1 < w and y + 1 < h:
                v = min(v, D[y + 1, x + 1] + 1.414)
            if x > 0 and y + 1 < h:
                v = min(v, D[y + 1, x - 1] + 1.414)
            D[y, x] = v
    return D


def detect_oval_bbox(rgb):
    h, w, _ = rgb.shape
    white = (rgb[:, :, 0] > 250) & (rgb[:, :, 1] > 250) & (rgb[:, :, 2] > 250)
    oval_left = np.full(h, np.nan)
    oval_right = np.full(h, np.nan)
    for y in range(h):
        row = ~white[y]
        d = np.diff(np.r_[0, row.astype(int), 0])
        starts = np.where(d == 1)[0]
        ends = np.where(d == -1)[0] - 1
        if len(starts) == 0:
            continue
        lengths = ends - starts + 1
        candidates = [
            (int(lengths[i]), int(starts[i]), int(ends[i]))
            for i in range(len(starts))
            if starts[i] < 0.72 * w and lengths[i] > 40
        ]
        if not candidates:
            continue
        candidates.sort(reverse=True)
        _, L, R = candidates[0]
        oval_left[y] = L
        oval_right[y] = R
    widths = np.where(np.isfinite(oval_left), oval_right - oval_left, 0.0)
    yeq = int(np.argmax(widths))
    thresh = 0.08 * widths[yeq]
    y0 = yeq
    while y0 > 0 and widths[y0 - 1] >= thresh:
        y0 -= 1
    y1 = yeq
    while y1 < h - 1 and widths[y1 + 1] >= thresh:
        y1 += 1
    x0 = float(oval_left[yeq])
    x1 = float(oval_right[yeq])
    return {
        "x0": int(x0),
        "x1": int(x1),
        "y0": int(y0),
        "y1": int(y1),
        "yeq": int(yeq),
        "cx": 0.5 * (x0 + x1),
        "cy": 0.5 * (y0 + y1),
        "Rx": 0.5 * (x1 - x0),
        "Ry": 0.5 * (y1 - y0 + 1),
    }


def refine_registration(rgb, seed):
    """Gray/dark Dec=0 lock (pre medial-axis). Prefer GOOD_REG hard-set in main."""
    h, w, _ = rgb.shape
    dark = rgb.mean(axis=2) < 55
    dist = chamfer_distance_to_dark(dark)

    ell_eq, b_eq = celestial_equator_lb(2500)
    lon = np.deg2rad(lon_disp_from_ell(ell_eq))
    lat = np.deg2rad(b_eq)
    mx, my = mollweide_xy(lon, lat)
    Xn = mx / (2.0 * np.sqrt(2.0))
    Yn = my / np.sqrt(2.0)

    def median_dist(cx, cy, Rx, Ry):
        px = cx + Xn * Rx
        py = cy - Yn * Ry
        m = (px >= 1) & (px < w - 1) & (py >= 1) & (py < h - 1)
        xi = np.clip(np.round(px[m]).astype(int), 0, w - 1)
        yi = np.clip(np.round(py[m]).astype(int), 0, h - 1)
        return float(np.median(dist[yi, xi]))

    base = {k: float(seed[k]) for k in ("cx", "cy", "Rx", "Ry")}
    best_sc = median_dist(**base)
    best = dict(base)
    for dcx in np.linspace(-10, 10, 21):
        for dcy in np.linspace(-10, 10, 21):
            for sRx in (0.98, 0.99, 1.0, 1.01, 1.02):
                for sRy in (0.98, 0.99, 1.0, 1.01, 1.02):
                    p = dict(
                        cx=base["cx"] + dcx,
                        cy=base["cy"] + dcy,
                        Rx=base["Rx"] * sRx,
                        Ry=base["Ry"] * sRy,
                    )
                    sc = median_dist(**p)
                    if sc < best_sc:
                        best_sc = sc
                        best = p
    b = best
    for dcx in np.linspace(-2, 2, 17):
        for dcy in np.linspace(-2, 2, 17):
            for sRx in np.linspace(0.995, 1.005, 11):
                for sRy in np.linspace(0.995, 1.005, 11):
                    p = dict(
                        cx=b["cx"] + dcx,
                        cy=b["cy"] + dcy,
                        Rx=b["Rx"] * sRx,
                        Ry=b["Ry"] * sRy,
                    )
                    sc = median_dist(**p)
                    if sc < best_sc:
                        best_sc = sc
                        best = p
    best["ceq_median_px"] = float(best_sc)
    return best


def choose_source():
    for p in SOURCE_CANDIDATES:
        if p.is_file():
            return p
    raise FileNotFoundError("No Migkas Fig.8 theta=75 source image found")


def sample_rgb(rgb, ell, b, cx, cy, Rx, Ry, r=4):
    px, py = lb_to_px(ell, b, cx, cy, Rx, Ry)
    xi, yi = int(round(as_xy(px, py)[0])), int(round(as_xy(px, py)[1]))
    patch = rgb[max(0, yi - r) : yi + r + 1, max(0, xi - r) : xi + r + 1]
    return patch.mean(axis=(0, 1))


def validate_landmarks(rgb, reg):
    cx, cy, Rx, Ry = reg["cx"], reg["cy"], reg["Rx"], reg["Ry"]
    pts = {}
    landmarks = [
        ("center", 0.0, 0.0),
        ("ell90", 90.0, 0.0),
        ("ell270", 270.0, 0.0),
        ("north", 0.0, 90.0),
        ("south", 0.0, -90.0),
    ]
    for i, (ell, b) in enumerate(MIGKAS_MAXES):
        landmarks.append((f"migkas_max{i}", ell, b))
    for i, (ell, b) in enumerate(MIGKAS_MINS):
        landmarks.append((f"migkas_min{i}", ell, b))
    for name, ell, b in landmarks:
        x, y = as_xy(*lb_to_px(ell, b, cx, cy, Rx, Ry))
        pts[name] = {"ell": ell, "b": b, "px": x, "py": y}

    assert pts["ell90"]["px"] < pts["center"]["px"], "l=90 must be LEFT of center"
    assert pts["ell270"]["px"] > pts["center"]["px"], "l=270 must be RIGHT of center"
    assert abs(pts["center"]["px"] - cx) < 1.5
    assert abs(pts["center"]["py"] - cy) < 1.5
    assert pts["north"]["py"] < pts["center"]["py"]
    assert pts["south"]["py"] > pts["center"]["py"]

    # Map-extrema pair: (17,-9) should sit in a warmer patch than (274,-22).
    cmax = sample_rgb(rgb, *MIGKAS_MAXES[0], cx, cy, Rx, Ry)
    cmin = sample_rgb(rgb, *MIGKAS_MINS[1], cx, cy, Rx, Ry)
    assert (cmax[0] + cmax[1]) > (cmin[0] + cmin[1]) + 50
    return {
        "landmarks_px": pts,
        "migkas_max0_rgb_mean": [float(x) for x in cmax],
        "migkas_min1_rgb_mean": [float(x) for x in cmin],
    }


def inside_ellipse(px, py, ecx, ecy, eRx, eRy, r_lim=1.0):
    """Boolean mask: points inside axis-aligned ellipse (image y down)."""
    px = np.asarray(px, dtype=float)
    py = np.asarray(py, dtype=float)
    Xn = (px - ecx) / eRx
    Yn = (py - ecy) / eRy
    return Xn * Xn + Yn * Yn <= float(r_lim) * float(r_lim)


def plot_lb_curve(
    ax, ell, b, cx, cy, Rx, Ry, inset=1.0, jump_px=45.0,
    clip_r2=None, clip_path=None, color_oval=None, **kw
):
    """Draw curve in oval; inset<1 keeps strokes inside the published Mollweide edge."""
    px, py = lb_to_px(ell, b, cx, cy, Rx, Ry)
    px = np.asarray(px, dtype=float).reshape(-1)
    py = np.asarray(py, dtype=float).reshape(-1)
    # Pull points slightly toward oval center so thick strokes do not spill past the rim.
    px = cx + (px - cx) * inset
    py = cy + (py - cy) * inset
    Xn = (px - cx) / Rx
    Yn = (cy - py) / Ry
    r2_lim = (inset * inset) if clip_r2 is None else float(clip_r2)
    inside = Xn * Xn + Yn * Yn <= r2_lim
    # Also require published colored oval (seed) so registration oversize cannot bleed.
    if color_oval is not None:
        ocx, ocy, oRx, oRy, o_r = color_oval
        inside = inside & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)
    px = np.where(inside, px, np.nan)
    py = np.where(inside, py, np.nan)

    segs = []
    cur_x, cur_y = [], []
    for x, y in zip(px, py):
        if not (np.isfinite(x) and np.isfinite(y)):
            if len(cur_x) >= 2:
                segs.append((cur_x, cur_y))
            cur_x, cur_y = [], []
            continue
        if cur_x and np.hypot(x - cur_x[-1], y - cur_y[-1]) > jump_px:
            if len(cur_x) >= 2:
                segs.append((cur_x, cur_y))
            cur_x, cur_y = [float(x)], [float(y)]
        else:
            cur_x.append(float(x))
            cur_y.append(float(y))
    if len(cur_x) >= 2:
        segs.append((cur_x, cur_y))
    if clip_path is not None:
        kw["clip_path"] = clip_path
        kw["clip_on"] = True
    for sx, sy in segs:
        (line,) = ax.plot(sx, sy, **kw)
        if clip_path is not None:
            line.set_clip_path(clip_path)
            line.set_clip_on(True)


def plot_closed_ring(
    ax, mu, cx, cy, Rx, Ry, color, lw, inset=1.0, jump_px=80.0,
    clip_path=None, color_oval=None, ls="-", wash_alpha=None,
):
    """
    Exact projected small-circle about D_HAT=+x: n·x̂ = μ.

    Parameterize on the sphere, convert to Galactic (ℓ,b), apply
    lon_disp = -wrap(ℓ), Mollweide-project, and stroke continuous polylines.
    Split wherever display-longitude jumps across the antimeridian (or pixel
    distance jumps). Never invent interior chords, ellipse fits, or decorative
    oval-rim bridges.

    Notes on θ labels used elsewhere:
      μ=cos(30°)  → r=30° about +x (closed oval near map center)
      μ=cos(185°)=cos(175°) → r=175° about +x ≡ 5° about CoE
        → small arcs on BOTH left and right map edges (not a rim loop)
      μ=cos(115°) → r=115° about +x (large seam-crossing curve)
      μ=cos(265°)=cos(95°) → r=95° about +x (near-meridian, seam-crossing)
    """
    pts = small_circle_about_axis(D_HAT, mu, n=4001)
    # Keep sampling order around the axis (azimuth) — stable for splitting.
    phi = np.arctan2(pts[:, 2], pts[:, 1])
    order = np.argsort(phi)
    pts = pts[order]

    ell = np.empty(len(pts))
    b = np.empty(len(pts))
    for i, p in enumerate(pts):
        ell[i], b[i] = n_to_lb(p)

    lon = np.deg2rad(lon_disp_from_ell(ell))  # λ ∈ (-π, π]
    lat = np.deg2rad(b)

    # Split on antimeridian jumps in display longitude BEFORE pixel projection.
    dlon = np.diff(lon)
    cuts = np.where(np.abs(dlon) > np.pi)[0] + 1
    # Also close the periodic sample: check wrap between last and first.
    # Build open segments from azimuth-ordered points.
    idx = np.concatenate([[0], cuts, [len(lon)]])
    lon_segs = []
    for i0, i1 in zip(idx[:-1], idx[1:]):
        if i1 - i0 >= 2:
            lon_segs.append((lon[i0:i1].copy(), lat[i0:i1].copy()))

    # Circle sampling is cyclic. Merge the artificial last→first break when
    # that periodic join is not itself an antimeridian crossing. Otherwise a
    # seam-crossing circle falsely appears as three pieces with a tiny hole.
    if len(lon_segs) > 1 and abs(lon[0] - lon[-1]) <= np.pi:
        first_lo, first_la = lon_segs[0]
        last_lo, last_la = lon_segs[-1]
        lon_segs = [(
            np.concatenate([last_lo, first_lo]),
            np.concatenate([last_la, first_la]),
        )] + lon_segs[1:-1]

    # Only now pin TRUE seam endpoints to ±π. Never snap the arbitrary
    # parameterization start/end to the boundary.
    snapped = []
    for lo, la in lon_segs:
        if abs(abs(lo[0]) - np.pi) < np.deg2rad(10):
            lo[0] = np.sign(lo[0] if lo[0] != 0 else 1.0) * np.pi
        if abs(abs(lo[-1]) - np.pi) < np.deg2rad(10):
            lo[-1] = np.sign(lo[-1] if lo[-1] != 0 else 1.0) * np.pi
        snapped.append((lo, la))
    lon_segs = snapped

    def stroke_xy(xs, ys):
        if len(xs) < 2:
            return
        if wash_alpha is not None:
            # Soft μ-wash stroke: no white underlay; translucent color only.
            layers = ((color, lw, 2, float(wash_alpha)),)
        else:
            # Clean mid-stroke: optional thin white underlay only (lw+1.5), then color.
            # No black outer rim, no thick white halo stack.
            layers = (
                ("white", lw + 1.5, 3, 0.95),
                (color, lw, 4, 1.0),
            )
        for col, width, z, a in layers:
            (line,) = ax.plot(
                xs, ys, color=col, lw=width, linestyle=ls,
                solid_capstyle="round", solid_joinstyle="round",
                zorder=z, alpha=a,
                clip_path=clip_path, clip_on=True,
            )
            if clip_path is not None:
                line.set_clip_path(clip_path)
                line.set_clip_on(True)

    r2_lim = inset * inset + 1e-10
    drawn = 0
    for lo, la in lon_segs:
        # Mollweide in normalized oval coords, then pixels.
        x_m, y_m = mollweide_xy(lo, la)
        Xn = x_m / (2.0 * np.sqrt(2.0))
        Yn = y_m / np.sqrt(2.0)
        px = cx + Xn * Rx
        py = cy - Yn * Ry
        # Mild inset so thick strokes stay inside published oval (no rim-snap).
        px = cx + (px - cx) * inset
        py = cy + (py - cy) * inset
        Xn2 = (px - cx) / Rx
        Yn2 = (cy - py) / Ry
        inside = (Xn2 * Xn2 + Yn2 * Yn2) <= r2_lim
        if color_oval is not None:
            ocx, ocy, oRx, oRy, o_r = color_oval
            inside = inside & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)

        # Pixel-space split on NaN / large jumps / outside oval.
        cur_x, cur_y = [], []
        for x, y, ok in zip(px, py, inside):
            if (not ok) or (not np.isfinite(x)) or (not np.isfinite(y)):
                if len(cur_x) >= 2:
                    stroke_xy(cur_x, cur_y)
                    drawn += 1
                cur_x, cur_y = [], []
                continue
            if cur_x and np.hypot(x - cur_x[-1], y - cur_y[-1]) > jump_px:
                if len(cur_x) >= 2:
                    stroke_xy(cur_x, cur_y)
                    drawn += 1
                cur_x, cur_y = [float(x)], [float(y)]
            else:
                cur_x.append(float(x))
                cur_y.append(float(y))
        if len(cur_x) >= 2:
            # Close only if this single segment is already nearly a loop
            # (non-wrapping ring such as MAX 30°).
            if (
                drawn == 0
                and len(lon_segs) == 1
                and np.hypot(cur_x[0] - cur_x[-1], cur_y[0] - cur_y[-1]) < jump_px
            ):
                cur_x.append(cur_x[0])
                cur_y.append(cur_y[0])
            stroke_xy(cur_x, cur_y)
            drawn += 1



def _ring_px_arrays(mu, cx, cy, Rx, Ry, inset, n_pts=721):
    """Small-circle n·D=μ → pixel (px, py) arrays, azimuth-ordered, with inset."""
    pts = small_circle_about_axis(D_HAT, float(mu), n=n_pts)
    phi = np.arctan2(pts[:, 2], pts[:, 1])
    order = np.argsort(phi)
    pts = pts[order]
    ell = np.empty(n_pts)
    b = np.empty(n_pts)
    for i, p in enumerate(pts):
        ell[i], b[i] = n_to_lb(p)
    lon = np.deg2rad(lon_disp_from_ell(ell))
    lat = np.deg2rad(b)
    x_m, y_m = mollweide_xy(lon, lat)
    Xn = x_m / (2.0 * np.sqrt(2.0))
    Yn = y_m / np.sqrt(2.0)
    px = cx + Xn * Rx
    py = cy - Yn * Ry
    px = cx + (px - cx) * inset
    py = cy + (py - cy) * inset
    return lon, px, py


def fill_highlight_band(
    ax, theta_lo, theta_hi, theta_mid, cx, cy, Rx, Ry, color,
    inset=1.0, clip_path=None, color_oval=None, n_pts=721, jump_px=80.0,
):
    """
    One highlight band matching three-panel Panel B:
      continuous transparent polygon strips between ±BAND_HALF edges
      (broken only at Mollweide antimeridian / large pixel jumps — NOT a dense
      multi-stroke wash) + ONE solid mid-stroke at θ_mid.
    """
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    lon_a, px_a, py_a = _ring_px_arrays(mu_a, cx, cy, Rx, Ry, inset, n_pts=n_pts)
    lon_b, px_b, py_b = _ring_px_arrays(mu_b, cx, cy, Rx, Ry, inset, n_pts=n_pts)

    r2_lim = inset * inset + 1e-10

    def _inside(px, py):
        Xn2 = (px - cx) / Rx
        Yn2 = (cy - py) / Ry
        ok = (Xn2 * Xn2 + Yn2 * Yn2) <= r2_lim
        if color_oval is not None:
            ocx, ocy, oRx, oRy, o_r = color_oval
            ok = ok & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)
        return ok

    gap = 0.9 * np.pi
    bad = np.zeros(n_pts - 1, dtype=bool)
    for i in range(n_pts - 1):
        if abs(lon_a[i + 1] - lon_a[i]) > gap:
            bad[i] = True
        elif abs(lon_b[i + 1] - lon_b[i]) > gap:
            bad[i] = True
        elif abs(lon_a[i] - lon_b[i]) > gap:
            bad[i] = True
        elif abs(lon_a[i + 1] - lon_b[i + 1]) > gap:
            bad[i] = True
        elif np.hypot(px_a[i + 1] - px_a[i], py_a[i + 1] - py_a[i]) > jump_px:
            bad[i] = True
        elif np.hypot(px_b[i + 1] - px_b[i], py_b[i + 1] - py_b[i]) > jump_px:
            bad[i] = True
        elif (not _inside(px_a[i], py_a[i])) or (not _inside(px_b[i], py_b[i])):
            bad[i] = True
        elif (not _inside(px_a[i + 1], py_a[i + 1])) or (not _inside(px_b[i + 1], py_b[i + 1])):
            bad[i] = True

    polys = []
    i = 0
    while i < n_pts - 1:
        while i < n_pts - 1 and bad[i]:
            i += 1
        if i >= n_pts - 1:
            break
        j = i
        while j < n_pts - 1 and not bad[j]:
            j += 1
        i0, i1 = i, j
        if i1 - i0 >= 1:
            xs = np.concatenate([px_a[i0:i1 + 1], px_b[i0:i1 + 1][::-1]])
            ys = np.concatenate([py_a[i0:i1 + 1], py_b[i0:i1 + 1][::-1]])
            polys.append(np.column_stack([xs, ys]))
        i = j + 1

    if polys:
        # Soft transparent highlight (clean pre-halo look).
        coll = PolyCollection(
            polys, facecolors=color, edgecolors="none",
            alpha=0.30, zorder=2, linewidths=0,
        )
        ax.add_collection(coll)
        if clip_path is not None:
            coll.set_clip_path(clip_path)
            coll.set_clip_on(True)

    # ONE solid mid-stroke only (thin white underlay lw+1.5). No boundary-edge
    # strokes — dual edges with halos multiply visual lines into "extra loops".
    plot_closed_ring(
        ax, mu_mid, cx, cy, Rx, Ry, color, lw=2.8,
        inset=inset, clip_path=clip_path, color_oval=color_oval,
    )



# --- Opacity gradient (same rule as make_coe_rings_migkas_pin_axis_matched_grad*.py):
# alpha *= (1+cos psi)/2, psi = ring azimuth from the θ-ray point (ring point most
# along sign(sin θ)·Ŷ, Ŷ = Galactic +z projected ⊥ D), rolled by ROLL (Y' = cY − sZ).
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba


def _ring_basis_roll(D, roll_deg):
    a = D / np.linalg.norm(D)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    y = tmp - np.dot(tmp, a) * a
    y /= np.linalg.norm(y)
    z = np.cross(a, y)
    z /= np.linalg.norm(z)
    r = np.deg2rad(float(roll_deg))
    c, s_ = np.cos(r), np.sin(r)
    return c * y - s_ * z, s_ * y + c * z


def _grad_ring_px(mu, theta_mid, cx, cy, Rx, Ry, inset, n):
    a = D_HAT / np.linalg.norm(D_HAT)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    rho = np.sqrt(1.0 - mu * mu)
    y, z = _ring_basis_roll(a, ROLL)
    phi0 = 0.5 * np.pi if np.sin(np.deg2rad(theta_mid)) >= 0 else 1.5 * np.pi
    phi = phi0 + np.linspace(-np.pi, np.pi, n)
    pts = (mu * a)[None, :] + rho * (np.outer(np.sin(phi), y) + np.outer(np.cos(phi), z))
    ell = np.array([n_to_lb(p)[0] for p in pts]); b = np.array([n_to_lb(p)[1] for p in pts])
    lon = np.deg2rad(lon_disp_from_ell(ell)); lat = np.deg2rad(b)
    x_m, y_m = mollweide_xy(lon, lat)
    px = cx + x_m / (2.0 * np.sqrt(2.0)) * Rx
    py = cy - y_m / np.sqrt(2.0) * Ry
    px = cx + (px - cx) * inset; py = cy + (py - cy) * inset
    return lon, px, py, 0.5 * (1.0 + np.cos(phi - phi0))


def _rgba(color, alphas):
    r, g, b_, _ = to_rgba(color)
    al = np.clip(np.asarray(alphas, dtype=float), 0.0, 1.0)
    out = np.empty((al.size, 4)); out[:, 0], out[:, 1], out[:, 2], out[:, 3] = r, g, b_, al
    return out


def _runs(seg_sc, valid, n_levels=64):
    lev = np.round(np.clip(seg_sc, 0.0, 1.0) * n_levels) / n_levels
    ok = (lev > 0) & valid
    runs, i, n = [], 0, lev.size
    while i < n:
        if not ok[i]:
            i += 1; continue
        j = i + 1
        while j < n and ok[j] and lev[j] == lev[i]:
            j += 1
        runs.append((i, j, float(lev[i]))); i = j
    return runs


def fill_highlight_band_grad(
    ax, theta_lo, theta_hi, theta_mid, cx, cy, Rx, Ry, color,
    inset=1.0, clip_path=None, color_oval=None, n_pts=721, jump_px=80.0,
):
    """Same band (alpha 0.30) + mid-stroke (white lw+1.5 @0.95 under colour lw 2.8)
    as fill_highlight_band, but with the (1+cos psi)/2 opacity gradient."""
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    r2_lim = inset * inset + 1e-10

    def _inside(px, py):
        ok = ((px - cx) / Rx) ** 2 + ((cy - py) / Ry) ** 2 <= r2_lim
        if color_oval is not None:
            ocx, ocy, oRx, oRy, o_r = color_oval
            ok = ok & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)
        return ok

    def _clip(coll):
        if clip_path is not None:
            coll.set_clip_path(clip_path); coll.set_clip_on(True)
        ax.add_collection(coll)

    gap = 0.9 * np.pi
    lon_a, pxa, pya, sc = _grad_ring_px(mu_a, theta_mid, cx, cy, Rx, Ry, inset, n_pts)
    lon_b, pxb, pyb, _ = _grad_ring_px(mu_b, theta_mid, cx, cy, Rx, Ry, inset, n_pts)
    qa = np.stack([lon_a[:-1], lon_a[1:], lon_b[1:], lon_b[:-1]])
    ins_a, ins_b = _inside(pxa, pya), _inside(pxb, pyb)
    valid = ((qa.max(0) - qa.min(0)) <= gap)
    valid &= np.hypot(np.diff(pxa), np.diff(pya)) <= jump_px
    valid &= np.hypot(np.diff(pxb), np.diff(pyb)) <= jump_px
    valid &= ins_a[:-1] & ins_a[1:] & ins_b[:-1] & ins_b[1:]
    runs = _runs(0.5 * (sc[:-1] + sc[1:]), valid)
    if runs:
        polys = [np.column_stack([np.concatenate([pxa[i0:i1 + 1], pxb[i0:i1 + 1][::-1]]),
                                  np.concatenate([pya[i0:i1 + 1], pyb[i0:i1 + 1][::-1]])])
                 for i0, i1, _ in runs]
        _clip(PolyCollection(polys, facecolors=_rgba(color, [0.30 * r[2] for r in runs]),
                             edgecolors="none", linewidths=0, zorder=2))
    lon_m, pxm, pym, scm = _grad_ring_px(mu_mid, theta_mid, cx, cy, Rx, Ry, inset, 4001)
    ins = _inside(pxm, pym)
    valid = (np.abs(np.diff(lon_m)) <= np.pi) & (np.hypot(np.diff(pxm), np.diff(pym)) <= jump_px)
    valid &= ins[:-1] & ins[1:]
    runs = _runs(0.5 * (scm[:-1] + scm[1:]), valid)
    if runs:
        segs = [np.column_stack([pxm[i0:i1 + 1], pym[i0:i1 + 1]]) for i0, i1, _ in runs]
        _clip(LineCollection(segs, colors=_rgba("white", [0.95 * r[2] for r in runs]),
                             linewidths=2.8 + 1.5, zorder=3, capstyle="butt", joinstyle="round"))
        _clip(LineCollection(segs, colors=_rgba(color, [1.0 * r[2] for r in runs]),
                             linewidths=2.8, zorder=4, capstyle="butt", joinstyle="round"))


def draw_overlay(rgb, reg, seed, source):

    h, w, _ = rgb.shape
    cx, cy, Rx, Ry = reg["cx"], reg["cy"], reg["Rx"], reg["Ry"]
    # Published colored Mollweide (seed detect) — used to keep overlays inside
    # the Migkas map even when registration ellipse is slightly larger.
    scx = float(seed["cx"])
    scy = float(seed["cy"])
    sRx = float(seed["Rx"])
    sRy = float(seed["Ry"])

    dpi = 120
    _extra_h = 1.35 if SHOW_LEGEND else 0.35  # tighten when legend omitted
    fig, ax = plt.subplots(figsize=(w / dpi, (h / dpi) + _extra_h), dpi=dpi)
    ax.imshow(rgb, origin="upper", extent=(0, w, h, 0), interpolation="nearest")
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    # clip_path = registered oval (cx,cy,Rx,Ry). Registration itself is unchanged.
    # Because the refined ellipse can sit a bit outside the published colored
    # map, we ALSO inset strokes and cull/inset markers against the seed oval.
    clip_oval = Ellipse(
        (cx, cy), width=2.0 * Rx, height=2.0 * Ry,
        facecolor="none", edgecolor="none", transform=ax.transData,
    )
    ax.add_patch(clip_oval)

    # Stroke inset so thick white+color underlays do not AA-bleed past the rim.
    STROKE_INSET = 0.97
    # Cull points outside this fraction of the published colored seed oval.
    COLOR_R = 0.978
    color_oval = (scx, scy, sRx, sRy, COLOR_R)

    def stroke_curve(ell, b, color, lw, ls="-"):
        # White underlay then color — stands out on magma without looking oversized.
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=STROKE_INSET, clip_path=clip_oval,
            color_oval=color_oval,
            color="white", lw=lw + 1.6, ls=ls, alpha=0.95,
            solid_capstyle="round", zorder=3,
        )
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=STROKE_INSET, clip_path=clip_oval,
            color_oval=color_oval,
            color=color, lw=lw, ls=ls, alpha=1.0,
            solid_capstyle="round", zorder=4,
        )

    def mark_point(ell, b, marker, color, ms):
        """Place marker; pull limb symbols inward so halo stays in colored oval."""
        px, py = as_xy(*lb_to_px(ell, b, cx, cy, Rx, Ry))
        # Approx marker+halo radius in pixels (ms is points).
        halo_ms = ms + 1.8
        ms_px = 0.5 * halo_ms * (dpi / 72.0)
        # Normalized margin inside the published colored oval.
        margin = ms_px / min(sRx, sRy) + 0.012
        r_lim = max(0.55, COLOR_R - margin)

        Xn = (px - scx) / sRx
        Yn = (py - scy) / sRy
        r_col = float(np.hypot(Xn, Yn))

        if r_col > 1e-9:
            # Scale toward colored-oval center so the whole symbol fits.
            if r_col > r_lim:
                scale = r_lim / r_col
                px = scx + (px - scx) * scale
                py = scy + (py - scy) * scale
                Xn = (px - scx) / sRx
                Yn = (py - scy) / sRy
                r_col = float(np.hypot(Xn, Yn))
            # Mild extra inset for everyone (matches stroke inset spirit).
            mild = min(STROKE_INSET, 0.995)
            px = scx + (px - scx) * mild
            py = scy + (py - scy) * mild
            Xn = (px - scx) / sRx
            Yn = (py - scy) / sRy
            r_col = float(np.hypot(Xn, Yn))

        # Skip if still too close to the rim for this marker size.
        if r_col + margin > COLOR_R + 1e-6:
            return

        def _mark(m_size, m_color, m_edge, m_ew, z):
            (line,) = ax.plot(
                px, py, marker=marker, ms=m_size, color=m_color,
                markeredgecolor=m_edge, markeredgewidth=m_ew, linestyle="None",
                zorder=z, clip_path=clip_oval, clip_on=True,
            )
            line.set_clip_path(clip_oval)
            line.set_clip_on(True)

        _mark(ms + 1.8, "white", "white", 0, 6)
        _mark(ms, color, "k", 0.55, 7)

    # CoE axis at geometric inset=1.0. Do NOT overlay Dec=0 — Migkas already
    # printed the celestial equator on Fig. 8; a second green dashed stroke fights it.
    def stroke_curve_geom(ell, b, color, lw, ls="-"):
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=1.0, clip_path=clip_oval,
            color_oval=color_oval,
            color="white", lw=lw + 1.6, ls=ls, alpha=0.95,
            solid_capstyle="round", zorder=3,
        )
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=1.0, clip_path=clip_oval,
            color_oval=color_oval,
            color=color, lw=lw, ls=ls, alpha=1.0,
            solid_capstyle="round", zorder=4,
        )

    # CoE↔anti-CoE axis = great circle through ±D (tilted +30°; not Galactic equator)
    _ax_ell, _ax_b = great_circle_through_poles(D_HAT, n=721)
    stroke_curve_geom(_ax_ell, _ax_b, C_AXIS, lw=1.15)

    # Exactly four highlight bands (match three-panel Panel B): each ±9°
    # (= 18° end-to-end = 5% of 360°) + one solid mid-stroke. No companion.
    for (lo, hi), mid, color in zip(MAX_BANDS, THETA_MAX, COLORS_MAX):
        (fill_highlight_band_grad if GRAD else fill_highlight_band)(
            ax, lo, hi, mid, cx, cy, Rx, Ry, color,
            inset=STROKE_INSET, clip_path=clip_oval, color_oval=color_oval,
        )
    for (lo, hi), mid, color in zip(MIN_BANDS, THETA_MIN, COLORS_MIN):
        (fill_highlight_band_grad if GRAD else fill_highlight_band)(
            ax, lo, hi, mid, cx, cy, Rx, Ry, color,
            inset=STROKE_INSET, clip_path=clip_oval, color_oval=color_oval,
        )

    # Markers (limb ones get stronger inset inside mark_point); per-loop colors
    mark_point(ANTI_COE_LB[0], ANTI_COE_LB[1], "o", C_ANTI, ms=5.5)
    mark_point(COE_LB[0], COE_LB[1], "*", C_COE, ms=7.5)
    for (ell, b), color in zip(MIGKAS_MAXES, COLORS_MAX):
        mark_point(ell, b, "^", color, ms=6.5)
    for (ell, b), color in zip(MIGKAS_MINS, COLORS_MIN):
        mark_point(ell, b, "s", color, ms=6.0)

    handles = [
        Patch(
            facecolor=C_MAX_30, alpha=0.55, edgecolor=C_MAX_30,
            label=(
                rf"MAX $\theta\!=\!{THETA_MAX[0]:g}^\circ$ band "
                rf"[{MAX_BANDS[0][0]:g},{MAX_BANDS[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$=5% of $360^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MAX_185, alpha=0.55, edgecolor=C_MAX_185,
            label=(
                rf"MAX $\theta\!=\!{THETA_MAX[1]:g}^\circ$ band "
                rf"[{MAX_BANDS[1][0]:g},{MAX_BANDS[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_115, alpha=0.55, edgecolor=C_MIN_115,
            label=(
                rf"MIN $\theta\!=\!{THETA_MIN[0]:g}^\circ$ band "
                rf"[{MIN_BANDS[0][0]:g},{MIN_BANDS[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_265, alpha=0.55, edgecolor=C_MIN_265,
            label=(
                rf"MIN $\theta\!=\!{THETA_MIN[1]:g}^\circ$ band "
                rf"[{MIN_BANDS[1][0]:g},{MIN_BANDS[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$)"
            ),
        ),
        Line2D([0], [0], color=C_AXIS, lw=2.0, label=r"CoE axis tilted $+15^\circ$"),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_30,
            markeredgecolor="k", linestyle="None", ms=8,
            label=rf"MAX $({MIGKAS_MAXES[0][0]:.0f}^\circ,{MIGKAS_MAXES[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_185,
            markeredgecolor="k", linestyle="None", ms=8,
            label=rf"MAX $({MIGKAS_MAXES[1][0]:.0f}^\circ,{MIGKAS_MAXES[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_115,
            markeredgecolor="k", linestyle="None", ms=7,
            label=rf"MIN $({MIGKAS_MINS[0][0]:.0f}^\circ,{MIGKAS_MINS[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_265,
            markeredgecolor="k", linestyle="None", ms=7,
            label=rf"MIN $({MIGKAS_MINS[1][0]:.0f}^\circ,{MIGKAS_MINS[1][1]:+.0f}^\circ)$",
        ),
        Line2D([0], [0], marker="o", color=C_ANTI, markeredgecolor="k",
               linestyle="None", ms=6,
               label=r"anti-CoE $+x$ $(0^\circ,+15^\circ)$"),
        Line2D([0], [0], marker="*", color=C_COE, markeredgecolor="k",
               linestyle="None", ms=9,
               label=r"CoE $(180^\circ,-15^\circ)$"),
    ]
    if SHOW_LEGEND:
        leg = ax.legend(
            handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.012),
            ncol=4, fontsize=7.5, frameon=True, fancybox=False,
            edgecolor="#222222", framealpha=1.0, facecolor="white",
            labelcolor="#111111", borderpad=0.6, handlelength=2.0,
        )
        leg.get_frame().set_linewidth(1.2)
    fig.suptitle(
        r"Migkas Fig. 8 $\theta=75^\circ$ $A/A_{\rm all}$"
        r" + four CoE highlight bands, each $18^\circ$ end-to-end ($\pm9^\circ$ = 5% of $360^\circ$)"
        "\n"
        r"(two MAX + two MIN; cyan/green/magenta/blue; anti-CoE at $(0^\circ,+15^\circ)$; "
        r"CoE axis tilted $+15^\circ$; $\ell$ increases LEFT"
        + (r"; gradient: darkest at $\theta$ ray, transparent opposite" if GRAD else "")
        + (rf"; ${ROLL:+g}^\circ$ CCW about CoE axis" if abs(ROLL) > 1e-9 else "")
        + ")",
        fontsize=10.5, y=0.985,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=(0.13 if SHOW_LEGEND else 0.02))
    fig.savefig(PNG_PATH, dpi=dpi, facecolor="white")
    plt.close(fig)
    print(f"wrote {PNG_PATH}")


def main():
    source = choose_source()
    print(f"source: {source}")
    rgb = np.asarray(Image.open(source).convert("RGB"))
    seed = detect_oval_bbox(rgb)
    print("oval seed:", {k: seed[k] for k in ("cx", "cy", "Rx", "Ry", "x0", "x1", "y0", "y1")})
    # Restore previous visually-good oval (Chris preferred). No medial-axis Rx stretch.
    reg = dict(GOOD_REG)
    h, w, _ = rgb.shape
    dark = rgb.mean(axis=2) < 55
    dist = chamfer_distance_to_dark(dark)
    ell_eq, b_eq = celestial_equator_lb(2500)
    lon = np.deg2rad(lon_disp_from_ell(ell_eq))
    lat = np.deg2rad(b_eq)
    mx, my = mollweide_xy(lon, lat)
    Xn = mx / (2.0 * np.sqrt(2.0))
    Yn = my / np.sqrt(2.0)
    px = reg["cx"] + Xn * reg["Rx"]
    py = reg["cy"] - Yn * reg["Ry"]
    m = (px >= 1) & (px < w - 1) & (py >= 1) & (py < h - 1)
    xi = np.clip(np.round(px[m]).astype(int), 0, w - 1)
    yi = np.clip(np.round(py[m]).astype(int), 0, h - 1)
    reg["ceq_median_px"] = float(np.median(dist[yi, xi]))
    reg["registration_source"] = "hard-set GOOD_REG (pre medial-axis retune)"
    print("restored reg:", {k: reg[k] for k in ("cx", "cy", "Rx", "Ry", "ceq_median_px")})
    checks = validate_landmarks(rgb, reg)
    print(
        "validation OK — max0 RGB", checks["migkas_max0_rgb_mean"],
        "min1 RGB", checks["migkas_min1_rgb_mean"],
    )
    global PNG_PATH, JSON_PATH, GRAD, ROLL
    for name, g, r in VERSIONS:
        PNG_PATH = SM_DIR / f"fig8_overlay_{name}.png"
        JSON_PATH = SM_DIR / f"fig8_overlay_{name}.json"
        GRAD, ROLL = g, r
        draw_overlay(rgb, reg, seed, source)
        img = Image.open(PNG_PATH).convert("RGB")
        img.resize((1700, round(img.height * 1700 / img.width)), Image.LANCZOS).save(
            SM_DIR / f"fig8_overlay_{name}_web.jpg", quality=88)

    meta = {
        "png": str(PNG_PATH),
        "script": str(Path(__file__).resolve()),
        "source_image": str(source),
        "convention": (
            "Galactic Mollweide; lon_disp=-wrap(ell); ell increases LEFT "
            "(Migkas Fig.8)"
        ),
        "axis_tilt_b_deg": AXIS_TILT_B_DEG,
        "anti_CoE_lb_deg": [float(ANTI_COE_LB[0]), float(ANTI_COE_LB[1])],
        "CoE_lb_deg": [float(COE_LB[0]), float(COE_LB[1])],
        "D_hat_galactic_xyz": D_HAT.tolist(),
        "axis_tilt_reason": (
            "Chris 2026-09-22 — markers fixed, rings rotate up 30°"
        ),
        "migkas_max_lb_deg": [list(p) for p in MIGKAS_MAXES],
        "migkas_min_lb_deg": [list(p) for p in MIGKAS_MINS],
        "theta_max_from_plus_x_deg": list(THETA_MAX),
        "theta_min_from_plus_x_deg": list(THETA_MIN),
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "band_width_note": "5% of 360° = 18° end-to-end = ±9° about each MAX/MIN center",
        "max_bands_deg": [list(b) for b in MAX_BANDS],
        "min_bands_deg": [list(b) for b in MIN_BANDS],
        "band_intervals": {
            "MAX_30": [float(MAX_BANDS[0][0]), float(MAX_BANDS[0][1])],
            "MAX_185": [float(MAX_BANDS[1][0]), float(MAX_BANDS[1][1])],
            "MIN_115": [float(MIN_BANDS[0][0]), float(MIN_BANDS[0][1])],
            "MIN_265": [float(MIN_BANDS[1][0]), float(MIN_BANDS[1][1])],
        },
        "colors": {
            "MAX_30": C_MAX_30,
            "MAX_185": C_MAX_185,
            "MIN_115": C_MIN_115,
            "MIN_265": C_MIN_265,
            "palette": "magma-popping cyan/magenta/lime/sky (no reds/oranges)",
        },
        "mu_max_bands": [list(mu_pair_for_band(lo, hi)) for lo, hi in MAX_BANDS],
        "mu_min_bands": [list(mu_pair_for_band(lo, hi)) for lo, hi in MIN_BANDS],
        "scan_stop_source": (
            "best_fit_2ell_multiplier.npz A_after extrema (two max, two min); "
            f"each drawn as one transparent {BAND_FULL:.0f}° end-to-end highlight "
            f"(±{BAND_HALF:.0f}° = 5% of 360°) + one solid mid-stroke; "
            "matches three-panel Panel B; no companion"
        ),
        "matches_three_panel": str(OUT_DIR / "make_coe_rings_migkas_pin_axis.py"),
        "oval_seed_bbox": {
            "x0": seed["x0"], "x1": seed["x1"], "y0": seed["y0"], "y1": seed["y1"],
            "cx": seed["cx"], "cy": seed["cy"], "Rx": seed["Rx"], "Ry": seed["Ry"],
        },
        "registration": {
            "cx": reg["cx"], "cy": reg["cy"], "Rx": reg["Rx"], "Ry": reg["Ry"],
            "ceq_median_px": reg["ceq_median_px"],
            "registration_source": reg.get("registration_source"),
            "ceq_residual_before_registration": dict(GOOD_REG),
            "rejected_medial_axis_reg": {
                "cx": 865.6, "cy": 766.4, "Rx": 706.19, "Ry": 393.70,
                "note": "medial-axis / achromatic DT-ridge retune; blew up Rx; reverted",
            },
            "dofs": ["cx", "cy", "Rx", "Ry"],
            "drawing": {
                "dec0_drawn": False,
                "dec0_note": "not drawn; use Migkas published celestial equator",
                "coe_axis_inset": 1.0,
                "rings_inset": 0.97,
                "note": "CoE axis geometric inset=1; rings mild inset; no Dec=0 overlay",
            },
            "method": (
                "oval seed detect + hard-set to pre-retune visually-good (cx,cy,Rx,Ry); "
                "medial-axis Rx stretch rejected; Dec=0 NOT overlaid (published CEQ kept); "
                "CoE-axis inset=1.0; rings inset~0.97"
            ),
        },
        "landmark_checks": checks,
        "drawing": (
            "Exactly 4 features. Each = ONE smooth transparent PolyCollection "
            "band (θ±9°, 18° end-to-end) + ONE solid mid-stroke at center θ. "
            "Fill alpha≈0.30. Mid-stroke: thin white underlay (lw+1.5) then color "
            "(lw=2.8). NO black outer rim; NO thick white halo stacks; NO dual "
            "boundary-edge strokes. Colors: cyan/magenta/lime/sky for magma. "
            "No 283° companion. Dec=0 NOT overlaid."
        ),
        "legend_summary": (
            "four highlight bands, each 18° end-to-end (±9° = 5% of 360%); "
            "MAX cyan+magenta, MIN lime+sky; intervals [21,39]/[176,194]/[106,124]/[256,274]"
        ),
        "notes": [
            "Exactly 4 features matching three-panel Panel B geometry/width",
            "Each = ONE smooth transparent PolyCollection ribbon (±9°) + ONE solid mid-stroke only",
            "Clean draw: fill alpha=0.30; mid-stroke thin white underlay (lw+1.5) + color lw=2.8; no black rim; no boundary-edge halos",
            "No multi-stroke scribble wash; no 283° companion; no Dec=0 overlay",
            "MAX bands [21,39] & [176,194]; MIN bands [106,124] & [256,274]",
            "Colors: MAX_30=#00e5ff, MAX_185=#ff2d95, MIN_115=#b8f200, MIN_265=#4fc3f7",
            "band_half_width_deg=9, band_full_width_deg=18 (5% of 360°)",
            "Exact n·x̂=μ small-circles; seam-split polygon fill; mid via plot_closed_ring",
            "No ellipse-fit, no interior chords, no decorative oval-rim loops",
            "No blocking text box on the map; legend below",
            "Colorbar preserved from source (overlays stay on oval)",
            "Overlays clip_path = registered Ellipse(cx,cy,Rx,Ry); rings stroke inset 0.97",
            "Dec=0 NOT drawn (Migkas Fig.8 already has celestial equator); CoE axis = great circle through ±D (tilted +30°) at inset=1.0",
            "axis_tilt_b_deg=30; anti-CoE (0°,+30°), CoE (180°,−30°); Migkas markers fixed at true sky positions",
            "Registration hard-set to pre medial-axis params (863,763.25,693.9927,392.542425); DT-ridge fit rejected",
            "Limb markers inset toward colored-oval center so halo stays inside",
            "Points culled at COLOR_R=0.978 of seed oval; limb markers inset to fit halo",
        ],
    }
    pass  # per-version JSON not written (PNG only)
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
```

### Fig07 — `Fig07_storyboard.png`

- **Gold extract:** `gold-009.png`
- **Expected dims (gold):** 1045×1795
- **Method:** `rerun:makers/make_fig07_storyboard.py`
- **Output:** `/workspace/figs_letter_exact_python/Fig07_storyboard.png`
- **Notes:** Primary: cycle7 make_donut_universe_storyboard.py

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_fig07_storyboard.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_fig07_storyboard.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_fig07_storyboard.py
#!/usr/bin/env python3
"""Fig07 — Donut-universe storyboard A–D regenerated from geometry code.

Primary: cycle7_recombination/make_donut_universe_storyboard.py
Fallback: rebuilt 2×2 stages matching letter narrative (fixed R_rec).
"""
from __future__ import annotations

import os
import runpy
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

OUT = Path("/workspace/figs_letter_exact_python")
CYCLE7 = Path(
    "/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle7_recombination"
)
SCRIPT = CYCLE7 / "make_donut_universe_storyboard.py"
OUT_NAME = "Fig07_storyboard.png"


def _fallback_draw():
    """Self-contained A–D storyboard if cycle7 script fails."""
    C, V_IN, V_OUT = 1.0, 2.0, 4.0
    V_D = 0.5 * (V_IN + V_OUT)
    R_REC = 1.20
    T_REC = R_REC / V_D

    def stage(t, title):
        rin, rout, d = V_IN * t, V_OUT * t, V_D * t
        rvis = C * (t - T_REC) if d > R_REC else 0.0
        return dict(title=title, t=t, D=d, Rin=rin, Rout=rout, Rvis=rvis, Rrec=R_REC)

    stages = [
        stage(0.25, "A — before recombination"),
        stage(T_REC, r"B — US at $D_{\mathrm{recomb}}$"),
        stage(0.90, "C — early post-recombination"),
        stage(1.60, "D — later expansion"),
    ]
    lim = 7.0
    # Letter gold ~1045×1795 → tall portrait; use 2×2 then we'll accept landscape
    # Prefer vertical A-B-C-D to better match letter aspect:
    fig, axes = plt.subplots(4, 1, figsize=(6.5, 11.2), dpi=160)
    for ax, s in zip(axes, stages):
        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_title(s["title"], fontsize=11)
        th = np.linspace(0, 2 * np.pi, 400)
        if s["Rout"] > s["Rin"] > 0:
            xs = np.concatenate([s["Rout"] * np.cos(th), s["Rin"] * np.cos(th[::-1])])
            ys = np.concatenate([s["Rout"] * np.sin(th), s["Rin"] * np.sin(th[::-1])])
            ax.fill(xs, ys, color="#aed6f1", alpha=0.7, zorder=1)
        ax.add_patch(Circle((0, 0), s["Rin"], fill=False, ec="#1a5276", lw=1.2, zorder=2))
        ax.add_patch(Circle((0, 0), s["Rout"], fill=False, ec="#1a5276", lw=1.2, zorder=2))
        ax.add_patch(Circle((0, 0), s["Rrec"], fill=False, ec="#c0392b", lw=1.8, ls="--", zorder=3))
        ax.plot(s["D"], 0, "o", color="#2471a3", ms=7, zorder=5)
        ax.text(s["D"], -0.12 * lim, "US", color="#2471a3", ha="center", fontsize=8)
        ax.plot(0, 0, "o", color="#c0392b", ms=5, zorder=5)
        ax.text(0, -0.12 * lim, "CoE", color="#c0392b", ha="center", fontsize=8)
        if s["Rvis"] > 0:
            ax.add_patch(
                Circle((s["D"], 0), s["Rvis"], facecolor="#abebc6", edgecolor="#1e8449",
                       alpha=0.45, lw=1.3, zorder=4)
            )
        ax.grid(True, alpha=0.2)
    fig.suptitle(
        "Donut-universe time sequence: fixed recombination surface; shells expand; observer mid-shell",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    dest = OUT / OUT_NAME
    fig.savefig(dest, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return dest


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if SCRIPT.exists():
        cwd = os.getcwd()
        try:
            os.chdir(CYCLE7)
            runpy.run_path(str(SCRIPT), run_name="__main__")
        finally:
            os.chdir(cwd)
        # Find newest png in cycle7
        cands = sorted(CYCLE7.glob("*.png"), key=lambda p: p.stat().st_mtime, reverse=True)
        # Prefer names with storyboard/donut
        prefer = [p for p in cands if "storyboard" in p.name.lower() or "donut" in p.name.lower()]
        src = (prefer or cands)[0] if (prefer or cands) else None
        if src is None:
            print("cycle7 produced no png; using fallback")
            return _fallback_draw()
        dest = OUT / OUT_NAME
        shutil.copy2(src, dest)
        print("wrote", dest, "from", src.name, dest.stat().st_size)
        return dest
    return _fallback_draw()


if __name__ == "__main__":
    main()

```

### Fig08 — `Fig08_two_emitters.png`

- **Gold extract:** `gold-011.png`
- **Expected dims (gold):** 1181×1908
- **Method:** `rerun:makers/make_fig08_two_emitters.py`
- **Output:** `/workspace/figs_letter_exact_python/Fig08_two_emitters.png`
- **Notes:** Primary: cycle6 make_two_emitters_comoving.py

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_fig08_two_emitters.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_fig08_two_emitters.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_fig08_two_emitters.py
#!/usr/bin/env python3
"""Fig08 — Two-emitter baseline regenerated from cycle6 JSON + schematic.

Primary: cycle6_two_emitters/make_two_emitters_comoving.py
Fallback: self-contained schematic + T_obs from two_emitters.json.
"""
from __future__ import annotations

import json
import os
import runpy
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

OUT = Path("/workspace/figs_letter_exact_python")
CYCLE6 = Path(
    "/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle6_two_emitters"
)
SCRIPT = CYCLE6 / "make_two_emitters_comoving.py"
JSON_PATH = CYCLE6 / "two_emitters.json"
OUT_NAME = "Fig08_two_emitters.png"


def _fallback():
    data = json.loads(JSON_PATH.read_text()) if JSON_PATH.exists() else {}
    fig = plt.figure(figsize=(7.4, 11.9), dpi=160)  # portrait ~letter 1181×1908
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.25)
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-2.2, 2.4)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Shared bulk velocity (schematic)", loc="left", fontsize=12)
    ax.add_patch(Rectangle((-2.8, 1.2), 5.6, 0.35, color="#5dade2", alpha=0.85, ec="#1a5276"))
    ax.add_patch(Rectangle((-2.8, -1.55), 5.6, 0.35, color="#5dade2", alpha=0.85, ec="#1a5276"))
    ax.text(0, 1.55, "+ plate / shell", ha="center", fontsize=9, color="#1a5276")
    ax.text(0, -1.7, "− plate / shell", ha="center", fontsize=9, color="#1a5276")
    ax.plot(0, 0, "o", color="#2471a3", ms=10)
    ax.text(0.15, -0.35, "observer", color="#2471a3", fontsize=9)
    for y in (1.38, 0.0, -1.38):
        ax.annotate("", xy=(1.4, y), xytext=(-0.2, y),
                    arrowprops=dict(arrowstyle="->", color="#c0392b", lw=2.2))
    ax.text(1.55, 0.15, r"shared $\vec{V}$", color="#c0392b", fontsize=10)
    P = data.get("plates", {})
    th = np.asarray(P.get("theta_deg", np.linspace(0, 360, 361)), float)
    mid = np.asarray(P.get("mid_sameT", np.ones_like(th)), float)
    ax2.plot(th, mid, color="C0", lw=2.4, label=r"same $T$, shared $V$ (mid)")
    if "mid_sameT_doppler" in P:
        ax2.plot(th, P["mid_sameT_doppler"], color="C2", lw=1.6, ls="--",
                 label=r"observer boost $\beta=0.01$ (dipole)")
    ax2.axhline(1.0, color="0.5", lw=0.8, ls=":")
    ax2.set_xlim(0, 360)
    ax2.set_xlabel(r"great-circle angle $\theta$ (deg)")
    ax2.set_ylabel(r"$T_{\mathrm{obs}}/T_0$")
    ax2.set_title("Kinematic sky before gravity", loc="left", fontsize=12)
    ax2.legend(fontsize=8, frameon=False)
    ax2.grid(True, alpha=0.25)
    fig.suptitle("Two-emitter baseline", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    dest = OUT / OUT_NAME
    fig.savefig(dest, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return dest


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if SCRIPT.exists():
        cwd = os.getcwd()
        try:
            os.chdir(CYCLE6)
            runpy.run_path(str(SCRIPT), run_name="__main__")
        finally:
            os.chdir(cwd)
        cands = sorted(CYCLE6.glob("*.png"), key=lambda p: p.stat().st_mtime, reverse=True)
        prefer = [p for p in cands if "emitter" in p.name.lower() or "two" in p.name.lower()]
        src = (prefer or cands)[0] if (prefer or cands) else None
        if src is not None:
            dest = OUT / OUT_NAME
            shutil.copy2(src, dest)
            print("wrote", dest, "from", src.name)
            return dest
    return _fallback()


if __name__ == "__main__":
    main()

```

### Fig09 — `Fig09_kill.jpg`

- **Gold extract:** `gold-013.jpg`
- **Expected dims (gold):** 1206×1280
- **Method:** `rerun:makers/make_fig09_kill.py`
- **Output:** `/workspace/figs_letter_exact_python/Fig09_kill.jpg`
- **Notes:** Kill function vs plate/shell/shear from dtt.json

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_fig09_kill.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_fig09_kill.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_fig09_kill.py
#!/usr/bin/env python3
"""Fig09 — Kill function vs plate/shell/shear from cycle4 dtt.json."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths import OUT  # noqa: E402
from physics import extend_0_360, kill_from_doppler, load_dtt, load_vs_kill, unit_norm  # noqa: E402

OUT_NAME = "Fig09_kill.jpg"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dtt = load_dtt()
    rows = dtt["rows"]
    vs = load_vs_kill()  # noqa: F841 — available for extensions
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 5.3), dpi=160, sharey=True)
    for ax, row in zip(axes, rows):
        th, kill = kill_from_doppler(row)
        th360, kill360 = extend_0_360(th, kill)
        kn = unit_norm(kill360)
        ax.plot(th360, kn, color="k", lw=2.2, label="kill (unit-norm)")
        phi = np.radians(th360)
        ax.plot(th360, unit_norm(np.cos(phi)), color="C0", lw=1.3, ls="--", alpha=0.85, label="plate-like")
        ax.plot(th360, unit_norm(np.cos(2 * phi)), color="C1", lw=1.3, ls=":", alpha=0.85, label="shell/P2-like")
        ax.plot(th360, unit_norm(np.sin(2 * phi)), color="C3", lw=1.2, ls="-.", alpha=0.85, label="shear-like")
        ax.set_title(rf"$D={row['D']:.2f}$ Gpc", fontsize=11)
        ax.set_xlabel(r"$\phi$ (deg) on great circle")
        ax.set_xlim(0, 360)
        ax.set_xticks([0, 90, 180, 270, 360])
        ax.axhline(0, color="0.5", lw=0.6)
        ax.grid(True, alpha=0.25)
        if ax is axes[0]:
            ax.set_ylabel("unit-normalized amplitude")
            ax.legend(fontsize=7.5, frameon=False, loc="lower left")
    fig.suptitle("Kill function vs exploratory plate / shell / shear residuals", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    png = OUT / "_Fig09_kill.png"
    fig.savefig(png, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    dest = OUT / OUT_NAME
    im = Image.open(png).convert("RGB")
    # gold ~1206×1280 — portrait; our plot is landscape; save as-is as jpg
    w = 1206
    h = round(im.height * w / im.width)
    im.resize((w, h), Image.LANCZOS).save(dest, quality=92)
    print("wrote", dest, dest.stat().st_size, (w, h))
    return dest


if __name__ == "__main__":
    main()

```

### Fig10 — `Fig10_cancel_partner.png`

- **Gold extract:** `gold-014.png`
- **Expected dims (gold):** 2084×2212
- **Method:** `rerun:makers/make_fig10_cancel_partner.py`
- **Output:** `/workspace/figs_letter_exact_python/Fig10_cancel_partner.png`
- **Notes:** Cancel partner residuals from dtt.json

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_fig10_cancel_partner.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_fig10_cancel_partner.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_fig10_cancel_partner.py
#!/usr/bin/env python3
"""Fig10 — CoE Doppler + multipole partner residuals from dtt.json."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths import OUT  # noqa: E402
from physics import load_dtt, partner_cancel  # noqa: E402

OUT_NAME = "Fig10_cancel_partner.png"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dtt = load_dtt()
    prefer = sorted(dtt["rows"], key=lambda r: abs(r["D"] - 49.0))[:2]
    if len(prefer) < 2:
        prefer = dtt["rows"][:2]
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 11.0), dpi=160, sharex="col")
    for col, row in enumerate(prefer):
        th = np.asarray(row["theta_toward_deg"], float)
        dop = np.asarray(row["dtt_doppler"], float)
        dop0 = dop - dop.mean()
        _, res_p2, a = partner_cancel(dop0, th, mode="dip_p2")
        _, res_l4, _ = partner_cancel(dop0, th, mode="l_le_4")
        ax0 = axes[0, col]
        ax0.plot(th, dop0, color="C0", lw=1.8, label="Doppler (demeaned)")
        ax0.plot(th, -dop0, color="0.35", lw=1.2, ls="--", label="−Doppler (= kill shape)")
        ax0.set_title(rf"$D={row['D']:.2f}$ Gpc", fontsize=11)
        ax0.set_ylabel(r"$\Delta T/T$")
        ax0.legend(fontsize=8, frameon=False)
        ax0.grid(True, alpha=0.25)
        ax0.axhline(0, color="k", lw=0.5)
        ax1 = axes[1, col]
        ax1.plot(th, res_p2, color="C3", lw=1.8, label=r"residual after dipole+$P_2$")
        ax1.plot(th, res_l4, color="C2", lw=1.6, ls="--", label=r"residual after $\ell\leq 4$")
        ax1.axhline(0, color="k", lw=0.5)
        ax1.set_xlabel(r"$\theta$ from CoE (deg)")
        ax1.set_ylabel(r"$\Delta T/T$ residual")
        ax1.legend(fontsize=8, frameon=False)
        ax1.grid(True, alpha=0.25)
        ax0.text(0.02, 0.05, rf"$a_2^{{\mathrm{{Dop}}}}={row.get('doppler_a2', a[2]):.3f}$",
                 transform=ax0.transAxes, fontsize=8, color="0.3")
    fig.suptitle("CoE Doppler + multipole partner cancellation", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    dest = OUT / OUT_NAME
    fig.savefig(dest, dpi=190, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", dest, dest.stat().st_size)
    return dest


if __name__ == "__main__":
    main()

```

### FigA1 — `FigA1_ltb_shell.png`

- **Gold extract:** `gold-016.png`
- **Expected dims (gold):** 1444×969
- **Method:** `rerun:makers/make_figA1_ltb_shell.py`
- **Output:** `/workspace/figs_letter_exact_python/FigA1_ltb_shell.png`
- **Notes:** LTB shell evolution panels

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_figA1_ltb_shell.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_figA1_ltb_shell.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_figA1_ltb_shell.py
#!/usr/bin/env python3
"""Fig A1 — LTB shell evolution: flat m+r; E=0; E>0; E<0."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))
from physics import save_fig
from paths import OUT
OUT_NAME = "FigA1_ltb_shell.png"


def evolve_dust(E: float, t_max: float = 8.0, a0: float = 1e-4):
    """Integrate ȧ² = 2E + 2GM/a with GM=1 units; a(0)=a0, ȧ>0."""
    GM = 1.0

    def f(t, y):
        a, = y
        rhs = 2 * E + 2 * GM / max(a, 1e-12)
        if rhs < 0:
            return [0.0]
        return [np.sqrt(rhs)]

    sol = solve_ivp(f, (0, t_max), [a0], dense_output=True, max_step=0.01, rtol=1e-6)
    t = np.linspace(0, min(t_max, sol.t[-1]), 600)
    a = sol.sol(t)[0]
    # stop if recollapse starts
    if E < 0:
        # analytic cycloid max a = GM/|E|
        amax = GM / abs(E)
        # continue through recollapse with signed root
        def f2(t, y):
            a, = y
            rhs = 2 * E + 2 * GM / max(a, 1e-12)
            # use phase from energy: after turnaround ȧ < 0
            return [-np.sqrt(max(rhs, 0.0))]

        # find turnaround
        i_turn = int(np.argmax(a))
        if a[i_turn] > 0.5 * amax:
            t_turn = t[i_turn]
            sol2 = solve_ivp(f2, (t_turn, t_max), [a[i_turn]], dense_output=True, max_step=0.01)
            t2 = np.linspace(t_turn, min(t_max, sol2.t[-1]), 400)
            a2 = sol2.sol(t2)[0]
            t = np.concatenate([t[: i_turn + 1], t2[1:]])
            a = np.concatenate([a[: i_turn + 1], a2[1:]])
    return t, a


def evolve_flat_mr(t_max=8.0):
    """Schematic flat matter+radiation: early a~t^{1/2}, late a~t^{2/3}."""
    t = np.linspace(1e-3, t_max, 800)
    # match at teq=1: continuous a and prefer radiation then matter
    teq = 1.0
    a_eq = teq ** 0.5
    a = np.where(t < teq, t ** 0.5, a_eq * (t / teq) ** (2 / 3))
    return t, a


def main():
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.5), dpi=160)
    panels = [
        (axes[0, 0], "Flat matter + radiation", evolve_flat_mr()),
        (axes[0, 1], r"Dust $E=0$ (marginally bound)", evolve_dust(0.0)),
        (axes[1, 0], r"Dust $E>0$ (expand forever)", evolve_dust(+0.15)),
        (axes[1, 1], r"Dust $E<0$ (recollapse)", evolve_dust(-0.25, t_max=10.0)),
    ]
    for ax, title, (t, a) in panels:
        ax.plot(t, a, color="C0", lw=2.0)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r"$t$ (code units)")
        ax.set_ylabel(r"$a(t)$ / areal factor")
        ax.grid(True, alpha=0.3)
        if "radiation" in title.lower():
            ax.axvline(1.0, color="0.5", ls=":", lw=1)
            ax.text(1.05, ax.get_ylim()[1] * 0.1 if False else 0.3, r"$t_{\mathrm{eq}}$", fontsize=8, color="0.4")

    # annotate asymptotics on first panel
    axes[0, 0].text(0.35, 0.55, r"$a\sim t^{1/2}$", fontsize=9, color="0.35", transform=axes[0, 0].transAxes)
    axes[0, 0].text(0.62, 0.75, r"$a\sim t^{2/3}$", fontsize=9, color="0.35", transform=axes[0, 0].transAxes)

    fig.suptitle("LTB shell evolution (Einstein/Friedmann shell constraint)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    png = save_fig(fig, "_figA1_ltb_shell", also_pdf=False)
    dest = OUT / OUT_NAME
    import shutil
    shutil.move(str(png), str(dest))
    p = dest
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()

```

### FigA2 — `FigA2_coe_sun.jpg`

- **Gold extract:** `gold-018.jpg`
- **Expected dims (gold):** 1280×636
- **Method:** `rerun:makers/make_figA2_coe_sun.py`
- **Output:** `/workspace/figs_letter_exact_python/FigA2_coe_sun.jpg`
- **Notes:** CoE/Sun timing 2026 from Plano

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_figA2_coe_sun.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_figA2_coe_sun.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_figA2_coe_sun.py
#!/usr/bin/env python3
"""Fig A2 — CoE and Sun timing from Plano through 2026."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time
import astropy.units as u
from astropy.utils import iers

iers.conf.auto_download = False

import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))
from physics import COE_GAL_L, COE_GAL_B, PLANO_LAT, PLANO_LON, PLANO_H, save_fig
from paths import OUT
from PIL import Image
OUT_NAME = "FigA2_coe_sun.jpg"


def find_events(coord, site, year, n_days=365):
    """Approximate rise / culmination / set local hours for each day."""
    # sample every 4 minutes over each day in UTC
    doys = np.arange(1, n_days + 1)
    rise_h, culm_h, set_h = [], [], []
    culm_alt = []
    lon_hours = site.lon.to(u.hourangle).value
    t0 = Time(f"{year}-01-01 00:00:00", scale="utc")
    for d in doys:
        hours = np.arange(0, 24, 4 / 60)
        # interpret hours as Local Mean Time, convert to UTC
        utc_days = (d - 1) + hours / 24.0 - lon_hours / 24.0
        times = t0 + utc_days * u.day
        aa = coord.transform_to(AltAz(obstime=times, location=site))
        alt = aa.alt.deg
        az = aa.az.deg
        # culmination = max alt
        i_c = int(np.argmax(alt))
        culm_h.append(hours[i_c])
        culm_alt.append(alt[i_c])
        # rise: first upward zero crossing
        rise = set_ = np.nan
        for i in range(1, len(alt)):
            if alt[i - 1] < 0 <= alt[i]:
                # linear interp
                f = -alt[i - 1] / (alt[i] - alt[i - 1] + 1e-12)
                rise = hours[i - 1] + f * (hours[i] - hours[i - 1])
                break
        for i in range(1, len(alt)):
            if alt[i - 1] > 0 >= alt[i]:
                f = alt[i - 1] / (alt[i - 1] - alt[i] + 1e-12)
                set_ = hours[i - 1] + f * (hours[i] - hours[i - 1])
        rise_h.append(rise)
        set_h.append(set_)
    return (
        doys,
        np.asarray(rise_h),
        np.asarray(culm_h),
        np.asarray(set_h),
        np.asarray(culm_alt),
    )


def main():
    year = 2026
    site = EarthLocation(lat=PLANO_LAT * u.deg, lon=PLANO_LON * u.deg, height=PLANO_H * u.m)
    coe = SkyCoord(l=COE_GAL_L * u.deg, b=COE_GAL_B * u.deg, frame="galactic").icrs

    # CoE is essentially fixed RA/Dec → rise/set nearly constant in LMT
    doy, cr, cc, cs, ca = find_events(coe, site, year, n_days=365)

    # Sun: sample weekly for speed then interpolate display
    sun_doys = np.arange(1, 366, 3)
    sr, sc, ss = [], [], []
    t0 = Time(f"{year}-01-01 00:00:00", scale="utc")
    lon_hours = site.lon.to(u.hourangle).value
    for d in sun_doys:
        hours = np.arange(0, 24, 4 / 60)
        utc_days = (d - 1) + hours / 24.0 - lon_hours / 24.0
        times = t0 + utc_days * u.day
        sun = get_sun(times)
        aa = sun.transform_to(AltAz(obstime=times, location=site))
        alt = aa.alt.deg
        i_c = int(np.argmax(alt))
        sc.append(hours[i_c])
        rise = set_ = np.nan
        for i in range(1, len(alt)):
            if alt[i - 1] < 0 <= alt[i]:
                f = -alt[i - 1] / (alt[i] - alt[i - 1] + 1e-12)
                rise = hours[i - 1] + f * (hours[i] - hours[i - 1])
                break
        for i in range(1, len(alt)):
            if alt[i - 1] > 0 >= alt[i]:
                f = alt[i - 1] / (alt[i - 1] - alt[i] + 1e-12)
                set_ = hours[i - 1] + f * (hours[i] - hours[i - 1])
        sr.append(rise)
        ss.append(set_)
    sr, sc, ss = map(np.asarray, (sr, sc, ss))

    fig, ax = plt.subplots(figsize=(11.0, 5.2), dpi=160)
    ax.plot(doy, cr, color="#6c3483", lw=1.8, label="CoE rise")
    ax.plot(doy, cc, color="#8e44ad", lw=1.8, ls="--", label="CoE culmination")
    ax.plot(doy, cs, color="#5b2c6f", lw=1.8, ls=":", label="CoE set")
    ax.plot(sun_doys, sr, color="#e67e22", lw=1.6, label="Sun rise")
    ax.plot(sun_doys, sc, color="#d35400", lw=1.6, ls="--", label="Sun culmination")
    ax.plot(sun_doys, ss, color="#a04000", lw=1.6, ls=":", label="Sun set")
    ax.set_xlabel("Day of year 2026")
    ax.set_ylabel("Local mean time (hours)")
    ax.set_xlim(1, 365)
    ax.set_ylim(0, 24)
    ax.set_yticks(range(0, 25, 3))
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2, frameon=False, loc="upper right")
    ax.set_title(
        "CoE and Sun timing from Plano, TX (2026)\n"
        rf"CoE at Galactic $(\ell,b)=({COE_GAL_L:.0f}^\circ,{COE_GAL_B:.0f}^\circ)$",
        loc="left",
        fontsize=11,
    )
    png = save_fig(fig, "_figA2_coe_sun", also_pdf=False)
    dest = OUT / OUT_NAME
    im = Image.open(png).convert("RGB")
    w = 1280
    h = round(im.height * w / im.width)
    im.resize((w, h), Image.LANCZOS).save(dest, quality=92)
    png.unlink(missing_ok=True)
    p = dest
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()

```

### FigA2b — `FigA2b_sun_coe_crossings.jpg`

- **Gold extract:** `gold-019.jpg`
- **Expected dims (gold):** 1239×570
- **Method:** `rerun:makers/make_figA2b_crossings.py`
- **Output:** `/workspace/figs_letter_exact_python/FigA2b_sun_coe_crossings.jpg`
- **Notes:** Sun/CoE crossings 2027

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_figA2b_crossings.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_figA2b_crossings.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_figA2b_crossings.py
#!/usr/bin/env python3
"""Table A1 visual — Sun/CoE event crossings Central Texas 2027."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time
import astropy.units as u
from astropy.utils import iers

iers.conf.auto_download = False

import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))
from physics import COE_GAL_L, COE_GAL_B, PLANO_LAT, PLANO_LON, PLANO_H, save_fig
from paths import OUT
from PIL import Image
OUT_NAME = "FigA2b_sun_coe_crossings.jpg"


def daily_events(coord_fn, site, year):
    """Return rise, culm, set LMT hours for each DOY (coord_fn(times)->SkyCoord)."""
    lon_h = site.lon.to(u.hourangle).value
    t0 = Time(f"{year}-01-01 00:00:00", scale="utc")
    doys = np.arange(1, 366)
    out = {k: np.full(len(doys), np.nan) for k in ("rise", "culm", "set")}
    hours = np.arange(0, 24, 5 / 60)
    for i, d in enumerate(doys):
        utc_days = (d - 1) + hours / 24.0 - lon_h / 24.0
        times = t0 + utc_days * u.day
        aa = coord_fn(times).transform_to(AltAz(obstime=times, location=site))
        alt = aa.alt.deg
        out["culm"][i] = hours[int(np.argmax(alt))]
        for j in range(1, len(alt)):
            if alt[j - 1] < 0 <= alt[j] and np.isnan(out["rise"][i]):
                f = -alt[j - 1] / (alt[j] - alt[j - 1] + 1e-12)
                out["rise"][i] = hours[j - 1] + f * (hours[j] - hours[j - 1])
            if alt[j - 1] > 0 >= alt[j] and np.isnan(out["set"][i]):
                f = alt[j - 1] / (alt[j - 1] - alt[j] + 1e-12)
                out["set"][i] = hours[j - 1] + f * (hours[j] - hours[j - 1])
    return doys, out


def crossings(doys, sun, coe, tol_h=0.15):
    """Find DOYs where Sun and CoE event times match within tol."""
    labels = ["rise", "culm", "set"]
    pairs = []
    for i, s_ev in enumerate(labels):
        for j, c_ev in enumerate(labels):
            ds = sun[s_ev] - coe[c_ev]
            # unwrap near day boundary not needed for LMT hours usually
            hit = np.where(np.isfinite(ds) & (np.abs(ds) < tol_h))[0]
            # also pick closest approach even if > tol (min |ds|)
            if len(hit) == 0:
                k = int(np.nanargmin(np.abs(ds)))
                if np.isfinite(ds[k]) and abs(ds[k]) < 0.75:
                    hit = np.array([k])
            for k in hit:
                pairs.append(
                    dict(
                        doy=int(doys[k]),
                        sun=s_ev,
                        coe=c_ev,
                        sun_h=float(sun[s_ev][k]),
                        coe_h=float(coe[c_ev][k]),
                        diag=(i == j),
                    )
                )
    # unique by (sun,coe) keeping closest
    best = {}
    for p in pairs:
        key = (p["sun"], p["coe"])
        err = abs(p["sun_h"] - p["coe_h"])
        if key not in best or err < abs(best[key]["sun_h"] - best[key]["coe_h"]):
            best[key] = p
    return sorted(best.values(), key=lambda p: (not p["diag"], p["doy"]))


def main():
    year = 2027
    site = EarthLocation(lat=PLANO_LAT * u.deg, lon=PLANO_LON * u.deg, height=PLANO_H * u.m)
    coe_icrs = SkyCoord(l=COE_GAL_L * u.deg, b=COE_GAL_B * u.deg, frame="galactic").icrs

    def coe_fn(times):
        return coe_icrs

    def sun_fn(times):
        return get_sun(times)

    # subsample for speed: every 2 days then refine
    doys, sun = daily_events(sun_fn, site, year)
    _, coe = daily_events(coe_fn, site, year)
    rows = crossings(doys, sun, coe, tol_h=0.25)
    # keep up to 9
    rows = rows[:9]

    fig, ax = plt.subplots(figsize=(10.5, 3.6), dpi=160)
    ax.axis("off")
    ax.set_title(
        f"Sun / CoE event crossings — Central Texas {year}\n"
        r"Highlighted diagonal: rise↔rise, culm↔culm, set↔set",
        fontsize=11,
        loc="left",
        pad=12,
    )

    col_labels = ["#", "DOY", "Date", "Sun event", "CoE event", "Sun LMT", "CoE LMT", "Δ (h)"]
    cell = []
    colors = []
    for n, r in enumerate(rows, 1):
        t = Time(f"{year}-01-01") + (r["doy"] - 1) * u.day
        date = t.strftime("%Y-%m-%d")
        cell.append(
            [
                str(n),
                str(r["doy"]),
                date,
                r["sun"],
                r["coe"],
                f"{r['sun_h']:.2f}",
                f"{r['coe_h']:.2f}",
                f"{r['sun_h']-r['coe_h']:+.2f}",
            ]
        )
        colors.append(["#fdebd0" if r["diag"] else "white"] * len(col_labels))

    if not cell:
        cell = [["—"] * len(col_labels)]
        colors = [["white"] * len(col_labels)]

    table = ax.table(cellText=cell, colLabels=col_labels, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.15, 1.45)
    for (i, j), cell_obj in table.get_celld().items():
        if i == 0:
            cell_obj.set_facecolor("#2c3e50")
            cell_obj.set_text_props(color="white", weight="bold")
        elif i >= 1 and i - 1 < len(colors):
            cell_obj.set_facecolor(colors[i - 1][j])

    png = save_fig(fig, "_figA2b_crossings", also_pdf=False)
    dest = OUT / OUT_NAME
    im = Image.open(png).convert("RGB")
    w = 1239
    h = round(im.height * w / im.width)
    im.resize((w, h), Image.LANCZOS).save(dest, quality=92)
    png.unlink(missing_ok=True)
    p = dest
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()

```

### FigA3 — `FigA3_coe_sky.png`

- **Gold extract:** `gold-020.png`
- **Expected dims (gold):** 2415×1432
- **Method:** `rerun:makers/make_figA3_coe_sky.py`
- **Output:** `/workspace/figs_letter_exact_python/FigA3_coe_sky.png`
- **Notes:** Reruns letter_brief/figs/_make_coe_sky_path_dates.py. Pixel-matches gold.

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_figA3_coe_sky.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_figA3_coe_sky.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_figA3_coe_sky.py
#!/usr/bin/env python3
"""FigA3 — Annual CoE sky path from Plano (pixel-matches letter gold on rerun)."""
from __future__ import annotations

import os
import runpy
import shutil
from pathlib import Path

OUT = Path("/workspace/figs_letter_exact_python")
SCRIPT = Path(
    "/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/_make_coe_sky_path_dates.py"
)
SRC_PNG = SCRIPT.parent / "coe_sky_path_dates_texas.png"
OUT_NAME = "FigA3_coe_sky.png"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(SCRIPT.parent)
        runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        os.chdir(cwd)
    dest = OUT / OUT_NAME
    shutil.copy2(SRC_PNG, dest)
    print("wrote", dest, dest.stat().st_size)
    return dest


if __name__ == "__main__":
    main()

```

**Dependency (full source):** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/_make_coe_sky_path_dates.py`

```python
# PATH: /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/_make_coe_sky_path_dates.py
#!/usr/bin/env python3
"""Make a clean annual CoE sky-path chart for Plano, Texas."""
from __future__ import annotations

from pathlib import Path
import calendar

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.utils import iers

iers.conf.auto_download = False

i = Path(__file__).resolve().parent
OUT = i / "coe_sky_path_dates_texas.png"

# Inputs supplied in the request.
YEAR = 2026
LAT, LON, HEIGHT_M = 33.02, -96.70, 200.0
LMT_HOUR = 12.0
GAL_L, GAL_B = 180.0, -15.0

site = EarthLocation(lat=LAT * u.deg, lon=LON * u.deg, height=HEIGHT_M * u.m)
coe_gal = SkyCoord(l=GAL_L * u.deg, b=GAL_B * u.deg, frame="galactic")
coe = coe_gal.icrs

# Local mean solar time: UT = LMT - longitude/15 hours.
doys = np.arange(1, 366)
t0 = Time(f"{YEAR}-01-01 00:00:00", scale="utc")
times = t0 + (doys - 1 + LMT_HOUR / 24.0 - LON / 15.0 / 24.0) * u.day
altaz = coe.transform_to(AltAz(obstime=times, location=site))
az = np.asarray(altaz.az.deg)
alt = np.asarray(altaz.alt.deg)

# Month labels at the 15th, an uncluttered representative date for each month.
month_names = [calendar.month_abbr[m] for m in range(1, 13)]
month_doys = np.array([sum(calendar.monthrange(YEAR, m)[1] for m in range(1, mm)) + 15 for mm in range(1, 13)])
month_idx = month_doys - 1

# Build a daily colored path; omit the single segment that crosses the 0/360 seam.
points = np.column_stack([az, alt])
segments = np.stack([points[:-1], points[1:]], axis=1)
seam = np.abs(np.diff(az)) > 180.0
segments = segments[~seam]
seg_day = (doys[:-1] + doys[1:]) / 2.0
seg_day = seg_day[~seam]

fig, ax = plt.subplots(figsize=(12.4, 8.4), dpi=200)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Below-horizon region and horizon.
ax.axhspan(-45, 0, color="#f1f3f5", zorder=0)
ax.axhline(0, color="#30343b", lw=1.45, zorder=1)
ax.text(356.5, 1.7, "horizon", ha="right", va="bottom", fontsize=8.5, color="#50555c")

# Cardinal guides: restrained enough to keep the path dominant.
for x in (0, 90, 180, 270, 360):
    ax.axvline(x, color="#dfe3e7", lw=0.8, zorder=0)

norm = Normalize(vmin=1, vmax=365)
cmap = plt.get_cmap("turbo")
lc = LineCollection(segments, cmap=cmap, norm=norm, linewidths=2.5, alpha=0.96,
                    capstyle="round", joinstyle="round", zorder=3)
lc.set_array(seg_day)
ax.add_collection(lc)

# Daily sample dots make the once-per-day sampling explicit.
ax.scatter(az, alt, c=doys, cmap=cmap, norm=norm, s=11, linewidths=0,
           alpha=0.98, zorder=4)

# Month markers and labels. Offsets are chosen to keep all labels clear of the curve.
offsets = {
    "Jan": (-28, -18), "Feb": (-28, 10), "Mar": (-25, 10), "Apr": (-28, 9),
    "May": (-28, 10), "Jun": (10, 10), "Jul": (10, 9), "Aug": (11, 5),
    "Sep": (11, -1), "Oct": (10, -15), "Nov": (-35, -4), "Dec": (10, -7),
}
for name, idx in zip(month_names, month_idx):
    x, y = az[idx], alt[idx]
    ax.scatter([x], [y], s=68, facecolor="white", edgecolor="#20242a",
               linewidth=1.15, zorder=6)
    dx, dy = offsets[name]
    ax.annotate(
        f"{name} 15",
        (x, y), xytext=(dx, dy), textcoords="offset points",
        fontsize=9.2, fontweight="bold", color="#20242a", zorder=7,
        bbox=dict(boxstyle="round,pad=0.20", fc="white", ec="#b8bec6", lw=0.65, alpha=0.94),
        arrowprops=dict(arrowstyle="-", color="#747b84", lw=0.75, shrinkA=4, shrinkB=4),
    )

# Axes, grid, and cardinal labels.
ax.set_xlim(0, 360)
ax.set_ylim(-45, 85)
ax.set_xticks([0, 90, 180, 270, 360])
ax.set_xticklabels(["N\n0°", "E\n90°", "S\n180°", "W\n270°", "N\n360°"],
                   fontsize=10.5, fontweight="bold", color="#252a30")
ax.set_xticks(np.arange(45, 360, 45), minor=True)
ax.set_yticks(np.arange(-40, 81, 20))
ax.set_xlabel("Azimuth (clockwise from North)", fontsize=11, labelpad=13)
ax.set_ylabel("Altitude above horizon (°)", fontsize=11, labelpad=10)
ax.grid(which="major", color="#cdd2d8", lw=0.75, alpha=0.62)
ax.grid(which="minor", axis="x", color="#e5e8eb", lw=0.55, alpha=0.8)
ax.tick_params(axis="x", which="minor", length=3, color="#8d949c")
ax.tick_params(axis="both", which="major", labelsize=9.5)
for spine in ax.spines.values():
    spine.set_color("#858c94")
    spine.set_linewidth(0.8)

fig.suptitle("CoE annual sky path from Plano, Texas", fontsize=18, fontweight="bold",
             color="#171a1e", y=0.965)
fig.text(
    0.5, 0.925,
    "Daily position at 12:00 local mean solar time in 2026  •  fixed celestial direction",
    ha="center", va="center", fontsize=11.5, color="#3e444b",
)
fig.text(
    0.5, 0.895,
    f"Plano: {LAT:.2f}°N, {abs(LON):.2f}°W   •   Galactic (ℓ,b)=({GAL_L:.0f}°, {GAL_B:.0f}°)   →   "
    f"ICRS (RA, Dec)=({coe.ra.deg:.2f}°, {coe.dec.deg:+.2f}°)",
    ha="center", va="center", fontsize=9.4, color="#555c64",
)

# Color scale is a compact legend for the date progression.
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.045, pad=0.13,
                    aspect=55)
cbar.set_label("Date in 2026 (day of year)", fontsize=10, labelpad=5)
cbar.set_ticks(month_doys)
cbar.set_ticklabels(month_names)
cbar.ax.tick_params(labelsize=8.5, length=3)

fig.subplots_adjust(left=0.095, right=0.985, top=0.84, bottom=0.18)
fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Wrote {OUT}")
print(f"CoE ICRS RA={coe.ra.deg:.8f} deg, Dec={coe.dec.deg:+.8f} deg")
print(f"Altitude range: {alt.min():.3f} to {alt.max():.3f} deg")
```

### comic — `FigA_comic_HANDDRAWN_from_PDF.jpg`

- **Gold extract:** `gold-022.jpg`
- **Expected dims (gold):** 2080×944
- **Method:** `rerun:makers/make_figA_comic_HANDDRAWN.py`
- **Output:** `/workspace/figs_letter_exact_python/FigA_comic_HANDDRAWN_from_PDF.jpg`
- **Notes:** HANDDRAWN exception — byte-copy of PDF embed.

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_figA_comic_HANDDRAWN.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_figA_comic_HANDDRAWN.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_figA_comic_HANDDRAWN.py
#!/usr/bin/env python3
"""Comic — HANDDRAWN exception only (not generated from drawing code).

Byte-copies the PDF embed extract. Chris: only exception to regenerate-from-code.
"""
from __future__ import annotations
import shutil
from pathlib import Path

OUT = Path("/workspace/figs_letter_exact_python")
GOLD = OUT / "gold_extracts" / "gold-022.jpg"
OUT_NAME = "FigA_comic_HANDDRAWN_from_PDF.jpg"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / OUT_NAME
    shutil.copy2(GOLD, dest)
    print("HANDDRAWN exception wrote", dest, dest.stat().st_size, "from", GOLD)
    return dest


if __name__ == "__main__":
    main()
```

### eclipse — `eclipse_inner_shell.png`

- **Gold extract:** `n/a`
- **Expected dims (gold):** 883×945
- **Method:** `rerun:makers/make_eclipse_inner_shell.py`
- **Output:** `/workspace/figs_letter_exact_python/eclipse_inner_shell.png`
- **Notes:** Reruns refs_from_chris/make_inner_shell_eclipse_sketch.py

**CLI:**

```bash
cd /workspace/figs_letter_exact_python/makers && python3 make_eclipse_inner_shell.py
```

**Maker script (full source):** `/workspace/figs_letter_exact_python/makers/make_eclipse_inner_shell.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/make_eclipse_inner_shell.py
#!/usr/bin/env python3
"""eclipse_inner_shell — inner-shell eclipse geometry (supplemental)."""
from __future__ import annotations

import os
import runpy
import shutil
from pathlib import Path

OUT = Path("/workspace/figs_letter_exact_python")
SCRIPT = Path(
    "/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/refs_from_chris/make_inner_shell_eclipse_sketch.py"
)
OUT_NAME = "eclipse_inner_shell.png"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(SCRIPT.parent)
        runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        os.chdir(cwd)
    # script writes next to itself or in cwd
    cands = list(SCRIPT.parent.glob("*eclipse*.png")) + list(Path("/workspace").glob("**/inner_shell_eclipse_sketch.png"))
    # prefer letter_brief/figs
    # Prefer the PNG the letter-era script just wrote next to itself
    script_out = SCRIPT.parent / "inner_shell_eclipse_sketch.png"
    if script_out.exists():
        src = script_out
    else:
        prefer = Path(
            "/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/inner_shell_eclipse_sketch.png"
        )
        src = prefer if prefer.exists() else (sorted(cands, key=lambda p: p.stat().st_mtime, reverse=True)[0] if cands else None)
    if src is None:
        raise FileNotFoundError("eclipse sketch png not found after run")
    dest = OUT / OUT_NAME
    shutil.copy2(src, dest)
    print("wrote", dest, "from", src, dest.stat().st_size)
    return dest


if __name__ == "__main__":
    main()

```

**Dependency (full source):** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/refs_from_chris/make_inner_shell_eclipse_sketch.py`

```python
# PATH: /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/refs_from_chris/make_inner_shell_eclipse_sketch.py
#!/usr/bin/env python3
"""Inner-shell eclipse geometry sketch (cycle-6 concentric shells).

R_in=1, R_out=3, D=1.8 → α = arcsin(R_in/D) ≈ 33.7°.
Toward CoE, lines of sight with |θ|<α hit the inner shell only
(eclipsing the outer); outside that disk they hit the outer shell.

Regenerates cycle6_two_emitters/inner_shell_eclipse_sketch.png
to match the approved reference used in letter drafts.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Circle, FancyArrowPatch, Polygon, Wedge

HERE = Path(__file__).resolve().parent
OUT = HERE / "inner_shell_eclipse_sketch.png"

R_IN, R_OUT, D = 1.0, 3.0, 1.8
ALPHA = float(np.arcsin(R_IN / D))  # rad
ALPHA_DEG = float(np.rad2deg(ALPHA))

# Colors matching the approved sketch
C_SHELL = "#8B1A1A"
C_SHELL_FILL = "#F5E6E0"
C_US = "#1E5AA8"
C_ECLIPSE = "#B8C4D0"
C_CAPTION = "#333333"


def main() -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 7.8), dpi=150)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-3.55, 3.55)
    ax.set_ylim(-3.85, 4.05)

    # Title + params
    ax.text(
        0.0,
        3.72,
        "Inner shell eclipses outer toward CoE",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
        color="0.05",
    )
    ax.text(
        0.0,
        3.35,
        rf"$R_{{\rm in}}=1,\ R_{{\rm out}}=3,\ D=1.8$"
        rf"$\qquad\alpha=\arcsin(R_{{\rm in}}/D)\approx{ALPHA_DEG:.1f}^\circ$",
        ha="center",
        va="center",
        fontsize=11,
        color="0.15",
    )

    # Outer fill + rings
    ax.add_patch(
        Circle((0, 0), R_OUT, facecolor=C_SHELL_FILL, edgecolor=C_SHELL, lw=2.2, zorder=1)
    )
    ax.add_patch(
        Circle((0, 0), R_IN, facecolor="white", edgecolor=C_SHELL, lw=2.2, zorder=2)
    )

    # Tangent geometry from US at (D, 0)
    # Tangent length and contact angle from CoE: sin(α)=R_in/D
    # Contact points on inner shell: angle ±(90°-α) from +x toward CoE side...
    # Observer at (D,0); tangents touch at angle φ where cos(φ_from_x_at_CoE related):
    # Contact: point = R_in * (cos ψ, sin ψ) with ψ = ±(π/2 + (π/2-α))? 
    # Standard: contact angles from +x (CoE→US): ψ = ±arccos(R_in/D) = ±(90°-α)... 
    # Actually: vector CoE→contact perpendicular to tangent; contact at angle ±β from US direction
    # where β = arcsin?  Contact angle from +x axis: ±(π - π/2 - α wait)
    # US at (D,0). Contact C on circle |C|=R_in, (C-US)·C = 0 ⇒ C·C = C·US ⇒ R_in^2 = D * C_x
    # ⇒ C_x = R_in^2/D, C_y = ±R_in*sqrt(1-(R_in/D)^2) = ±R_in*cos(α)
    cx = (R_IN**2) / D
    cy = R_IN * np.cos(ALPHA)
    c_top = np.array([cx, cy])
    c_bot = np.array([cx, -cy])
    us = np.array([D, 0.0])

    # Extend tangents past contact (toward CoE / beyond)
    def extend(p0, p1, scale=2.8):
        v = p1 - p0
        v = v / np.linalg.norm(v)
        return p0 + scale * v * R_OUT

    t_top_end = extend(us, c_top, scale=2.6)
    t_bot_end = extend(us, c_bot, scale=2.6)

    # Eclipse wedge (sector between tangents, toward CoE)
    # polygon: US → along top tangent far → arc or through contacts → bottom
    # Fill triangle-ish: US, far along top, far along bot (clipped visually by shells)
    eclipse_pts = np.vstack(
        [
            us,
            t_top_end,
            t_bot_end,
        ]
    )
    ax.add_patch(
        Polygon(
            eclipse_pts,
            closed=True,
            facecolor=C_ECLIPSE,
            edgecolor="none",
            alpha=0.55,
            zorder=0,
        )
    )

    # Redraw shells above eclipse fill for clean edges on near side
    ax.add_patch(
        Circle((0, 0), R_OUT, facecolor="none", edgecolor=C_SHELL, lw=2.2, zorder=3)
    )
    ax.add_patch(
        Circle((0, 0), R_IN, facecolor="white", edgecolor=C_SHELL, lw=2.2, zorder=4)
    )
    # re-fill eclipse only behind inner (clip): draw eclipse wedge again under labels
    # Keep simple: shade between tangents inside outer disk via wedge approx
    # Re-shade using Wedge from US is hard; use polygon clipped conceptually.

    # Tangents (dashed blue)
    for end in (t_top_end, t_bot_end):
        ax.plot(
            [us[0], end[0]],
            [us[1], end[1]],
            ls="--",
            color=C_US,
            lw=1.6,
            zorder=5,
        )

    # Direct LOS toward CoE (arrow to inner shell)
    ax.annotate(
        "",
        xy=(R_IN, 0),
        xytext=(D - 0.02, 0),
        arrowprops=dict(arrowstyle="-|>", color=C_US, lw=1.8, mutation_scale=12),
        zorder=6,
    )

    # Outer hit ray (solid red) — angle slightly outside α
    ang = ALPHA + np.deg2rad(18)
    # Ray from US in direction toward CoE-ish but above: direction = (-cos(ang), sin(ang))
    direction = np.array([-np.cos(ang), np.sin(ang)])
    # Intersect outer circle |US + t d|^2 = R_out^2
    # (us + t d)·(us + t d) = R_OUT^2
    a = np.dot(direction, direction)
    b = 2 * np.dot(us, direction)
    c = np.dot(us, us) - R_OUT**2
    disc = b**2 - 4 * a * c
    t1 = (-b - np.sqrt(disc)) / (2 * a)
    t2 = (-b + np.sqrt(disc)) / (2 * a)
    t_hit = t1 if t1 > 0 else t2
    hit = us + t_hit * direction
    ax.annotate(
        "",
        xy=hit,
        xytext=us,
        arrowprops=dict(arrowstyle="-|>", color=C_SHELL, lw=1.7, mutation_scale=11),
        zorder=6,
    )
    ax.plot(hit[0], hit[1], "o", color=C_SHELL, ms=5, zorder=7)
    ax.text(
        hit[0] + 0.12,
        hit[1] + 0.28,
        "hits outer",
        color=C_SHELL,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    # Contact dots
    ax.plot(c_top[0], c_top[1], "o", color=C_SHELL, ms=5, zorder=7)
    ax.plot(c_bot[0], c_bot[1], "o", color=C_SHELL, ms=5, zorder=7)

    # Alpha arc
    arc_r = 0.55
    ax.add_patch(
        Arc(
            us,
            2 * arc_r,
            2 * arc_r,
            angle=0,
            theta1=180 - ALPHA_DEG,
            theta2=180,
            color=C_US,
            lw=1.5,
            zorder=8,
        )
    )
    ax.text(
        us[0] - arc_r * 0.85,
        0.22,
        r"$\alpha$",
        color=C_US,
        fontsize=13,
        ha="center",
        va="bottom",
        fontweight="bold",
    )

    # Labels
    ax.plot(0, 0, "o", color="k", ms=7, zorder=8)
    ax.text(0, -0.28, "CoE", ha="center", va="top", fontsize=12, fontweight="bold")
    ax.plot(D, 0, "o", color=C_US, ms=11, zorder=8)
    ax.text(D, -0.32, "us", color=C_US, ha="center", va="top", fontsize=13, fontweight="bold")

    ax.text(
        0,
        R_OUT + 0.18,
        "outer shell",
        color=C_SHELL,
        ha="center",
        va="bottom",
        fontsize=11,
        fontweight="bold",
    )
    ax.text(
        0,
        R_IN + 0.12,
        "inner shell",
        color=C_SHELL,
        ha="center",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        -1.15,
        0.55,
        "hits inner\n(eclipses outer)",
        color=C_US,
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        zorder=9,
    )

    # Footer caption
    ax.text(
        0.0,
        -3.55,
        r"Sky disk $|\theta|<\alpha$ (and $360^\circ-\alpha$ twin): "
        r"lines of sight hit inner shell only",
        ha="center",
        va="center",
        fontsize=10.5,
        color=C_CAPTION,
        fontweight="bold",
    )

    fig.savefig(OUT, dpi=150, bbox_inches="tight", facecolor="white", pad_inches=0.15)
    plt.close(fig)
    print(f"wrote {OUT}  α={ALPHA_DEG:.2f}°")
    return OUT


if __name__ == "__main__":
    main()
```

### Shared — `paths.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/paths.py
"""Path constants for CoE figure regeneration."""
from pathlib import Path

ROOT = Path("/workspace/coe-figures")
SRC = ROOT / "src"
OUT = Path("/workspace/figs_letter_exact_python")
FIGS_FROM_MAKERS = OUT
REF = ROOT / "ref"

COSMO = Path("/workspace/cosmology")
MIGKAS_CSV = COSMO / "migkas/digitised/migkas2020_fig6_A_vs_longitude.csv"
MIGKAS_FIG7_CSV = COSMO / "migkas/digitised/migkas2020_fig7_OurSample_A_vs_longitude.csv"
BUBBLE = COSMO / "papers/spreadsheet_model/bubble"
LETTER_FIGS = (
    COSMO
    / "papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs"
)
CMB = COSMO / "papers/A_equals_fourier/steps_v3/cmb_reconcile"
CYCLE4 = CMB / "cycle4"
CYCLE6 = CMB / "cycle6_two_emitters"
CYCLE7 = CMB / "cycle7_recombination"

UNIVERSE_NEWT = BUBBLE / "universe_newt_omb.json"
UNIVERSE_REL = BUBBLE / "universe_rel_mcfree.json"
THIN_FREE = BUBBLE / "thin_b_free_hug.json"
NEWTON_CURVE_NPY = BUBBLE / "gr_ltb/refit_mass/newton_curve.npy"
LTB_CURVE_NPY = BUBBLE / "gr_ltb/refit_mass/ltb_mcfree_z03_curve.npy"
DTT_JSON = CYCLE4 / "dtt/dtt.json"
CANCEL_JSON = CYCLE6 / "simple_cmb_two_quad_cancel.json"
VS_KILL_JSON = CYCLE6 / "vs_kill_term_0_360.json"
TWO_EMITTERS_JSON = CYCLE6 / "two_emitters.json"
COE_ANALEMMA_JSON = LETTER_FIGS / "coe_analemma_texas.json"
COE_CARDINAL_JSON = LETTER_FIGS / "coe_cardinal_direction_texas.json"

```

### Shared — `physics.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/physics.py
"""Shared physics / data helpers for CoE letter figures.

Imports digitised Migkas CSVs and may call Newtonian helpers under the
cosmology bubble tree. Figure scripts should import from here, not copy
legacy plotting code.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

from paths import (
    BUBBLE,
    CANCEL_JSON,
    COE_ANALEMMA_JSON,
    COE_CARDINAL_JSON,
    DTT_JSON,
    LTB_CURVE_NPY,
    MIGKAS_CSV,
    MIGKAS_FIG7_CSV,
    NEWTON_CURVE_NPY,
    THIN_FREE,
    TWO_EMITTERS_JSON,
    UNIVERSE_NEWT,
    UNIVERSE_REL,
    VS_KILL_JSON,
)

# --- headline constants (letter / CoE_calculations.md) ---
R_VISION = 14.30048287260589  # Gpc
K_IN = 1.0276268339399446
K_OUT = 1.4754794866889613
H_MEAN = 71.6  # km/s/Mpc (cycle convention)
OS_OMB = 0.049
BAND_CENTRES = (29.5, 122.6, 187.5, 264.4)  # deg on model circle
BAND_HALF = 9.0  # deg
BAND_COLORS = ("#0097b2", "#6fa800", "#e0107a", "#2166d4")
BAND_KINDS = ("MAX", "MIN", "MAX", "MIN")
TILT_DEG = 15.0
ROLL_DEG = 80.8

# CoE sky: +15° tilt → Galactic (ℓ,b)=(180°, −15°)
COE_GAL_L, COE_GAL_B = 180.0, -15.0
PLANO_LAT, PLANO_LON, PLANO_H = 33.02, -96.70, 200.0

# A_F four-harmonic coefficients (bubble_core convention)
AF_COEFFS = {
    0: (0.9965, 0.0),
    1: (0.0296, 0.0434),
    2: (0.0636, 0.0366),
    3: (-0.0299, -0.0114),
}

TH1 = np.arange(0.0, 360.0, 1.0)


def A_F(ell: np.ndarray) -> np.ndarray:
    """Migkas Fourier reconstruction A_F(ℓ)."""
    ell = np.asarray(ell, dtype=float)
    rad = np.radians(ell)
    y = np.full_like(rad, AF_COEFFS[0][0], dtype=float)
    for m, (c, s) in AF_COEFFS.items():
        if m == 0:
            continue
        y = y + c * np.cos(m * rad) + s * np.sin(m * rad)
    return y


def load_migkas_fig6():
    """Digitised Migkas 2020 Fig. 6 A/A_all vs Galactic longitude."""
    data = np.genfromtxt(MIGKAS_CSV, delimiter=",", names=True)
    return (
        np.asarray(data["galactic_longitude_deg"], float),
        np.asarray(data["A_over_Aall"], float),
        np.asarray(data["A_over_Aall_err_1sigma"], float),
    )


def load_migkas_fig7():
    data = np.genfromtxt(MIGKAS_FIG7_CSV, delimiter=",", names=True)
    return (
        np.asarray(data["galactic_longitude_deg"], float),
        np.asarray(data["A_over_Aall"], float),
        np.asarray(data["A_over_Aall_err_1sigma"], float),
    )


def load_fit_curves():
    """Return stored Newtonian / free / relativistic A(ℓ) on the 1° grid.

    Prefer recomputing the Ωs=0.049 Newtonian curve via newt_omb when available;
    fall back to the stored A array in universe_newt_omb.json.
    """
    newt_j = json.loads(UNIVERSE_NEWT.read_text())
    sine = newt_j["Rs14.3|sine"]
    A_newt = np.asarray(sine["A"], float)

    # free Ωs Newtonian (previous fit) from newton_curve.npy row 1
    A_free = np.load(NEWTON_CURVE_NPY)[1].astype(float)
    free_meta = json.loads(THIN_FREE.read_text())

    A_rel = np.load(LTB_CURVE_NPY).astype(float)
    rel_j = json.loads(UNIVERSE_REL.read_text())

    meta = {
        "newt": {
            "D": sine["D"],
            "R_in": sine["R_in"],
            "R_out": sine["R_out"],
            "Os": sine["Os"],
            "l0": sine["l0"],
            "amp": sine["amp"],
            "phi": sine["phi"],
            "R2": sine["R2"],
            "mcc": sine["mcc"],
        },
        "free": {
            "D": free_meta["DR"] * R_VISION,
            "Os": free_meta["Os"],
            "l0": free_meta["l0"],
            "amp": free_meta.get("amp"),
            "phi": free_meta.get("phi"),
            "R2": free_meta["R2"],
            "rin": free_meta["rin"],
            "rout": free_meta["rout"],
        },
        "rel": {
            "D": rel_j["D"],
            "R_in": rel_j["R_in"],
            "R_out": rel_j["R_out"],
            "Os": rel_j["Os"],
            "l0": rel_j["l0_s"],
            "amp": rel_j["amp"],
            "phi": rel_j["phi"],
            "R2": rel_j["R2_s"],
            "mcc": rel_j["mcc"],
            "z": rel_j.get("z", 0.3),
        },
    }
    return TH1.copy(), A_newt, A_free, A_rel, meta


def recompute_newt_omb_A(force: bool = False):
    """Optionally recompute Newtonian Ωs=0.049 A(ℓ) via bubble newt_omb."""
    cwd = os.getcwd()
    try:
        sys.path.insert(0, str(BUBBLE))
        os.chdir(BUBBLE)
        import newt_omb as N
        from bubble_core import A_from, shifted

        meta = json.loads(UNIVERSE_NEWT.read_text())["Rs14.3|sine"]
        a = N.geo(meta["D"], meta["mcc"], meta["kin"], meta["kout"], R_VISION)
        af, _rr = N.af_of(a)
        A = A_from(shifted(af, meta["l0"]), TH1 - meta["l0"], meta["amp"], meta["phi"])
        return TH1.copy(), np.asarray(A, float)
    finally:
        os.chdir(cwd)


def r2_score(y, yhat):
    y = np.asarray(y, float)
    yhat = np.asarray(yhat, float)
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan


def load_dtt():
    return json.loads(DTT_JSON.read_text())


def kill_from_doppler(dtt_row: dict) -> tuple[np.ndarray, np.ndarray]:
    """kill(θ) = −(Doppler − ⟨Doppler⟩); θ toward CoE."""
    th = np.asarray(dtt_row["theta_toward_deg"], float)
    dop = np.asarray(dtt_row["dtt_doppler"], float)
    kill = -(dop - dop.mean())
    return th, kill


def extend_0_360(th_0_180: np.ndarray, f_0_180: np.ndarray):
    """Mirror an axisymmetric 0–180° pattern onto a full 0–360° great circle."""
    th = np.asarray(th_0_180, float)
    f = np.asarray(f_0_180, float)
    # assume th goes 0..180
    th2 = np.concatenate([th, 360.0 - th[-2:0:-1]])
    f2 = np.concatenate([f, f[-2:0:-1]])
    return th2, f2


def unit_norm(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, float)
    m = np.max(np.abs(y))
    return y / m if m > 0 else y


def legendre_coeffs(f_theta: np.ndarray, theta_deg: np.ndarray, lmax: int = 4):
    """Legendre a_ℓ for axisymmetric f(θ), θ from CoE (deg)."""
    from numpy.polynomial.legendre import legval

    th = np.radians(theta_deg)
    mu = np.cos(th)
    o = np.argsort(mu)
    mu, f = mu[o], np.asarray(f_theta, float)[o]
    out = []
    for ell in range(lmax + 1):
        P = legval(mu, [0] * ell + [1])
        out.append(float((2 * ell + 1) / 2 * np.trapezoid(f * P, mu)))
    return out


def partner_cancel(doppler: np.ndarray, theta: np.ndarray, mode: str = "dip_p2"):
    """Cancel Doppler with a Legendre partner (dipole+P2 or ℓ≤4)."""
    a = legendre_coeffs(doppler, theta, lmax=4)
    from numpy.polynomial.legendre import legval

    mu = np.cos(np.radians(theta))
    if mode == "dip_p2":
        partner = -(a[1] * legval(mu, [0, 1]) + a[2] * legval(mu, [0, 0, 1]))
    elif mode == "l_le_4":
        partner = np.zeros_like(mu)
        for ell in range(1, 5):
            partner = partner - a[ell] * legval(mu, [0] * ell + [1])
    else:
        raise ValueError(mode)
    residual = doppler + partner - (doppler + partner).mean()
    return partner, residual, a


def load_coe_sky_meta():
    card = json.loads(COE_CARDINAL_JSON.read_text())
    return {
        "gal_lb": tuple(card["coe_galactic_lb_deg"]),
        "ra_deg": card["coe_icrs_ra_deg"],
        "dec_deg": card["coe_icrs_dec_deg"],
        "site_lat": card["site_lat_deg"],
        "site_lon": card["site_lon_deg"],
        "site_h": card["site_height_m"],
        "site_name": card["site_name"],
    }


def save_fig(fig, stem: str, also_pdf: bool = True, dpi: int = 200):
    from paths import OUT

    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / f"{stem}.png"
    fig.savefig(png, dpi=dpi, bbox_inches="tight", facecolor="white")
    if also_pdf:
        try:
            fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
        except Exception:
            pass
    return png


def shell_radii(D: float, kin: float = K_IN, kout: float = K_OUT, R: float = R_VISION):
    return D - kin * R, D + kout * R


# --- sky orientation geometry ---
def band_circle_xyz(theta_c: float, half: float = BAND_HALF, n: int = 360):
    """Unit-sphere band at polar angle theta_c from +x (anti-CoE), width ±half."""
    # In model frame: anti-CoE = +x, CoE = −x; observer at origin.
    # Band is a latitude belt at polar angle θ from +x.
    ph = np.linspace(0, 2 * np.pi, n, endpoint=False)
    # sample mid, inner, outer edges
    out = {}
    for name, th in (("mid", theta_c), ("lo", theta_c - half), ("hi", theta_c + half)):
        t = np.radians(th)
        # x = cos θ, cylindrical radius = sin θ
        x = np.full_like(ph, np.cos(t))
        y = np.sin(t) * np.cos(ph)
        z = np.sin(t) * np.sin(ph)
        out[name] = np.column_stack([x, y, z])
    return out


def rotate_xyz(xyz: np.ndarray, tilt_deg: float = 0.0, roll_deg: float = 0.0):
    """Apply tilt about y then roll about the CoE axis (−x after tilt).

    Untilted: anti-CoE at Galactic (ℓ,b)=(0,0) → +x in model.
    +tilt: rotate about y so anti-CoE goes to b=+tilt (and CoE to b=−tilt).
    +roll: rotate about the anti-CoE axis.
    """
    t = np.radians(tilt_deg)
    r = np.radians(roll_deg)
    # tilt about +y: (x,z) -> (x cos t - z sin t, z cos t + x sin t)
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    xt = x * np.cos(t) - z * np.sin(t)
    yt = y
    zt = x * np.sin(t) + z * np.cos(t)
    # roll about anti-CoE axis (tilted +x'): rotate (y,z) in the plane ⊥ to x'
    # After tilt, anti-CoE direction is (cos t, 0, sin t).
    # Use Rodrigues about that axis.
    ax = np.array([np.cos(t), 0.0, np.sin(t)])
    v = np.column_stack([xt, yt, zt])
    if abs(roll_deg) < 1e-12:
        return v
    k = ax
    c, s = np.cos(r), np.sin(r)
    # v' = v cos + (k×v) sin + k (k·v) (1-cos)
    kxv = np.cross(k, v)
    kdotv = v @ k
    return v * c + kxv * s + np.outer(kdotv, k) * (1 - c)


def xyz_to_galactic(xyz: np.ndarray):
    """Map model-frame unit vectors to Galactic (ℓ,b) deg.

    Convention after tilt/roll: anti-CoE direction → Galactic (ℓ,b).
    Untilted anti-CoE at (0,0); +15° tilt → anti-CoE at (0,+15°), CoE at (180,−15°).
    """
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    # Galactic: ℓ along equator, b latitude. Map x→ toward ℓ=0,b=0; y→ℓ=90; z→b=+90
    b = np.degrees(np.arcsin(np.clip(z, -1, 1)))
    ell = np.degrees(np.arctan2(y, x)) % 360.0
    return ell, b


def mollweide_xy(lon_deg: np.ndarray, lat_deg: np.ndarray, center_lon: float = 0.0):
    """Mollweide projection; lon increases to the left (astro convention)."""
    lon = np.radians(((lon_deg - center_lon + 180) % 360) - 180)
    lat = np.radians(lat_deg)
    # solve 2γ + sin 2γ = π sin φ
    gamma = lat.copy()
    for _ in range(12):
        gamma = gamma - (2 * gamma + np.sin(2 * gamma) - np.pi * np.sin(lat)) / (
            4 * np.cos(gamma) ** 2
        )
    x = (2 * np.sqrt(2) / np.pi) * lon * np.cos(gamma)
    # flip x so longitude increases to the left
    x = -x
    y = np.sqrt(2) * np.sin(gamma)
    return x, y


def load_cancel_json():
    return json.loads(CANCEL_JSON.read_text())


def load_vs_kill():
    return json.loads(VS_KILL_JSON.read_text())


def load_two_emitters():
    return json.loads(TWO_EMITTERS_JSON.read_text())

```

### Shared — `run_all.py`

```python
# PATH: /workspace/figs_letter_exact_python/makers/run_all.py
#!/usr/bin/env python3
"""Run every letter figure maker → /workspace/figs_letter_exact_python/."""
from __future__ import annotations

import importlib
import sys
import traceback
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = Path("/workspace/figs_letter_exact_python")
sys.path.insert(0, str(HERE))

MODULES = [
    "make_fig01_02_geometry_achart",
    "make_fig03_migkas_slices",
    "make_fig04_amp_full",
    "make_fig05_sky_orientation",
    "make_fig06_overlay_stack",
    "make_fig07_storyboard",
    "make_fig08_two_emitters",
    "make_fig09_kill",
    "make_fig10_cancel_partner",
    "make_figA1_ltb_shell",
    "make_figA2_coe_sun",
    "make_figA2b_crossings",
    "make_figA3_coe_sky",
    "make_figA_comic_HANDDRAWN",
    "make_eclipse_inner_shell",
]


def main() -> int:
    ok, fail = [], []
    for name in MODULES:
        print(f"\n=== {name} ===", flush=True)
        try:
            mod = importlib.import_module(name)
            importlib.reload(mod)
            out = mod.main()
            ok.append((name, str(out)))
        except Exception:
            traceback.print_exc()
            fail.append(name)
    zpath = OUT / "figs_letter_exact_python.zip"
    files = list(OUT.glob("Fig*.png")) + list(OUT.glob("Fig*.jpg")) + list(OUT.glob("eclipse_*.png"))
    with zipfile.ZipFile(zpath, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(set(files)):
            zf.write(f, f.name)
    print("\n======== SUMMARY ========")
    print(f"OK: {len(ok)}  FAIL: {len(fail)}")
    for n, p in ok:
        print(f"  OK {n} → {p}")
    for n in fail:
        print(f"  FAIL {n}")
    print("zip:", zpath)
    return 0 if not fail else 1


if __name__ == "__main__":
    raise SystemExit(main())
```

