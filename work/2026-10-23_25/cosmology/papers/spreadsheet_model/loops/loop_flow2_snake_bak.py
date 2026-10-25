import json, sys, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.patches import Ellipse, FancyArrowPatch, Patch
from matplotlib.lines import Line2D
from PIL import Image
sys.path.insert(0, "/workspace/cosmology/figures"); sys.path.insert(0, "../bubble")
from coe_half_loop_geometry import lb_to_n, n_to_lb, ring_basis, MIGKAS_MAXES, MIGKAS_MINS, lon_disp_from_ell, great_circle_through_poles
_here = os.getcwd(); os.chdir("../bubble")
import analysis_curve as AC
from bubble_core import TH1, ML, MA, MLO, MHI, AFfour
os.chdir(_here)
F = json.load(open("loop_fit.json")); d = json.load(open("../bubble/thin_b_free_hug.json"))
A = AC.curve("ptshell", d); l0 = F["l0"]
FIT = F["fits"]["b0_only"]; UNT = F["fits"]["untilted"]
COLS = {"MAX1": "#0097b2", "MAX2": "#e0107a", "MIN1": "#6fa800", "MIN2": "#2166d4"}
pairs = F["pairs"]
keys = []
cm = {"MAX": 0, "MIN": 0}
for p in pairs: cm[p["kind"]] += 1; keys.append("%s%d" % (p["kind"], cm[p["kind"]]))
BH = 9.0
def moll(l, b):
    lon = np.radians(lon_disp_from_ell(np.asarray(l, float))); lat = np.radians(np.asarray(b, float))
    th = lat.copy()
    for _ in range(40):
        f = 2 * th + np.sin(2 * th) - np.pi * np.sin(lat); th = th - f / (2 + 2 * np.cos(2 * th) + 1e-12)
    return 2 * np.sqrt(2) / np.pi * lon * np.cos(th), np.sqrt(2) * np.sin(th)
def pts_lb(P):
    b = np.degrees(np.arcsin(np.clip(P[:, 2], -1, 1))); l = np.degrees(np.arctan2(P[:, 1], P[:, 0])) % 360; return l, b
def ring(D, Y, Z, t, phi):
    t = np.radians(t); return np.cos(t) * D[None] + np.sin(t) * (np.outer(np.sin(phi), Y) + np.outer(np.cos(phi), Z))
def sky_axes(ax):
    tt = np.linspace(0, 2 * np.pi, 400); ax.plot(2 * np.sqrt(2) * np.cos(tt), np.sqrt(2) * np.sin(tt), color="k", lw=0.9)
    for bb in (-60, -30, 0, 30, 60):
        x, y = moll(np.linspace(0.01, 359.99, 300), np.full(300, bb)); o = np.argsort(x); ax.plot(x[o], y[o], color="0.85", lw=0.5, zorder=0)
    for ll in range(0, 360, 30):
        x, y = moll(np.full(200, ll if ll else 0.001), np.linspace(-89.9, 89.9, 200)); ax.plot(x, y, color="0.85", lw=0.5, zorder=0)
    for ll in (90, 0, 270):
        x, y = moll(ll, 0); ax.text(x, y - 0.1, "%d°" % ll, fontsize=7, ha="center", va="top", color="0.35")
    for bb in (-60, -30, 30, 60):
        x, y = moll(180.01, bb); ax.text(x - 0.08, y, "%+d°" % bb, fontsize=6.5, ha="right", va="center", color="0.35")
    ax.set_xlim(-2.95, 2.95); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal"); ax.axis("off")
def draw_loops(ax, lat_b, lon_a, roll, gradient, to_xy=None, lw=1.6, under=None):
    D = lb_to_n(lon_a, lat_b); Y, Z = ring_basis(D, roll)
    to_xy = to_xy or (lambda l, b: moll(l, b))
    phi0 = np.pi / 2; phi = phi0 + np.linspace(-np.pi, np.pi, 721)
    alpha = 1 - 0.95 * np.abs(np.angle(np.exp(1j * (phi - phi0)))) / np.pi if gradient else np.ones_like(phi)
    for k, p in zip(keys, pairs):
        c = COLS[k]; t = p["theta"]
        lo, hi = ring(D, Y, Z, t - BH, phi), ring(D, Y, Z, t + BH, phi); mid = ring(D, Y, Z, t, phi)
        x1, y1 = to_xy(*pts_lb(lo)); x2, y2 = to_xy(*pts_lb(hi)); xm, ym = to_xy(*pts_lb(mid))
        quads, qa, segs, sa = [], [], [], []
        span = 0.8 if to_xy is None or True else 0.8
        scale = max(np.ptp(xm), 1e-9)
        for i in range(len(phi) - 1):
            q = np.array([[x1[i], y1[i]], [x1[i + 1], y1[i + 1]], [x2[i + 1], y2[i + 1]], [x2[i], y2[i]]])
            if np.max(np.ptp(q, axis=0)) > JUMP: continue
            quads.append(q); qa.append(0.5 * (alpha[i] + alpha[i + 1]))
            s = np.array([[xm[i], ym[i]], [xm[i + 1], ym[i + 1]]])
            if np.max(np.ptp(s, axis=0)) < JUMP: segs.append(s); sa.append(0.5 * (alpha[i] + alpha[i + 1]))
        from matplotlib.colors import to_rgba
        r, g, b_, _ = to_rgba(c)
        ax.add_collection(PolyCollection(quads, facecolors=[(r, g, b_, 0.34 * a) for a in qa], edgecolors="none", zorder=2))
        if under: ax.add_collection(LineCollection(segs, colors=[(1, 1, 1, 0.9 * a) for a in sa], linewidths=lw + 1.8, zorder=3, capstyle="butt"))
        ax.add_collection(LineCollection(segs, colors=[(r, g, b_, a) for a in sa], linewidths=lw, zorder=4, capstyle="butt"))
    return D, Y
def markers(ax, D, to_xy=None, ms=7.5, anti=True):
    to_xy = to_xy or (lambda l, b: moll(l, b))
    for k, p in zip(keys, pairs):
        x, y = to_xy(*p["migkas"]); ax.plot(x, y, "^" if p["kind"] == "MAX" else "s", ms=ms, color=COLS[k], mec="k", mew=0.7, zorder=8)
    if anti:
        la, ba = n_to_lb(D); x, y = to_xy(la, ba); ax.plot(x, y, "o", ms=6, color="#8e1b1b", mec="k", zorder=8)
        lc, bc = n_to_lb(-D); x, y = to_xy(lc, bc); ax.plot(x, y, "*", ms=9, color="#6b3fa0", mec="k", zorder=8)
def axis_line(ax, D, Y, to_xy=None, color="#5b2c6f"):
    to_xy = to_xy or (lambda l, b: moll(l, b))
    t = np.linspace(0, 2 * np.pi, 721); P = np.outer(np.cos(t), D) + np.outer(np.sin(t), Y)
    x, y = to_xy(*pts_lb(P)); segs = [np.array([[x[i], y[i]], [x[i + 1], y[i + 1]]]) for i in range(720) if abs(x[i + 1] - x[i]) < JUMP and abs(y[i + 1] - y[i]) < JUMP]
    ax.add_collection(LineCollection(segs, colors=color, linewidths=1.1, linestyles="--", zorder=5))
JUMP = 0.8
FW, FH = 19.0, 20.0
fig = plt.figure(figsize=(FW, FH))
def R(x_in, top_in, w_in, h_in): return [x_in / FW, 1 - (top_in + h_in) / FH, w_in / FW, h_in / FH]
BLAST = "/home/box/agent-data/agents/196ab413-0708-4f81-99fc-58a32c5366c5/attachments/6b94c3bbca538c4948a97d607ab242f375b79826f33c12618e471cd7b0fe047a.png"
bimg = np.asarray(Image.open(BLAST).convert("RGB")); bw_in = 12.2; bh_in = bw_in * bimg.shape[0] / bimg.shape[1]
MW, MH = 8.6, 4.3
rects = {0: R(0.15, 0.75, bw_in, bh_in), 1: R(13.45, 1.25, 5.3, 4.6),
         3: R(0.6, 8.4, MW, MH), 2: R(10.15, 8.4, MW, MH), 4: R(0.6, 13.95, MW, MH), 5: R(10.15, 13.95, MW, MH)}
geo = {k: (r[0], r[1], r[2]) for k, r in rects.items()}
axs = {k: fig.add_axes(r, projection=("mollweide" if k in (2, 3, 4) else None)) for k, r in rects.items()}
# 0. Chris's homologous-blast diagram, embedded as sent (both vector drawings + full proof box)
a = axs[0]; a.imshow(bimg, interpolation="lanczos"); a.axis("off")
a.set_title("0. Homologous blast $\\vec V=H\\vec R$ from the CoE (your diagram)", fontsize=12, pad=4, loc="left")
fig.text(R(0.3, 0, 1, 0)[0], 1 - (0.75 + bh_in + 0.12) / FH,
         "Link to step 1: pure homologous expansion $V=HR$ gives an isotropic Hubble flow, so the CoE drops out. The observed anisotropy comes from gravity decelerating\n"
         "different positions by different amounts (a departure from homology), which is what your model A captures.",
         fontsize=10.5, color="#8e1b1b", style="italic", va="top", ha="left")
# 1. A(l)
a = axs[1]; th2 = np.append(TH1, 360); A2 = np.append(A, A[0])
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=2.8, color="0.5", ecolor="0.82", lw=0.7, label="Migkas Fig. 6")
a.plot(th2, AFfour(th2), color="0.25", lw=1.0, ls=":", label="$A_F$")
a.plot(th2, A2, color="k", lw=1.8, label="your model A, thin shell ($R^2$=%.3f)" % F["R2"])
for k, p in zip(keys, pairs):
    L = p["l_model"]; a.axvspan(L - BH, L + BH, color=COLS[k], alpha=0.22, lw=0)
    a.axvline(L, color=COLS[k], lw=1.2)
    a.text(L, 1.215, "%s\n$\\ell$=%.1f°" % (p["kind"], L), color=COLS[k], ha="center", va="bottom", fontsize=8, fontweight="bold")
a.set_xlim(0, 360); a.set_xticks(range(0, 361, 45)); a.set_ylim(0.76, 1.21)
a.set_xlabel("Galactic longitude $\\ell$ (deg)", fontsize=9); a.set_ylabel("$A/A_{\\rm all}$ (higher = faster)", fontsize=9)
a.legend(fontsize=7.5, loc="lower left", frameon=False)
a.set_title("1. $A(\\ell)$, latest fit: $\\ell_0$=%.2f°, Amp=%.3f, $\\varphi$=%.3f rad" % (l0, F["amp"], F["phi"]), fontsize=11, pad=34)
# ======== steps 2-5 in Chris's original format (make_coe_rings_gradient_loops.py panel B, make_migkas_fig8_gradient_loops_overlay.py) ========
import coe_gradient_loop_geometry as CG
import make_coe_rings_gradient_loops as RB
import make_migkas_fig8_gradient_loops_overlay as OV
from coe_half_loop_geometry import celestial_equator_lb, great_circle_through_poles, mu_pair_for_band
from coe_gradient_loop_geometry import small_circle_full_gradient, segment_scale, quantized_runs, rgba_array
MIN_OFF = 2.0   # display-only offset (deg) applied to BOTH MIN loops so the near-coincident pair is seen as two lines
# loop list in Chris's order MAX-1, MAX-2, MIN-1, MIN-2 (pairs are already in this order)
LOOPS = []
for k, p in zip(keys, pairs):
    t = p["theta"]; td = t + (MIN_OFF if p["kind"] == "MIN" else 0.0)
    lo, hi = td - BH, td + BH
    if hi > 180 and lo < 180 and abs(td - 180) < BH: hi = 179.999   # small ring round the CoE: band is a disc of radius |180-t|+9
    LOOPS.append(dict(key=k, kind=p["kind"], theta=t, tdisp=td, lo=lo, hi=hi, l=p["l_model"], mig=p["migkas"]))
CB = dict(COLS)   # loop colours = step-1 A(l) band colours (MAX1 cyan, MAX2 magenta, MIN1 green, MIN2 blue)
MKF = "#2b2b2b"   # Migkas marker fill (neutral); edge = paired loop colour
CO = dict(COLS)
LWB = {"MAX1": 2.4, "MAX2": 3.2, "MIN1": 1.9, "MIN2": 1.9}
LWO = {"MAX1": 3.0, "MAX2": 3.4, "MIN1": 2.6, "MIN2": 2.6}
OV.FILL_ALPHA = 0.42

def ribbon_B(ax, D, L, roll, lw, underlay=False, n_pts=721):
    """= RB.fill_mollweide_gradient_ribbon with a per-loop line width (+ optional white underlay)."""
    mu_a, mu_b = mu_pair_for_band(L["lo"], L["hi"]); mu_mid = float(np.cos(np.radians(L["tdisp"]))); color = CB[L["key"]]
    lon_a, lat_a, sc = RB._full_ring_lonlat_rad(D, mu_a, L["tdisp"], roll, n_pts=n_pts)
    lon_b, lat_b, _ = RB._full_ring_lonlat_rad(D, mu_b, L["tdisp"], roll, n_pts=n_pts)
    gap = 0.9 * np.pi; ql = np.stack([lon_a[:-1], lon_a[1:], lon_b[1:], lon_b[:-1]])
    runs = quantized_runs(segment_scale(sc), (ql.max(0) - ql.min(0)) <= gap)
    if runs: ax.add_collection(PolyCollection(RB._strip_polys_2d(lon_a, lat_a, lon_b, lat_b, runs), facecolors=rgba_array(color, [RB.FILL_ALPHA_B * r[2] for r in runs]), edgecolors="none", zorder=2, linewidths=0))
    lon_m, lat_m, sc_m = RB._full_ring_lonlat_rad(D, mu_mid, L["tdisp"], roll, n_pts=1441)
    runs = quantized_runs(segment_scale(sc_m), np.abs(np.diff(lon_m)) <= gap)
    if runs:
        pl = RB._polylines_2d(lon_m, lat_m, runs)
        if underlay: ax.add_collection(LineCollection(pl, colors=rgba_array("white", [0.95 * r[2] for r in runs]), linewidths=lw + 1.6, zorder=3.9, capstyle="butt", joinstyle="round"))
        ax.add_collection(LineCollection(pl, colors=rgba_array(color, [RB.LINE_ALPHA * r[2] for r in runs]), linewidths=lw, zorder=4, capstyle="butt", joinstyle="round"))

def panel_B(ax, b0, roll, gradient, title):
    CG.ALPHA_FLOOR = 0.12 if gradient else 1.0
    D = lb_to_n(l0, b0)
    for b in [-60, -30, 0, 30, 60]:
        ell = np.linspace(0, 360, 361); RB.plot_lb_curve_mollweide(ax, ell, np.full_like(ell, b), color="#bdc3c7", lw=0.5, zorder=0)
    for e0 in range(0, 360, 30):
        bb = np.linspace(-90, 90, 181); RB.plot_lb_curve_mollweide(ax, np.full_like(bb, e0), bb, color="#bdc3c7", lw=0.5, zorder=0)
    ls, bs = great_circle_through_poles(D, n=721); RB.plot_lb_curve_mollweide(ax, ls, bs, color=RB.C_AXIS, lw=2.0, zorder=2)
    for L in LOOPS: ribbon_B(ax, D, L, roll, LWB[L["key"]], underlay=(L["key"] == "MAX2"))
    le, be = celestial_equator_lb(721); RB.plot_lb_curve_mollweide(ax, le, be, color=RB.C_EQ, lw=1.8, ls="--", zorder=2)
    la, ba = n_to_lb(D); lc, bc = n_to_lb(-D)
    ax.plot(np.radians(lon_disp_from_ell(la)), np.radians(ba), "o", ms=8, color=RB.C_ANTI, mec="k", mew=0.6, zorder=6)
    ax.plot(np.radians(lon_disp_from_ell(lc)), np.radians(bc), "*", ms=14, color=RB.C_COE, mec="k", mew=0.6, zorder=6)
    for L in LOOPS:
        ax.plot(np.radians(lon_disp_from_ell(L["mig"][0])), np.radians(L["mig"][1]), "^" if L["kind"] == "MAX" else "s", ms=10, color=MKF, mec=CB[L["key"]], mew=2.2, zorder=7)
    te = np.array([135, 90, 45, 0, 315, 270, 225]); ax.set_xticks(np.radians(lon_disp_from_ell(te))); ax.set_xticklabels([r"$%d^\circ$" % e for e in te], fontsize=7)
    ax.set_yticks(np.radians([-60, -30, 0, 30, 60])); ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$+30^\circ$", r"$+60^\circ$"], fontsize=7.5)
    ax.grid(False); ax.set_xlabel(r"Galactic $\ell$ (increases LEFT)", fontsize=8); ax.set_ylabel(r"Galactic $b$", fontsize=8)
    ax.set_title(title, fontsize=10, pad=6)
    CG.ALPHA_FLOOR = 0.12
    return D

def coe_inset(rect, D, roll, gradient, cols, lws, bg="white", ring_under=True):
    """Zoom on the small MAX-2 ring round the CoE (tangent-plane view, ell increases LEFT)."""
    CG.ALPHA_FLOOR = 0.12 if gradient else 1.0
    ai = fig.add_axes(rect); ai.set_facecolor(bg); c = -np.asarray(D, float); lc, bc = n_to_lb(c)
    lr, br = np.radians(lc), np.radians(bc); el = np.array([-np.sin(lr), np.cos(lr), 0]); eb = np.array([-np.sin(br) * np.cos(lr), -np.sin(br) * np.sin(lr), np.cos(br)])
    proj = lambda P: (-np.degrees(P @ el / (P @ c)), np.degrees(P @ eb / (P @ c)))
    L = [q for q in LOOPS if q["key"] == "MAX2"][0]; col = cols["MAX2"]
    mu_a, mu_b = mu_pair_for_band(L["lo"], L["hi"])
    pa, _, sc = small_circle_full_gradient(D, mu_a, L["tdisp"], roll_deg=roll, n=361); pb, _, _ = small_circle_full_gradient(D, mu_b, L["tdisp"], roll_deg=roll, n=361)
    xa, ya = proj(pa); xb, yb = proj(pb); runs = quantized_runs(segment_scale(sc))
    ai.add_collection(PolyCollection(RB._strip_polys_2d(xa, ya, xb, yb, runs), facecolors=rgba_array(col, [0.30 * r[2] for r in runs]), edgecolors="none", zorder=2))
    pm, _, scm = small_circle_full_gradient(D, float(np.cos(np.radians(L["tdisp"]))), L["tdisp"], roll_deg=roll, n=721)
    xm, ym = proj(pm); runs = quantized_runs(segment_scale(scm)); pl = RB._polylines_2d(xm, ym, runs)
    if ring_under: ai.add_collection(LineCollection(pl, colors=rgba_array("white", [0.95 * r[2] for r in runs]), linewidths=lws + 1.5, zorder=3, capstyle="butt"))
    ai.add_collection(LineCollection(pl, colors=rgba_array(col, [r[2] for r in runs]), linewidths=lws, zorder=4, capstyle="butt"))
    Y, _ = ring_basis(D, roll); tp = np.cos(np.radians(L["theta"])) * D + np.sin(np.radians(L["theta"])) * Y   # theta-ray point
    t = np.linspace(-0.5, 0.5, 50); gc = np.outer(np.cos(t), c) + np.outer(np.sin(t), np.cross(np.cross(c, [0, 0, 1]), c) / np.linalg.norm(np.cross(c, [0, 0, 1])))
    gx, gy = proj(gc); ai.plot(gx, gy, color=RB.C_AXIS, lw=1.2, zorder=1)
    ai.plot(0, 0, "*", ms=10, color=RB.C_COE, mec="k", mew=0.5, zorder=6)
    mx, my = proj(lb_to_n(*L["mig"])[None])[0][0], proj(lb_to_n(*L["mig"])[None])[1][0]
    ai.plot(mx, my, "^", ms=7, color=MKF, mec=col, mew=1.8, zorder=7)
    for s in (-1, 1):
        pass
    R = 17; ai.set_xlim(R, -R); ai.set_ylim(-R, R); ai.set_aspect("equal"); ai.set_xticks([]); ai.set_yticks([])
    for sp in ai.spines.values(): sp.set_edgecolor("0.3"); sp.set_linewidth(0.8)
    ai.text(0.5, 0.97, "zoom: CoE ±%d°" % R, transform=ai.transAxes, fontsize=6.3, ha="center", va="top", color=("0.2" if bg == "white" else "white"))
    ai.text(0.5, 0.03, "MAX-2 ring r=%.1f°" % (180 - L["theta"]), transform=ai.transAxes, fontsize=6.3, ha="center", va="bottom", color=("0.2" if bg == "white" else "white"))
    CG.ALPHA_FLOOR = 0.12
    return ai

tt = lambda F_: ", ".join("%+.1f" % v for v in F_["ring"])
# 2. untilted
D = panel_B(axs[2], 0.0, UNT["roll"], False, "2. Four full loops, untilted (anti-CoE at $(%.1f^\\circ,0^\\circ)$)\nloop radius $\\theta=\\ell_{\\rm ext}-\\ell_0$; ring rms miss %.1f°" % (l0, UNT["rms_ring"]))
# 3. tilt fit
D = panel_B(axs[3], FIT["b"], FIT["roll"], False, "3. Best-fit tilt $b_0$=%+.1f°: anti-CoE $(%.1f^\\circ,%+.1f^\\circ)$\nring misses %s°; all inside the ±9° bands" % (FIT["b"], l0, FIT["b"], tt(FIT)))
# 4. gradient
D = panel_B(axs[4], FIT["b"], FIT["roll"], True, "4. Full loops, opacity gradient darkest at $\\theta$ ray, 12%% opposite\n(roll %.1f° about CoE; MIN-1 dark on one side, MIN-2 on the other)" % FIT["roll"])
# 5. overlay (Chris's overlay functions, colours, white/purple CoE-axis stroke; Migkas' own black celestial equator)
rgb = np.asarray(Image.open("/workspace/cosmology/map_compare/FULL_migkas_theta75.png").convert("RGB"))
seed = OV.detect_oval_bbox(rgb); reg = dict(OV.GOOD_REG)
cx, cy, Rx, Ry = reg["cx"], reg["cy"], reg["Rx"], reg["Ry"]
a = axs[5]; a.set_facecolor("white"); a.imshow(rgb, origin="upper", extent=(0, rgb.shape[1], rgb.shape[0], 0), interpolation="nearest")
x0, x1 = cx - Rx - 150, cx + Rx + 520; y0, y1 = cy - Ry - 110, cy + Ry + 130
a.set_xlim(x0, x1); a.set_ylim(y1, y0); a.set_aspect("equal"); a.axis("off")
clip_oval = Ellipse((cx, cy), width=2 * Rx, height=2 * Ry, facecolor="none", edgecolor="none", transform=a.transData); a.add_patch(clip_oval)
color_oval = (float(seed["cx"]), float(seed["cy"]), float(seed["Rx"]), float(seed["Ry"]), 0.978)
D = lb_to_n(l0, FIT["b"]); CG.ALPHA_FLOOR = 0.12
gl, gb = great_circle_through_poles(D, n=721)
OV.plot_lb_curve(a, gl, gb, cx, cy, Rx, Ry, clip_path=clip_oval, color_oval=color_oval, color="white", lw=1.15 + 1.6, alpha=0.95, zorder=3)
OV.plot_lb_curve(a, gl, gb, cx, cy, Rx, Ry, clip_path=clip_oval, color_oval=color_oval, color=OV.C_AXIS, lw=1.15, zorder=4)
_orig_ring = OV.plot_full_ring_gradient
for L in LOOPS:
    OV.plot_full_ring_gradient = (lambda *aa, lw=2.8, _w=LWO[L["key"]], **kk: _orig_ring(*aa, lw=_w, **kk))
    OV.fill_full_gradient_band(a, D, L["lo"], L["hi"], L["tdisp"], FIT["roll"], cx, cy, Rx, Ry, CO[L["key"]], inset=0.97, clip_path=clip_oval, color_oval=color_oval)
OV.plot_full_ring_gradient = _orig_ring
def mk(l, b, m, c, ms, edge="k", mew=0.55):
    px, py = OV.lb_to_px(l, b, cx, cy, Rx, Ry); px, py = cx + (px - cx) * 0.985, cy + (py - cy) * 0.985
    a.plot(px, py, m, ms=ms + 1.8, color="white", mec="white", zorder=6); a.plot(px, py, m, ms=ms, color=c, mec=edge, mew=mew, zorder=7)
la, ba = n_to_lb(D); lc, bc = n_to_lb(-D)
mk(la, ba, "o", OV.C_ANTI, 5.5); mk(lc, bc, "*", OV.C_COE, 8)
for L in LOOPS: mk(L["mig"][0], L["mig"][1], "^" if L["kind"] == "MAX" else "s", MKF, 6.5, edge=CO[L["key"]], mew=1.8)
a.set_title("5. Overlay on Migkas Fig. 8, $\\theta=75^\\circ$ ($A\\propto H_0^2$): four full loops,\nopacity gradient darkest at $\\theta$ ray (CoE axis white; celestial equator black)", fontsize=10, pad=4)
# insets on the small MAX-2 ring (steps 2-5)
IW, IH = 1.15 / FW, 1.15 / FH
for k, gr, cols, bb in ((2, False, CB, UNT["b"]), (3, False, CB, FIT["b"]), (4, True, CB, FIT["b"])):
    x, y, w = geo[k]; coe_inset([x + w - IW + 0.01, y - 0.004, IW, IH], lb_to_n(l0, bb), (UNT if k == 2 else FIT)["roll"], gr, cols, 2.6)
x, y, w = geo[5]; coe_inset([x + w - IW + 0.004, y + 0.004, IW, IH], lb_to_n(l0, FIT["b"]), FIT["roll"], True, CO, 2.6)
# arrows between panels (snake: 0 -> 1 right; 1 -> 2 down; 2 -> 3 left; 3 -> 4 down; 4 -> 5 right)
def arrow(p0, p1, rad=0.0):
    fig.add_artist(FancyArrowPatch((p0[0] / FW, 1 - p0[1] / FH), (p1[0] / FW, 1 - p1[1] / FH), transform=fig.transFigure,
                   arrowstyle="simple,head_width=12,head_length=14,tail_width=5", color="0.35", connectionstyle="arc3,rad=%s" % rad))
arrow((12.45, 3.6), (13.0, 3.6)); arrow((16.1, 6.8), (16.1, 7.75)); arrow((10.0, 10.55), (9.35, 10.55))
arrow((4.9, 13.0), (4.9, 13.45)); arrow((9.35, 16.1), (10.0, 16.1))
# legend (panel-B style) listing all four loops, plus the overlay colour key
def lab(L): return "%s-%s loop $\\theta$=%.1f° ($\\ell$=%.1f°), band ±9°%s" % (L["kind"], L["key"][-1], L["theta"], L["l"], (" (drawn +%.0f°)" % MIN_OFF) if L["kind"] == "MIN" else "")
h = [Line2D([0], [0], color=RB.C_AXIS, lw=2.5, label="CoE$\\leftrightarrow$anti-CoE axis (great circle)")]
h += [Patch(facecolor=CB[L["key"]], alpha=0.55, edgecolor=CB[L["key"]], label=lab(L)) for L in LOOPS]
h += [Line2D([0], [0], color=RB.C_EQ, lw=1.8, ls="--", label="celestial equator (Dec=0)")]
h += [Line2D([0], [0], marker="^" if L["kind"] == "MAX" else "s", color="w", markerfacecolor=MKF, markeredgecolor=CB[L["key"]], markeredgewidth=2.0, ms=8, ls="None",
             label="Migkas %s (%.0f°, %+.0f°)" % (L["kind"], L["mig"][0], L["mig"][1])) for L in LOOPS]
h += [Line2D([0], [0], marker="o", color="w", markerfacecolor=RB.C_ANTI, markeredgecolor="k", ms=8, ls="None", label="anti-CoE"),
      Line2D([0], [0], marker="*", color="w", markerfacecolor=RB.C_COE, markeredgecolor="k", ms=11, ls="None", label="CoE")]
fig.legend(handles=h, loc="lower center", ncol=4, fontsize=8.2, frameon=True, fancybox=True, bbox_to_anchor=(0.5, 0.013))
fig.text(0.5, 0.004, "Display offset: the two MIN loops (true cone radii %.1f° and %.1f° about anti-CoE) are each drawn %.0f° further along $\\theta$ so both lines are seen; all numbers use the true positions. "
         "Insets zoom on the %.1f°-radius MAX-2 ring round the CoE." % (LOOPS[2]["theta"], 360 - LOOPS[3]["theta"], MIN_OFF, 180 - LOOPS[1]["theta"]), ha="center", va="bottom", fontsize=8, style="italic", color="0.25")
fig.suptitle("From the homologous blast and your model A($\\ell$) to loops on the sky: thin-shell best fit, full loops, fitted tilt, gradients and the Migkas overlay"
             "   (Galactic Mollweide; $\\ell$ increases LEFT)", fontsize=14, y=0.992)
fig.savefig("../loop_flow_thinshell.png", dpi=150)
im = Image.open("../loop_flow_thinshell.png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../loop_flow_thinshell_web.jpg", quality=88)
print("ok", im.size)
