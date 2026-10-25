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
fig = plt.figure(figsize=(19, 11.6))
Hh = 0.33
geo = {0: (0.015, 0.55, 0.30), 1: (0.395, 0.55, 0.255), 2: (0.69, 0.55, 0.295), 3: (0.015, 0.10, 0.295), 4: (0.34, 0.10, 0.295), 5: (0.665, 0.10, 0.325)}
axs = {k: fig.add_axes([x, y, w, Hh]) for k, (x, y, w) in geo.items()}
# 0. homologous blast
a = axs[0]; a.set_aspect("equal"); a.axis("off"); a.set_xlim(-0.6, 7.6); a.set_ylim(-3.0, 2.3)
def vec(p, q, c, lw=1.8, ls="-", head=True):
    a.annotate("", xy=q, xytext=p, arrowprops=dict(arrowstyle="-|>" if head else "-", color=c, lw=lw, ls=ls, shrinkA=0, shrinkB=0))
CoE = np.array([0.0, 0.0]); O = np.array([2.0, 0.0]); Pp = np.array([3.05, 1.25]); Hs = 0.55
tt = np.linspace(0, 2 * np.pi, 200); a.plot(O[0] + 1.62 * np.cos(tt), O[1] + 1.62 * np.sin(tt), color="0.6", lw=0.8, ls="--")
vec(CoE, O, "#c44e52"); vec(CoE, Pp, "#c44e52", lw=1.2, ls="--"); vec(O, Pp, "0.2", lw=1.0, ls=":")
vec(O, O + Hs * O, "#4c72b0", lw=2.0); vec(Pp, Pp + Hs * Pp, "#55a868", lw=2.0); vec(Pp, Pp + Hs * (Pp - O), "#dd4444", lw=2.0)
for p_, c in ((CoE, "#c44e52"), (O, "#4c72b0"), (Pp, "#55a868")): a.plot(*p_, "o", color=c, ms=6, zorder=5)
a.text(-0.15, -0.25, "CoE", color="#c44e52", fontsize=8, ha="right", va="top"); a.text(2.0, -0.28, "us (O)", color="#4c72b0", fontsize=8, ha="center", va="top")
a.text(3.15, 1.1, "galaxy P", color="#55a868", fontsize=7.5, ha="left", va="top")
a.text(1.0, -0.12, "$\\vec R_O$", color="#c44e52", fontsize=8.5, ha="center", va="top"); a.text(1.2, 0.85, "$\\vec R_P$", color="#c44e52", fontsize=8.5)
a.text(2.62, 0.5, "$\\vec r$", fontsize=8.5); a.text(3.15, 0.2, "$\\vec V_O=H\\vec R_O$", color="#4c72b0", fontsize=7.5)
a.text(4.35, 2.05, "$\\vec V_P=H\\vec R_P$", color="#55a868", fontsize=7.5, ha="center"); a.text(3.05, 2.1, "$\\vec V_H$", color="#dd4444", fontsize=8, ha="right")
# velocity triangle
T0 = np.array([5.6, 0.2]); VP = Hs * Pp; VO = Hs * O
vec(T0, T0 + VP, "#55a868", lw=1.8); vec(T0 + VP, T0 + VP - VO, "#4c72b0", lw=1.8); vec(T0, T0 + VP - VO, "#dd4444", lw=1.8)
a.text(T0[0] + VP[0] / 2 + 0.25, T0[1] + VP[1] / 2 - 0.2, "$\\vec V_P$", color="#55a868", fontsize=8); a.text(T0[0] + VP[0] - VO[0] / 2, T0[1] + VP[1] + 0.12, "$-\\vec V_O$", color="#4c72b0", fontsize=8, ha="center")
a.text(T0[0] - 0.1, T0[1] + 0.45, "$\\vec V_H$", color="#dd4444", fontsize=8, ha="right")
a.text(6.1, -0.35, "$\\vec V_H=\\vec V_P-\\vec V_O=H\\vec r$", fontsize=8, ha="center", va="top")
a.text(-0.5, -1.95, "Proof: every point moves radially, $\\vec V=H\\vec R$. So $\\vec V_O=H\\vec R_O$, $\\vec V_P=H\\vec R_P$, and\n"
       "$\\vec V_H=\\vec V_P-\\vec V_O=H(\\vec R_P-\\vec R_O)=H\\vec r$: it depends only on $\\vec r$, not on where the CoE is.\n"
       "Every observer sees $|\\vec V_H|=H r$ in every direction, i.e. an isotropic Hubble flow.",
       fontsize=7.6, va="top", ha="left", bbox=dict(fc="#f5f5f5", ec="0.7", pad=4))
a.text(-0.5, -3.05, "Pure homology leaves no CoE signature. The anisotropy comes from gravity slowing different\n"
       "positions by different amounts (a departure from homology), which your model A captures (step 1).",
       fontsize=7.6, va="top", ha="left", color="#8e1b1b", style="italic")
a.set_title("0. Homologous blast $\\vec V=H\\vec R$ from the CoE", fontsize=10.5, pad=30)
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
a.set_title("1. $A(\\ell)$, latest fit: $\\ell_0$=%.2f°, Amp=%.3f, $\\varphi$=%.3f rad" % (l0, F["amp"], F["phi"]), fontsize=10.5, pad=30)
# 2. untilted
a = axs[2]; sky_axes(a); D, Y = draw_loops(a, 0.0, l0, UNT["roll"], False); markers(a, D)
x, y = moll(l0, 0); a.text(x - 0.36, y + 0.22, "anti-CoE\n(%.1f°, 0°)" % l0, fontsize=7, color="#8e1b1b", ha="right")
a.set_title("2. Full loops, untilted (anti-CoE on the plane)\nloop radius $\\theta=\\ell_{\\rm ext}-\\ell_0$; ring rms miss %.1f°" % UNT["rms_ring"], fontsize=10.5)
# 3. tilt fit
a = axs[3]; sky_axes(a); D, Y = draw_loops(a, FIT["b"], l0, FIT["roll"], False); markers(a, D); axis_line(a, D, Y)
x, y = moll(l0, FIT["b"]); a.text(x - 0.36, y + 0.22, "anti-CoE\n(%.1f°, %+.1f°)" % (l0, FIT["b"]), fontsize=7, color="#8e1b1b", ha="right")
a.set_title("3. Best-fit tilt $b_0$=%+.1f° puts Migkas extrema on the loops\nring misses %s°; all inside the ±9° bands" % (FIT["b"], ", ".join("%+.1f" % v for v in FIT["ring"])), fontsize=10.5)
# 4. gradient
a = axs[4]; sky_axes(a); D, Y = draw_loops(a, FIT["b"], l0, FIT["roll"], True); markers(a, D); axis_line(a, D, Y)
a.set_title("4. Loop gradients: darkest at the $\\theta$ ray, fading opposite\n(fitted roll %.1f° about the CoE axis; the two MIN loops nearly coincide)" % FIT["roll"], fontsize=10.5)
# 5. overlay
sys.path.insert(0, "/workspace/cosmology/figures")
from make_migkas_fig8_gradient_loops_overlay import GOOD_REG, lb_to_px
rgb = np.asarray(Image.open("/workspace/cosmology/map_compare/FULL_migkas_theta75.png").convert("RGB"))
cx, cy, Rx, Ry = GOOD_REG["cx"], GOOD_REG["cy"], GOOD_REG["Rx"], GOOD_REG["Ry"]
a = axs[5]; a.imshow(rgb, origin="upper", extent=(0, rgb.shape[1], rgb.shape[0], 0))
x0, x1 = cx - Rx - 150, min(rgb.shape[1], cx + Rx + 330); y0, y1 = cy - Ry - 110, cy + Ry + 130
a.set_xlim(x0, x1); a.set_ylim(y1, y0); a.set_aspect("equal"); a.axis("off")
to_px = lambda l, b: tuple(np.asarray(v, float) for v in lb_to_px(l, b, cx, cy, Rx, Ry))
clip = Ellipse((cx, cy), 2 * Rx * 0.985, 2 * Ry * 0.985, transform=a.transData, fc="none", ec="none"); a.add_patch(clip)
JUMP = 0.8 * Rx / 2.83
nb = len(a.collections)
D, Y = draw_loops(a, FIT["b"], l0, FIT["roll"], True, to_xy=to_px, lw=2.0, under=True); axis_line(a, D, Y, to_xy=to_px, color="white")
for c in a.collections[nb:]: c.set_clip_path(clip)
markers(a, D, to_xy=to_px, ms=8)
a.set_title("5. Overlay on Migkas Fig. 8, $\\theta=75^\\circ$ ($A\\propto H_0^2$)", fontsize=10.5)
JUMP = 0.8
# arrows between panels
def arrow(p0, p1, rad=0.0):
    fig.add_artist(FancyArrowPatch(p0, p1, transform=fig.transFigure, arrowstyle="simple,head_width=10,head_length=12,tail_width=4", color="0.35", connectionstyle="arc3,rad=%s" % rad))
arrow((0.318, 0.715), (0.342, 0.715)); arrow((0.655, 0.715), (0.688, 0.715))
arrow((0.84, 0.535), (0.16, 0.50), rad=-0.05); arrow((0.312, 0.265), (0.338, 0.265)); arrow((0.637, 0.265), (0.663, 0.265))
# legend
h = [Patch(fc=COLS[k], alpha=0.6, label="%s loop, $\\ell$=%.1f° ($\\theta$=%.1f°, band ±9°)" % (p["kind"], p["l_model"], p["theta"])) for k, p in zip(keys, pairs)]
h += [Line2D([0], [0], marker="^" if p["kind"] == "MAX" else "s", color="w", mfc=COLS[k], mec="k", ms=8, ls="None",
             label="Migkas %s (%.0f°, %+.0f°)" % (p["kind"], p["migkas"][0], p["migkas"][1])) for k, p in zip(keys, pairs)]
h += [Line2D([0], [0], marker="o", color="w", mfc="#8e1b1b", mec="k", ms=7, ls="None", label="anti-CoE (+x)"),
      Line2D([0], [0], marker="*", color="w", mfc="#6b3fa0", mec="k", ms=10, ls="None", label="CoE"),
      Line2D([0], [0], color="#5b2c6f", ls="--", lw=1.2, label="$\\theta$-ray great circle (orbit plane)")]
fig.legend(handles=h, loc="lower center", ncol=4, fontsize=8.5, frameon=True, bbox_to_anchor=(0.5, 0.0))
fig.suptitle("From the homologous blast and your model A($\\ell$) to loops on the sky: thin-shell best fit, full loops, fitted tilt, gradients and the Migkas overlay"
             "   (Galactic Mollweide; $\\ell$ increases LEFT)", fontsize=13, y=0.985)
fig.savefig("../loop_flow_thinshell.png", dpi=150)
im = Image.open("../loop_flow_thinshell.png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../loop_flow_thinshell_web.jpg", quality=88)
print("ok", im.size)
