"""Journal page: A(l) + Chris's own 3-panel figure (make_coe_rings_gradient_loops.main) for three transformations
+ his Fig. 8 overlay (make_migkas_fig8_gradient_loops_overlay.main). Only colours, theta angles and axis orientation
(and the numbers in titles/legends) are changed; figures, layout and styling are his, unchanged."""
import json, sys, os, re, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from pathlib import Path
from PIL import Image
sys.path.insert(0, "/workspace/cosmology/figures"); sys.path.insert(0, "../bubble")
import importlib.util
_sp = importlib.util.spec_from_file_location("coe_gradient_loop_geometry", "/workspace/cosmology/figures/coe_gradient_loop_geometry.cos_bak.py")
CG = importlib.util.module_from_spec(_sp); sys.modules["coe_gradient_loop_geometry"] = CG; _sp.loader.exec_module(CG)
import make_coe_rings_gradient_loops as RB
import make_migkas_fig8_gradient_loops_overlay as OV
from coe_half_loop_geometry import lb_to_n, n_to_lb, MIGKAS_MAXES, MIGKAS_MINS, band_edges, mu_pair_for_band, BAND_HALF
_here = os.getcwd(); os.chdir("../bubble")
import analysis_curve as AC
from bubble_core import TH1, ML, MA, MLO, MHI, AFfour
os.chdir(_here)
OUT = Path("stack_rows").resolve(); OUT.mkdir(exist_ok=True)
F = json.load(open("loop_fit.json")); d = json.load(open("../bubble/thin_b_free_hug.json"))
A = AC.curve("ptshell", d); l0 = F["l0"]; FIT = F["fits"]["b0_only"]; UNT = F["fits"]["untilted"]
P = {(p["kind"], round(p["migkas"][0])): p for p in F["pairs"]}
pMAX1, pMAX2, pMIN1, pMIN2 = P[("MAX", 17)], P[("MAX", 185)], P[("MIN", 115)], P[("MIN", 274)]
C = {"MAX1": "#0097b2", "MIN1": "#6fa800", "MAX2": "#e0107a", "MIN2": "#2166d4"}
THMAX = [pMAX1["theta"], pMAX2["theta"]]; THMIN = [pMIN1["theta"], pMIN2["theta"]]
# ---- colours and angles (only changes) ----
for M in (RB, OV):
    M.C_MAX_30, M.C_MAX_185, M.C_MIN_115, M.C_MIN_265 = C["MAX1"], C["MAX2"], C["MIN1"], C["MIN2"]
    M.COLORS_MAX = [C["MAX1"], C["MAX2"]]; M.COLORS_MIN = [C["MIN1"], C["MIN2"]]
OV.THETA_MAX, OV.THETA_MIN = list(THMAX), list(THMIN); OV.MAX_BANDS, OV.MIN_BANDS = band_edges(THMAX), band_edges(THMIN)
RB.STRIP_THETA_MAX, RB.STRIP_THETA_MIN = list(THMAX), list(THMIN)
RB.OUT_DIR = OUT; OV.OUT_DIR = OUT
def rays(ax, geo):   # his draw_panelA_theta_rays with the emphasised ray/arc at the MAX-1 angle instead of the fixed 30 deg
    th1 = geo["theta_max_from_plus_x_deg"][0]
    for mid, color in list(zip(geo["theta_max_from_plus_x_deg"], RB.COLORS_MAX)) + list(zip(geo["theta_min_from_plus_x_deg"], RB.COLORS_MIN)):
        t = np.radians(mid); x1, y1 = np.cos(t), np.sin(t); main = abs(mid - th1) < 0.5
        ax.plot([0.0, x1], [0.0, y1], color=color, lw=3.2 if main else 1.8, alpha=0.70 if main else 0.45, zorder=3, solid_capstyle="round")
        ax.plot([x1], [y1], "o", color=color, ms=4.5 if main else 3.5, alpha=0.9, zorder=5, markeredgecolor="none")
    phi = np.linspace(0.0, np.radians(th1), 60); ax.plot(0.38 * np.cos(phi), 0.38 * np.sin(phi), color=RB.C_MAX_30, lw=1.6, alpha=0.85, zorder=3)
    ax.text(0.38 * np.cos(np.radians(th1 / 2)) + 0.06, 0.38 * np.sin(np.radians(th1 / 2)) + 0.04, r"$\theta$", color=RB.C_MAX_30, fontsize=10, fontweight="bold", ha="left", va="bottom", zorder=6)
RB.draw_panelA_theta_rays = rays
# ---- axis orientation per transformation ----
STATE = {}
def pose_for(b0):
    D = lb_to_n(l0, b0); return {"anti_lb": n_to_lb(D), "coe_lb": n_to_lb(-D), "n_anti": D, "n_coe": -D, "D": D, "tilt_b_deg": float(b0)}
def geometry(tilt_b, roll_deg):
    pose = pose_for(tilt_b); mb, nb = band_edges(THMAX), band_edges(THMIN)
    n_max = [lb_to_n(*lb) for lb in MIGKAS_MAXES]; n_min = [lb_to_n(*lb) for lb in MIGKAS_MINS]
    return {**pose, "n_maxes": n_max, "n_mins": n_min, "theta_max_from_plus_x_deg": list(THMAX), "theta_min_from_plus_x_deg": list(THMIN),
            "max_bands_deg": mb, "min_bands_deg": nb, "mu_max_bands": [mu_pair_for_band(*b) for b in mb], "mu_min_bands": [mu_pair_for_band(*b) for b in nb],
            "sep_deg": float(np.degrees(np.arccos(np.clip(np.dot(n_max[0], n_min[1]), -1, 1)))), "band_half_width_deg": BAND_HALF,
            "roll_about_coe_deg": float(roll_deg), "full_gradient_loops": True, "gradient_rule": CG.GRADIENT_RULE}
RB.compute_geometry = geometry; OV.axis_pose = pose_for
# ---- update the numbers in his titles/legends to the latest parameters (text only) ----
exact = {int(round(t)): t for t in THMAX + THMIN}
def fix(s):
    b0 = STATE["b0"]; roll = STATE["roll"]; an, cn = n_to_lb(lb_to_n(l0, b0)), n_to_lb(-lb_to_n(l0, b0))
    s = re.sub(r"(?<!anti-)CoE (\$?)\(\$?\s*-?\d+\^\\circ,\s*[+-−]?\s*\d+\^\\circ\s*\$?\)(\$?)", lambda m: "CoE %s(%.1f^\\circ,%+.1f^\\circ)%s" % (m.group(1), cn[0], cn[1], m.group(2)), s)
    _r = "CoE $(%.1f^\\circ,%+.1f^\\circ)$" % cn; s = re.sub(r"(?<!anti-)CoE \$\((\d+)\^\\circ,\s*[+-−]\d+\^\\circ\)\$", lambda m: _r, s)
    s = re.sub(r"anti-CoE( at| \$\+x\$| \+x)? \$\(\s*-?\d+\^\\circ,\s*[+-]?\d+\^\\circ\)\$", lambda m: "anti-CoE%s $(%.1f^\\circ,%+.1f^\\circ)$" % (m.group(1) or "", an[0], an[1]), s)
    _r2 = "tilted $%+.1f^\\circ$" % b0; s = re.sub(r"tilted \$\+\d+\^\\circ\$", lambda m: _r2, s)
    _r3 = "$+%.1f^\\circ$ CCW" % roll; s = re.sub(r"\$\+\d+\^\\circ\$ CCW", lambda m: _r3, s)
    def th(m): return re.sub(r"(\d+)\^\\circ", lambda n: ("%.1f^\\circ" % exact[int(n.group(1))]) if int(n.group(1)) in exact else n.group(0), m.group(0))
    s = re.sub(r"\\theta\\!=\\![^$]*", th, s)
    s = s.replace("fading to 12% opacity opposite", "transparent opposite").replace(r"$\alpha\propto1-0.88|\psi|/\pi$", r"$\alpha\propto(1+\cos\psi)/2$").replace(r"$\alpha = 1 - 0.88|\psi|/\pi$", r"$\alpha\propto(1+\cos\psi)/2$")
    s = s.replace("A. Model frame", r"A. Galactic $\ell$ frame").replace(r"axisymmetric about $x$ = CoE$\leftrightarrow$anti-CoE", r"$\ell$ CCW from $+x$; CoE axis at $\ell_0$=%.1f°" % l0)
    if not STATE["grad"]:
        s = s.replace(r"opacity gradient darkest at $\theta$ ray, fading to 12% opacity opposite", "uniform opacity (no gradient)")
        s = s.replace(r"$\alpha\propto1-0.88|\psi|/\pi$, $\psi$ = ring azimuth from $\theta$-ray point  ·  ", "")
        s = s.replace(r" ($\alpha\propto1-0.88|\psi|/\pi$)", "")
        s = s.replace(r"darkest at $\theta$ ray, fading to 12% opacity opposite; ", "").replace("full loops, opacity gradient", "full loops, uniform opacity").replace("full loop, gradient", "full loop, uniform").replace("full loops, gradient", "full loops, uniform")
    s = s.replace("-0.0^\\circ", "0.0^\\circ").replace("+0.0^\\circ", "0.0^\\circ")
    return s
def fix_fig(fig):
    objs = list(fig.texts) + ([fig._suptitle] if fig._suptitle else [])
    for ax in fig.axes:
        objs += list(ax.texts) + [ax.title, ax._left_title, ax._right_title]
        if ax.get_legend(): objs += list(ax.get_legend().get_texts())
    for lg in fig.legends: objs += list(lg.get_texts())
    for t in objs:
        s = t.get_text()
        if s: t.set_text(fix(s))
_sf = Figure.savefig
from mpl_toolkits.mplot3d import proj3d
SIZES = []
def match_sphere(fig):
    """only size change: zoom panel C so the drawn sphere diameter equals panel A's circle diameter"""
    ax3 = [a for a in fig.axes if a.name == "3d"]; ax2 = [a for a in fig.axes if a.name == "rectilinear"]
    if not ax3 or not ax2: return
    aC, aA = ax3[0], ax2[0]; z = 1.0
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 90), np.linspace(-np.pi / 2, np.pi / 2, 45))
    X, Y, Z = (np.cos(u) * np.cos(v)).ravel(), (np.sin(u) * np.cos(v)).ravel(), np.sin(v).ravel()
    for _ in range(6):
        aC.set_box_aspect((1, 1, 1), zoom=z); fig.canvas.draw()
        dA = np.ptp(aA.transData.transform([(-1, 0), (1, 0)])[:, 0])
        x2, y2, _ = proj3d.proj_transform(X, Y, Z, aC.get_proj()); P_ = aC.transData.transform(np.c_[x2, y2])
        dC = max(np.ptp(P_[:, 0]), np.ptp(P_[:, 1]))
        if abs(dC / dA - 1) < 0.01: break
        z *= dA / dC
    SIZES.append((round(float(dA), 1), round(float(dC), 1), round(z, 3)))
def savefig_fixed(self, *a, **k):
    fix_fig(self); match_sphere(self); return _sf(self, *a, **k)
Figure.savefig = savefig_fixed


# ---- panel A in Galactic l: rotate the model-frame drawing by l0 (l=0 at +x, CCW); loops land at their own l ----
import matplotlib.transforms as mtransforms
_panelA = RB.panel_A_model_frame
LTXT = []
def panelA_gal(ax, geo, roll_deg):
    _panelA(ax, geo, roll_deg)
    rot = mtransforms.Affine2D().rotate_deg(l0); c, s_ = np.cos(np.radians(l0)), np.sin(np.radians(l0))
    R2 = lambda x, y: (c * x - s_ * y, s_ * x + c * y)
    for ln in ax.lines: ln.set_transform(rot + ax.transData)
    for co in ax.collections: co.set_transform(rot + ax.transData)
    for t in list(ax.texts):
        txt = t.get_text()
        if isinstance(t, matplotlib.text.Annotation):
            t.xy = R2(*t.xy); t.set_position(R2(*t.get_position())); continue
        if "CoE" in txt or "theta\\!=" in txt: t.remove(); continue
        x, y = t.get_position()
        if y > -1.3: t.set_position(R2(x, y))      # keep the caption block under the circle where it was
    # thin l reference diameters + l ticks/labels
    for a_ in (0, 90):
        u = np.array([np.cos(np.radians(a_)), np.sin(np.radians(a_))]); ax.plot([-u[0], u[0]], [-u[1], u[1]], color="0.55", lw=0.6, ls=(0, (2, 2)), zorder=1.5)
    for a_, ha, va, off in ((0, "left", "top", (0.02, -0.03)), (90, "center", "bottom", (0, 0.02)), (180, "right", "bottom", (-0.02, 0.03)), (270, "center", "top", (0, -0.02))):
        u = np.array([np.cos(np.radians(a_)), np.sin(np.radians(a_))]); ax.plot(*zip(u * 1.0, u * 1.05), color="0.35", lw=0.9, zorder=2)
        ax.text(*(u * 1.06 + off), r"$\ell$=%d°" % a_, color="0.35", fontsize=7, ha=ha, va=va, zorder=6)
    ea = np.array([c, s_])
    ax.text(*(ea * 1.20), "anti-CoE\n$\\ell$=%.1f°" % l0, color=RB.C_ANTI, fontsize=8, ha="left", va="center")
    ax.text(*(-ea * 1.20), "CoE\n$\\ell$=%.1f°" % ((l0 + 180) % 360), color=RB.C_COE, fontsize=8, ha="right", va="center")
    for mid, col in list(zip(geo["theta_max_from_plus_x_deg"], RB.COLORS_MAX)) + list(zip(geo["theta_min_from_plus_x_deg"], RB.COLORS_MIN)):
        L = (mid + l0) % 360; u = np.array([np.cos(np.radians(L)), np.sin(np.radians(L))]); pos = u * 1.30; va = "center"
        if abs(L - 187.5) < 1: pos = pos + [0, 0.05]; va = "bottom"
        LTXT.append(ax.text(*pos, r"$\ell\!=\!%.1f^\circ$" % L, color=col, fontsize=7.5, ha="center", va=va, fontweight="bold"))
    ax.figure.canvas.draw(); rr_ = ax.figure.canvas.get_renderer()
    bbs = [(t.get_text()[:16], t.get_window_extent(rr_)) for t in ax.texts if t.get_text()]
    from itertools import combinations
    for (n1, b1), (n2, b2) in combinations(bbs, 2):
        if b1.overlaps(b2): print("A-OVERLAP", repr(n1), repr(n2))
RB.panel_A_model_frame = panelA_gal
import io, contextlib
RUNS = []
for ROLL, TAG in ((45.0, "roll45"), (FIT["roll"], "roll80p8")):
    STATE.update(b0=FIT["b"], roll=ROLL, grad=True)
    png, pdf, js = RB.main(tilt_b=FIT["b"], roll_deg=ROLL)
    three = OUT / ("threepanel_%s.png" % TAG); os.replace(png, three)
    for f in (pdf, js): os.remove(f)
    with contextlib.redirect_stdout(io.StringIO()): OV.main(tilt_b=FIT["b"], roll_deg=ROLL)
    ov = sorted(OUT.glob("migkas_fig8_gradient_loops_*.png"), key=os.path.getmtime)[-1]; ovp = OUT / ("overlay_%s.png" % TAG); os.replace(ov, ovp)
    for f in OUT.glob("migkas_fig8_gradient_loops_*.json"): os.remove(f)
    RUNS.append((TAG, three, ovp))
Figure.savefig = _sf
# ---- A vs angle, Figure-2 style ----
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
U = json.load(open("../bubble/universe_thin.json"))["free"]
th2 = np.append(TH1, 360); Af2 = AFfour(th2); A2 = np.append(A, A[0])
fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, 5.25), sharex=True, gridspec_kw={"height_ratios": [2.0, 1], "hspace": 0.06})
for k, p in (("MAX1", pMAX1), ("MIN1", pMIN1), ("MAX2", pMAX2), ("MIN2", pMIN2)):
    L = p["l_model"]
    for ax_ in (a, b): ax_.axvspan(L - 9, L + 9, color=C[k], alpha=0.22, lw=0, zorder=0)
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
a.plot(th2, Af2, color="C0", lw=1.5, label="$A_F$")
a.plot(th2, A2, color="C3", lw=1.6, label="thin shell + core, $\\Omega_s$=%.2f, $\\ell_0$=%.1f° ($R^2$=%.3f; vs $A_F$ %.3f)" % (U["Os"], U["l0"], U["R2"], U["R2F"]))
b.plot(th2, 100 * (A2 - Af2) / Af2, color="C3", lw=1.4, label="model")
a.set_ylabel("$A/A_{\\rm all}$"); a.set_ylim(0.62, 1.22); a.legend(fontsize=7.2, loc="lower left", frameon=False)
b.axhline(0, color="k", lw=0.6); b.plot(ML, 100 * (MA - AFfour(ML)) / AFfour(ML), "o", ms=2.5, color="0.45", label="data")
b.set_ylabel("$100\\,(\\cdot-A_F)/A_F$ [%]", fontsize=8); b.set_xlabel("Galactic longitude $\\ell$ (deg)")
b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45)); b.legend(fontsize=7, frameon=False, loc="upper right", ncol=2)
agr = OUT / "Agraph.png"; fig.savefig(agr, dpi=200, bbox_inches="tight", pad_inches=0.04); plt.close(fig)
# ---- stack at natural sizes ----
def trim(im, pad=12):
    g = np.asarray(im.convert("L")) < 250; ys, xs = np.nonzero(g)
    return im.crop((max(xs.min() - pad, 0), max(ys.min() - pad, 0), min(xs.max() + pad + 1, im.width), min(ys.max() + pad + 1, im.height)))
for TAG, three, ovp in RUNS:
    parts = [trim(Image.open(p).convert("RGB")) for p in (agr, three, ovp)]
    W = max(p.width for p in parts); H = sum(p.height for p in parts)
    page = Image.new("RGB", (W, H), "white"); y = 0
    for p in parts: page.paste(p, ((W - p.width) // 2, y)); y += p.height
    page.save("../loop_stack_%s.png" % TAG); w = page.copy(); w.thumbnail((1700, 1700)); w.save("../loop_stack_%s_web.jpg" % TAG, quality=88)
    for n, p in (("loop_threepanel_" + TAG, parts[1]), ("loop_overlay_" + TAG, parts[2])):
        p.save("../%s.png" % n); q = p.copy(); q.thumbnail((1700, 1700)); q.save("../%s_web.jpg" % n, quality=88)
    print(TAG, [p.size for p in parts], page.size)
q = parts[0].copy(); parts[0].save("../loop_Agraph.png"); q.thumbnail((1700, 1700)); q.save("../loop_Agraph_web.jpg", quality=88)
print("sizes", SIZES)
