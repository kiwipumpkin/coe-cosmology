"""Journal page: A(l) + Chris's own 3-panel figure (make_coe_rings_gradient_loops.main) for three transformations
+ his Fig. 8 overlay (make_migkas_fig8_gradient_loops_overlay.main). Only colours, theta angles and axis orientation
(and the numbers in titles/legends) are changed; figures, layout and styling are his, unchanged."""
import json, sys, os, re, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from pathlib import Path
from PIL import Image
sys.path.insert(0, "/workspace/cosmology/figures"); sys.path.insert(0, "../bubble")
import coe_gradient_loop_geometry as CG
import make_coe_rings_gradient_loops as RB
import make_migkas_fig8_gradient_loops_overlay as OV
from coe_half_loop_geometry import lb_to_n, n_to_lb, MIGKAS_MAXES, MIGKAS_MINS, band_edges, mu_pair_for_band, BAND_HALF
_here = os.getcwd(); os.chdir("../bubble")
import analysis_curve as AC
from bubble_core import TH1, ML, MA, MLO, MHI
os.chdir(_here)
OUT = Path("journal_rows").resolve()
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
# ---- rows 1-3: his figure, unchanged apart from the above ----
rows = []
for name, b0, roll, grad in (("row1_untilted", 0.0, 0.0, False), ("row2_tilt", FIT["b"], 0.0, False), ("row3_gradient", FIT["b"], FIT["roll"], True)):
    STATE.update(b0=b0, roll=roll, grad=grad); CG.ALPHA_FLOOR = 0.12 if grad else 1.0
    png, pdf, js = RB.main(tilt_b=b0, roll_deg=roll)
    tgt = OUT / (name + ".png"); os.replace(png, tgt); rows.append(tgt)
    for f in (pdf, js): os.remove(f)
CG.ALPHA_FLOOR = 0.12
# ---- row 4: his overlay ----
STATE.update(b0=FIT["b"], roll=FIT["roll"], grad=True)
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()): OV.main(tilt_b=FIT["b"], roll_deg=FIT["roll"])
ov = sorted(OUT.glob("migkas_fig8_gradient_loops_*.png"), key=os.path.getmtime)[-1]; tgt = OUT / "row4_overlay.png"; os.replace(ov, tgt)
for f in OUT.glob("migkas_fig8_gradient_loops_*.json"): os.remove(f)
Figure.savefig = _sf
# ---- row 0: A vs angle ----
W0 = Image.open(rows[0]).width
fig, a = plt.subplots(figsize=(W0 / 160, 3.3), dpi=160)
th2 = np.append(TH1, 360); A2 = np.append(A, A[0])
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.8", lw=0.8, label="Migkas Fig. 6")
a.plot(th2, A2, color="k", lw=1.8, label="your model A, thin shell")
for k, p in (("MAX1", pMAX1), ("MIN1", pMIN1), ("MAX2", pMAX2), ("MIN2", pMIN2)):
    L = p["l_model"]; a.axvspan(L - 9, L + 9, color=C[k], alpha=0.22, lw=0); a.axvline(L, color=C[k], lw=1.3)
    a.text(L, 1.215, "%s  $\\ell$=%.1f°" % (p["kind"], L), color=C[k], ha="center", va="bottom", fontsize=10, fontweight="bold")
a.set_xlim(0, 360); a.set_xticks(range(0, 361, 45)); a.set_ylim(0.76, 1.21)
a.set_xlabel("Galactic longitude $\\ell$ (deg)"); a.set_ylabel("$A/A_{\\rm all}$ (higher = faster)")
a.legend(fontsize=9, loc="lower left", frameon=False)
a.set_title("Your model A($\\ell$), thin-shell best fit: $\\ell_0$ = %.2f°, Amp = %.3f, $\\varphi$ = %.3f rad, $R^2$ = %.3f (bands ±9°)" % (l0, F["amp"], F["phi"], F["R2"]), fontsize=12, fontweight="bold", pad=26)
fig.subplots_adjust(left=0.06, right=0.985, top=0.80, bottom=0.17); fig.savefig(OUT / "row0_A.png", dpi=160, facecolor="white"); plt.close(fig)
# ---- caption strips (label + one line) with a small down arrow ----
def strip(label, cap, W):
    fig = plt.figure(figsize=(W / 160, 0.62), dpi=160)
    fig.patches.append(matplotlib.patches.FancyArrow(0.5, 0.97, 0, -0.30, width=0.0025, head_width=0.012, head_length=0.18, color="0.45", transform=fig.transFigure, length_includes_head=True))
    fig.text(0.012, 0.22, label, fontsize=12, fontweight="bold", va="center"); fig.text(0.09, 0.22, cap, fontsize=11.5, va="center", style="italic")
    b = io.BytesIO(); fig.savefig(b, dpi=160, facecolor="white"); plt.close(fig); b.seek(0); return Image.open(b).convert("RGB")
caps = [("Transf. 1", "Transformation 1: each extremum becomes a full loop of radius $\\theta=\\ell_{\\rm ext}-\\ell_0$ about anti-CoE (%.1f°, 0°), untilted; ring RMS miss %.1f°" % (l0, UNT["rms_ring"])),
        ("Transf. 2", "Transformation 2: tilt CoE axis $b_0$ = %+.1f° (anti-CoE (%.1f°, %+.1f°)); ring RMS miss %.1f° to %.1f°" % (FIT["b"], l0, FIT["b"], UNT["rms_ring"], FIT["rms_ring"])),
        ("Transf. 3", "Transformation 3: roll %.1f° about the CoE axis; loop gradient darkest at the $\\theta$ ray, fading to 12%% opposite; $\\theta$-ray RMS miss %.1f°" % (FIT["roll"], FIT["rms_ray"])),
        ("Transf. 4", "Transformation 4: overlay the tilted, rolled gradient loops on Migkas Fig. 8, $\\theta=75^\\circ$ ($A\\propto H_0^2$)")]
parts = [Image.open(OUT / "row0_A.png").convert("RGB")]
for (lab, cap), r in zip(caps, rows + [tgt]):
    im = Image.open(r).convert("RGB")
    if im.width != W0: im = im.resize((W0, int(im.height * W0 / im.width)), Image.LANCZOS)
    parts += [strip(lab, cap, W0), im]
Hs = sum(p.height for p in parts); page = Image.new("RGB", (W0, Hs), "white"); y = 0
for p in parts: page.paste(p, (0, y)); y += p.height
page.save("../loop_journal_page.png", dpi=(160, 160))
w = page.copy(); w.thumbnail((1700, 1700)); w.save("../loop_journal_page_web.jpg", quality=88)
cuts = [0, sum(p.height for p in parts[:5]), sum(p.height for p in parts[:7]), page.height]   # legible web parts: rows 0-2 | row 3 | overlay
for i in range(3):
    q = page.crop((0, cuts[i], page.width, cuts[i + 1])); q.thumbnail((1700, 1700)); q.save("../loop_journal_page_web_part%d.jpg" % (i + 1), quality=88)
print("A vs C diameters (px at 160 dpi), zoom:", SIZES); print("page", page.size, "web", w.size, [p.size for p in parts])
