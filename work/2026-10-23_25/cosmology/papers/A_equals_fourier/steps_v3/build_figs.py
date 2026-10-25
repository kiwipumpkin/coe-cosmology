"""Appends the loop / map figure steps (Part B) to the notebook cell list C (called from build_nb.py)."""
import os, time
FIGS = "/workspace/cosmology/papers/A_equals_fourier/steps_v3/figs"
FD = "/workspace/cosmology/figures"; SM = "/workspace/cosmology/papers/spreadsheet_model"; LP = SM + "/loops"
def src(path, repl=()):
    s = open(path).read()
    for a, b in repl:
        assert a in s, (path, a); s = s.replace(a, b)
    return s
def mt(p): return time.strftime("%Y-%m-%d %H:%M CT", time.localtime(os.path.getmtime(p)))

def add(md, py):
    md(r"""
---
# Part B — Loops and maps

Part B redraws every loop and map figure made earlier in the project, using **the latest script for each figure, run verbatim**. The only edit to each script is where it writes its output: the images go to `steps_v3/figs/` so the approved originals on disk are not overwritten. That edit is listed in each step. After each script runs, the new image is shown next to the approved reference image and compared pixel by pixel.

**Shared geometry.** With a pure CoE flow, $A$ is axisymmetric about the CoE axis, so any value read at angle $\theta$ from the anti-CoE direction holds on the whole small circle (cone) of half-angle $\theta$ about that direction. The model $A(\ell)$ after the $2\ell$ step has two maxima ($\theta\approx30^\circ$, $185^\circ$) and two minima ($\theta\approx115^\circ$, $265^\circ$). Each one becomes a loop: a transparent band of $\pm9^\circ$ ($18^\circ$ end to end, 5% of $360^\circ$) about the centre $\theta$, plus one solid mid-line. The axis can be tilted to latitude $b_0$ and rolled about itself. Points on a loop are $\hat n=\cos\theta\,\hat a+\sin\theta(\cos\psi\,\hat e_1+\sin\psi\,\hat e_2)$, with $\hat a$ the anti-CoE unit vector and $\hat e_{1,2}$ perpendicular to it. They are converted to Galactic $(\ell,b)$ and drawn on a Mollweide projection (equal-area, Snyder 1987, *Map Projections — A Working Manual*, USGS PP 1395) with $\ell$ increasing to the left, as in Migkas Fig. 8.

**Data.** Migkas et al. 2020 (A&A 636, A15; arXiv:2004.03305): Table 2 directions of the $A$ extrema, and the Fig. 8 $\theta=75^\circ$ $A/A_{\rm all}$ map (background image loaded from `/workspace/cosmology/map_compare/FULL_migkas_theta75.png`). Markers: MAX at $(17^\circ,-9^\circ)$ and $(185^\circ,0^\circ)$; MIN at $(115^\circ,0^\circ)$ and $(274^\circ,-22^\circ)$.
""")
    py(r'''
# Part B setup - run each figure script verbatim in its own namespace; outputs go to steps_v3/figs
import os, sys, shlex, matplotlib, numpy as np
from pathlib import Path
from PIL import Image
from IPython.core.magic import register_cell_magic
from IPython.display import display
FIGS = Path("/workspace/cosmology/papers/A_equals_fourier/steps_v3/figs"); FIGS.mkdir(exist_ok=True)
for d in ("journal_rows", "stack_rows"): (FIGS / d).mkdir(exist_ok=True)
FIGDIR = "/workspace/cosmology/figures"
if FIGDIR not in sys.path: sys.path.insert(0, FIGDIR)

@register_cell_magic
def runscript(line, cell):
    """%%runscript <script name> [--cwd DIR] [--argv ARGS...]: exec the cell as that script (__name__ == '__main__'),
    __file__ placed in steps_v3/figs (so OUT_DIR = figs), matplotlib defaults, then restore cwd/argv/rc/modules."""
    tok = shlex.split(line); name = tok[0]; cwd = FIGS; argv = []
    if "--cwd" in tok: cwd = Path(tok[tok.index("--cwd") + 1])
    if "--argv" in tok: argv = tok[tok.index("--argv") + 1:]
    old = (os.getcwd(), sys.argv[:], matplotlib.rcParams.copy(), dict(sys.modules), sys.path[:])
    import matplotlib.pyplot as plt
    try:
        os.chdir(cwd); sys.argv = [name] + argv; matplotlib.rcdefaults(); matplotlib.use("Agg")
        g = {"__name__": "__main__", "__file__": str(FIGS / name)}
        exec(compile(cell, str(FIGS / name), "exec"), g)
    finally:
        plt.close("all"); os.chdir(old[0]); sys.argv = old[1]; matplotlib.rcParams.update(old[2])
        for k in list(sys.modules):
            if k not in old[3] and (k.startswith("coe_") or k.startswith("make_")): del sys.modules[k]
        for k, v in old[3].items():
            if k.startswith("coe_") or k.startswith("make_"): sys.modules[k] = v
        sys.path[:] = old[4]
        if FIGDIR not in sys.path: sys.path.insert(0, FIGDIR)

MATCH = []
def compare(new, ref, label=None):
    """approved reference (left) | notebook re-run (right), side by side; prints sizes and pixel differences."""
    from IPython.display import Image as IImage
    a = Image.open(ref).convert("RGB"); b = Image.open(new).convert("RGB")
    same = a.size == b.size; bb = b if same else b.resize(a.size, Image.LANCZOS)
    d = np.abs(np.asarray(a, float) - np.asarray(bb, float)); mad = d.mean(); frac = (d.max(axis=2) > 40).mean()
    verdict = "identical" if d.max() == 0 else ("match (tiny rendering differences)" if frac < 0.005 else "DIFFERS")
    MATCH.append(dict(figure=label or Path(new).name, reference=str(ref), size_ref=a.size, size_new=b.size, mean_abs_diff=round(mad, 3), frac_px_diff=round(frac, 5), verdict=verdict))
    print(f"{Path(new).name}: new {b.size} vs ref {a.size} | mean |diff| = {mad:.3f}/255, pixels differing >40: {100*frac:.3f}% -> {verdict}")
    print("left: approved reference " + str(ref).replace("/workspace/cosmology/", "") + "   |   right: notebook re-run steps_v3/figs/" + Path(new).name)
    W = 1100; h = round(a.size[1] * W / a.size[0]); pair = Image.new("RGB", (2 * W + 24, h), "white")
    pair.paste(a.resize((W, h), Image.LANCZOS), (0, 0)); pair.paste(b.resize((W, h), Image.LANCZOS), (W + 24, 0))
    out = FIGS / ("cmp_" + Path(new).stem + ".jpg"); pair.save(out, quality=85); display(IImage(filename=str(out)))
''')
    P = FIGS
    steps = []
    # (title, markdown, script, runscript line args, replacements, extra-code, [(new, ref)])
    steps.append((r"Step 9 — Mollweide maps of Migkas Table 2 $A$ extrema, and Fig. 8 with the published markers",
      r"""**What it shows.** (i) A clean 2D Galactic Mollweide map of the Migkas et al. 2020 Table 2 $A$ maximum and minimum directions for every sample and cone size, with no CoE loops and no model curves. (ii) The same published MAX/MIN markers placed on the real Migkas Fig. 8 ($\theta=75^\circ$) $A/A_{\rm all}$ map.
**Geometry.** Mollweide: $x=\tfrac{2\sqrt2}{\pi}\lambda\cos\vartheta$, $y=\sqrt2\sin\vartheta$, with $2\vartheta+\sin2\vartheta=\pi\sin b$ and $\lambda=-\mathrm{wrap}(\ell)$, so $\ell=0$ is at the centre and $\ell$ increases to the left. The Fig. 8 image is registered by its oval: centre and semi-axes $(c_x,c_y,R_x,R_y)=(863.0,763.25,693.99,392.54)$ px.""",
      FD + "/make_migkas_A_extrema_maps.py", "", [], "",
      [("migkas_A_extrema_table2_mollweide.png", FD + "/migkas_A_extrema_table2_mollweide.png"), ("migkas_fig8_theta75_A_extrema_only.png", FD + "/migkas_fig8_theta75_A_extrema_only.png")]))
    steps.append((r"Step 10 — First three-panel loop figure (axis tilted $+30^\circ$)",
      r"""**What it shows.** Panel A: the flat model-frame angle circle with the four extrema bands (two reds for MAX, two blues for MIN) and Galactic $\ell$ ticks. Panel B: the four loops on a Galactic Mollweide map with the Migkas markers. Panel C: the same loops on a 3D sphere. In this first version the CoE axis is tilted $+30^\circ$: anti-CoE at $(0^\circ,+30^\circ)$, CoE at $(180^\circ,-30^\circ)$.
**Maths.** Loop centres $\theta=30^\circ,185^\circ$ (MAX) and $115^\circ,265^\circ$ (MIN) from the anti-CoE direction. Bands $[\theta-9^\circ,\theta+9^\circ]$. A small circle of radius $\theta$ about $\hat a$ is traced by $\psi\in[0,2\pi)$.""",
      FD + "/make_coe_rings_migkas_pin_axis.py", "", [], "", [("coe_rings_migkas_pin_axis.png", FD + "/coe_rings_migkas_pin_axis.png")]))
    steps.append((r"Step 11 — The Step 10 loops overlaid on Migkas Fig. 8",
      r"""**What it shows.** The four MAX/MIN bands of Step 10 (axis $+30^\circ$) drawn on the real Migkas Fig. 8 $\theta=75^\circ$ map, in high-contrast colours against magma (cyan, hot pink, lime, sky blue), with the published MAX triangles and MIN squares. Drawing is clipped to the coloured oval.""",
      FD + "/make_migkas_fig8_coe_rings_overlay.py", "", [], "", [("migkas_fig8_theta75_with_coe_rings.png", FD + "/migkas_fig8_theta75_with_coe_rings.png")]))
    steps.append((r"Step 12 — Half-loops: CoE axis $+15^\circ$, and the $+45^\circ$ CCW roll pair",
      r"""**What it shows.** The same three panels, but with the axis tilted $+15^\circ$ (anti-CoE $(0^\circ,+15^\circ)$) and each loop kept only on its $\theta$-ray side, an open $180^\circ$ arc. First run: roll 0. Second run: roll $45^\circ$ counter-clockwise about the CoE axis (looking along $+D$).
**Maths.** Half-cut: keep ring points with $(\hat n\cdot\hat Y)\,\sin\theta_{\rm mid}\ge0$ (top half when $\sin\theta>0$, bottom when $\sin\theta<0$). A roll by $\rho$ rotates $\hat e_{1,2}$ about $\hat a$: $\hat e_1'=\cos\rho\,\hat e_1+\sin\rho\,\hat e_2$.""",
      FD + "/make_coe_rings_half_loops.py", "--argv --roll 0", [], "main(tilt_b=15.0, roll_deg=45.0)   # the +45 deg CCW roll version",
      [("coe_rings_half_loops_tilt15.png", FD + "/coe_rings_half_loops_tilt15.png"), ("coe_rings_half_loops_tilt15_rot45ccw.png", FD + "/coe_rings_half_loops_tilt15_rot45ccw.png")]))
    steps.append((r"Step 13 — Half-loops on Migkas Fig. 8 (roll 0 and $+45^\circ$ CCW)",
      r"""**What it shows.** The Step 12 half-arcs drawn on the Migkas Fig. 8 $\theta=75^\circ$ map, with the same registration and marker positions as Step 11. The Migkas markers stay fixed; no Dec $=0$ line is drawn.""",
      FD + "/make_migkas_fig8_half_loops_overlay.py", "--argv --roll 0", [], "main(tilt_b=15.0, roll_deg=45.0)",
      [("migkas_fig8_half_loops_tilt15.png", FD + "/migkas_fig8_half_loops_tilt15.png"), ("migkas_fig8_half_loops_tilt15_rot45ccw.png", FD + "/migkas_fig8_half_loops_tilt15_rot45ccw.png")]))
    steps.append((r"Step 14 — Full loops with an opacity gradient (roll 0 and $+45^\circ$ CCW)",
      r"""**What it shows.** Full loops instead of half-loops. Each loop is drawn as short segments whose opacity fades with ring azimuth $\psi$ from the $\theta$-ray point: darkest there, fully transparent opposite. Panel A is always drawn at roll 0; the roll applies to Panels B and C only.
**Maths.** $\alpha(\psi)=\alpha_0\,\tfrac12(1+\cos\psi)$, where $\psi$ is measured from the ring point furthest along $\mathrm{sign}(\sin\theta)\,\hat Y$ ($\hat Y$ rolled).""",
      FD + "/make_coe_rings_gradient_loops.py", "--argv --roll 0", [], "main(tilt_b=15.0, roll_deg=45.0)",
      [("coe_rings_gradient_loops_tilt15.png", FD + "/coe_rings_gradient_loops_tilt15.png"), ("coe_rings_gradient_loops_tilt15_rot45ccw.png", FD + "/coe_rings_gradient_loops_tilt15_rot45ccw.png")]))
    steps.append((r"Step 15 — Gradient loops on Migkas Fig. 8 (roll 0 and $+45^\circ$ CCW)",
      r"""**What it shows.** The Step 14 gradient loops on the Migkas Fig. 8 map, with the same registration, markers and gradient rule.""",
      FD + "/make_migkas_fig8_gradient_loops_overlay.py", "--argv --roll 0", [], "main(tilt_b=15.0, roll_deg=45.0)",
      [("migkas_fig8_gradient_loops_tilt15.png", FD + "/migkas_fig8_gradient_loops_tilt15.png"), ("migkas_fig8_gradient_loops_tilt15_rot45ccw.png", FD + "/migkas_fig8_gradient_loops_tilt15_rot45ccw.png")]))
    for sc, stem, title, extra in (
        ("make_coe_rings_migkas_pin_axis_matched.py", "loop_threepanel_matched", "Step 16a — Matched three-panel: 2 MAX + 2 MIN loops, $\\pm9^\\circ$ bands, axis $+15^\\circ$",
         r"""**What it shows.** The final three-panel layout: two MAX loops (two reds) and two MIN loops (two blues), each a clean transparent $\pm9^\circ$ band with a mid-line, and no orange loop. Axis tilted $+15^\circ$. Panel A is the flat model-frame angle circle with Galactic $\ell$ ticks; Panel B the Mollweide map; Panel C the 3D sphere, drawn **at the same diameter as Panel A** (the canvas was widened rather than anything shrunk)."""),
        ("make_coe_rings_migkas_pin_axis_matched_grad.py", "loop_threepanel_matched_grad", "Step 16b — Matched three-panel with the opacity gradient",
         r"""**What it shows.** Step 16a with each loop fading from its $\theta$-ray point, $\alpha\propto\tfrac12(1+\cos\psi)$, in Panels B and C."""),
        ("make_coe_rings_migkas_pin_axis_matched_grad_rot45.py", "loop_threepanel_matched_grad_rot45", "Step 16c — Matched three-panel, gradient, $45^\\circ$ roll",
         r"""**What it shows.** Step 16b with the loops rolled $45^\circ$ about the CoE axis. This is the $45^\circ$ roll version of the three-panel figure."""),
        ("make_coe_rings_migkas_pin_axis_matched_grad_rot_fit.py", "loop_threepanel_matched_grad_rotfit_rot81", "Step 16d — Matched three-panel, gradient, fitted $80.8^\\circ$ roll",
         r"""**What it shows.** Step 16b with the roll set to the fitted $80.8^\circ$: the roll that brings the loops' $\theta$-ray points closest to the Migkas markers."""),
        ("make_coe_rings_notilt.py", "loop_threepanel_notilt", "Step 17 — Untilted, unrolled three-panel (anti-CoE $(0^\\circ,0^\\circ)$, CoE $(180^\\circ,0^\\circ)$)",
         r"""**What it shows.** The matched three-panel with no tilt and no roll. The CoE axis lies in the Galactic plane, with anti-CoE at $(\ell,b)=(0^\circ,0^\circ)$ and CoE at $(180^\circ,0^\circ)$, so each loop is a small circle centred on the plane. Band centres (from the strip) are $\theta_{\max}=29.5^\circ,187.5^\circ$ and $\theta_{\min}=122.6^\circ,264.4^\circ$, each $\pm9^\circ$. Panel C is at the same diameter as Panel A.""")):
        steps.append((title, extra, FD + "/" + sc, "", [('"/workspace/cosmology/papers/spreadsheet_model/%s' % stem, '"%s/%s' % (P, stem))], "",
                      [(stem + ".png", SM + "/" + stem + ".png")]))
    steps.append((r"Step 18 — The four-loop overlay on Migkas Fig. 8, recoloured for contrast against magma (4 versions)",
      r"""**What it shows.** The matched loops of Step 16 drawn on the Migkas Fig. 8 $\theta=75^\circ$ map, recoloured to stand out against the magma colour map. Four versions from one script: matched; with gradient; gradient with $45^\circ$ roll; gradient with the fitted $80.8^\circ$ roll.""",
      FD + "/make_migkas_fig8_overlay_matched.py", "", [('SM_DIR = Path("/workspace/cosmology/papers/spreadsheet_model")', 'SM_DIR = Path("%s")' % P)], "",
      [("fig8_overlay_%s.png" % n, SM + "/fig8_overlay_%s.png" % n) for n in ("matched", "matched_grad", "matched_grad_rot45", "matched_grad_rot81")]))
    lrep = [('"../loop_', '"%s/loop_' % P), ('"../%s.png"', '"%s/%%s.png"' % P), ('"../%s_web.jpg"', '"%s/%%s_web.jpg"' % P)]
    steps.append((r"Step 19 — Loop-flow page (the four transformation steps, earlier 26 Sep version)",
      r"""**What it shows.** The step-by-step loop-flow figure: the $A(\ell)$ graph with its extrema bands, then the loops untilted, tilted, and with the gradient and roll, ending on the Migkas Fig. 8 overlay. The explanatory text comes from `stack_text.py`. Run from the `loops/` folder, whose `loop_fit.json` holds the fitted tilt and roll.""",
      LP + "/loop_flow2.py", "--cwd " + LP, [l for l in lrep if l[0] == '"../loop_'], "",
      [("loop_flow_thinshell.png", SM + "/loop_flow_thinshell.png")]))
    steps.append((r"Step 20 — Journal page of the loop transformations (26 Sep)",
      r"""**What it shows.** The journal-style page: the $A$ graph (row 0), then Transformation 1 (each extremum becomes a full loop of radius $\theta=\ell_{\rm ext}-\ell_0$ about anti-CoE, untilted), Transformation 2 (tilt), Transformation 3 (roll, with the gradient darkest at the $\theta$ ray), and the Fig. 8 overlay. Each caption gives the RMS miss.""",
      LP + "/loop_journal.py", "--cwd " + LP, [('"../loop_', '"%s/loop_' % P), ('Path("journal_rows")', 'Path("%s/journal_rows")' % P)], "",
      [("loop_journal_page.png", SM + "/loop_journal_page.png")]))
    steps.append((r"Step 21 — Loop stacks for the $45^\circ$ and $80.8^\circ$ rolls (26 Sep)",
      r"""**What it shows.** For each roll ($45^\circ$ and the fitted $80.8^\circ$): the $A$ graph, the three-panel loops in the Galactic-$\ell$ frame, and the Fig. 8 overlay, stacked on one page. This script loads the gradient module version it was written with (`coe_gradient_loop_geometry.cos_bak.py`).""",
      LP + "/loop_stack.py", "--cwd " + LP, lrep + [('Path("stack_rows")', 'Path("%s/stack_rows")' % P)], "",
      [("loop_stack_roll45.png", SM + "/loop_stack_roll45.png"), ("loop_stack_roll80p8.png", SM + "/loop_stack_roll80p8.png")]))
    steps.append((r"Step 22 — Stacked figures: 3-column and 4-column (stack_4col) versions",
      r"""**What it shows.** The summary stack: the model $A(\ell)$ graph (Newtonian $\Omega_s=0.049$ plus the relativistic curve) at the top, then four loop rows (untilted/no roll; $+15^\circ$ tilt; $+15^\circ$ tilt with gradient; $+15^\circ$ tilt with gradient and $+80.8^\circ$ roll), plus the Migkas Fig. 8 overlays (matched, gradient, $80.8^\circ$ roll). The 3-column version puts the overlays in a row at the bottom; the 4-column version puts each overlay beside its loop row at the original scale. This script composes existing images: the approved row and overlay images on disk, including three row images and two overlay images kept in the chat attachments folder.""",
      FD + "/make_stack_clean.py", "", [('P = "/workspace/cosmology/papers/spreadsheet_model/"', 'P = "/workspace/cosmology/papers/spreadsheet_model/"\nOUTP = "%s/"   # notebook: outputs redirected' % P),
                                        ('im.save(P + name + ".jpg"', 'im.save(OUTP + name + ".jpg"'), ('.save(P + name + "_web.jpg"', '.save(OUTP + name + "_web.jpg"')], "",
      [("stack_3col_clean.jpg", SM + "/stack_3col_clean.jpg"), ("stack_4col_clean.jpg", SM + "/stack_4col_clean.jpg")]))
    BUBD = SM + "/bubble"
    steps.append((r"Step 23 — Thin-shell article, Figure 1 (model panels) and Figure 2 ($A$(angle) fits incl. relativistic 0.885)",
      r"""**What it shows.** The two figures of `universe_bestfit_thinshell.pdf`, from its `.tex` source (`\includegraphics{universe_thin_drawing.png}` and `{universe_thin_Achart.png}`).
*Figure 1*: three scales at the relativistic best fit. (a) Hubble velocity vectors: in a Sedov–Taylor blast wave (Taylor 1950; Sedov 1959) $\vec V=k\,\vec r$, so $\vec V_P-\vec V_O=k\vec r$ is Hubble's law and the CoE drops out. (b) The Hubble bubble to scale in Gpc: ordinary-matter shell ($\Omega_s=0.049$) from $R_{\rm in}=18.1$ to $R_{\rm out}=49.7$ Gpc around an empty interior and the CoE core, the $D=32.9$ Gpc circle through Earth, and the $R=14.3$ Gpc vision sphere. (c) The cosmic shell with the core's Schwarzschild radius $r_s=2GM_c/c^2$.
*Figure 2*: Migkas points with error bars, $A_F$ and the three model curves, normalised to the sky mean: Newtonian thin shell ($R^2=0.909$), Newtonian $\Omega_s=0.049$ ($0.893$), and relativistic LTB ($0.885$). The lower panel shows percent differences from $A_F$.
The script also rewrites the table file `universe_thin.json`; that copy is redirected to `steps_v3/figs/` too, so the data file on disk is untouched. Run from the `bubble/` folder, which holds its fit inputs.""",
      BUBD + "/universe_thin.py", "--cwd " + BUBD,
      [('open("universe_thin.json", "w")', 'open("%s/universe_thin.json", "w")' % P), ('"../universe_thin_', '"%s/universe_thin_' % P), ('open("../" + f + ".png")', 'open("%s/" + f + ".png")' % P), ('im.save("../" + f + "_web.jpg"', 'im.save("%s/" + f + "_web.jpg"' % P)], "",
      [("universe_thin_drawing.png", SM + "/universe_thin_drawing.png"), ("universe_thin_Achart.png", SM + "/universe_thin_Achart.png")]))
    for title, text, script, args, repl, extra, pairs in steps:
        changes = "; ".join("`%s` → `%s`" % (a.strip('"')[:70], b.strip('"')[:70]) for a, b in repl) or "none. The script's `OUT_DIR` follows `__file__`, which is set to `steps_v3/figs/`"
        md("## %s\n\n%s\n\n*Source script:* `%s` (last modified %s), run verbatim. *Output-path edit:* %s. *Reference image:* %s." %
           (title, text, script, mt(script), changes, ", ".join("`%s`" % r.replace("/workspace/cosmology/", "") for _, r in pairs)))
        py(("%%%%runscript %s %s\n" % (os.path.basename(script), args)).replace("  ", " ") + src(script, repl))
        if extra:
            py("%%%%runscript %s_roll45 --argv --roll 45\n# second run of the same script with the +45 deg CCW roll: load its definitions (everything above its __main__ block), then call main()\n%s\nmain(tilt_b=15.0, roll_deg=45.0)" % (os.path.basename(script), "exec(open(%r).read().split('if __name__ == \"__main__\":')[0])" % script))
        py("\n".join("compare(FIGS / %r, %r)" % (n, r) for n, r in pairs))
    py(r'''
# Part B summary - match check for every figure
from IPython.display import Markdown
t = "| figure | reference | size ref | size new | mean abs diff (/255) | pixels differing >40 | verdict |\n|---|---|---|---|---|---|---|\n"
t += "\n".join(f"| {m['figure']} | {m['reference'].replace('/workspace/cosmology/', '')} | {m['size_ref'][0]}x{m['size_ref'][1]} | {m['size_new'][0]}x{m['size_new'][1]} | {m['mean_abs_diff']:.3f} | {100*m['frac_px_diff']:.3f}% | {m['verdict']} |" for m in MATCH)
display(Markdown(t))
''')
