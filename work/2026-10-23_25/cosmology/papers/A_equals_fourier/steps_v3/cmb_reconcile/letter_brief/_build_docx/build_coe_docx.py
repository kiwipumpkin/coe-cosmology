#!/usr/bin/env python3
"""Build CoE_letter_for_Word.docx: full draft text + all letter.tex figures + equation PNGs."""
from __future__ import annotations

import io
import os
import re
import subprocess
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from docx import Document
from docx.shared import Inches, Pt, RGBColor, Twips
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

ROOT = Path("/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief")
FIGS = ROOT / "figs"
EQDIR = ROOT / "_eq_pngs"
OUT = ROOT / "CoE_letter_for_Word.docx"
EQDIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Equation renderer (matplotlib mathtext; pdflatex fallback)
# ---------------------------------------------------------------------------

def render_eq_mpl(latex: str, out_path: Path, fontsize: float = 18, dpi: int = 200) -> Path:
    """Render a display equation with matplotlib mathtext."""
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.patch.set_facecolor("white")
    # Wrap in $...$ if not already
    s = latex.strip()
    if not (s.startswith("$") or s.startswith(r"\(") or s.startswith(r"\[")):
        s = f"${s}$"
    # Convert \[ \] / \( \) to $
    s = s.replace(r"\[", "$").replace(r"\]", "$").replace(r"\(", "$").replace(r"\)", "$")
    try:
        fig.text(0.5, 0.5, s, fontsize=fontsize, ha="center", va="center")
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight", pad_inches=0.15,
                    facecolor="white", edgecolor="none")
        plt.close(fig)
        return out_path
    except Exception as e:
        plt.close(fig)
        raise e


def render_eq_pdflatex(latex: str, out_path: Path, fontsize: str = "12pt") -> Path:
    """Render via pdflatex + pdftoppm + whitespace trim (no preview.sty needed)."""
    from PIL import Image, ImageChops
    body = latex.strip()
    if body.startswith("$"):
        body = body.strip("$")
    # Wide page for long displays
    width = "30cm" if len(body) > 80 else "16cm"
    tex = rf"""\documentclass[12pt]{{article}}
\usepackage{{amsmath,amssymb}}
\usepackage[paperwidth={width},paperheight=3.5cm,margin=2mm]{{geometry}}
\pagestyle{{empty}}
\begin{{document}}
\[
{body}
\]
\end{{document}}
"""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "eq.tex").write_text(tex)
        r = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "eq.tex"],
            cwd=td, capture_output=True, text=True,
        )
        if r.returncode != 0 or not (td / "eq.pdf").exists():
            raise RuntimeError(f"pdflatex failed:\n{r.stdout[-1500:]}")
        subprocess.run(
            ["pdftoppm", "-png", "-r", "240", "-singlefile", str(td / "eq.pdf"), str(td / "eq")],
            check=True, capture_output=True,
        )
        im = Image.open(td / "eq.png").convert("RGB")
        bg = Image.new("RGB", im.size, (255, 255, 255))
        diff = ImageChops.difference(im, bg)
        bbox = diff.getbbox()
        if bbox:
            pad = 10
            x0, y0, x1, y1 = bbox
            x0 = max(0, x0 - pad); y0 = max(0, y0 - pad)
            x1 = min(im.width, x1 + pad); y1 = min(im.height, y1 + pad)
            im = im.crop((x0, y0, x1, y1))
        im.save(out_path)
    return out_path


def render_eq(name: str, latex: str, display: bool = True) -> Path:
    out = EQDIR / f"{name}.png"
    # Prefer pdflatex for quality / complex macros
    try:
        return render_eq_pdflatex(latex, out)
    except Exception as e1:
        try:
            # matplotlib fallback: simplify common macros
            s = latex
            s = s.replace(r"\mathrm", r"\mathrm")
            s = s.replace(r"\vec{", r"\mathbf{")
            return render_eq_mpl(s, out, fontsize=16 if display else 14)
        except Exception as e2:
            print(f"WARN: eq {name} failed: {e1} | {e2}")
            # Last resort: write placeholder text image
            fig = plt.figure(figsize=(6, 0.6))
            fig.text(0.5, 0.5, latex[:120], fontsize=10, ha="center", va="center",
                     family="monospace", wrap=True)
            fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            return out


# ---------------------------------------------------------------------------
# Document helpers
# ---------------------------------------------------------------------------

def set_run_font(run, size=11, bold=False, italic=False, name="Calibri"):
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.name = name
    r = run._element
    rPr = r.get_or_add_rPr()
    rFonts = rPr.find(qn("w:rFonts"))
    if rFonts is None:
        rFonts = OxmlElement("w:rFonts")
        rPr.insert(0, rFonts)
    rFonts.set(qn("w:ascii"), name)
    rFonts.set(qn("w:hAnsi"), name)


def add_para(doc, text, *, style=None, size=11, bold=False, italic=False,
             space_after=8, space_before=0, align=None):
    p = doc.add_paragraph(style=style) if style else doc.add_paragraph()
    if align is not None:
        p.alignment = align
    pf = p.paragraph_format
    pf.space_after = Pt(space_after)
    pf.space_before = Pt(space_before)
    pf.line_spacing = 1.15
    if text:
        run = p.add_run(text)
        set_run_font(run, size=size, bold=bold, italic=italic)
    return p


def add_heading_custom(doc, text, level=1):
    sizes = {0: 18, 1: 14, 2: 12, 3: 11}
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(14 if level <= 1 else 10)
    p.paragraph_format.space_after = Pt(6)
    run = p.add_run(text)
    set_run_font(run, size=sizes.get(level, 11), bold=True)
    return p


def add_equation_image(doc, png_path: Path, width_in: float = 5.5):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(6)
    run = p.add_run()
    # Cap width by image aspect — don't blow up tiny eqs
    from PIL import Image
    with Image.open(png_path) as im:
        w_px, h_px = im.size
    # target ~ 0.35–0.5 in height for single-line, wider for multi
    natural_w = w_px / 220.0  # inches at 220 dpi
    use_w = min(width_in, max(1.2, natural_w * 1.05))
    run.add_picture(str(png_path), width=Inches(use_w))
    return p


def add_figure(doc, img_path: Path, caption: str, width_in: float = 6.3, fig_label: str | None = None):
    if not img_path.exists():
        add_para(doc, f"[MISSING FIGURE: {img_path.name}]", italic=True, size=10)
        return
    from PIL import Image
    with Image.open(img_path) as im:
        w_px, h_px = im.size
    aspect = h_px / max(w_px, 1)
    # Cap height roughly ~ 8 inches for very tall stacks
    max_h = 8.0
    use_w = width_in
    if aspect * use_w > max_h:
        use_w = max_h / aspect

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(10)
    p.paragraph_format.space_after = Pt(2)
    run = p.add_run()
    run.add_picture(str(img_path), width=Inches(use_w))

    cap = caption
    if fig_label:
        cap = f"{fig_label}. {caption}" if not caption.startswith(fig_label) else caption
    cp = doc.add_paragraph()
    cp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cp.paragraph_format.space_before = Pt(2)
    cp.paragraph_format.space_after = Pt(10)
    run = cp.add_run(cap)
    set_run_font(run, size=9, italic=True)


def clean_inline_math_to_unicode(text: str) -> str:
    """Replace simple $...$ with unicode-ish plain; leave complex for eq images."""
    # Already handled via equation images for displays; for remaining inline,
    # strip $ and keep latex-ish readable text, with common replacements.
    repl = [
        (r"\\sim", "∼"),
        (r"\\approx", "≈"),
        (r"\\lesssim", "≲"),
        (r"\\gtrsim", "≳"),
        (r"\\times", "×"),
        (r"\\cdot", "·"),
        (r"\\ell", "ℓ"),
        (r"\\theta", "θ"),
        (r"\\phi", "φ"),
        (r"\\Omega", "Ω"),
        (r"\\omega", "ω"),
        (r"\\beta", "β"),
        (r"\\nu", "ν"),
        (r"\\mu", "μ"),
        (r"\\Delta", "Δ"),
        (r"\\delta", "δ"),
        (r"\\circ", "°"),
        (r"\\mathrm\{([^}]*)\}", r"\1"),
        (r"\\mathbf\{([^}]*)\}", r"\1"),
        (r"\\textit\{([^}]*)\}", r"\1"),
        (r"\\emph\{([^}]*)\}", r"\1"),
        (r"\\vec\{([^}]*)\}", r"\1⃗"),
        (r"_\{([^}]*)\}", r"_\1"),
        (r"\^\{([^}]*)\}", r"^\1"),
        (r"\\,", " "),
        (r"\\;", " "),
        (r"\\ ", " "),
        (r"\\%", "%"),
        (r"\\&", "&"),
    ]
    def sub_math(m):
        s = m.group(1)
        for a, b in repl:
            s = re.sub(a, b, s)
        s = s.replace("{", "").replace("}", "")
        s = s.replace("\\", "")
        return s

    text = re.sub(r"\$([^$]+)\$", sub_math, text)
    text = re.sub(r"\\\((.+?)\\\)", sub_math, text)
    return text


def add_body(doc, text, **kw):
    text = clean_inline_math_to_unicode(text)
    # Remove draft editorial notes that look like TODOs at end of sentences
    return add_para(doc, text, **kw)


# ---------------------------------------------------------------------------
# Render all equations first
# ---------------------------------------------------------------------------

EQUATIONS = {
    "vh": r"\vec{V}_H = \vec{V}_g - \vec{V}_O",
    "einstein": r"G_{\mu\nu} = 8\pi G\, T_{\mu\nu}",
    "ltb_metric": (
        r"\mathrm{d}s^{2} = -\mathrm{d}t^{2}"
        r" + \frac{(R')^{2}}{1+2E(r)}\,\mathrm{d}r^{2}"
        r" + R(t,r)^{2}\,\mathrm{d}\Omega^{2}"
    ),
    "friedmann": r"\dot{R}^{2} = 2E(r) + \frac{2GM(r)}{R}",
    "mass": r"M(r) = \frac{4\pi}{3}\, R^{3}\, \nu\, \beta(r)",
    "ha": (
        r"H_{a}^{2} \propto \frac{\nu_{m}\beta_{m}(r)}{a^{3}}"
        r" + \frac{\nu_{r}\beta_{r}(r)}{a^{4}}"
        r" - \frac{k(r)}{a^{2}}"
    ),
    "newt_geom": (
        r"D \approx 49.17\,\mathrm{Gpc},\quad"
        r"R_{\mathrm{in}} \approx 34.48\,\mathrm{Gpc},\quad"
        r"R_{\mathrm{out}} \approx 70.27\,\mathrm{Gpc},\quad"
        r"\Omega_{s} = 0.049,\quad"
        r"\ell_{0} = 13^{\circ}"
    ),
    "pref_rel": (
        r"D \approx 32.9\,\mathrm{Gpc},\quad"
        r"R_{\mathrm{in}} \approx 18.1\,\mathrm{Gpc},\quad"
        r"R_{\mathrm{out}} \approx 49.7\,\mathrm{Gpc}"
    ),
    "shell_rule": r"D \gtrsim 14.7\,\mathrm{Gpc}\quad (R \approx 14.3\,\mathrm{Gpc})",
    "vhb": r"V_{\mathrm{HB}} = 0 \quad\text{until}\quad D_{\mathrm{recomb}}",
}

print("Rendering equations...")
for name, latex in EQUATIONS.items():
    p = render_eq(name, latex)
    print(f"  {name}: {p.name} ({p.stat().st_size} bytes)")

# Also need Pillow
try:
    from PIL import Image
except ImportError:
    subprocess.check_call(["/tmp/docx_venv/bin/pip", "install", "pillow", "-q"])
    from PIL import Image

# ---------------------------------------------------------------------------
# Build document
# ---------------------------------------------------------------------------

doc = Document()

# Narrow margins a bit
for section in doc.sections:
    section.top_margin = Inches(0.75)
    section.bottom_margin = Inches(0.75)
    section.left_margin = Inches(0.85)
    section.right_margin = Inches(0.85)

# Title block
add_heading_custom(doc, "Centre-of-Expansion shell cosmology: from Migkas anisotropy to a flat CMB", level=0)
add_para(doc, "Chris Alameddin", size=12, align=WD_ALIGN_PARAGRAPH.CENTER, space_after=2)
add_para(doc, "1 October 2026", size=11, italic=True, align=WD_ALIGN_PARAGRAPH.CENTER, space_after=12)

# Abstract
add_heading_custom(doc, "Abstract", level=1)
add_body(
    doc,
    "Directional Hubble measurements from Migkas et al. fix an off-centre Sedov–Taylor / "
    "Lemaître–Tolman–Bondi (LTB) shell geometry about a Centre of Expansion (CoE). The same "
    "geometry is then subjected to an independent consistency check against the cosmic microwave "
    "background (CMB). This letter treats both the Migkas sky-dependent cluster Hubble rate and "
    "the Planck-level CMB isotropy of order $10^{-5}$ as correct empirical facts. The cosmological "
    "principle alone does not account for the observed Hubble anisotropy; the scientific goal is "
    "coexistence of that anisotropy with a flat last-scattering sky, not replacement of either "
    "dataset. We define the CoE as the geometric origin of the radial shell structure; the observer "
    "(US) sits off-centre inside a thick ordinary-matter shell—the donut-universe configuration; "
    "and we ask whether the gravitational multipole implied by that placement can be cancelled by "
    "a kill function: an electromagnetic temperature multipole ($\\Delta T/T$ of kinematic–Doppler "
    "character) opposite in sign to the gravitational last-scattering multipole and sharing its "
    "angular shape. In this geometry the apparent accelerated expansion of the late universe is an "
    "illusion—no dark energy and no actual outward acceleration of the shell are required—arising "
    "instead as a geometric, off-centre sampling artifact. A contracting universe would likewise "
    "present, on average, the appearance of accelerated expansion, but with a different angular "
    "signature that can distinguish the two cases.",
)

# 1. Introduction
add_heading_custom(doc, "1. Introduction", level=1)
add_body(
    doc,
    "Migkas et al. report a directional dependence of the cluster-derived Hubble rate across the sky. "
    "Standard cosmology assumes large-scale isotropy, and the CMB is isotropic to roughly one part in "
    "$10^5$ according to Planck. Those two statements appear to pull in opposite directions: a "
    "sky-dependent low-redshift expansion rate sits uneasily next to a nearly featureless last-scattering "
    "surface if one insists that the cosmological principle alone must govern both.",
)
add_body(
    doc,
    "This letter refuses that forced choice. We treat the standard CMB analysis and the CoE shell "
    "geometry as both correct within their domains. The Migkas anisotropy is interpreted as free-fall "
    "relative to an off-centre CoE, fitted with a thick shell whose Newtonian and relativistic LTB "
    "realisations are compared in Table 1 and Figures 2–4. The LTB origin is the CoE, not the observer: "
    "we sit at radial coordinate $r = D$, and the observed anisotropy is an off-centre sampling of the "
    "radial mass-multiplier profile $\\beta(r)$, not a direction-dependent scale factor $a(t,\\theta)$. "
    "The CMB is then used only to stress-test that fixed geometry. If the same off-centre placement "
    "produces a gravitational multipole on the last-scattering sky, we ask what electromagnetic "
    "temperature pattern—the kill function—would cancel it, and whether that pattern morphologically "
    "matches a pure kinematic Doppler contribution from emitting surfaces at the same radii. The open "
    "problem is coexistence—preserving the Migkas $A(\\ell)$ fit while keeping CMB multipoles of order "
    "$10^{-5}$—not discarding either dataset.",
)
add_body(
    doc,
    "A further consequence of the same geometry is that the late-time appearance of accelerated expansion "
    "need not signal a dark-energy component or an actual outward acceleration of the material shell. "
    "Off-centre sampling of the expansion field can produce that appearance as an artefact. A contracting "
    "CoE geometry would also show, on average, the appearance of accelerated expansion, but with a "
    "different angular pattern on the sky—a natural observational distinguisher for future work.",
)

# 2. Methods
add_heading_custom(doc, "2. Methods and model", level=1)
add_heading_custom(doc, "2.1 Geometry and kinematic picture", level=2)
add_body(doc, "Local Hubble velocity is free-fall toward the CoE minus our own motion,")
add_equation_image(doc, EQDIR / "vh.png", width_in=3.5)
add_body(
    doc,
    "The observer (US) sits inside a thick shell at comoving distance $D$ from the CoE. Together, the "
    "CoE and the thick shell with US inside constitute the donut-universe configuration. The vision "
    "radius used for today’s-field diagnostics is $R \\approx 14.3\\,\\mathrm{Gpc}$. The shell rule with "
    "that vision radius requires $D \\gtrsim 14.7\\,\\mathrm{Gpc}$, so that the observer remains mid-shell "
    "between the inner and outer edges rather than sitting outside the shell or inside an empty cavity "
    "in a way that would violate the fitted anisotropy.",
)
add_body(doc, "Figure 1 shows the geometry: the CoE, the observer (US) inside the thick shell, and the vision sphere of radius $R$.")
add_figure(
    doc, FIGS / "att14_universe_thin_drawing.png",
    "Geometry: CoE, observer (US) inside a thick shell, and the vision sphere of radius R.",
    width_in=5.5, fig_label="Figure 1",
)

# Preferred relativistic fit callout
add_body(doc, "Preferred relativistic fit:")
add_equation_image(doc, EQDIR / "pref_rel.png", width_in=5.8)

add_heading_custom(doc, "2.2 Fitting procedure", level=2)
add_body(doc, "The analysis proceeds in five steps.")
for i, step in enumerate([
    "Compare Migkas $A/A_{\\mathrm{all}}$ versus Galactic longitude $\\ell$ (their Fig. 6; our Fig. 2), and fit by redshift slice (our Fig. 3), with orientation against Migkas Fig. 8.",
    "Fit a CoE plus thick shell in both Newtonian gravity and relativistic LTB.",
    "Orient the model on the Galactic sky: untilted, then $+15^\\circ$ tilt, then a radial gradient, then an $+80.8^\\circ$ roll about the CoE axis, matching Migkas Fig. 8.",
    "Evolve the donut-universe picture—shells expanding outward faster than light in a coordinate sense, while the sphere of events visible to us grows at $c$, with recombination as a spherical surface of fixed physical radius around the CoE.",
    "Stress-test the CMB: compute the gravitational multipole on the last-scattering sky, form the kill function (the electromagnetic temperature multipole opposite in sign to that gravitational multipole), and compare exploratory two-emitter and multipole-partner constructions.",
], 1):
    add_body(doc, f"{i}. {step}", space_after=4)

# Table 1
add_heading_custom(doc, "Table 1. Newtonian and relativistic shell+core fits to Migkas A(ℓ)", level=3)
table = doc.add_table(rows=5, cols=4)
table.style = "Table Grid"
headers = ["", "Newt. Ω_s free", "Newt. Ω_s=0.049", "Rel. LTB Ω_s=0.049"]
data = [
    ["D (Gpc)", "93.7", "49.2", "32.9"],
    ["Ω_s", "1.01", "0.049", "0.049"],
    ["ℓ_0 (deg)", "12.8", "13.0", "13.5"],
    ["R² vs Migkas", "0.909", "0.893", "0.885"],
    ["R² vs A_F", "0.987", "0.977", "0.963"],
]
# Actually rows=6 with header
# rebuild
tbl = doc.tables[-1]
# Clear and redo properly
doc._body._body.remove(tbl._tbl)

table = doc.add_table(rows=6, cols=4)
table.style = "Table Grid"
row0 = ["", "Newt. Ωₛ free", "Newt. Ωₛ=0.049", "Rel. LTB Ωₛ=0.049"]
rows = [
    row0,
    ["D (Gpc)", "93.7", "49.2", "32.9"],
    ["Ωₛ", "1.01", "0.049", "0.049"],
    ["ℓ₀ (deg)", "12.8", "13.0", "13.5"],
    ["R² vs Migkas", "0.909", "0.893", "0.885"],
    ["R² vs A_F", "0.987", "0.977", "0.963"],
]
for i, rowdata in enumerate(rows):
    for j, val in enumerate(rowdata):
        cell = table.rows[i].cells[j]
        cell.text = ""
        p = cell.paragraphs[0]
        run = p.add_run(val)
        set_run_font(run, size=9, bold=(i == 0 or j == 0))
add_para(doc, "", space_after=4)

add_body(
    doc,
    "Figure 2 shows Migkas $A/A_{\\mathrm{all}}$ versus Galactic longitude together with the Newtonian "
    "and relativistic shell-plus-core model curves. Figure 3 shows the same $A(\\ell)$ pattern in three "
    "redshift slices ($0.03$–$0.07$, $z > 0.07$, and all-$z$). Fit parameters are held fixed; the model "
    "is recomputed at each slice’s median sampling distance, so that the $A(\\ell)$ maxima and minima "
    "move with that distance. Figure 3 also displays the $D$–core-mass degeneracy. Figure 4 shows the "
    "full Migkas $A/A_{\\mathrm{all}}$ profile and residual panels for the Newtonian and relativistic fits.",
)

add_figure(
    doc, FIGS / "att14_Achart.png",
    "Migkas A/A_all versus Galactic longitude, with Newtonian and relativistic shell+core model curves.",
    width_in=6.3, fig_label="Figure 2",
)
add_figure(
    doc, FIGS / "att14_migkas_slices.png",
    "Migkas A(ℓ) in three redshift slices (0.03–0.07, z>0.07, and all-z), with the CoE model "
    "recomputed at each slice’s median distance (fit parameters fixed), plus the D–core-mass degeneracy.",
    width_in=6.3, fig_label="Figure 3",
)

add_heading_custom(doc, "2.3 Einstein and Friedmann equations for the CoE+shell", level=2)
add_body(
    doc,
    "The LTB metric is written about the CoE, not about us. The observer (US) sits at comoving distance "
    "$D$ from that centre. Sky anisotropy is therefore an off-centre sampling of radial structure "
    "$\\beta(r)$, not a direction-dependent scale factor $a(t,\\theta)$. Einstein’s equation is",
)
add_equation_image(doc, EQDIR / "einstein.png", width_in=3.2)
add_body(doc, "For concentric dust shells about the CoE the Lemaître–Tolman–Bondi line element is")
add_equation_image(doc, EQDIR / "ltb_metric.png", width_in=5.8)
add_body(
    doc,
    "with areal radius $R(t,r)$, $R' \\equiv \\partial_r R$, and $\\dot{R} \\equiv \\partial_t R$. "
    "Each shell obeys the Friedmann-like constraint",
)
add_equation_image(doc, EQDIR / "friedmann.png", width_in=3.8)
add_body(
    doc,
    "The interior mass uses a background density scale $\\nu$ and a radial mass multiplier $\\beta(r)$ "
    "(not a scale factor),",
)
add_equation_image(doc, EQDIR / "mass.png", width_in=3.8)
add_body(
    doc,
    "so that $\\beta = 1$ recovers critical-density normalisation. In the CoE+shell geometry $\\beta(r)$ "
    "is small in the cavity, rises through the thick shell, and returns toward the exterior background "
    "beyond $R_{\\mathrm{out}}$.",
)
add_body(doc, "Including radiation on the same radial scaffolding,")
add_equation_image(doc, EQDIR / "ha.png", width_in=5.5)
add_body(
    doc,
    "with $H_a = \\dot{a}/a$ and $a(t,r)$ the local expansion factor of each shell. Early times are "
    "radiation-dominated ($a \\sim t^{1/2}$); late times are dust-dominated ($a \\sim t^{2/3}$). "
    "Curvature $k(r)$ follows from $E(r)$ and introduces no angular dependence into the metric coefficients.",
)
add_body(doc, "For the Newtonian $\\Omega_s = 0.049$ CoE+shell fit in Table 1, the solved geometry matching Migkas $A(\\ell)$ is")
add_equation_image(doc, EQDIR / "newt_geom.png", width_in=6.2)
add_body(doc, "The shell rule with vision radius:")
add_equation_image(doc, EQDIR / "shell_rule.png", width_in=4.5)
add_body(
    doc,
    "The fitted $D$ lies above that bound with the observer mid-shell between $R_{\\mathrm{in}}$ and "
    "$R_{\\mathrm{out}}$.",
)

add_figure(
    doc, FIGS / "att14_amp_full.png",
    "Full Migkas A/A_all profile and residual panel for the Newtonian and relativistic fits.",
    width_in=6.3, fig_label="Figure 4",
)

add_heading_custom(doc, "2.4 Sky orientation", level=2)
add_body(
    doc,
    "The model must be placed on the Galactic sky. Starting from an untilted, no-roll configuration, "
    "we apply a $+15^\\circ$ tilt, then a radial gradient, then an $+80.8^\\circ$ roll about the CoE axis. "
    "Figure 5 shows that sequence: each row presents the model frame, the Mollweide sky map, and a "
    "three-dimensional view. Figure 6 overlays the predicted MAX/MIN bands on Migkas Fig. 8 for tilt "
    "alone, tilt with gradient, and tilt with gradient and roll.",
)
add_figure(
    doc, FIGS / "att14_four_row.jpg",
    "Sky orientation sequence: untilted; +15° tilt; tilt with gradient; tilt, gradient, and +80.8° roll. "
    "Each row shows the model frame, the Mollweide sky map, and a 3D view.",
    width_in=6.3, fig_label="Figure 5",
)
add_body(
    doc,
    "The untilted-to-tilted-to-rolled sequence places the CoE axis on the Galactic sky. With a $+15^\\circ$ "
    "tilt alone, the model MAX/MIN bands already track Migkas Fig. 8. Adding a radial gradient sharpens "
    "contrast; the $+80.8^\\circ$ roll aligns the MAX/MIN pattern. Fit parameters stay fixed; only sky "
    "placement varies.",
)
add_figure(
    doc, FIGS / "att14_overlay_stack.png",
    "Model MAX/MIN bands overlaid on Migkas Fig. 8 for +15° tilt; tilt with gradient; and tilt, "
    "gradient, and +80.8° roll.",
    width_in=5.5, fig_label="Figure 6",
)

add_heading_custom(doc, "2.5 Time evolution, donut universe, and CMB stress test", level=2)
add_body(
    doc,
    "In the donut-universe picture the CoE is fixed while thick shells expand outward; the sphere of "
    "events visible to us grows at $c$. Recombination is a spherical surface of fixed physical radius "
    "around the CoE. The observer remains mid-shell between the inner and outer edges (S1 and S2). "
    "Critically, the Hubble-bubble volume about US satisfies",
)
add_equation_image(doc, EQDIR / "vhb.png", width_in=4.5)
add_body(
    doc,
    "(the CoE distance at which vision turns on). After that crossing, the visible sphere grows at $c$. "
    "Scale discussion may be framed in terms of the time-dependent volumes $V_{\\mathrm{HB}}(t)$, "
    "$V_{\\mathrm{S1}}(t)$, and $V_{\\mathrm{S2}}(t)$ and the corresponding volume fractions of the "
    "Hubble bubble relative to the shell annulus and the outer sphere. Figure 7 shows the time sequence "
    "(frames A–D): shells cross the recombination surface while the visible sphere grows more slowly, "
    "with S2 filled in every frame.",
)
add_figure(
    doc, FIGS / "att14_storyboard.png",
    "Donut-universe time sequence (A–D): fixed recombination surface; shells expand through it; "
    "S2 is filled in every frame; the observer sits mid-shell.",
    width_in=4.8, fig_label="Figure 7",
)
add_body(
    doc,
    "Figure 8 shows the two-emitter baseline: if opposite plates (or shells) share the observer’s bulk "
    "velocity, the kinematic sky is flat before gravitational multipoles are included. The same "
    "off-centre geometry still produces a gravitational multipole on the last-scattering sky. The kill "
    "function is the electromagnetic temperature multipole—$\\Delta T/T$ of kinematic–Doppler photon "
    "character—required to cancel that gravitational multipole: opposite in sign, same angular shape. "
    "Figure 9 compares the unit-normalized kill function with exploratory plate, shell, and shear "
    "residuals at three distances $D$. Figure 10 shows residuals after adding a multipole partner: "
    "after a dipole-plus-$P_2$ cancellation, and after an $\\ell \\le 4$ partner.",
)
add_figure(
    doc, FIGS / "att14_two_emitters.png",
    "Two-emitter baseline. If opposite plates (or shells) share the observer’s bulk velocity, the "
    "kinematic sky is flat before gravitational multipoles are included.",
    width_in=4.5, fig_label="Figure 8",
)
add_figure(
    doc, FIGS / "att14_kill.jpg",
    "Unit-normalized kill function compared with exploratory plate, shell, and shear residuals at "
    "three distances D.",
    width_in=6.3, fig_label="Figure 9",
)
add_figure(
    doc, FIGS / "att14_cancel_partner.png",
    "CoE Doppler term plus a multipole partner: residuals after a dipole+P₂ cancellation and after "
    "an ℓ≤4 partner.",
    width_in=6.3, fig_label="Figure 10",
)

# 3. Results
add_heading_custom(doc, "3. Results", level=1)
add_body(
    doc,
    "Table 1 summarises the Newtonian and relativistic shell-plus-core fits to Migkas $A(\\ell)$:",
)
add_body(
    doc,
    "Newtonian, $\\Omega_s$ free: $D \\approx 93.7\\,\\mathrm{Gpc}$, $\\Omega_s \\approx 1.01$, "
    "$\\ell_0 \\approx 12.8^\\circ$, $R^2$ versus Migkas $0.909$, $R^2$ versus $A_F$ $0.987$.",
)
add_body(
    doc,
    "Newtonian, $\\Omega_s = 0.049$: $D \\approx 49.2\\,\\mathrm{Gpc}$, $\\Omega_s = 0.049$, "
    "$\\ell_0 \\approx 13.0^\\circ$, $R^2$ versus Migkas $0.893$, $R^2$ versus $A_F$ $0.977$.",
)
add_body(
    doc,
    "Relativistic LTB, $\\Omega_s = 0.049$: $D \\approx 32.9\\,\\mathrm{Gpc}$, $\\Omega_s = 0.049$, "
    "$\\ell_0 \\approx 13.5^\\circ$, $R^2$ versus Migkas $0.885$, $R^2$ versus $A_F$ $0.963$.",
)
add_body(doc, "For the Newtonian $\\Omega_s = 0.049$ CoE+shell fit, the solved geometry matching Migkas $A(\\ell)$ is")
add_equation_image(doc, EQDIR / "newt_geom.png", width_in=6.2)
add_body(
    doc,
    "The preferred relativistic headline quoted in the geometry section is",
)
add_equation_image(doc, EQDIR / "pref_rel.png", width_in=5.8)
add_body(
    doc,
    "In both cases the fitted $D$ lies above the shell-rule bound $D \\gtrsim 14.7\\,\\mathrm{Gpc}$ "
    "(with vision radius $R \\approx 14.3\\,\\mathrm{Gpc}$), and the observer sits mid-shell between "
    "$R_{\\mathrm{in}}$ and $R_{\\mathrm{out}}$.",
)
add_body(
    doc,
    "The redshift-slice tests (Fig. 3) show that the CoE model, recomputed at each slice’s median "
    "distance with fit parameters held fixed, tracks the Migkas $A(\\ell)$ pattern across the low-$z$, "
    "high-$z$, and all-$z$ samples; the loci of the $A(\\ell)$ maxima and minima therefore move with "
    "sampling distance. The $D$–core-mass degeneracy is made explicit. The full amplitude profile and "
    "residuals (Fig. 4) confirm that both the Newtonian $\\Omega_s = 0.049$ and relativistic LTB curves "
    "follow the Migkas longitude dependence without being tuned to the CMB.",
)
add_body(
    doc,
    "On the sky, the untilted-to-tilted-to-rolled sequence (Fig. 5) places the model MAX/MIN loci into "
    "agreement with Migkas Fig. 8 once the $+15^\\circ$ tilt, radial gradient, and $+80.8^\\circ$ roll are "
    "applied (Fig. 6).",
)
add_body(
    doc,
    "For the CMB stress test, exploratory two-emitter and multipole-partner constructions follow the "
    "kill-function morphology at large $D$ and remove most of the demeaned power through $\\ell \\le 4$ "
    "(Figs. 8–10), but they do not yet reach the CMB bound of order $10^{-5}$. The morphological match "
    "between the gravitational kill and the electromagnetic Doppler pattern is already present; closing "
    "the residual to $10^{-5}$ on the full light cone remains the open calculation.",
)

# 4. Discussion
add_heading_custom(doc, "4. Discussion", level=1)
add_body(
    doc,
    "The Migkas anisotropy is interpreted here as a gravitational effect: free-fall relative to an "
    "off-centre CoE. On the last-scattering sky that same gravity produces a temperature multipole. "
    "Separately, a pure electromagnetic Doppler pattern—from radial motion of emitting surfaces—has "
    "nearly the same angular shape. The kill function is precisely that electromagnetic temperature "
    "multipole written opposite in sign to the gravitational residual. The result of the present stress "
    "test is that this gravitational kill tracks the electromagnetic Doppler morphology across the "
    "distances scanned. That morphological agreement is unexpected if gravity and the photon bath are "
    "unrelated, and it suggests a geometric route to CMB flatness that does not discard the Migkas fit. "
    "The open light-cone problem is whether that morphological match can be driven to the $10^{-5}$ "
    "amplitude required by Planck once gravitational and kinematic contributions are integrated together "
    "along null geodesics.",
)
add_body(
    doc,
    "The apparent accelerated expansion of the late universe is, in this geometry, an illusion. No "
    "dark-energy component and no actual outward acceleration of the material shell are required. What "
    "looks like acceleration is a geometric artefact of off-centre sampling of the expansion field by "
    "an observer at distance $D$ from the CoE. A contracting universe under analogous off-centre and "
    "light-cone sampling would also present, on average, the appearance of accelerated expansion, but "
    "with a different angular signature on the sky. That difference is a natural distinguisher for "
    "future comparison against directional Hubble and supernova data.",
)
add_body(
    doc,
    "In the donut-universe time evolution, $V_{\\mathrm{HB}} = 0$ until US passes $D_{\\mathrm{recomb}}$; "
    "afterward the visible sphere grows at $c$ while the shells continue to expand faster than light in "
    "the coordinate sense used here. Volume fractions $V_{\\mathrm{HB}}/V_{\\mathrm{ann}}$, "
    "$V_{\\mathrm{HB}}/V_{\\mathrm{S1}}$, and $V_{\\mathrm{S2}}/V_{\\mathrm{S1}}$ provide a compact "
    "language for tracking how the Hubble bubble fills the shell annulus relative to the outer and "
    "inner spheres as the geometry evolves.",
)
add_body(
    doc,
    "The stance remains coexistence. The cosmological principle does not explain the directional Hubble "
    "measurements; those measurements fix the CoE+shell. The Planck CMB isotropy remains a hard external "
    "constraint. Nothing in the present letter claims that the $10^{-5}$ bound has been met. The open "
    "problem is whether a full time-dependent light-cone calculation—evolving $D(t)$, $R_{\\mathrm{in}}(t)$, "
    "and $R_{\\mathrm{out}}(t)$ from the present back through recombination, integrating null geodesics, "
    "and combining gravitational and kinematic contributions—can keep $\\ell \\ge 2$ below $\\sim 10^{-5}$ "
    "while the low-redshift $A(\\ell)$ pattern remains Migkas-compatible. If a residual remains, the next "
    "step is to identify the smallest additional degree of freedom (timing of recombination, "
    "shell-thickness history, or a mild partner multipole) that removes it without reintroducing a "
    "fixed-curve conflict across redshift slices.",
)
add_body(
    doc,
    "Appendix material supports the dynamical and observational reading of the same geometry without "
    "changing the fit. Figure A1 collects LTB shell-evolution solutions for a flat matter-plus-radiation "
    "track together with dust branches at $E = 0$, $E > 0$ (expand forever), and $E < 0$ (recollapse). "
    "Because the CoE is a fixed direction on the sky, it is a place one can point to: from Central Texas "
    "that direction rises in the northeast and sets in the northwest on a sidereal schedule fixed by the "
    "fitted CoE coordinates. Figure A2 compares CoE rise, culmination, and set with the Sun through 2026; "
    "Table A1 lists the Sun–CoE event crossings for 2027. Figure A3 gives the annual local-sky path of "
    "the CoE from Plano at noon local mean time through 2026—a large oval locus of a fixed celestial "
    "direction in horizon coordinates, not a solar figure-eight.",
)

# 5. Conclusions
add_heading_custom(doc, "5. Conclusions", level=1)
add_body(
    doc,
    "Migkas et al. fix an off-centre CoE+shell—the donut-universe geometry with US at $r = D$—that fits "
    "directional Hubble across redshift slices. Key Newtonian $\\Omega_s = 0.049$ parameters are "
    "$D \\approx 49.17\\,\\mathrm{Gpc}$, $R_{\\mathrm{in}} \\approx 34.48\\,\\mathrm{Gpc}$, "
    "$R_{\\mathrm{out}} \\approx 70.27\\,\\mathrm{Gpc}$, and $\\ell_0 = 13^\\circ$; the preferred "
    "relativistic headline is $D \\approx 32.9\\,\\mathrm{Gpc}$. Gravity on that geometry predicts a "
    "last-scattering multipole whose kill function—an electromagnetic temperature multipole opposite in "
    "sign and matching in angular shape—tracks the kinematic Doppler pattern from the same radii. "
    "Apparent late-time acceleration is an illusion of off-centre sampling, requiring neither dark energy "
    "nor outward shell acceleration; a contracting geometry would share the average appearance of "
    "acceleration but differ in angular signature. The cosmological principle does not account for the "
    "Migkas anisotropy; the CMB isotropy remains an independent constraint. The open problem is no longer "
    "whether the gravitational and electromagnetic shapes talk to each other—they do—but whether a full "
    "time-dependent light-cone calculation can make them cancel at the $10^{-5}$ level without spoiling "
    "the Migkas anisotropy. The goal is coexistence of both datasets under one off-centre shell geometry.",
)

# What we do next (from letter.tex)
add_heading_custom(doc, "What we do next", level=1)
for i, step in enumerate([
    "Evolve $D(t)$, $R_{\\mathrm{in}}(t)$, and $R_{\\mathrm{out}}(t)$ from the present back through recombination, with the observer mid-shell.",
    "Integrate null geodesics on that history to obtain $T(\\theta)$ on our sky, including gravitational and kinematic contributions together.",
    "Ask whether a shared bulk velocity of the emitting surfaces, plus the gravitational multipole, can keep $\\ell\\ge 2$ below $\\sim 10^{-5}$ while the low-redshift $A(\\ell)$ pattern remains Migkas-compatible.",
    "If a residual remains, identify the smallest additional degree of freedom (timing of recombination, shell thickness history, or a mild partner multipole) that removes it without reintroducing a fixed-curve conflict across redshift slices.",
], 1):
    add_body(doc, f"{i}. {step}", space_after=4)

# References
add_heading_custom(doc, "References", level=1)
refs = [
    "1. K. Migkas, G. Schellenberger, T. H. Reiprich, F. Pacaud, M. E. Ramos-Ceja, and L. Lovisari, Astron. Astrophys. 636, A15 (2020); arXiv:2004.03305.",
    "2. L. I. Sedov, J. Appl. Math. Mech. 10, 241 (1946).",
    "3. G. I. Taylor, Proc. R. Soc. Lond. A 201, 159 (1950).",
    "4. G. Lemaître, Ann. Soc. Sci. Bruxelles A53, 51 (1933).",
    "5. R. C. Tolman, Proc. Natl. Acad. Sci. USA 20, 169 (1934).",
    "6. H. Bondi, Mon. Not. R. Astron. Soc. 107, 410 (1947).",
    "7. Planck Collaboration, Astron. Astrophys. 641, A7 (2020).",
]
for r in refs:
    add_para(doc, r, size=10, space_after=3)

# Acknowledgments
add_heading_custom(doc, "Acknowledgments", level=1)
add_body(
    doc,
    "I thank Luke for the important question that sharpened this work, and Dr. Rendler for his support "
    "and guidance. Computational assistance is gratefully acknowledged.",
)
add_para(doc, "Respectfully,  Chris Alameddin", bold=True, space_before=8)

# Appendix
add_heading_custom(doc, "Appendix", level=0)
add_heading_custom(doc, "Shell expansion and pointing at the CoE", level=1)
add_body(
    doc,
    "The LTB shell equation integrates to the usual early radiation and late matter asymptotics, with "
    "the sign of $E(r)$ fixing whether a dust shell expands forever or recollapses. Figure A1 collects "
    "those solutions for a flat matter+radiation track, dust with $E=0$, $E>0$, and $E<0$.",
)
add_figure(
    doc, FIGS / "att14_ltb_shell.png",
    "LTB shell evolution. The flat matter+radiation track is shown with the E=0 dust solution and the "
    "E>0 (expand-forever) and E<0 (recollapsing) branches.",
    width_in=6.3, fig_label="Figure A1",
)
add_body(
    doc,
    "Because the CoE is a fixed direction on the sky in our geometry, it is a place one can point to. "
    "From Central Texas (where this work was done), that direction rises in the northeast and sets in "
    "the northwest on a sidereal schedule fixed by the fitted CoE coordinates. Figure A2 compares CoE "
    "rise, culmination, and set with the Sun through 2026. Table A1 lists the nine Sun–CoE event "
    "crossings for 2027; the highlighted diagonal is rise with rise, highest with highest, and set with set.",
)
add_figure(
    doc, FIGS / "att14_coe_sun.jpg",
    "CoE and Sun timing from Plano through 2026: rise, culmination, and set compared with the Sun.",
    width_in=6.3, fig_label="Figure A2",
)
add_figure(
    doc, FIGS / "att14_coe_crossings.jpg",
    "Sun/CoE event crossings for Central Texas in 2027.",
    width_in=6.0, fig_label="Table A1",
)
add_body(
    doc,
    "Figure A3 gives the annual track of the Centre of Expansion (CoE) on the local sky as seen from "
    "Plano at 12:00 local mean time on each day of 2026. The month labels make the calendar progression "
    "explicit: as the Earth moves through its orbit, one fixed celestial direction traces a large oval "
    "in local horizon coordinates. This is a locus on the sky, not the Sun’s familiar figure-eight "
    "(analemma).",
)
add_figure(
    doc, FIGS / "att14_coe_sky.png",
    "Annual local-sky path of the Centre of Expansion (CoE) from Plano, Texas, evaluated at 12:00 local "
    "mean time on each day of 2026. The colour gradient and labelled dates mark the progression through "
    "the year; the large oval is the locus of the fixed celestial direction in local horizon coordinates, "
    "not a solar figure-eight.",
    width_in=6.3, fig_label="Figure A3",
)
add_figure(
    doc, FIGS / "att14_coe_comic.jpg",
    "CoE pointing comic (appendix).",
    width_in=4.5, fig_label="Figure A4",
)

doc.save(OUT)
print(f"\nWrote {OUT}")
print(f"Size: {OUT.stat().st_size / 1e6:.2f} MB")

# Word / page estimate
from docx.oxml.ns import qn as _qn
paras = len(doc.paragraphs)
words = 0
for p in doc.paragraphs:
    words += len(p.text.split())
print(f"Paragraphs: {paras}, approx words (text only): {words}")
print("Figures embedded: 15 att14_* + equation PNGs")
