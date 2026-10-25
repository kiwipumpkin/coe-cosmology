#!/usr/bin/env python3
"""
CoE geometry schematic: allowed CMB dipole along D.

Chris markup (2026-09-22):
  - Orange lens blob on far left (between R_CMB and outer shells)
  - Thick purple lensed ray: from CMB side, skims below lens, bends down to us
  - Red squiggly radiation arrows from upper-left / left toward interior
  - Vision shell radius = 2r (was r)
  - Remove prior near-CoE kink-as-lens illustration

Outputs PNG + PDF under figures/, plus a brief JSON change note.
Also refreshes drive_upload_payloads copy when present.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

OUT = Path(__file__).resolve().parent
PNG = OUT / "coe_geometry_cmb_dipole_schematic.png"
PDF = OUT / "coe_geometry_cmb_dipole_schematic.pdf"
JSON_NOTE = OUT / "coe_geometry_cmb_dipole_schematic_changes.json"

DRIVE_COPY = (
    OUT.parent
    / "drive_upload_payloads"
    / "papers"
    / "stage_setting_2026-09-20"
    / "figures"
    / "coe_geometry_cmb_dipole_schematic.png"
)

# Wong / Tol palette (colorblind-safe)
C_BLUE = "#0072B2"
C_TEAL = "#009E73"
C_ORANGE = "#E69F00"
C_VERM = "#D55E00"
C_PURPLE = "#CC79A7"
C_SKY = "#56B4E9"
C_GRAY = "#666666"
C_BLACK = "#222222"
C_GREEN = "#009E73"


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 220,
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def squiggle_arrow(ax, x0, y0, x1, y1, *, n_waves=4, amp=0.045, color=C_VERM, lw=1.6, z=5):
    """Wavy polyline with arrowhead (radiation)."""
    t = np.linspace(0.0, 1.0, 120)
    x = x0 + t * (x1 - x0)
    y = y0 + t * (y1 - y0)
    # perpendicular unit
    dx, dy = x1 - x0, y1 - y0
    L = np.hypot(dx, dy) + 1e-12
    px, py = -dy / L, dx / L
    y = y + amp * np.sin(2 * np.pi * n_waves * t) * np.linspace(1.0, 0.35, t.size)
    x = x + amp * 0.15 * np.sin(2 * np.pi * n_waves * t + 0.4) * px
    ax.plot(x, y, color=color, lw=lw, solid_capstyle="round", zorder=z)
    ax.annotate(
        "",
        xy=(x1, y1),
        xytext=(x[-8], y[-8]),
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            lw=lw,
            mutation_scale=11,
            shrinkA=0,
            shrinkB=0,
        ),
        zorder=z + 1,
    )


def irregular_blob(ax, cx, cy, rx, ry, *, facecolor, edgecolor, alpha=0.85, z=6, seed=7):
    """Soft irregular ellipse (lens object), not a perfect circle."""
    rng = np.random.default_rng(seed)
    th = np.linspace(0, 2 * np.pi, 60, endpoint=True)
    jitter = 1.0 + 0.12 * np.sin(3 * th + 0.5) + 0.08 * np.sin(5 * th + 1.2)
    jitter += 0.04 * rng.normal(size=th.size)
    jitter[-1] = jitter[0]
    xs = cx + rx * jitter * np.cos(th)
    ys = cy + ry * jitter * np.sin(th)
    ax.fill(xs, ys, facecolor=facecolor, edgecolor=edgecolor, lw=1.8, alpha=alpha, zorder=z)
    # inner highlight
    ax.fill(
        cx + 0.45 * (xs - cx),
        cy + 0.45 * (ys - cy) + 0.02,
        facecolor="#FFE6A8",
        edgecolor="none",
        alpha=0.45,
        zorder=z + 1,
    )


def lensed_ray_path(blob_xy, us_xy, r_cmb):
    """
    Thick purple path: from CMB (left), skim just below orange lens, bend down to us.
    Returns Nx2 array.
    """
    bx, by = blob_xy
    ux, uy = us_xy
    # start on / just inside R_CMB on the left, slightly below blob height
    start = np.array([-0.92 * r_cmb, by - 0.02])
    # skim point just below blob
    skim = np.array([bx + 0.05, by - 0.16])
    # continue a bit past the lens before bending down
    mid = np.array([bx + 0.65, by - 0.28])
    # approach observer from upper-left
    approach = np.array([ux - 0.40, uy + 0.28])
    end = np.array([ux, uy])

    # piecewise cubic via Catmull-Rom / dense lerp of Bezier segments
    ctrl = np.vstack([start, skim, mid, approach, end])
    segs = []
    for i in range(len(ctrl) - 1):
        p0 = ctrl[max(i - 1, 0)]
        p1 = ctrl[i]
        p2 = ctrl[i + 1]
        p3 = ctrl[min(i + 2, len(ctrl) - 1)]
        tt = np.linspace(0, 1, 40, endpoint=False)
        # Catmull-Rom
        t = tt[:, None]
        a0 = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3
        a1 = p0 - 2.5 * p1 + 2.0 * p2 - 0.5 * p3
        a2 = -0.5 * p0 + 0.5 * p2
        a3 = p1
        segs.append(a0 * t**3 + a1 * t**2 + a2 * t + a3)
    segs.append(end[None, :])
    return np.vstack(segs)


def draw(ax: plt.Axes) -> dict:
    # Geometry in units of D (observer offset)
    D = 1.0
    r = 0.36          # vision-shell radius
    r_vision = r      # Chris: vision is r
    # Sn to Sn+1 spacing is 2r
    R_S3 = 0.85
    R_S2 = R_S3 + 2.0 * r
    R_S1 = R_S2 + 2.0 * r
    R_CMB = R_S1 + 0.55
    shell_lw = 1.6

    coe = np.array([0.0, 0.0])
    us = np.array([D, 0.0])

    ax.set_aspect("equal")
    ax.set_xlim(-3.15, 2.25)
    ax.set_ylim(-2.35, 2.55)
    ax.axis("off")

    # --- R_CMB dashed ---
    ax.add_patch(
        Circle(
            coe,
            R_CMB,
            fill=False,
            ec=C_GRAY,
            lw=1.4,
            ls=(0, (5, 4)),
            zorder=1,
        )
    )
    ax.text(
        0.35,
        R_CMB + 0.08,
        r"$R_{\mathrm{CMB}}$ (last-scattering shell about CoE)",
        fontsize=10,
        color=C_GRAY,
        ha="center",
        va="bottom",
    )

    # --- concentric CoE shells ---
    for R, color, label, lxy in [
        (R_S1, C_BLUE, r"$S_1$", (0.55, R_S1 + 0.08)),
        (R_S2, C_SKY, r"$S_2$", (0.35, R_S2 + 0.08)),
        (R_S3, C_ORANGE, r"$S_3$", (0.15, R_S3 + 0.08)),
    ]:
        ax.add_patch(Circle(coe, R, fill=False, ec=color, lw=shell_lw, zorder=2))
        # small label chip on upper-right arc
        ang = np.deg2rad(55)
        lx, ly = R * np.cos(ang), R * np.sin(ang)
        ax.text(
            lx,
            ly,
            label,
            fontsize=9,
            color=color,
            ha="center",
            va="center",
            bbox=dict(
                boxstyle="round,pad=0.18",
                fc="white",
                ec=color,
                lw=1.0,
                alpha=0.95,
            ),
            zorder=8,
        )

    # --- CoE marker (glow is CoE, NOT the lens) ---
    ax.add_patch(
        Circle(coe, 0.16, facecolor="#F2B8B8", edgecolor=C_VERM, lw=1.35, alpha=0.92, zorder=4)
    )
    ax.plot(*coe, "o", color=C_VERM, ms=8, zorder=7, markeredgecolor="white", markeredgewidth=0.7)
    ax.text(0.0, -0.28, "CoE", fontsize=11, color=C_VERM, ha="center", va="top", fontweight="bold")

    # --- D axis ---
    ax.plot([0, D], [0, 0], ls="--", color=C_VERM, lw=1.2, zorder=3)
    ax.annotate(
        "",
        xy=(D - 0.02, -0.42),
        xytext=(0.02, -0.42),
        arrowprops=dict(arrowstyle="<->", color=C_BLACK, lw=1.2, mutation_scale=10),
        zorder=5,
    )
    ax.text(D / 2, -0.52, r"$D$", fontsize=12, color=C_BLACK, ha="center", va="top")

    # --- observer + vision shell (radius r) ---
    ax.add_patch(
        Circle(
            us,
            r_vision,
            facecolor="#C8E6C9",
            edgecolor=C_GREEN,
            lw=1.6,
            alpha=0.45,
            zorder=3,
        )
    )
    ax.plot(*us, "o", color="#004C2E", ms=8, zorder=9, markeredgecolor="white", markeredgewidth=0.7)
    ax.text(D, -0.18, "us", fontsize=11, color="#004C2E", ha="center", va="top", fontweight="bold")
    ax.text(
        D + r_vision * 0.15,
        r_vision + 0.08,
        r"vision shell $r$",
        fontsize=10,
        color=C_GREEN,
        ha="left",
        va="bottom",
    )

    # sight-lines from us (keep prior schematic feel)
    for ang_deg in (25, 55, -20):
        ang = np.deg2rad(ang_deg)
        ax.plot(
            [D, D + 1.35 * np.cos(ang)],
            [0, 1.35 * np.sin(ang)],
            color="#888888",
            lw=0.9,
            zorder=2,
        )

    # anti-CoE shorter path
    ax.annotate(
        "",
        xy=(D + 0.85, 0.0),
        xytext=(D + 0.08, 0.0),
        arrowprops=dict(arrowstyle="-|>", color=C_GREEN, lw=2.0, mutation_scale=12),
        zorder=6,
    )
    ax.text(
        D + 0.95,
        0.12,
        "anti-CoE\n(shorter path)",
        fontsize=9,
        color=C_GREEN,
        ha="left",
        va="bottom",
    )

    # --- orange lens blob (far left, outside/near shells, toward CoE side) ---
    blob_x = -0.5 * (R_S1 + R_CMB)  # between S1 and R_CMB
    blob_y = 0.22
    irregular_blob(
        ax,
        blob_x,
        blob_y,
        rx=0.22,
        ry=0.17,
        facecolor=C_ORANGE,
        edgecolor="#B36B00",
        alpha=0.9,
        z=6,
    )
    # --- thick purple lensed ray: skim below blob, bend down to us ---
    ray = lensed_ray_path((blob_x, blob_y), us, R_CMB)
    ax.plot(
        ray[:, 0],
        ray[:, 1],
        color=C_PURPLE,
        lw=3.2,
        solid_capstyle="round",
        solid_joinstyle="round",
        zorder=7,
        alpha=0.95,
    )
    # arrowheads: inbound from CMB, and into observer
    ax.annotate(
        "",
        xy=tuple(ray[12]),
        xytext=tuple(ray[2]),
        arrowprops=dict(arrowstyle="-|>", color=C_PURPLE, lw=2.4, mutation_scale=14),
        zorder=8,
    )
    ax.annotate(
        "",
        xy=tuple(us),
        xytext=tuple(ray[-18]),
        arrowprops=dict(
            arrowstyle="-|>",
            color=C_PURPLE,
            lw=2.6,
            mutation_scale=14,
            shrinkA=0,
            shrinkB=4,
        ),
        zorder=8,
    )
    ax.text(
        -1.85,
        blob_y + 0.05,
        "lensed ray\n(skim / eclipse-like)",
        fontsize=9.5,
        color=C_PURPLE,
        ha="center",
        va="bottom",
        fontweight="bold",
    )
    ax.text(
        0.35,
        0.38,
        "toward CoE\n(longer path)",
        fontsize=9.5,
        color=C_PURPLE,
        ha="center",
        va="bottom",
        fontweight="bold",
    )

    # --- red squiggly radiation from upper-left / left toward interior ---
    squiggle_arrow(ax, -2.75, 0.55, -2.15, 0.45, n_waves=3, amp=0.04, color="#C41E3A", lw=1.8)
    squiggle_arrow(ax, -2.75, 0.20, -2.15, 0.18, n_waves=3, amp=0.035, color="#C41E3A", lw=1.8)
    squiggle_arrow(ax, -2.75, -0.15, -2.15, -0.05, n_waves=3, amp=0.035, color="#C41E3A", lw=1.8)
    # diagonal upper-left toward shells
    squiggle_arrow(ax, -2.35, 1.55, -1.55, 1.05, n_waves=3, amp=0.05, color="#C41E3A", lw=1.8)
    squiggle_arrow(ax, -2.05, 1.85, -1.25, 1.25, n_waves=3, amp=0.05, color="#C41E3A", lw=1.8)
    ax.text(
        -2.78,
        0.85,
        r"radiation from" + "\n" + r"$R_{\mathrm{CMB}}$",
        fontsize=9,
        color="#C41E3A",
        ha="left",
        va="bottom",
        fontweight="bold",
    )

    # title
    ax.set_title(
        r"CoE geometry: allowed CMB dipole along $\mathbf{D}$",
        fontsize=14,
        pad=10,
        fontweight="bold",
        color=C_BLACK,
    )

    # caption box
    caption = (
        r"Off-center CoE geometry (schematic). CMB from $R_{\mathrm{CMB}}\gg D$; "
        r"vision shell radius $r$; $S_n\!\to\!S_{n+1}$ spacing $2r$; offset $D$. Opposite sides of $\mathbf{D}$ "
        r"$\Rightarrow$ geometric CMB dipole; lensing (orange) near the CoE-facing "
        r"sky warps paths (purple skim / eclipse-like ray). Late Hubble quadrupole "
        r"lives on cluster /$r$ scales (not shown as CMB quad)."
    )
    ax.text(
        0.5,
        -0.02,
        caption,
        transform=ax.transAxes,
        fontsize=8.5,
        color="#444444",
        ha="center",
        va="top",
        wrap=True,
        bbox=dict(boxstyle="round,pad=0.45", fc="#F4F4F4", ec="#CCCCCC", lw=0.8),
    )

    return {
        "D": D,
        "r": r,
        "r_vision": r_vision,
        "R_S1": R_S1,
        "R_S2": R_S2,
        "R_S3": R_S3,
        "R_CMB": R_CMB,
        "lens_blob_xy": [float(blob_x), float(blob_y)],
        "shell_lw": shell_lw,
    }


def main() -> None:
    apply_style()
    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    meta = draw(ax)
    fig.tight_layout(rect=(0.01, 0.08, 0.99, 0.98))
    fig.savefig(PNG, dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(PDF, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    note = {
        "updated": "2026-09-22",
        "outputs": [str(PNG), str(PDF)],
        "changes": [
            "Replaced prior near-CoE kink-as-lens illustration with orange irregular lens blob on far left (between R_CMB and S1).",
            "Thick purple lensed ray now comes from CMB side, skims just below the orange lens, then bends down to observer (us).",
            "Added red squiggly radiation arrows from left and upper-left toward interior shells.",
            "Vision shell radius = r; concentric shell spacing Sn to Sn+1 = 2r.",
            "Kept CoE, observer at D, concentric S1/S2/S3, R_CMB dashed, toward-CoE / anti-CoE path labels, title about allowed CMB dipole along D.",
        ],
        "geometry": meta,
        "note_cmb_dipole_estimate": (
            "cmb_dipole_estimate_coe_geometry.png is a quantitative dipole-budget "
            "plot (not this schematic); left unchanged."
        ),
        "script": str(OUT / "make_coe_geometry_cmb_dipole_schematic.py"),
    }
    JSON_NOTE.write_text(json.dumps(note, indent=2) + "\n")

    if DRIVE_COPY.parent.is_dir():
        shutil.copy2(PNG, DRIVE_COPY)
        note["drive_copy"] = str(DRIVE_COPY)
        JSON_NOTE.write_text(json.dumps(note, indent=2) + "\n")

    print(f"Wrote {PNG}")
    print(f"Wrote {PDF}")
    print(f"Wrote {JSON_NOTE}")
    if DRIVE_COPY.parent.is_dir():
        print(f"Updated {DRIVE_COPY}")


if __name__ == "__main__":
    main()
