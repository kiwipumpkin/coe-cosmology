#!/usr/bin/env python3
"""
8-case Mollweide small-circle validation matching Chris's reference centers
(radius 20°). Exact equation n·n0 = cos r; seam-crossing cases drawn as
boundary segments on BOTH edges — no interior ellipse, no rim bridges.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from spherical_small_circle_geometry import (
    lb_to_n,
    lon_disp_from_ell,
    mollweide_xy,
    n_to_lb,
    project_small_circle_mollweide,
    small_circle_points,
    wrap_lon_deg,
)

OUT = Path(__file__).resolve().parent / "mollweide_small_circle_validation.png"

# Match the user's reference infographic as closely as possible
CASES = [
    (1, 0.0, +60.0, "#e74c3c"),
    (2, 120.0, +45.0, "#3498db"),
    (3, 240.0, +30.0, "#27ae60"),
    (4, 60.0, 0.0, "#e67e22"),
    (5, 180.0, 0.0, "#8e44ad"),  # seam: must split L/R
    (6, 300.0, -30.0, "#a04000"),
    (7, 90.0, -45.0, "#16a085"),
    (8, 200.0, -60.0, "#c0392b"),  # near seam
]
R = 20.0


def draw_wire_sphere(ax, n_mer=12, n_par=6):
    u = np.linspace(0, 2 * np.pi, 48)
    v = np.linspace(-np.pi / 2, np.pi / 2, 24)
    xs = np.outer(np.cos(u), np.cos(v))
    ys = np.outer(np.sin(u), np.cos(v))
    zs = np.outer(np.ones_like(u), np.sin(v))
    ax.plot_surface(xs, ys, zs, color="#dfe6e9", alpha=0.18, linewidth=0, shade=False)
    kw = dict(color="#95a5a6", lw=0.35, alpha=0.5)
    for i in range(n_mer):
        phi = 2 * np.pi * i / n_mer
        t = np.linspace(-np.pi / 2, np.pi / 2, 60)
        ax.plot(np.cos(t) * np.cos(phi), np.cos(t) * np.sin(phi), np.sin(t), **kw)
    for j in range(1, n_par):
        b = -np.pi / 2 + np.pi * j / n_par
        t = np.linspace(0, 2 * np.pi, 60)
        r = np.cos(b)
        ax.plot(r * np.cos(t), r * np.sin(t), np.full_like(t, np.sin(b)), **kw)


def draw_mollweide_grid(ax):
    # Boundary
    t = np.linspace(-np.pi, np.pi, 361)
    xb = 2 * np.sqrt(2) * np.cos(t)
    # Mollweide boundary is the ellipse x^2/(2√2)^2 + y^2/(√2)^2 = 1
    # Parametrize via lon=±π at all aux θ:
    th = np.linspace(-np.pi / 2, np.pi / 2, 181)
    x_r = (2 * np.sqrt(2) / np.pi) * np.pi * np.cos(th)
    y_r = np.sqrt(2) * np.sin(th)
    ax.plot(x_r, y_r, "k-", lw=0.9)
    ax.plot(-x_r, y_r, "k-", lw=0.9)
    for b in range(-75, 76, 15):
        ell = np.linspace(0, 360, 361)
        lon = np.deg2rad(lon_disp_from_ell(ell))
        lat = np.full_like(lon, np.deg2rad(b))
        # split at seam
        d = np.diff(lon)
        cuts = np.where(np.abs(d) > np.pi)[0] + 1
        idx = np.r_[0, cuts, len(lon)]
        for i0, i1 in zip(idx[:-1], idx[1:]):
            if i1 - i0 < 2:
                continue
            x, y = mollweide_xy(lon[i0:i1], lat[i0:i1])
            ax.plot(x, y, color="#bdc3c7", lw=0.4, zorder=0)
    for ell0 in range(0, 360, 30):
        b = np.linspace(-90, 90, 181)
        lon = np.full_like(b, np.deg2rad(lon_disp_from_ell(ell0)))
        lat = np.deg2rad(b)
        x, y = mollweide_xy(lon, lat)
        ax.plot(x, y, color="#bdc3c7", lw=0.4, zorder=0)
    ax.set_aspect("equal")
    ax.set_xlim(-2.9, 2.9)
    ax.set_ylim(-1.55, 1.55)
    ax.axis("off")
    # tick labels
    for ell, lab in [(180, "180°"), (90, "90°"), (0, "0°"), (270, "270°"), (180, "180°")]:
        pass
    for ell, lab, ha in [
        (180, r"$180^\circ$", "right"),
        (90, r"$90^\circ$", "center"),
        (0, r"$0^\circ$", "center"),
        (270, r"$270^\circ$", "center"),
    ]:
        x, y = mollweide_xy(np.deg2rad(lon_disp_from_ell(ell)), np.array([0.0]))
        # place 180 twice
        if ell == 180:
            ax.text(-2.82, -0.08, lab, fontsize=6, ha="left", va="top", color="#333")
            ax.text(2.82, -0.08, lab, fontsize=6, ha="right", va="top", color="#333")
        else:
            ax.text(float(np.asarray(x).ravel()[0]), -0.08, lab, fontsize=6, ha=ha, va="top", color="#333")
    for b, lab in [(90, r"$+90^\circ$"), (45, r"$+45^\circ$"), (0, r"$0^\circ$"), (-45, r"$-45^\circ$"), (-90, r"$-90^\circ$")]:
        x, y = mollweide_xy(np.array([0.0]), np.deg2rad([b]))
        ax.text(0.05, float(np.asarray(y).ravel()[0]), lab, fontsize=6, ha="left", va="center", color="#333")


def main():
    fig = plt.figure(figsize=(14.5, 17.5), dpi=140)
    fig.suptitle(
        r"Polar-type circles (small circles) at different locations and their Mollweide projections"
        "\n"
        r"Exact locus $\mathbf{n}\cdot\mathbf{n}_0=\cos r$ with $r=20^\circ$; "
        r"Galactic Mollweide, $\ell$ increases LEFT ($\lambda=-\mathrm{wrap}(\ell)$)"
        "\n"
        r"Seam-crossing circles = multiple boundary segments (no interior chords / rim bridges)",
        fontsize=11,
        fontweight="bold",
        y=0.995,
    )

    # layout: 8 rows × (sphere | mollweide)
    for i, (num, l0, b0, color) in enumerate(CASES):
        row = i
        ax3 = fig.add_subplot(8, 2, 2 * row + 1, projection="3d")
        axm = fig.add_subplot(8, 2, 2 * row + 2)

        # --- sphere ---
        draw_wire_sphere(ax3)
        pts = small_circle_points(l0, b0, R, n=721)
        # close for plot
        pts_c = np.vstack([pts, pts[:1]])
        ax3.plot(pts_c[:, 0], pts_c[:, 1], pts_c[:, 2], color=color, lw=2.0, zorder=5)
        c = lb_to_n(l0, b0)
        ax3.scatter([c[0]], [c[1]], [c[2]], c=[color], s=28, depthshade=False, zorder=6)
        ax3.set_xlim(-1.15, 1.15)
        ax3.set_ylim(-1.15, 1.15)
        ax3.set_zlim(-1.15, 1.15)
        ax3.set_box_aspect((1, 1, 1))
        ax3.view_init(elev=18, azim=-60)
        ax3.set_axis_off()
        ax3.set_title(
            rf"{num}.  center $(\ell,b)=({l0:.0f}^\circ,{b0:+.0f}^\circ)$,  $r={R:.0f}^\circ$",
            fontsize=9,
            color=color,
            pad=2,
        )

        # --- Mollweide ---
        draw_mollweide_grid(axm)
        segs = project_small_circle_mollweide(l0, b0, R, n=2881)
        for x, y, lo, la in segs:
            axm.plot(x, y, color=color, lw=2.0, zorder=3, solid_capstyle="round")
        # center marker (may sit on seam for case 5 — plot both edges)
        if abs(wrap_lon_deg(l0)) == 180.0 or abs(wrap_lon_deg(l0) - (-180)) < 1e-9 or abs(l0 - 180) < 1e-9:
            for sgn in (-1.0, 1.0):
                x, y = mollweide_xy(np.array([sgn * np.pi]), np.deg2rad([b0]))
                axm.plot(x, y, "o", color=color, ms=4, zorder=4)
        else:
            x, y = mollweide_xy(
                np.deg2rad([lon_disp_from_ell(l0)]), np.deg2rad([b0])
            )
            axm.plot(x, y, "o", color=color, ms=4, zorder=4)

        nseg = len(segs)
        note = ""
        if num == 5:
            note = "\nSEAM: center on antimeridian → arcs on BOTH map edges"
        elif num == 8:
            note = "\nnear seam: split L/R if circle crosses λ=±π"
        axm.set_title(
            rf"Mollweide  ({nseg} segment{'s' if nseg != 1 else ''}){note}",
            fontsize=8,
            color="#222",
            pad=2,
        )

    fig.text(
        0.5,
        0.005,
        "Circles away from poles → shifted/tilted oval-like curves. "
        "Near poles → compressed arcs. "
        "If a circle crosses the chosen cut (λ=±π ↔ ℓ=180°), it appears on both left and right "
        "boundaries — topologically one circle, not a closed interior ellipse. "
        "Reference infographic case 5 drawing a closed oval near ℓ=90° is incorrect for center (180°,0°).",
        ha="center",
        va="bottom",
        fontsize=8,
        color="#222",
        wrap=True,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(OUT, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
