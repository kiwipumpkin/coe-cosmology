#!/usr/bin/env python3
"""
Diagnostic for CoE MAX/MIN scan-stop angles 30/185/115/265° from +x.

Shows: raw 3D locus, exact Mollweide projection, and for θ>180 the reduction
via μ=cos(θ)=cos(r) with r=arccos(μ)∈[0,180] (equiv. circle about −x).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from spherical_small_circle_geometry import (
    lb_to_n,
    lon_disp_from_ell,
    mollweide_xy,
    n_to_lb,
    project_small_circle_mollweide,
    small_circle_about_axis,
    split_antimeridian_lonlat,
    theta_scan_to_spherical,
)

OUT = Path(__file__).resolve().parent / "coe_ring_angle_diagnostic.png"

D_HAT = np.array([1.0, 0.0, 0.0])  # anti-CoE = +x at (ℓ,b)=(0,0)
THETAS = [
    ("MAX", 30.0, "#e74c3c"),
    ("MAX", 185.0, "#c0392b"),
    ("MIN", 115.0, "#2980b9"),
    ("MIN", 265.0, "#1a5276"),
]


def draw_wire(ax):
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(-np.pi / 2, np.pi / 2, 20)
    ax.plot_surface(
        np.outer(np.cos(u), np.cos(v)),
        np.outer(np.sin(u), np.cos(v)),
        np.outer(np.ones_like(u), np.sin(v)),
        color="#dfe6e9",
        alpha=0.15,
        linewidth=0,
        shade=False,
    )


def mollweide_axes(ax):
    th = np.linspace(-np.pi / 2, np.pi / 2, 181)
    xr = (2 * np.sqrt(2) / np.pi) * np.pi * np.cos(th)
    yr = np.sqrt(2) * np.sin(th)
    ax.plot(xr, yr, "k-", lw=0.8)
    ax.plot(-xr, yr, "k-", lw=0.8)
    for b in (-60, -30, 0, 30, 60):
        ell = np.linspace(0, 360, 361)
        lon = np.deg2rad(lon_disp_from_ell(ell))
        lat = np.full_like(lon, np.deg2rad(b))
        for lo, la in split_antimeridian_lonlat(lon, lat):
            x, y = mollweide_xy(lo, la)
            ax.plot(x, y, color="#bdc3c7", lw=0.35, zorder=0)
    for ell0 in (0, 90, 180, 270):
        b = np.linspace(-90, 90, 181)
        lon = np.full_like(b, np.deg2rad(lon_disp_from_ell(ell0)))
        x, y = mollweide_xy(lon, np.deg2rad(b))
        ax.plot(x, y, color="#bdc3c7", lw=0.35, zorder=0)
    ax.set_aspect("equal")
    ax.set_xlim(-2.9, 2.9)
    ax.set_ylim(-1.55, 1.55)
    ax.axis("off")


def main():
    fig = plt.figure(figsize=(14.0, 12.5), dpi=140)
    fig.suptitle(
        r"CoE ring angle diagnostic: scan stops $\theta$ from $+x$ (anti-CoE)"
        "\n"
        r"locus $\mathbf{n}\cdot\hat{\mathbf{x}}=\mu=\cos\theta$; "
        r"spherical radius $r=\arccos\mu\in[0,180^\circ]$ only",
        fontsize=12,
        fontweight="bold",
        y=0.995,
    )

    # Top summary table as text panel
    axT = fig.add_subplot(3, 1, 1)
    axT.axis("off")
    lines = [
        r"Scan angle θ (0–360° along a great-circle strip) vs spherical small-circle:",
        "",
    ]
    for kind, th, col in THETAS:
        r, r_coe, mu = theta_scan_to_spherical(th)
        extra = ""
        if th > 180:
            extra = (
                rf"  → θ>180° is NOT a spherical distance; "
                rf"$\mu=\cos{th:.0f}^\circ=\cos{r:.0f}^\circ$ "
                rf"so same cone as $r={r:.0f}^\circ$ from $+x$ "
                rf"(= ${r_coe:.0f}^\circ$ from CoE/$-x$)."
            )
        lines.append(
            rf"{kind}: θ={th:.0f}°   μ={mu:+.6f}   r_from_+x={r:.1f}°   "
            rf"r_from_CoE={r_coe:.1f}°{extra}"
        )
    lines += [
        "",
        "Hence MAX 185° ring ≡ 5° circle about CoE (ℓ=180°,b=0) — appears as TWO small "
        "arcs at the LEFT and RIGHT Mollweide edges near the equator.",
        "MIN 265° ring ≡ 95° circle about +x (same μ as cos 95°) — large near-meridian "
        "ovals, often split by the antimeridian cut.",
        "Drawing a decorative oval-rim ‘outer loop’ for μ≲−0.95 is geometrically false.",
    ]
    axT.text(
        0.02,
        0.98,
        "\n".join(lines),
        transform=axT.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        family="monospace",
        wrap=True,
        bbox=dict(boxstyle="round", facecolor="#f8f9f9", edgecolor="#555"),
    )

    # Row of 4: 3D
    for i, (kind, th, col) in enumerate(THETAS):
        ax = fig.add_subplot(3, 4, 4 + i + 1, projection="3d")
        draw_wire(ax)
        r, r_coe, mu = theta_scan_to_spherical(th)
        pts = small_circle_about_axis(D_HAT, mu, n=800)
        pts_c = np.vstack([pts, pts[:1]])
        ax.plot(pts_c[:, 0], pts_c[:, 1], pts_c[:, 2], color=col, lw=2.0)
        ax.plot([-1.2, 1.2], [0, 0], [0, 0], color="#5b2c6f", lw=1.5)
        ax.scatter([1], [0], [0], c=["#c0392b"], s=40, depthshade=False)
        ax.scatter([-1], [0], [0], c=["#6b3fa0"], s=50, marker="*", depthshade=False)
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.2, 1.2)
        ax.set_zlim(-1.2, 1.2)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=20, azim=-55)
        ax.set_axis_off()
        ax.set_title(
            rf"{kind} θ={th:.0f}° → r={r:.0f}° (+x)"
            + (rf" ≡ {r_coe:.0f}° (CoE)" if th > 180 else ""),
            fontsize=8,
            color=col,
        )

    # Row of 4: Mollweide
    for i, (kind, th, col) in enumerate(THETAS):
        ax = fig.add_subplot(3, 4, 8 + i + 1)
        mollweide_axes(ax)
        r, r_coe, mu = theta_scan_to_spherical(th)
        # Project via exact radius from +x center (0,0)
        segs = project_small_circle_mollweide(0.0, 0.0, r, n=2881)
        for x, y, lo, la in segs:
            ax.plot(x, y, color=col, lw=1.8, zorder=3)
        # anti-CoE / CoE markers
        x0, y0 = mollweide_xy(np.array([0.0]), np.array([0.0]))
        ax.plot(x0, y0, "o", color="#c0392b", ms=5, zorder=5)
        for sgn in (-1.0, 1.0):
            xc, yc = mollweide_xy(np.array([sgn * np.pi]), np.array([0.0]))
            ax.plot(xc, yc, "*", color="#6b3fa0", ms=8, zorder=5)
        ax.set_title(
            rf"Mollweide: {len(segs)} seg; μ={mu:+.3f}",
            fontsize=8,
            color=col,
        )

    handles = [
        Line2D([0], [0], color="#c0392b", marker="o", ls="None", label="anti-CoE +x (0,0)"),
        Line2D([0], [0], color="#6b3fa0", marker="*", ls="None", label="CoE (180,0) both edges"),
        Line2D([0], [0], color="#5b2c6f", lw=2, label="CoE↔anti-CoE axis"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9, frameon=True)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    fig.savefig(OUT, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT}")
    for kind, th, col in THETAS:
        r, r_coe, mu = theta_scan_to_spherical(th)
        print(f"{kind} θ={th}: mu={mu:.8f} r={r:.4f} r_CoE={r_coe:.4f}")


if __name__ == "__main__":
    main()
