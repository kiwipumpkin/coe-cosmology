#!/usr/bin/env python3
"""
Demo: on Mollweide, a Euclidean page-circle is round, but a true sky small-circle
about the map center projects to an oval/curve. Area-preserving ≠ shape-preserving.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from make_migkas_fig8_coe_rings_overlay import (
    D_HAT,
    lon_disp_from_ell,
    mollweide_xy,
    n_to_lb,
    small_circle_about_axis,
)

OUT = Path(__file__).resolve().parent / "mollweide_circle_vs_oval_demo.png"


def main():
    # Normalized Mollweide plane: Xn in [-1,1], Yn in [-1,1] (same as overlay)
    # Page coords: x = Xn * 2*sqrt(2), y = Yn * sqrt(2) for native; we plot Xn,Yn.
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), dpi=140)

    # ---- shared blank Mollweide grid (lon_disp = -wrap(ell); ell increases LEFT)
    def draw_grid(ax):
        # oval boundary
        t = np.linspace(0, 2 * np.pi, 721)
        ax.plot(np.cos(t), np.sin(t), color="#222", lw=1.4, zorder=2)
        # meridians (ell = 0, ±60, +/-120, 180) via lon_disp
        for ell in (0, 60, 120, 180, 240, 300):
            b = np.linspace(-89.5, 89.5, 400)
            lon = np.deg2rad(lon_disp_from_ell(np.full_like(b, ell)))
            lat = np.deg2rad(b)
            x, y = mollweide_xy(lon, lat)
            Xn = x / (2.0 * np.sqrt(2.0))
            Yn = y / np.sqrt(2.0)
            ax.plot(Xn, Yn, color="#bbbbbb", lw=0.6, zorder=1)
        # parallels
        for bdeg in (-60, -30, 0, 30, 60):
            ell = np.linspace(0, 360, 721, endpoint=False)
            lon = np.deg2rad(lon_disp_from_ell(ell))
            lat = np.deg2rad(np.full_like(ell, bdeg))
            x, y = mollweide_xy(lon, lat)
            Xn = x / (2.0 * np.sqrt(2.0))
            Yn = y / np.sqrt(2.0)
            # break at antimeridian jump
            d = np.hypot(np.diff(Xn), np.diff(Yn))
            brk = np.where(d > 0.15)[0]
            idx = np.r_[0, brk + 1, len(Xn)]
            for a, b in zip(idx[:-1], idx[1:]):
                if b - a >= 2:
                    ax.plot(Xn[a:b], Yn[a:b], color="#bbbbbb", lw=0.6, zorder=1)
        ax.set_aspect("equal")
        ax.set_xlim(-1.08, 1.08)
        ax.set_ylim(-1.08, 1.08)
        ax.axis("off")
        ax.axhline(0, color="#ddd", lw=0.5, zorder=0)
        ax.axvline(0, color="#ddd", lw=0.5, zorder=0)

    # Panel A: Euclidean circle on the page
    ax = axes[0]
    draw_grid(ax)
    for r, ls in ((0.35, "-"), (0.65, "-"), (0.92, "--")):
        t = np.linspace(0, 2 * np.pi, 400)
        ax.plot(
            r * np.cos(t), r * np.sin(t),
            color="#c0392b", lw=2.0, ls=ls, zorder=3,
            label="circle on the page" if r == 0.65 else None,
        )
    ax.set_title("A. Euclidean circles in the Mollweide plane", fontsize=11, pad=8)
    ax.legend(loc="lower center", fontsize=9, frameon=True, fancybox=False)
    ax.text(
        0, -1.02, r"perfectly round in $x$–$y$ (page geometry)",
        ha="center", va="top", fontsize=8.5, color="#444",
    )

    # Panel B: true sky small-circles about +x (map center / anti-CoE)
    ax = axes[1]
    draw_grid(ax)
    colors = {30: "#e74c3c", 60: "#2980b9", 90: "#8e44ad"}
    for theta in (30.0, 60.0, 90.0):
        mu = float(np.cos(np.deg2rad(theta)))
        pts = small_circle_about_axis(D_HAT, mu, n=2001)
        ell = np.array([n_to_lb(p)[0] for p in pts])
        b = np.array([n_to_lb(p)[1] for p in pts])
        lon = np.deg2rad(lon_disp_from_ell(ell))
        lat = np.deg2rad(b)
        x, y = mollweide_xy(lon, lat)
        Xn = x / (2.0 * np.sqrt(2.0))
        Yn = y / np.sqrt(2.0)
        phi = np.arctan2(pts[:, 2], pts[:, 1])
        order = np.argsort(phi)
        Xn, Yn = Xn[order], Yn[order]
        # segment on jumps
        segs_x, segs_y = [Xn[0]], [Yn[0]]
        drawn_label = False
        for i in range(1, len(Xn)):
            if np.hypot(Xn[i] - segs_x[-1], Yn[i] - segs_y[-1]) > 0.25:
                if len(segs_x) >= 2:
                    ax.plot(
                        segs_x, segs_y, color=colors[int(theta)], lw=2.0, zorder=3,
                        label=(
                            "constant angle from center (sky circles)"
                            if theta == 60 and not drawn_label else None
                        ),
                    )
                    drawn_label = True
                segs_x, segs_y = [Xn[i]], [Yn[i]]
            else:
                segs_x.append(Xn[i])
                segs_y.append(Yn[i])
        if len(segs_x) >= 2:
            # close if single loop
            close = np.hypot(segs_x[0] - segs_x[-1], segs_y[0] - segs_y[-1]) < 0.25
            xs = list(segs_x) + ([segs_x[0]] if close else [])
            ys = list(segs_y) + ([segs_y[0]] if close else [])
            ax.plot(
                xs, ys, color=colors[int(theta)], lw=2.0, zorder=3,
                label=(
                    "constant angle from center (sky circles)"
                    if theta == 60 and not drawn_label else None
                ),
            )
        # labels placed off the equator so they do not collide
        label_xy = {30: (0.22, 0.22), 60: (0.48, 0.42), 90: (0.72, 0.68)}
        lx, ly = label_xy[int(theta)]
        ax.text(
            lx, ly, rf"$\theta={theta:.0f}^\circ$",
            color=colors[int(theta)], fontsize=9, va="center", ha="left",
        )

    ax.set_title("B. True sky small-circles about map center (+x)", fontsize=11, pad=8)
    ax.legend(loc="lower center", fontsize=9, frameon=True, fancybox=False)
    ax.text(
        0, -1.02, "sky circles → ovals / curves on the page",
        ha="center", va="top", fontsize=8.5, color="#444",
    )

    fig.suptitle(
        "Mollweide preserves area, not shape; sky circles become ovals/curves",
        fontsize=12.5, y=0.98, fontweight="normal",
    )
    fig.text(
        0.5, 0.01,
        r"Convention: lon_disp $=-\mathrm{wrap}(\ell)$; $\ell$ increases LEFT; center $=$ anti-CoE $(0^\circ,0^\circ)$",
        ha="center", fontsize=8, color="#555",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(OUT, dpi=140, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
