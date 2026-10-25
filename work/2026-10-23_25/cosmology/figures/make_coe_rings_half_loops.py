#!/usr/bin/env python3
"""
CoE rings — theta-ray-side half-loops only (Chris 2026-09-23 whiteout, CORRECTED).

Flags (CLI or module-level defaults):
  AXIS_TILT_B_DEG = 15   -> anti-CoE (0 deg, +15 deg), CoE (180 deg, -15 deg)
  ROLL_ABOUT_COE_DEG = 0 or 45  -> CCW about +D looking along +D
  HALF_LOOPS = True       -> open 180 deg arcs only (never full loops)

Half-cut (theta-ray side — NOT outer limb):
  foreshortened ex=mu+0.045*rho*cos(phi), ey=rho*sin(phi);
  keep ey*sin(theta_mid) >= 0  (TOP if sin theta>0, BOTTOM if sin theta<0).
  3D: keep (n·Y_hat)*sin(theta_mid) >= 0. Discard sign(mu)*cos(phi)>=0.

Outputs depend on roll:
  roll=0  -> coe_rings_half_loops_tilt15.{png,pdf,json}
  roll=45 -> coe_rings_half_loops_tilt15_rot45ccw.{png,pdf,json}
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle, Patch, Polygon
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from coe_half_loop_geometry import (
    BAND_FULL,
    BAND_HALF,
    HALF_LOOP_RULE,
    MIGKAS_MAXES,
    MIGKAS_MINS,
    STRIP_THETA_MAX,
    STRIP_THETA_MIN,
    axis_pose,
    band_edges,
    celestial_equator_lb,
    foreshortened_half_band_polygon,
    foreshortened_half_ellipse_xy,
    great_circle_through_poles,
    lb_to_n,
    lon_disp_from_ell,
    mu_pair_for_band,
    n_to_lb,
    small_circle_theta_ray_half,
)

OUT_DIR = Path(__file__).resolve().parent

# Defaults (overridden by CLI)
AXIS_TILT_B_DEG = 15.0
ROLL_ABOUT_COE_DEG = 0.0
HALF_LOOPS = True

C_COE = "#6b3fa0"
C_ANTI = "#8e1b1b"
C_MAX_30 = "#c0392b"
C_MAX_185 = "#ec7063"
C_MIN_115 = "#1a5276"
C_MIN_265 = "#5dade2"
C_AXIS = "#5b2c6f"
C_EQ = "#1e8449"
COLORS_MAX = [C_MAX_30, C_MAX_185]
COLORS_MIN = [C_MIN_115, C_MIN_265]


def output_stem(tilt_b: float, roll: float) -> str:
    base = f"coe_rings_half_loops_tilt{int(round(tilt_b))}"
    if abs(roll) > 0.5:
        base += f"_rot{int(round(roll))}ccw"
    return base


def split_lon_segments(lon_rad, lat_rad, gap=np.pi):
    d = np.diff(lon_rad)
    cuts = np.where(np.abs(d) > gap)[0] + 1
    idx = np.concatenate([[0], cuts, [len(lon_rad)]])
    segs = []
    for i0, i1 in zip(idx[:-1], idx[1:]):
        if i1 - i0 >= 2:
            segs.append((lon_rad[i0:i1], lat_rad[i0:i1]))
    return segs


def plot_lb_curve_mollweide(ax, l_deg, b_deg, **kwargs):
    lon = np.deg2rad(lon_disp_from_ell(l_deg))
    lat = np.deg2rad(b_deg)
    for xr, yr in split_lon_segments(lon, lat):
        ax.plot(xr, yr, **kwargs)


def _half_ring_lb(D, mu, theta_mid, roll_deg, n_pts=721):
    pts, _ = small_circle_theta_ray_half(
        D, mu, theta_mid, roll_deg=roll_deg, n=n_pts
    )
    rl = np.array([n_to_lb(p)[0] for p in pts])
    rb = np.array([n_to_lb(p)[1] for p in pts])
    return rl, rb


def _half_ring_lonlat_rad(D, mu, theta_mid, roll_deg, n_pts=361):
    pts, _ = small_circle_theta_ray_half(
        D, float(mu), theta_mid, roll_deg=roll_deg, n=n_pts
    )
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    lon = np.deg2rad(lon_disp_from_ell(ls))
    lat = np.deg2rad(bs)
    return lon, lat


def fill_mollweide_half_ribbon(ax, D, theta_lo, theta_hi, theta_mid, color, roll_deg, n_pts=721):
    """Theta-ray-side Mollweide ribbon between +/-9 deg edges + solid mid half-arc."""
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    lon_a, lat_a = _half_ring_lonlat_rad(D, mu_a, theta_mid, roll_deg, n_pts=n_pts)
    lon_b, lat_b = _half_ring_lonlat_rad(D, mu_b, theta_mid, roll_deg, n_pts=n_pts)
    gap = 0.9 * np.pi
    bad = np.zeros(n_pts - 1, dtype=bool)
    for i in range(n_pts - 1):
        if abs(lon_a[i + 1] - lon_a[i]) > gap:
            bad[i] = True
        elif abs(lon_b[i + 1] - lon_b[i]) > gap:
            bad[i] = True
        elif abs(lon_a[i] - lon_b[i]) > gap:
            bad[i] = True
        elif abs(lon_a[i + 1] - lon_b[i + 1]) > gap:
            bad[i] = True
    polys = []
    i = 0
    while i < n_pts - 1:
        while i < n_pts - 1 and bad[i]:
            i += 1
        if i >= n_pts - 1:
            break
        j = i
        while j < n_pts - 1 and not bad[j]:
            j += 1
        i0, i1 = i, j
        if i1 - i0 >= 1:
            xs = np.concatenate([lon_a[i0 : i1 + 1], lon_b[i0 : i1 + 1][::-1]])
            ys = np.concatenate([lat_a[i0 : i1 + 1], lat_b[i0 : i1 + 1][::-1]])
            polys.append(np.column_stack([xs, ys]))
        i = j + 1
    if polys:
        ax.add_collection(
            PolyCollection(
                polys, facecolors=color, edgecolors="none",
                alpha=0.26, zorder=2, linewidths=0,
            )
        )
    rl, rb = _half_ring_lb(D, mu_mid, theta_mid, roll_deg, n_pts=max(n_pts, 721))
    plot_lb_curve_mollweide(ax, rl, rb, color=color, lw=2.6, alpha=0.95, zorder=4)


def draw_panelA_half_ribbon(ax, theta_lo, theta_hi, theta_mid, color, roll_deg):
    mu_lo = float(np.cos(np.deg2rad(theta_lo)))
    mu_hi = float(np.cos(np.deg2rad(theta_hi)))
    poly = foreshortened_half_band_polygon(
        mu_lo, mu_hi, theta_mid, roll_deg=roll_deg, n=160
    )
    patch = Polygon(
        poly, closed=True, facecolor=color, edgecolor="none",
        alpha=0.28, zorder=2, linewidth=0,
    )
    ax.add_patch(patch)
    clip_circ = Circle((0.0, 0.0), 1.0, transform=ax.transData)
    patch.set_clip_path(clip_circ)

    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    ex_m, ey_m, _ = foreshortened_half_ellipse_xy(
        mu_mid, theta_mid, roll_deg=roll_deg, n=256
    )
    ax.plot(ex_m, ey_m, color=color, lw=2.8, alpha=0.95, zorder=4, solid_capstyle="round")


def draw_panelA_theta_rays(ax, geo):
    items = list(zip(geo["theta_max_from_plus_x_deg"], COLORS_MAX)) + list(
        zip(geo["theta_min_from_plus_x_deg"], COLORS_MIN)
    )
    for mid, color in items:
        t = float(np.deg2rad(mid))
        x1, y1 = float(np.cos(t)), float(np.sin(t))
        is_max30 = abs(mid - 30.0) < 0.5
        lw = 3.2 if is_max30 else 1.8
        alpha = 0.70 if is_max30 else 0.45
        ax.plot(
            [0.0, x1], [0.0, y1],
            color=color, lw=lw, alpha=alpha, zorder=3, solid_capstyle="round",
        )
        ax.plot(
            [x1], [y1], "o", color=color, ms=4.5 if is_max30 else 3.5,
            alpha=0.9, zorder=5, markeredgecolor="none",
        )
    arc_r = 0.38
    phi = np.linspace(0.0, np.deg2rad(30.0), 60)
    ax.plot(
        arc_r * np.cos(phi), arc_r * np.sin(phi),
        color=C_MAX_30, lw=1.6, alpha=0.85, zorder=3,
    )
    ax.text(
        arc_r * np.cos(np.deg2rad(15.0)) + 0.06,
        arc_r * np.sin(np.deg2rad(15.0)) + 0.04,
        r"$\theta$",
        color=C_MAX_30, fontsize=10, fontweight="bold",
        ha="left", va="bottom", zorder=6,
    )


def panel_A_model_frame(ax, geo, roll_deg):
    roll_note = (
        rf"; $+{roll_deg:.0f}^\circ$ CCW roll about CoE"
        if abs(roll_deg) > 0.5
        else ""
    )
    ax.set_title(
        r"A. Model frame — $\theta$-ray-side half-loops only"
        "\n"
        r"axisymmetric about $x$ = CoE$\leftrightarrow$anti-CoE"
        f"{roll_note}\n"
        r"keep $(n\cdot\hat y)\sin\theta\geq0$ (same side as $\theta$-ray); "
        r"each $18^\circ$ end-to-end ($\pm9^\circ$)",
        fontsize=9.5,
        pad=8,
    )
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color="#7f8c8d", lw=1.2, zorder=1)
    ax.plot(0, 0, "ko", ms=6, zorder=5)
    ax.text(0.08, -0.12, "observer", fontsize=8, color="k", ha="left")
    ax.annotate(
        "",
        xy=(1.15, 0),
        xytext=(-1.15, 0),
        arrowprops=dict(arrowstyle="<->", color=C_AXIS, lw=2.5),
        zorder=4,
    )
    ax.text(1.18, 0.08, "anti-CoE\n(+x)", color=C_ANTI, fontsize=8, ha="left", va="bottom")
    ax.text(-1.18, 0.08, "CoE\n(−x)", color=C_COE, fontsize=8, ha="right", va="bottom")

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        draw_panelA_half_ribbon(ax, lo, hi, mid, color, roll_deg)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        draw_panelA_half_ribbon(ax, lo, hi, mid, color, roll_deg)

    draw_panelA_theta_rays(ax, geo)

    label_r = 1.30
    for mid, color in list(zip(geo["theta_max_from_plus_x_deg"], COLORS_MAX)) + list(
        zip(geo["theta_min_from_plus_x_deg"], COLORS_MIN)
    ):
        t = np.deg2rad(mid)
        ax.text(
            label_r * np.cos(t),
            label_r * np.sin(t),
            rf"$\theta\!=\!{mid:.0f}^\circ$",
            color=color, fontsize=7.5, ha="center", va="center", fontweight="bold",
        )

    ax.text(
        0.0,
        -1.38,
        (
            r"$\theta$-ray-side half-loops  ·  keep $y\,\mathrm{sgn}(\sin\theta)\geq0$  ·  "
            r"each $18^\circ$ end-to-end ($\pm9^\circ$)"
            + "\n"
            + r"(after-2$\ell$ $A(\ell)$: two MAX + two MIN; open $180^\circ$ arcs — no full loops)"
        ),
        fontsize=7.5,
        ha="center",
        va="top",
        style="italic",
        color="#34495e",
    )
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")


def panel_B_mollweide(ax, geo, roll_deg):
    tilt = geo["tilt_b_deg"]
    roll_note = (
        "\n" + rf"$+{roll_deg:.0f}^\circ$ CCW about CoE"
        if abs(roll_deg) > 0.5
        else ""
    )
    ax.set_title(
        rf"B. Migkas $A(\ell)$ half-loops on sky (anti-CoE at $(0^\circ,+{tilt:.0f}^\circ)$)"
        "\n"
        rf"(Galactic Mollweide; CoE axis tilted $+{tilt:.0f}^\circ$; "
        rf"$\theta$-ray-side half-loops only; $\ell$ increases LEFT){roll_note}",
        fontsize=9.5,
        pad=8,
    )
    for b in [-60, -30, 0, 30, 60]:
        ell = np.linspace(0, 360, 361)
        plot_lb_curve_mollweide(ax, ell, np.full_like(ell, b), color="#bdc3c7", lw=0.5, zorder=0)
    for ell0 in range(0, 360, 30):
        b = np.linspace(-90, 90, 181)
        plot_lb_curve_mollweide(ax, np.full_like(b, ell0), b, color="#bdc3c7", lw=0.5, zorder=0)

    D = geo["D"]
    ls, bs = great_circle_through_poles(D, n=721)
    plot_lb_curve_mollweide(ax, ls, bs, color=C_AXIS, lw=2.0, zorder=2)

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        fill_mollweide_half_ribbon(ax, D, lo, hi, mid, color, roll_deg)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        fill_mollweide_half_ribbon(ax, D, lo, hi, mid, color, roll_deg)

    ell_eq, b_eq = celestial_equator_lb(721)
    plot_lb_curve_mollweide(ax, ell_eq, b_eq, color=C_EQ, lw=1.8, ls="--", zorder=2)

    ax.plot(
        np.deg2rad(lon_disp_from_ell(geo["anti_lb"][0])),
        np.deg2rad(geo["anti_lb"][1]),
        marker="o", ms=8, color=C_ANTI,
        markerfacecolor=C_ANTI, markeredgecolor="k", markeredgewidth=0.6,
        zorder=6, linestyle="None",
    )
    ax.plot(
        np.deg2rad(lon_disp_from_ell(geo["coe_lb"][0])),
        np.deg2rad(geo["coe_lb"][1]),
        marker="*", ms=14, color=C_COE,
        markerfacecolor=C_COE, markeredgecolor="k", markeredgewidth=0.6,
        zorder=6, linestyle="None",
    )
    for (l, b), color in zip(MIGKAS_MAXES, COLORS_MAX):
        ax.plot(
            np.deg2rad(lon_disp_from_ell(l)), np.deg2rad(b),
            marker="^", ms=10, color=color,
            markeredgecolor="k", markeredgewidth=0.7, zorder=7, linestyle="None",
        )
    for (l, b), color in zip(MIGKAS_MINS, COLORS_MIN):
        ax.plot(
            np.deg2rad(lon_disp_from_ell(l)), np.deg2rad(b),
            marker="s", ms=10, color=color,
            markeredgecolor="k", markeredgewidth=0.7, zorder=7, linestyle="None",
        )

    xticks = np.deg2rad(lon_disp_from_ell(np.array([180, 90, 0, 270])))
    ax.set_xticks(xticks)
    ax.set_xticklabels([r"$180^\circ$", r"$90^\circ$", r"$0^\circ$", r"$270^\circ$"], fontsize=8)
    ax.set_yticks(np.deg2rad([-60, -30, 0, 30, 60]))
    ax.set_yticklabels(
        [r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$+30^\circ$", r"$+60^\circ$"], fontsize=8
    )
    ax.grid(False)
    ax.set_xlabel(r"Galactic $\ell$ (increases LEFT)", fontsize=8)
    ax.set_ylabel(r"Galactic $b$", fontsize=8)


def fill_sphere_half_ribbon(ax, D, theta_lo, theta_hi, theta_mid, color, roll_deg, n_pts=160):
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    p0, _ = small_circle_theta_ray_half(
        D, float(mu_a), theta_mid, roll_deg=roll_deg, n=n_pts
    )
    p1, _ = small_circle_theta_ray_half(
        D, float(mu_b), theta_mid, roll_deg=roll_deg, n=n_pts
    )
    verts = [[p0[j], p0[j + 1], p1[j + 1], p1[j]] for j in range(n_pts - 1)]
    coll = Poly3DCollection(
        verts, facecolors=color, edgecolors="none", linewidths=0.0,
        alpha=0.25, zorder=2,
    )
    ax.add_collection3d(coll)
    pts, _ = small_circle_theta_ray_half(
        D, mu_mid, theta_mid, roll_deg=roll_deg, n=n_pts
    )
    ax.plot(
        pts[:, 0], pts[:, 1], pts[:, 2],
        color=color, lw=2.6, alpha=0.95, zorder=5, solid_capstyle="round",
    )


def draw_wire_sphere(ax, n_mer=12, n_par=6, **kwargs):
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(-np.pi / 2, np.pi / 2, 40)
    xs = np.outer(np.cos(u), np.cos(v))
    ys = np.outer(np.sin(u), np.cos(v))
    zs = np.outer(np.ones_like(u), np.sin(v))
    ax.plot_surface(xs, ys, zs, color="#dfe6e9", alpha=0.15, linewidth=0, shade=False)
    kw = dict(color="#95a5a6", lw=0.4, alpha=0.45)
    kw.update(kwargs)
    for i in range(n_mer):
        phi = 2 * np.pi * i / n_mer
        t = np.linspace(-np.pi / 2, np.pi / 2, 80)
        ax.plot(np.cos(t) * np.cos(phi), np.cos(t) * np.sin(phi), np.sin(t), **kw)
    for j in range(1, n_par):
        b = -np.pi / 2 + np.pi * j / n_par
        t = np.linspace(0, 2 * np.pi, 80)
        r = np.cos(b)
        ax.plot(r * np.cos(t), r * np.sin(t), np.full_like(t, np.sin(b)), **kw)


def panel_C_sphere(ax, geo, roll_deg):
    th_max = geo["theta_max_from_plus_x_deg"]
    th_min = geo["theta_min_from_plus_x_deg"]
    roll_note = (
        rf"; $+{roll_deg:.0f}^\circ$ CCW about CoE"
        if abs(roll_deg) > 0.5
        else ""
    )
    ax.set_title(
        rf"C. Sphere 3D — $\theta$-ray-side half-loops{roll_note}"
        "\n"
        r"(open $180^\circ$ arcs; celestial equator dashed)",
        fontsize=9.5,
        pad=6,
    )
    draw_wire_sphere(ax)
    D = geo["D"]
    ell_eq, b_eq = celestial_equator_lb(400)
    eq_xyz = np.array([lb_to_n(l, b) for l, b in zip(ell_eq, b_eq)])
    ax.plot(eq_xyz[:, 0], eq_xyz[:, 1], eq_xyz[:, 2], color=C_EQ, lw=1.8, ls="--", zorder=3)
    ax.plot(
        [-D[0], D[0]], [-D[1], D[1]], [-D[2], D[2]],
        color=C_AXIS, lw=3.0, zorder=5,
    )
    ax.scatter(*(-D), color=C_COE, s=120, marker="*", edgecolors="k", linewidths=0.6,
               depthshade=False, zorder=8)
    ax.scatter(*(D), color=C_ANTI, s=80, edgecolors="k", linewidths=0.6, depthshade=False, zorder=8)
    ax.text(*(-D * 1.22), "CoE", color=C_COE, fontsize=8, fontweight="bold")
    ax.text(*(D * 1.18), "anti-CoE", color=C_ANTI, fontsize=8, fontweight="bold")
    ax.scatter([0], [0], [0], color="k", s=25, depthshade=False, zorder=9)

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        fill_sphere_half_ribbon(ax, D, lo, hi, mid, color, roll_deg)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        fill_sphere_half_ribbon(ax, D, lo, hi, mid, color, roll_deg)

    for n, color in zip(geo["n_maxes"], COLORS_MAX):
        ax.scatter(
            [n[0]], [n[1]], [n[2]],
            c=[color], marker="^", s=70,
            edgecolors="k", linewidths=0.7, depthshade=False, zorder=10,
        )
        ax.plot([0, n[0]], [0, n[1]], [0, n[2]], color=color, lw=0.8, ls=":", alpha=0.7)
    for n, color in zip(geo["n_mins"], COLORS_MIN):
        ax.scatter(
            [n[0]], [n[1]], [n[2]],
            c=[color], marker="s", s=60,
            edgecolors="k", linewidths=0.7, depthshade=False, zorder=10,
        )
        ax.plot([0, n[0]], [0, n[1]], [0, n[2]], color=color, lw=0.8, ls=":", alpha=0.7)

    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.set_zlim(-1.25, 1.25)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=18, azim=-55)
    ax.set_axis_off()
    ax.text2D(
        0.5, -0.02,
        rf"MAX $\theta\!=\!{th_max[0]:.0f}^\circ,{th_max[1]:.0f}^\circ$; "
        rf"MIN $\theta\!=\!{th_min[0]:.0f}^\circ,{th_min[1]:.0f}^\circ$  "
        rf"($\theta$-ray-side half-loops; $18^\circ$ end-to-end = $\pm9^\circ$)",
        transform=ax.transAxes, ha="center", fontsize=7.5, color="#34495e",
    )


def compute_geometry(tilt_b, roll_deg):
    pose = axis_pose(tilt_b)
    th_max = list(STRIP_THETA_MAX)
    th_min = list(STRIP_THETA_MIN)
    max_bands = band_edges(th_max)
    min_bands = band_edges(th_min)
    n_maxes = [lb_to_n(*lb) for lb in MIGKAS_MAXES]
    n_mins = [lb_to_n(*lb) for lb in MIGKAS_MINS]
    sep = float(
        np.rad2deg(np.arccos(np.clip(np.dot(n_maxes[0], n_mins[1]), -1.0, 1.0)))
    )
    return {
        **pose,
        "n_maxes": n_maxes,
        "n_mins": n_mins,
        "theta_max_from_plus_x_deg": th_max,
        "theta_min_from_plus_x_deg": th_min,
        "max_bands_deg": max_bands,
        "min_bands_deg": min_bands,
        "mu_max_bands": [mu_pair_for_band(lo, hi) for lo, hi in max_bands],
        "mu_min_bands": [mu_pair_for_band(lo, hi) for lo, hi in min_bands],
        "sep_deg": sep,
        "band_half_width_deg": BAND_HALF,
        "roll_about_coe_deg": float(roll_deg),
        "half_loops": True,
        "half_loop_rule": HALF_LOOP_RULE,
    }


def main(tilt_b=None, roll_deg=None):
    tilt_b = float(AXIS_TILT_B_DEG if tilt_b is None else tilt_b)
    roll_deg = float(ROLL_ABOUT_COE_DEG if roll_deg is None else roll_deg)
    assert HALF_LOOPS

    geo = compute_geometry(tilt_b, roll_deg)
    th_max = geo["theta_max_from_plus_x_deg"]
    th_min = geo["theta_min_from_plus_x_deg"]
    max_bands = geo["max_bands_deg"]
    min_bands = geo["min_bands_deg"]

    stem = output_stem(tilt_b, roll_deg)
    png_path = OUT_DIR / f"{stem}.png"
    pdf_path = OUT_DIR / f"{stem}.pdf"
    json_path = OUT_DIR / f"{stem}.json"

    roll_sup = (
        rf"; $+{roll_deg:.0f}^\circ$ CCW about CoE axis"
        if abs(roll_deg) > 0.5
        else ""
    )
    fig = plt.figure(figsize=(14.5, 5.6), dpi=160)
    fig.suptitle(
        r"Pure CoE axisymmetry: $\theta$-ray-side half-loops "
        r"(keep $(n\cdot\hat y)\sin\theta\geq0$)"
        "\n"
        rf"(each $18^\circ$ end-to-end = $\pm9^\circ$; "
        rf"MAX $\theta\!=\!{th_max[0]:.0f}^\circ/{th_max[1]:.0f}^\circ$; "
        rf"MIN $\theta\!=\!{th_min[0]:.0f}^\circ/{th_min[1]:.0f}^\circ$; "
        rf"anti-CoE $(0^\circ,+{tilt_b:.0f}^\circ)$, CoE axis tilted $+{tilt_b:.0f}^\circ$"
        f"{roll_sup})",
        fontsize=11,
        fontweight="bold",
        y=1.02,
    )

    axA = fig.add_subplot(1, 3, 1)
    panel_A_model_frame(axA, geo, 0.0)  # model frame never rolls
    axB = fig.add_subplot(1, 3, 2, projection="mollweide")
    panel_B_mollweide(axB, geo, roll_deg)
    axC = fig.add_subplot(1, 3, 3, projection="3d")
    panel_C_sphere(axC, geo, roll_deg)

    handles = [
        Line2D([0], [0], color=C_AXIS, lw=2.5, label=r"CoE$\leftrightarrow$anti-CoE axis"),
        Patch(
            facecolor=C_MAX_30, alpha=0.45, edgecolor=C_MAX_30,
            label=rf"MAX $\theta\!=\!{th_max[0]:.0f}^\circ$ $\theta$-ray half "
            rf"[{max_bands[0][0]:.0f},{max_bands[0][1]:.0f}]$^\circ$",
        ),
        Patch(
            facecolor=C_MAX_185, alpha=0.45, edgecolor=C_MAX_185,
            label=rf"MAX $\theta\!=\!{th_max[1]:.0f}^\circ$ $\theta$-ray half "
            rf"[{max_bands[1][0]:.0f},{max_bands[1][1]:.0f}]$^\circ$",
        ),
        Patch(
            facecolor=C_MIN_115, alpha=0.45, edgecolor=C_MIN_115,
            label=rf"MIN $\theta\!=\!{th_min[0]:.0f}^\circ$ $\theta$-ray half "
            rf"[{min_bands[0][0]:.0f},{min_bands[0][1]:.0f}]$^\circ$",
        ),
        Patch(
            facecolor=C_MIN_265, alpha=0.45, edgecolor=C_MIN_265,
            label=rf"MIN $\theta\!=\!{th_min[1]:.0f}^\circ$ $\theta$-ray half "
            rf"[{min_bands[1][0]:.0f},{min_bands[1][1]:.0f}]$^\circ$",
        ),
        Line2D([0], [0], color=C_EQ, lw=1.8, ls="--", label="celestial equator (Dec$=0$)"),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_30, markeredgecolor="k",
            ms=9, label=rf"MAX $({MIGKAS_MAXES[0][0]:.0f}^\circ,{MIGKAS_MAXES[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_185, markeredgecolor="k",
            ms=9, label=rf"MAX $({MIGKAS_MAXES[1][0]:.0f}^\circ,{MIGKAS_MAXES[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_115, markeredgecolor="k",
            ms=8, label=rf"MIN $({MIGKAS_MINS[0][0]:.0f}^\circ,{MIGKAS_MINS[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_265, markeredgecolor="k",
            ms=8, label=rf"MIN $({MIGKAS_MINS[1][0]:.0f}^\circ,{MIGKAS_MINS[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="o", color="w", markerfacecolor=C_ANTI, markeredgecolor="k",
            ms=8, label=rf"anti-CoE $({geo['anti_lb'][0]:.0f}^\circ,{geo['anti_lb'][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="*", color="w", markerfacecolor=C_COE, markeredgecolor="k",
            ms=11, label=rf"CoE $({geo['coe_lb'][0]:.0f}^\circ,{geo['coe_lb'][1]:+.0f}^\circ)$",
        ),
    ]
    fig.legend(
        handles=handles, loc="lower center", ncol=6, fontsize=7.0,
        frameon=True, fancybox=True, bbox_to_anchor=(0.5, -0.10),
    )
    fig.tight_layout()
    fig.savefig(png_path, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    payload = {
        "png": str(png_path),
        "pdf": str(pdf_path),
        "script": str(OUT_DIR / "make_coe_rings_half_loops.py"),
        "geometry_module": str(OUT_DIR / "coe_half_loop_geometry.py"),
        "HALF_LOOPS": True,
        "half_loop_rule": HALF_LOOP_RULE,
        "axis_tilt_b_deg": tilt_b,
        "roll_about_coe_deg": roll_deg,
        "migkas_max_lb_deg": [list(lb) for lb in MIGKAS_MAXES],
        "migkas_min_lb_deg": [list(lb) for lb in MIGKAS_MINS],
        "strip_theta_max_deg": list(STRIP_THETA_MAX),
        "strip_theta_min_deg": list(STRIP_THETA_MIN),
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "max_bands_deg": [list(b) for b in max_bands],
        "min_bands_deg": [list(b) for b in min_bands],
        "CoE_lb_deg": [round(geo["coe_lb"][0], 4), round(geo["coe_lb"][1], 4)],
        "anti_CoE_lb_deg": [round(geo["anti_lb"][0], 4), round(geo["anti_lb"][1], 4)],
        "D_hat_galactic_xyz": geo["D"].tolist(),
        "angular_sep_migkas_max0_min1_deg": round(geo["sep_deg"], 6),
        "orientation": (
            f"anti-CoE at (0°, +{tilt_b:.0f}°), CoE at (180°, −{tilt_b:.0f}°); "
            f"theta-ray-side half-loops (n·Y)*sin(theta)>=0; roll about CoE = {roll_deg:.0f}° CCW"
        ),
        "plot_convention": "Migkas Galactic Mollweide: lon_disp=-ell_wrapped (ell increases LEFT)",
        "colors": {
            "MAX_30": C_MAX_30,
            "MAX_185": C_MAX_185,
            "MIN_115": C_MIN_115,
            "MIN_265": C_MIN_265,
        },
    }
    json_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    print(f"\nWrote {png_path}")
    print(f"Wrote {pdf_path}")
    print(f"Wrote {json_path}")
    return png_path, pdf_path, json_path


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tilt", type=float, default=AXIS_TILT_B_DEG)
    p.add_argument("--roll", type=float, default=ROLL_ABOUT_COE_DEG,
                   help="CCW degrees about +D (0 or 45)")
    args = p.parse_args()
    main(tilt_b=args.tilt, roll_deg=args.roll)
