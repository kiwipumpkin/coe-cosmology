#!/usr/bin/env python3
"""
Migkas Fig. 8 overlay — FULL CoE loops with opacity gradient (Chris 2026-09-25).

Derived from make_migkas_fig8_half_loops_overlay.py (left unchanged).
Each FULL ring (±9° band fill + solid mid line with its white underlay) is drawn
with alpha multiplied by (1+cos psi)/2, psi = ring azimuth from the θ-ray point
(ring point most along sign(sinθ)·Ŷ, Ŷ from ring_basis rolled by roll_deg) —
same helper as make_coe_rings_gradient_loops.py (coe_gradient_loop_geometry).
Contiguous elements sharing a quantized alpha level (64 levels) are merged into
one strip polygon / polyline (no beading, few AA seams).

GOOD_REG, lon_disp = -wrap(ℓ), oval clipping, antimeridian splitting, colours,
markers all unchanged. Published map untouched (drawn underneath only).

Outputs:
  roll=0  → migkas_fig8_gradient_loops_tilt15.{png,json}
  roll=45 → migkas_fig8_gradient_loops_tilt15_rot45ccw.{png,json}
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
from PIL import Image

from coe_half_loop_geometry import (
    BAND_FULL,
    BAND_HALF,
    MIGKAS_MAXES,
    MIGKAS_MINS,
    STRIP_THETA_MAX,
    STRIP_THETA_MIN,
    axis_pose,
    band_edges,
    celestial_equator_lb,
    great_circle_through_poles,
    lon_disp_from_ell,
    mu_pair_for_band,
    n_to_lb,
)

from coe_gradient_loop_geometry import (
    GRADIENT_RULE,
    quantized_runs,
    rgba_array,
    segment_scale,
    small_circle_full_gradient,
)

OUT_DIR = Path(__file__).resolve().parent

AXIS_TILT_B_DEG = 15.0
ROLL_ABOUT_COE_DEG = 0.0
FULL_GRADIENT_LOOPS = True
FILL_ALPHA = 0.30
N_ALPHA_LEVELS = 64

SOURCE_CANDIDATES = [
    Path("/workspace/cosmology/map_compare/FULL_migkas_theta75.png"),
    Path("/workspace/cosmology/map_compare/migkas_fig8_theta75_crop.png"),
    Path("/workspace/cosmology/migkas/pub_maps/Fig8_full_hires.png"),
]

THETA_MAX = list(STRIP_THETA_MAX)
THETA_MIN = list(STRIP_THETA_MIN)
MAX_BANDS = band_edges(THETA_MAX)
MIN_BANDS = band_edges(THETA_MIN)

GOOD_REG = {
    "cx": 863.0,
    "cy": 763.25,
    "Rx": 693.9927,
    "Ry": 392.542425,
}

C_MAX_30 = "#00e5ff"
C_MAX_185 = "#ff2d95"
C_MIN_115 = "#b8f200"
C_MIN_265 = "#4fc3f7"
COLORS_MAX = [C_MAX_30, C_MAX_185]
COLORS_MIN = [C_MIN_115, C_MIN_265]
C_COE = "#6b3fa0"
C_ANTI = "#8e1b1b"
C_AXIS = "#5b2c6f"


def output_stem(tilt_b, roll):
    base = f"migkas_fig8_gradient_loops_tilt{int(round(tilt_b))}"
    if abs(roll) > 0.5:
        base += f"_rot{int(round(roll))}ccw"
    return base


def mollweide_xy(lon_rad, lat_rad, n_iter=64):
    lon = np.asarray(lon_rad, dtype=float)
    lat = np.asarray(lat_rad, dtype=float)
    target = np.pi * np.sin(np.clip(lat, -np.pi / 2, np.pi / 2))
    lo = np.full_like(lat, -np.pi / 2, dtype=float)
    hi = np.full_like(lat, np.pi / 2, dtype=float)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        g = 2.0 * mid + np.sin(2.0 * mid)
        lo = np.where(g < target, mid, lo)
        hi = np.where(g >= target, mid, hi)
    theta = 0.5 * (lo + hi)
    theta = np.where(lat >= np.pi / 2 - 1e-15, np.pi / 2, theta)
    theta = np.where(lat <= -np.pi / 2 + 1e-15, -np.pi / 2, theta)
    x = (2.0 * np.sqrt(2.0) / np.pi) * lon * np.cos(theta)
    y = np.sqrt(2.0) * np.sin(theta)
    return x, y


def lb_to_px(ell_deg, b_deg, cx, cy, Rx, Ry):
    lon = np.deg2rad(lon_disp_from_ell(np.asarray(ell_deg, dtype=float)))
    lat = np.deg2rad(np.asarray(b_deg, dtype=float))
    x, y = mollweide_xy(lon, lat)
    Xn = x / (2.0 * np.sqrt(2.0))
    Yn = y / np.sqrt(2.0)
    return cx + Xn * Rx, cy - Yn * Ry


def as_xy(px, py):
    return float(np.asarray(px).reshape(-1)[0]), float(np.asarray(py).reshape(-1)[0])


def chamfer_distance_to_dark(dark):
    h, w = dark.shape
    INF = 1e6
    D = np.where(dark, 0.0, INF)
    for y in range(h):
        for x in range(w):
            v = D[y, x]
            if x > 0:
                v = min(v, D[y, x - 1] + 1.0)
            if y > 0:
                v = min(v, D[y - 1, x] + 1.0)
            if x > 0 and y > 0:
                v = min(v, D[y - 1, x - 1] + 1.414)
            if x + 1 < w and y > 0:
                v = min(v, D[y - 1, x + 1] + 1.414)
            D[y, x] = v
    for y in range(h - 1, -1, -1):
        for x in range(w - 1, -1, -1):
            v = D[y, x]
            if x + 1 < w:
                v = min(v, D[y, x + 1] + 1.0)
            if y + 1 < h:
                v = min(v, D[y + 1, x] + 1.0)
            if x + 1 < w and y + 1 < h:
                v = min(v, D[y + 1, x + 1] + 1.414)
            if x > 0 and y + 1 < h:
                v = min(v, D[y + 1, x - 1] + 1.414)
            D[y, x] = v
    return D


def detect_oval_bbox(rgb):
    h, w, _ = rgb.shape
    white = (rgb[:, :, 0] > 250) & (rgb[:, :, 1] > 250) & (rgb[:, :, 2] > 250)
    oval_left = np.full(h, np.nan)
    oval_right = np.full(h, np.nan)
    for y in range(h):
        row = ~white[y]
        d = np.diff(np.r_[0, row.astype(int), 0])
        starts = np.where(d == 1)[0]
        ends = np.where(d == -1)[0] - 1
        if len(starts) == 0:
            continue
        lengths = ends - starts + 1
        candidates = [
            (int(lengths[i]), int(starts[i]), int(ends[i]))
            for i in range(len(starts))
            if starts[i] < 0.72 * w and lengths[i] > 40
        ]
        if not candidates:
            continue
        candidates.sort(reverse=True)
        _, L, R = candidates[0]
        oval_left[y] = L
        oval_right[y] = R
    widths = np.where(np.isfinite(oval_left), oval_right - oval_left, 0.0)
    yeq = int(np.argmax(widths))
    thresh = 0.08 * widths[yeq]
    y0 = yeq
    while y0 > 0 and widths[y0 - 1] >= thresh:
        y0 -= 1
    y1 = yeq
    while y1 < h - 1 and widths[y1 + 1] >= thresh:
        y1 += 1
    x0 = float(oval_left[yeq])
    x1 = float(oval_right[yeq])
    return {
        "x0": int(x0), "x1": int(x1), "y0": int(y0), "y1": int(y1),
        "yeq": int(yeq),
        "cx": 0.5 * (x0 + x1), "cy": 0.5 * (y0 + y1),
        "Rx": 0.5 * (x1 - x0), "Ry": 0.5 * (y1 - y0 + 1),
    }


def choose_source():
    for p in SOURCE_CANDIDATES:
        if p.is_file():
            return p
    raise FileNotFoundError("No Migkas Fig.8 theta=75 source image found")


def sample_rgb(rgb, ell, b, cx, cy, Rx, Ry, r=4):
    px, py = lb_to_px(ell, b, cx, cy, Rx, Ry)
    xi, yi = int(round(as_xy(px, py)[0])), int(round(as_xy(px, py)[1]))
    patch = rgb[max(0, yi - r) : yi + r + 1, max(0, xi - r) : xi + r + 1]
    return patch.mean(axis=(0, 1))


def validate_landmarks(rgb, reg):
    cx, cy, Rx, Ry = reg["cx"], reg["cy"], reg["Rx"], reg["Ry"]
    pts = {}
    landmarks = [
        ("center", 0.0, 0.0),
        ("ell90", 90.0, 0.0),
        ("ell270", 270.0, 0.0),
        ("north", 0.0, 90.0),
        ("south", 0.0, -90.0),
    ]
    for i, (ell, b) in enumerate(MIGKAS_MAXES):
        landmarks.append((f"migkas_max{i}", ell, b))
    for i, (ell, b) in enumerate(MIGKAS_MINS):
        landmarks.append((f"migkas_min{i}", ell, b))
    for name, ell, b in landmarks:
        x, y = as_xy(*lb_to_px(ell, b, cx, cy, Rx, Ry))
        pts[name] = {"ell": ell, "b": b, "px": x, "py": y}
    assert pts["ell90"]["px"] < pts["center"]["px"], "l=90 must be LEFT of center"
    assert pts["ell270"]["px"] > pts["center"]["px"], "l=270 must be RIGHT of center"
    assert abs(pts["center"]["px"] - cx) < 1.5
    assert abs(pts["center"]["py"] - cy) < 1.5
    assert pts["north"]["py"] < pts["center"]["py"]
    assert pts["south"]["py"] > pts["center"]["py"]
    cmax = sample_rgb(rgb, *MIGKAS_MAXES[0], cx, cy, Rx, Ry)
    cmin = sample_rgb(rgb, *MIGKAS_MINS[1], cx, cy, Rx, Ry)
    assert (cmax[0] + cmax[1]) > (cmin[0] + cmin[1]) + 50
    return {
        "landmarks_px": pts,
        "migkas_max0_rgb_mean": [float(x) for x in cmax],
        "migkas_min1_rgb_mean": [float(x) for x in cmin],
    }


def inside_ellipse(px, py, ecx, ecy, eRx, eRy, r_lim=1.0):
    px = np.asarray(px, dtype=float)
    py = np.asarray(py, dtype=float)
    Xn = (px - ecx) / eRx
    Yn = (py - ecy) / eRy
    return Xn * Xn + Yn * Yn <= float(r_lim) * float(r_lim)


def plot_lb_curve(
    ax, ell, b, cx, cy, Rx, Ry, inset=1.0, jump_px=45.0,
    clip_r2=None, clip_path=None, color_oval=None, **kw
):
    px, py = lb_to_px(ell, b, cx, cy, Rx, Ry)
    px = np.asarray(px, dtype=float).reshape(-1)
    py = np.asarray(py, dtype=float).reshape(-1)
    px = cx + (px - cx) * inset
    py = cy + (py - cy) * inset
    Xn = (px - cx) / Rx
    Yn = (cy - py) / Ry
    r2_lim = (inset * inset) if clip_r2 is None else float(clip_r2)
    inside = Xn * Xn + Yn * Yn <= r2_lim
    if color_oval is not None:
        ocx, ocy, oRx, oRy, o_r = color_oval
        inside = inside & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)
    px = np.where(inside, px, np.nan)
    py = np.where(inside, py, np.nan)
    segs = []
    cur_x, cur_y = [], []
    for x, y in zip(px, py):
        if not (np.isfinite(x) and np.isfinite(y)):
            if len(cur_x) >= 2:
                segs.append((cur_x, cur_y))
            cur_x, cur_y = [], []
            continue
        if cur_x and np.hypot(x - cur_x[-1], y - cur_y[-1]) > jump_px:
            if len(cur_x) >= 2:
                segs.append((cur_x, cur_y))
            cur_x, cur_y = [float(x)], [float(y)]
        else:
            cur_x.append(float(x))
            cur_y.append(float(y))
    if len(cur_x) >= 2:
        segs.append((cur_x, cur_y))
    if clip_path is not None:
        kw["clip_path"] = clip_path
        kw["clip_on"] = True
    for sx, sy in segs:
        (line,) = ax.plot(sx, sy, **kw)
        if clip_path is not None:
            line.set_clip_path(clip_path)
            line.set_clip_on(True)


def _full_ring_px_arrays(D, mu, theta_mid, roll_deg, cx, cy, Rx, Ry, inset, n_pts=721):
    """FULL ring → pixel arrays + per-point alpha scale (1+cos psi)/2.
    Sampled psi = -180°..+180° so the only closure point is the invisible one."""
    pts, _, scale = small_circle_full_gradient(
        D, float(mu), theta_mid, roll_deg=roll_deg, n=n_pts
    )
    ell = np.empty(n_pts)
    b = np.empty(n_pts)
    for i, p in enumerate(pts):
        ell[i], b[i] = n_to_lb(p)
    lon = np.deg2rad(lon_disp_from_ell(ell))
    lat = np.deg2rad(b)
    x_m, y_m = mollweide_xy(lon, lat)
    Xn = x_m / (2.0 * np.sqrt(2.0))
    Yn = y_m / np.sqrt(2.0)
    px = cx + Xn * Rx
    py = cy - Yn * Ry
    px = cx + (px - cx) * inset
    py = cy + (py - cy) * inset
    return lon, px, py, scale


def plot_full_ring_gradient(
    ax, D, mu, theta_mid, roll_deg, cx, cy, Rx, Ry, color, lw, inset=1.0, jump_px=80.0,
    clip_path=None, color_oval=None,
):
    """Stroke the FULL mid ring (white underlay + colour) with opacity gradient.
    Segments crossing the antimeridian / jumping / leaving the oval are dropped."""
    lon, px, py, sc = _full_ring_px_arrays(
        D, mu, theta_mid, roll_deg, cx, cy, Rx, Ry, inset, n_pts=2001
    )
    r2_lim = inset * inset + 1e-10
    gap = 0.9 * np.pi
    n = len(lon)
    ok_pt = ((px - cx) / Rx) ** 2 + ((cy - py) / Ry) ** 2 <= r2_lim
    if color_oval is not None:
        ocx, ocy, oRx, oRy, o_r = color_oval
        ok_pt = ok_pt & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)
    valid = (
        (np.abs(np.diff(lon)) <= gap)
        & (np.hypot(np.diff(px), np.diff(py)) <= jump_px)
        & ok_pt[:-1] & ok_pt[1:]
    )
    runs = quantized_runs(segment_scale(sc), valid, n_levels=N_ALPHA_LEVELS)
    if not runs:
        return
    lines = [np.column_stack([px[i0:i1 + 1], py[i0:i1 + 1]]) for i0, i1, _ in runs]
    levels = [r[2] for r in runs]
    for col, width, z, a in (("white", lw + 1.5, 3, 0.95), (color, lw, 4, 1.0)):
        lc = LineCollection(
            lines, colors=rgba_array(col, [a * L for L in levels]), linewidths=width,
            capstyle="butt", joinstyle="round", zorder=z,
        )
        ax.add_collection(lc)
        if clip_path is not None:
            lc.set_clip_path(clip_path)
            lc.set_clip_on(True)


def fill_full_gradient_band(
    ax, D, theta_lo, theta_hi, theta_mid, roll_deg, cx, cy, Rx, Ry, color,
    inset=1.0, clip_path=None, color_oval=None, n_pts=721, jump_px=80.0,
):
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    lon_a, px_a, py_a, sc = _full_ring_px_arrays(
        D, mu_a, theta_mid, roll_deg, cx, cy, Rx, Ry, inset, n_pts=n_pts
    )
    lon_b, px_b, py_b, _ = _full_ring_px_arrays(
        D, mu_b, theta_mid, roll_deg, cx, cy, Rx, Ry, inset, n_pts=n_pts
    )
    r2_lim = inset * inset + 1e-10

    def _inside(px, py):
        ok = ((px - cx) / Rx) ** 2 + ((cy - py) / Ry) ** 2 <= r2_lim
        if color_oval is not None:
            ocx, ocy, oRx, oRy, o_r = color_oval
            ok = ok & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)
        return ok

    gap = 0.9 * np.pi
    ia = _inside(px_a, py_a)
    ib = _inside(px_b, py_b)
    valid = (
        (np.abs(np.diff(lon_a)) <= gap)
        & (np.abs(np.diff(lon_b)) <= gap)
        & (np.abs(lon_a - lon_b)[:-1] <= gap)
        & (np.abs(lon_a - lon_b)[1:] <= gap)
        & (np.hypot(np.diff(px_a), np.diff(py_a)) <= jump_px)
        & (np.hypot(np.diff(px_b), np.diff(py_b)) <= jump_px)
        & ia[:-1] & ib[:-1] & ia[1:] & ib[1:]
    )
    runs = quantized_runs(segment_scale(sc), valid, n_levels=N_ALPHA_LEVELS)
    polys = [
        np.column_stack([
            np.concatenate([px_a[i0:i1 + 1], px_b[i0:i1 + 1][::-1]]),
            np.concatenate([py_a[i0:i1 + 1], py_b[i0:i1 + 1][::-1]]),
        ])
        for i0, i1, _ in runs
    ]
    if polys:
        coll = PolyCollection(
            polys, facecolors=rgba_array(color, [FILL_ALPHA * r[2] for r in runs]),
            edgecolors="none", zorder=2, linewidths=0,
        )
        ax.add_collection(coll)
        if clip_path is not None:
            coll.set_clip_path(clip_path)
            coll.set_clip_on(True)

    plot_full_ring_gradient(
        ax, D, mu_mid, theta_mid, roll_deg, cx, cy, Rx, Ry, color, lw=2.8,
        inset=inset, clip_path=clip_path, color_oval=color_oval,
    )


def draw_overlay(rgb, reg, seed, source, pose, roll_deg, png_path):
    h, w, _ = rgb.shape
    cx, cy, Rx, Ry = reg["cx"], reg["cy"], reg["Rx"], reg["Ry"]
    scx = float(seed["cx"])
    scy = float(seed["cy"])
    sRx = float(seed["Rx"])
    sRy = float(seed["Ry"])
    D = pose["D"]
    anti_lb = pose["anti_lb"]
    coe_lb = pose["coe_lb"]
    tilt = pose["tilt_b_deg"]

    dpi = 120
    fig, ax = plt.subplots(figsize=(w / dpi, (h / dpi) + 1.35), dpi=dpi)
    ax.imshow(rgb, origin="upper", extent=(0, w, h, 0), interpolation="nearest")
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    clip_oval = Ellipse(
        (cx, cy), width=2.0 * Rx, height=2.0 * Ry,
        facecolor="none", edgecolor="none", transform=ax.transData,
    )
    ax.add_patch(clip_oval)
    STROKE_INSET = 0.97
    COLOR_R = 0.978
    color_oval = (scx, scy, sRx, sRy, COLOR_R)

    def stroke_curve_geom(ell, b, color, lw, ls="-"):
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=1.0, clip_path=clip_oval,
            color_oval=color_oval,
            color="white", lw=lw + 1.6, ls=ls, alpha=0.95,
            solid_capstyle="round", zorder=3,
        )
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=1.0, clip_path=clip_oval,
            color_oval=color_oval,
            color=color, lw=lw, ls=ls, alpha=1.0,
            solid_capstyle="round", zorder=4,
        )

    def mark_point(ell, b, marker, color, ms):
        px, py = as_xy(*lb_to_px(ell, b, cx, cy, Rx, Ry))
        halo_ms = ms + 1.8
        ms_px = 0.5 * halo_ms * (dpi / 72.0)
        margin = ms_px / min(sRx, sRy) + 0.012
        r_lim = max(0.55, COLOR_R - margin)
        Xn = (px - scx) / sRx
        Yn = (py - scy) / sRy
        r_col = float(np.hypot(Xn, Yn))
        if r_col > 1e-9:
            if r_col > r_lim:
                scale = r_lim / r_col
                px = scx + (px - scx) * scale
                py = scy + (py - scy) * scale
                Xn = (px - scx) / sRx
                Yn = (py - scy) / sRy
                r_col = float(np.hypot(Xn, Yn))
            mild = min(STROKE_INSET, 0.995)
            px = scx + (px - scx) * mild
            py = scy + (py - scy) * mild
            Xn = (px - scx) / sRx
            Yn = (py - scy) / sRy
            r_col = float(np.hypot(Xn, Yn))
        if r_col + margin > COLOR_R + 1e-6:
            return

        def _mark(m_size, m_color, m_edge, m_ew, z):
            (line,) = ax.plot(
                px, py, marker=marker, ms=m_size, color=m_color,
                markeredgecolor=m_edge, markeredgewidth=m_ew, linestyle="None",
                zorder=z, clip_path=clip_oval, clip_on=True,
            )
            line.set_clip_path(clip_oval)
            line.set_clip_on(True)

        _mark(ms + 1.8, "white", "white", 0, 6)
        _mark(ms, color, "k", 0.55, 7)

    _ax_ell, _ax_b = great_circle_through_poles(D, n=721)
    stroke_curve_geom(_ax_ell, _ax_b, C_AXIS, lw=1.15)

    for (lo, hi), mid, color in zip(MAX_BANDS, THETA_MAX, COLORS_MAX):
        fill_full_gradient_band(
            ax, D, lo, hi, mid, roll_deg, cx, cy, Rx, Ry, color,
            inset=STROKE_INSET, clip_path=clip_oval, color_oval=color_oval,
        )
    for (lo, hi), mid, color in zip(MIN_BANDS, THETA_MIN, COLORS_MIN):
        fill_full_gradient_band(
            ax, D, lo, hi, mid, roll_deg, cx, cy, Rx, Ry, color,
            inset=STROKE_INSET, clip_path=clip_oval, color_oval=color_oval,
        )

    mark_point(anti_lb[0], anti_lb[1], "o", C_ANTI, ms=5.5)
    mark_point(coe_lb[0], coe_lb[1], "*", C_COE, ms=7.5)
    for (ell, b), color in zip(MIGKAS_MAXES, COLORS_MAX):
        mark_point(ell, b, "^", color, ms=6.5)
    for (ell, b), color in zip(MIGKAS_MINS, COLORS_MIN):
        mark_point(ell, b, "s", color, ms=6.0)

    roll_lab = (
        rf"; $+{roll_deg:.0f}^\circ$ CCW about CoE"
        if abs(roll_deg) > 0.5
        else ""
    )
    handles = [
        Patch(
            facecolor=C_MAX_30, alpha=0.55, edgecolor=C_MAX_30,
            label=rf"MAX $\theta\!=\!{THETA_MAX[0]:.0f}^\circ$ full loop, gradient "
            rf"[{MAX_BANDS[0][0]:.0f},{MAX_BANDS[0][1]:.0f}]$^\circ$",
        ),
        Patch(
            facecolor=C_MAX_185, alpha=0.55, edgecolor=C_MAX_185,
            label=rf"MAX $\theta\!=\!{THETA_MAX[1]:.0f}^\circ$ full loop, gradient "
            rf"[{MAX_BANDS[1][0]:.0f},{MAX_BANDS[1][1]:.0f}]$^\circ$",
        ),
        Patch(
            facecolor=C_MIN_115, alpha=0.55, edgecolor=C_MIN_115,
            label=rf"MIN $\theta\!=\!{THETA_MIN[0]:.0f}^\circ$ full loop, gradient "
            rf"[{MIN_BANDS[0][0]:.0f},{MIN_BANDS[0][1]:.0f}]$^\circ$",
        ),
        Patch(
            facecolor=C_MIN_265, alpha=0.55, edgecolor=C_MIN_265,
            label=rf"MIN $\theta\!=\!{THETA_MIN[1]:.0f}^\circ$ full loop, gradient "
            rf"[{MIN_BANDS[1][0]:.0f},{MIN_BANDS[1][1]:.0f}]$^\circ$",
        ),
        Line2D(
            [0], [0], color=C_AXIS, lw=2.0,
            label=rf"CoE axis tilted $+{tilt:.0f}^\circ$",
        ),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_30,
            markeredgecolor="k", linestyle="None", ms=8,
            label=rf"MAX $({MIGKAS_MAXES[0][0]:.0f}^\circ,{MIGKAS_MAXES[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_185,
            markeredgecolor="k", linestyle="None", ms=8,
            label=rf"MAX $({MIGKAS_MAXES[1][0]:.0f}^\circ,{MIGKAS_MAXES[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_115,
            markeredgecolor="k", linestyle="None", ms=7,
            label=rf"MIN $({MIGKAS_MINS[0][0]:.0f}^\circ,{MIGKAS_MINS[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_265,
            markeredgecolor="k", linestyle="None", ms=7,
            label=rf"MIN $({MIGKAS_MINS[1][0]:.0f}^\circ,{MIGKAS_MINS[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="o", color=C_ANTI, markeredgecolor="k",
            linestyle="None", ms=6,
            label=rf"anti-CoE $+x$ $(0^\circ,+{tilt:.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="*", color=C_COE, markeredgecolor="k",
            linestyle="None", ms=9,
            label=rf"CoE $(180^\circ,-{tilt:.0f}^\circ)$",
        ),
    ]
    leg = ax.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.012),
        ncol=4, fontsize=7.5, frameon=True, fancybox=False,
        edgecolor="#222222", framealpha=1.0, facecolor="white",
        labelcolor="#111111", borderpad=0.6, handlelength=2.0,
    )
    leg.get_frame().set_linewidth(1.2)
    fig.suptitle(
        r"Migkas Fig. 8 $\theta=75^\circ$ $A/A_{\rm all}$"
        r" + CoE full loops, opacity gradient darkest at $\theta$ ray, fading to 12% opacity opposite"
        "\n"
        rf"(anti-CoE at $(0^\circ,+{tilt:.0f}^\circ)$; CoE axis tilted $+{tilt:.0f}^\circ$; "
        rf"each $18^\circ$ end-to-end = $\pm9^\circ$; $\alpha\propto1-0.88|\psi|/\pi${roll_lab}; $\ell$ increases LEFT)",
        fontsize=10.5, y=0.985,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.13)
    fig.savefig(png_path, dpi=dpi, facecolor="white")
    plt.close(fig)
    print(f"wrote {png_path}")


def main(tilt_b=None, roll_deg=None):
    tilt_b = float(AXIS_TILT_B_DEG if tilt_b is None else tilt_b)
    roll_deg = float(ROLL_ABOUT_COE_DEG if roll_deg is None else roll_deg)
    assert FULL_GRADIENT_LOOPS

    pose = axis_pose(tilt_b)
    stem = output_stem(tilt_b, roll_deg)
    png_path = OUT_DIR / f"{stem}.png"
    json_path = OUT_DIR / f"{stem}.json"

    source = choose_source()
    print(f"source: {source}")
    rgb = np.asarray(Image.open(source).convert("RGB"))
    seed = detect_oval_bbox(rgb)
    reg = dict(GOOD_REG)
    h, w, _ = rgb.shape
    dark = rgb.mean(axis=2) < 55
    dist = chamfer_distance_to_dark(dark)
    ell_eq, b_eq = celestial_equator_lb(2500)
    lon = np.deg2rad(lon_disp_from_ell(ell_eq))
    lat = np.deg2rad(b_eq)
    mx, my = mollweide_xy(lon, lat)
    Xn = mx / (2.0 * np.sqrt(2.0))
    Yn = my / np.sqrt(2.0)
    px = reg["cx"] + Xn * reg["Rx"]
    py = reg["cy"] - Yn * reg["Ry"]
    m = (px >= 1) & (px < w - 1) & (py >= 1) & (py < h - 1)
    xi = np.clip(np.round(px[m]).astype(int), 0, w - 1)
    yi = np.clip(np.round(py[m]).astype(int), 0, h - 1)
    reg["ceq_median_px"] = float(np.median(dist[yi, xi]))
    reg["registration_source"] = "hard-set GOOD_REG (pre medial-axis retune)"
    checks = validate_landmarks(rgb, reg)
    print("validation OK")

    draw_overlay(rgb, reg, seed, source, pose, roll_deg, png_path)

    meta = {
        "png": str(png_path),
        "script": str(Path(__file__).resolve()),
        "geometry_module": str(OUT_DIR / "coe_half_loop_geometry.py"),
        "source_image": str(source),
        "gradient_module": str(OUT_DIR / "coe_gradient_loop_geometry.py"),
        "FULL_GRADIENT_LOOPS": True,
        "gradient_rule": GRADIENT_RULE,
        "alpha_levels": N_ALPHA_LEVELS,
        "base_alphas": {"fill": FILL_ALPHA, "white_underlay": 0.95, "mid_line": 1.0},
        "axis_tilt_b_deg": tilt_b,
        "roll_about_coe_deg": roll_deg,
        "anti_CoE_lb_deg": [float(pose["anti_lb"][0]), float(pose["anti_lb"][1])],
        "CoE_lb_deg": [float(pose["coe_lb"][0]), float(pose["coe_lb"][1])],
        "D_hat_galactic_xyz": pose["D"].tolist(),
        "migkas_max_lb_deg": [list(p) for p in MIGKAS_MAXES],
        "migkas_min_lb_deg": [list(p) for p in MIGKAS_MINS],
        "theta_max_from_plus_x_deg": list(THETA_MAX),
        "theta_min_from_plus_x_deg": list(THETA_MIN),
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "max_bands_deg": [list(b) for b in MAX_BANDS],
        "min_bands_deg": [list(b) for b in MIN_BANDS],
        "colors": {
            "MAX_30": C_MAX_30,
            "MAX_185": C_MAX_185,
            "MIN_115": C_MIN_115,
            "MIN_265": C_MIN_265,
        },
        "convention": (
            "Galactic Mollweide; lon_disp=-wrap(ell); ell increases LEFT "
            "(Migkas Fig.8)"
        ),
        "registration": {
            "cx": reg["cx"], "cy": reg["cy"], "Rx": reg["Rx"], "Ry": reg["Ry"],
            "ceq_median_px": reg["ceq_median_px"],
            "registration_source": reg.get("registration_source"),
        },
        "landmark_checks": checks,
        "matches_three_panel": str(OUT_DIR / "make_coe_rings_gradient_loops.py"),
        "drawing": (
            "FULL loops with opacity gradient. Each = PolyCollection strips between "
            "±9° edges + solid mid line (white underlay + colour LineCollections), "
            "alpha * (1+cos psi)/2 grouped into quantized levels; darkest at θ-ray "
            "point, alpha 0 at the opposite point. Antimeridian crossings dropped. "
            "Dec=0 NOT overlaid; published map untouched."
        ),
    }
    json_path.write_text(json.dumps(meta, indent=2) + "\n")
    print(f"wrote {json_path}")
    return png_path, json_path


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tilt", type=float, default=AXIS_TILT_B_DEG)
    p.add_argument("--roll", type=float, default=ROLL_ABOUT_COE_DEG)
    args = p.parse_args()
    main(tilt_b=args.tilt, roll_deg=args.roll)
