#!/usr/bin/env python3
"""
Superimpose CoE MAX/MIN highlight bands + Migkas markers onto the REAL Migkas
Fig. 8 theta=75 A/A_all Mollweide map.

Matches make_coe_rings_migkas_pin_axis.py (three-panel Panel B) geometry:
  Exactly FOUR features — each ONE transparent highlight band, 18° end-to-end
  (= ±9° = 5% of 360°) about the center θ + ONE solid mid-stroke:
    MAX 30°  → [21°, 39°]   color #00e5ff (bright cyan)
    MAX 185° → [176°, 194°] color #ff2d95 (magenta/hot pink)
    MIN 115° → [106°, 124°] color #b8f200 (lime)
    MIN 265° → [256°, 274°] color #4fc3f7 (sky blue)
  Drawing: clean PolyCollection fill + one mid-stroke (optional thin white
  underlay on mid only). NO black rim, NO thick white halo stacks, NO dual
  boundary-edge strokes. No 283° companion. No Dec=0 overlay.
  MAX triangles: (17°, −9°) and (185°, 0°)
  MIN squares:   (115°, 0°) and (274°, −22°)
  CoE axis tilted +30°: anti-CoE (0°,+30°), CoE (180°,−30°) (purple).

Registration (unchanged from working Fig. 8):
  GOOD_REG = (cx, cy, Rx, Ry) = (863.0, 763.25, 693.9927, 392.542425)
  Source: FULL_migkas_theta75.png
  CoE axis tilted +30°; clip/inset so overlays stay inside colored oval
  lon_disp = -wrap(ℓ)

Outputs:
  migkas_fig8_theta75_with_coe_rings.png
  migkas_fig8_theta75_with_coe_rings.json
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
from PIL import Image

OUT_DIR = Path(__file__).resolve().parent
SM_DIR = Path("/workspace/cosmology/papers/spreadsheet_model")
PNG_PATH = SM_DIR / "fig8_overlay_matched.png"
JSON_PATH = SM_DIR / "fig8_overlay_matched.json"

# Versions = panel B of each three-panel figure (grad: opacity gradient, roll: deg CCW about CoE axis)
VERSIONS = [
    ("matched", False, 0.0),
    ("matched_grad", True, 0.0),
    ("matched_grad_rot45", True, 45.0),
    ("matched_grad_rot81", True, 80.8),
    ("matched_rot81", False, 80.8),  # roll + all loops, no opacity gradient
]
GRAD = False
ROLL = 0.0
SHOW_LEGEND = False  # letter Fig.6: no legend box/strip

SOURCE_CANDIDATES = [
    Path("/workspace/cosmology/map_compare/FULL_migkas_theta75.png"),
    Path("/workspace/cosmology/map_compare/migkas_fig8_theta75_crop.png"),
    Path("/workspace/cosmology/migkas/pub_maps/Fig8_full_hires.png"),
]

# Strip A(ℓ) scan stops from +x (after-2ℓ A_after extrema) — four band centers.
# 5% of 360° = 18° end-to-end → half-width ±9° about each center (match three-panel).
BAND_HALF = 9.0
BAND_FULL = 2.0 * BAND_HALF  # 18°
THETA_MAX = [29.5, 187.5]
THETA_MIN = [122.6, 264.4]
MAX_BANDS = [(t - BAND_HALF, t + BAND_HALF) for t in THETA_MAX]  # [21,39], [176,194]
MIN_BANDS = [(t - BAND_HALF, t + BAND_HALF) for t in THETA_MIN]  # [106,124], [256,274]


def mu_pair_for_band(lo, hi):
    """Edge μ=cos(θ) for a scan-angle band (μ-space, not θ-space)."""
    return float(np.cos(np.deg2rad(lo))), float(np.cos(np.deg2rad(hi)))

# Visually-good registration (pre medial-axis retune). Hard-set as working default.
GOOD_REG = {
    "cx": 863.0,
    "cy": 763.25,
    "Rx": 693.9927,
    "Ry": 392.542425,
}

# Sky markers (ℓ, b): map extrema + strip stand-ins (same as 3-panel Panel B).
MIGKAS_MAXES = [(17.0, -9.0), (185.0, 0.0)]
MIGKAS_MINS = [(115.0, 0.0), (274.0, -22.0)]

# Galactic pose (Chris 2026-09-22): anti-CoE (+D/+x) at (0°, +30°); CoE at (180°, −30°)
ANTI_COE_LB = (0.0, 15.0)
COE_LB = (180.0, -15.0)
AXIS_TILT_B_DEG = 15.0


def lb_to_n(l_deg, b_deg):
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    return np.array(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], dtype=float
    )


D_HAT = lb_to_n(*ANTI_COE_LB)  # +D = anti-CoE = +x


def great_circle_through_poles(D, n=721):
    """Great circle through ±D (CoE axis meridian on the sky)."""
    a = D / np.linalg.norm(D)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([1.0, 0.0, 0.0])
    e2 = np.cross(a, tmp)
    e2 /= np.linalg.norm(e2)
    t = np.linspace(0.0, 2.0 * np.pi, n)
    pts = np.outer(np.cos(t), a) + np.outer(np.sin(t), e2)
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    return ls, bs

# High-contrast loop colors that pop on magma (avoid reds/oranges/corals)
C_MAX_30 = "#0097b2"    # cyan — MAX 29.5° (three-panel matched)
C_MAX_185 = "#e0107a"   # magenta — MAX 187.5°
C_MIN_115 = "#6fa800"   # green — MIN 122.6°
C_MIN_265 = "#2166d4"   # blue — MIN 264.4°
COLORS_MAX = [C_MAX_30, C_MAX_185]
COLORS_MIN = [C_MIN_115, C_MIN_265]
C_COE = "#6b3fa0"
C_ANTI = "#8e1b1b"  # dark red marker; distinct from loop reds (match three-panel)
C_AXIS = "#5b2c6f"
C_EQ = "#1e8449"


def wrap_lon_deg(lon):
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def lon_disp_from_ell(ell_deg):
    return -wrap_lon_deg(ell_deg)


def n_to_lb(n):
    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    b = float(np.rad2deg(np.arcsin(np.clip(n[2], -1.0, 1.0))))
    l = float(np.rad2deg(np.arctan2(n[1], n[0])) % 360.0)
    return l, b


def celestial_equator_lb(n=2500):
    ra_ngp = np.deg2rad(192.85948)
    dec_ngp = np.deg2rad(27.12825)
    l_ncp = np.deg2rad(122.93192)
    ra = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    dec = 0.0
    sin_b = np.clip(
        np.sin(dec) * np.sin(dec_ngp)
        + np.cos(dec) * np.cos(dec_ngp) * np.cos(ra - ra_ngp),
        -1.0,
        1.0,
    )
    b = np.arcsin(sin_b)
    y = np.cos(dec) * np.sin(ra - ra_ngp)
    x = (
        np.sin(dec) * np.cos(dec_ngp)
        - np.cos(dec) * np.sin(dec_ngp) * np.cos(ra - ra_ngp)
    )
    l = (l_ncp - np.arctan2(y, x)) % (2.0 * np.pi)
    return np.rad2deg(l), np.rad2deg(b)


def small_circle_about_axis(axis, mu, n=721):
    a = axis / np.linalg.norm(axis)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    tmp = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(a, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    t = np.linspace(0.0, 2.0 * np.pi, n)
    return (mu * a)[None, :] + rho * (
        np.outer(np.cos(t), e1) + np.outer(np.sin(t), e2)
    )


def mollweide_xy(lon_rad, lat_rad, n_iter=64):
    """Standard Mollweide; stable bracketed solve of 2θ+sin(2θ)=π sin(lat)."""
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
        "x0": int(x0),
        "x1": int(x1),
        "y0": int(y0),
        "y1": int(y1),
        "yeq": int(yeq),
        "cx": 0.5 * (x0 + x1),
        "cy": 0.5 * (y0 + y1),
        "Rx": 0.5 * (x1 - x0),
        "Ry": 0.5 * (y1 - y0 + 1),
    }


def refine_registration(rgb, seed):
    """Gray/dark Dec=0 lock (pre medial-axis). Prefer GOOD_REG hard-set in main."""
    h, w, _ = rgb.shape
    dark = rgb.mean(axis=2) < 55
    dist = chamfer_distance_to_dark(dark)

    ell_eq, b_eq = celestial_equator_lb(2500)
    lon = np.deg2rad(lon_disp_from_ell(ell_eq))
    lat = np.deg2rad(b_eq)
    mx, my = mollweide_xy(lon, lat)
    Xn = mx / (2.0 * np.sqrt(2.0))
    Yn = my / np.sqrt(2.0)

    def median_dist(cx, cy, Rx, Ry):
        px = cx + Xn * Rx
        py = cy - Yn * Ry
        m = (px >= 1) & (px < w - 1) & (py >= 1) & (py < h - 1)
        xi = np.clip(np.round(px[m]).astype(int), 0, w - 1)
        yi = np.clip(np.round(py[m]).astype(int), 0, h - 1)
        return float(np.median(dist[yi, xi]))

    base = {k: float(seed[k]) for k in ("cx", "cy", "Rx", "Ry")}
    best_sc = median_dist(**base)
    best = dict(base)
    for dcx in np.linspace(-10, 10, 21):
        for dcy in np.linspace(-10, 10, 21):
            for sRx in (0.98, 0.99, 1.0, 1.01, 1.02):
                for sRy in (0.98, 0.99, 1.0, 1.01, 1.02):
                    p = dict(
                        cx=base["cx"] + dcx,
                        cy=base["cy"] + dcy,
                        Rx=base["Rx"] * sRx,
                        Ry=base["Ry"] * sRy,
                    )
                    sc = median_dist(**p)
                    if sc < best_sc:
                        best_sc = sc
                        best = p
    b = best
    for dcx in np.linspace(-2, 2, 17):
        for dcy in np.linspace(-2, 2, 17):
            for sRx in np.linspace(0.995, 1.005, 11):
                for sRy in np.linspace(0.995, 1.005, 11):
                    p = dict(
                        cx=b["cx"] + dcx,
                        cy=b["cy"] + dcy,
                        Rx=b["Rx"] * sRx,
                        Ry=b["Ry"] * sRy,
                    )
                    sc = median_dist(**p)
                    if sc < best_sc:
                        best_sc = sc
                        best = p
    best["ceq_median_px"] = float(best_sc)
    return best


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

    # Map-extrema pair: (17,-9) should sit in a warmer patch than (274,-22).
    cmax = sample_rgb(rgb, *MIGKAS_MAXES[0], cx, cy, Rx, Ry)
    cmin = sample_rgb(rgb, *MIGKAS_MINS[1], cx, cy, Rx, Ry)
    assert (cmax[0] + cmax[1]) > (cmin[0] + cmin[1]) + 50
    return {
        "landmarks_px": pts,
        "migkas_max0_rgb_mean": [float(x) for x in cmax],
        "migkas_min1_rgb_mean": [float(x) for x in cmin],
    }


def inside_ellipse(px, py, ecx, ecy, eRx, eRy, r_lim=1.0):
    """Boolean mask: points inside axis-aligned ellipse (image y down)."""
    px = np.asarray(px, dtype=float)
    py = np.asarray(py, dtype=float)
    Xn = (px - ecx) / eRx
    Yn = (py - ecy) / eRy
    return Xn * Xn + Yn * Yn <= float(r_lim) * float(r_lim)


def plot_lb_curve(
    ax, ell, b, cx, cy, Rx, Ry, inset=1.0, jump_px=45.0,
    clip_r2=None, clip_path=None, color_oval=None, **kw
):
    """Draw curve in oval; inset<1 keeps strokes inside the published Mollweide edge."""
    px, py = lb_to_px(ell, b, cx, cy, Rx, Ry)
    px = np.asarray(px, dtype=float).reshape(-1)
    py = np.asarray(py, dtype=float).reshape(-1)
    # Pull points slightly toward oval center so thick strokes do not spill past the rim.
    px = cx + (px - cx) * inset
    py = cy + (py - cy) * inset
    Xn = (px - cx) / Rx
    Yn = (cy - py) / Ry
    r2_lim = (inset * inset) if clip_r2 is None else float(clip_r2)
    inside = Xn * Xn + Yn * Yn <= r2_lim
    # Also require published colored oval (seed) so registration oversize cannot bleed.
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


def plot_closed_ring(
    ax, mu, cx, cy, Rx, Ry, color, lw, inset=1.0, jump_px=80.0,
    clip_path=None, color_oval=None, ls="-", wash_alpha=None,
):
    """
    Exact projected small-circle about D_HAT=+x: n·x̂ = μ.

    Parameterize on the sphere, convert to Galactic (ℓ,b), apply
    lon_disp = -wrap(ℓ), Mollweide-project, and stroke continuous polylines.
    Split wherever display-longitude jumps across the antimeridian (or pixel
    distance jumps). Never invent interior chords, ellipse fits, or decorative
    oval-rim bridges.

    Notes on θ labels used elsewhere:
      μ=cos(30°)  → r=30° about +x (closed oval near map center)
      μ=cos(185°)=cos(175°) → r=175° about +x ≡ 5° about CoE
        → small arcs on BOTH left and right map edges (not a rim loop)
      μ=cos(115°) → r=115° about +x (large seam-crossing curve)
      μ=cos(265°)=cos(95°) → r=95° about +x (near-meridian, seam-crossing)
    """
    pts = small_circle_about_axis(D_HAT, mu, n=4001)
    # Keep sampling order around the axis (azimuth) — stable for splitting.
    phi = np.arctan2(pts[:, 2], pts[:, 1])
    order = np.argsort(phi)
    pts = pts[order]

    ell = np.empty(len(pts))
    b = np.empty(len(pts))
    for i, p in enumerate(pts):
        ell[i], b[i] = n_to_lb(p)

    lon = np.deg2rad(lon_disp_from_ell(ell))  # λ ∈ (-π, π]
    lat = np.deg2rad(b)

    # Split on antimeridian jumps in display longitude BEFORE pixel projection.
    dlon = np.diff(lon)
    cuts = np.where(np.abs(dlon) > np.pi)[0] + 1
    # Also close the periodic sample: check wrap between last and first.
    # Build open segments from azimuth-ordered points.
    idx = np.concatenate([[0], cuts, [len(lon)]])
    lon_segs = []
    for i0, i1 in zip(idx[:-1], idx[1:]):
        if i1 - i0 >= 2:
            lon_segs.append((lon[i0:i1].copy(), lat[i0:i1].copy()))

    # Circle sampling is cyclic. Merge the artificial last→first break when
    # that periodic join is not itself an antimeridian crossing. Otherwise a
    # seam-crossing circle falsely appears as three pieces with a tiny hole.
    if len(lon_segs) > 1 and abs(lon[0] - lon[-1]) <= np.pi:
        first_lo, first_la = lon_segs[0]
        last_lo, last_la = lon_segs[-1]
        lon_segs = [(
            np.concatenate([last_lo, first_lo]),
            np.concatenate([last_la, first_la]),
        )] + lon_segs[1:-1]

    # Only now pin TRUE seam endpoints to ±π. Never snap the arbitrary
    # parameterization start/end to the boundary.
    snapped = []
    for lo, la in lon_segs:
        if abs(abs(lo[0]) - np.pi) < np.deg2rad(10):
            lo[0] = np.sign(lo[0] if lo[0] != 0 else 1.0) * np.pi
        if abs(abs(lo[-1]) - np.pi) < np.deg2rad(10):
            lo[-1] = np.sign(lo[-1] if lo[-1] != 0 else 1.0) * np.pi
        snapped.append((lo, la))
    lon_segs = snapped

    def stroke_xy(xs, ys):
        if len(xs) < 2:
            return
        if wash_alpha is not None:
            # Soft μ-wash stroke: no white underlay; translucent color only.
            layers = ((color, lw, 2, float(wash_alpha)),)
        else:
            # Clean mid-stroke: optional thin white underlay only (lw+1.5), then color.
            # No black outer rim, no thick white halo stack.
            layers = (
                ("white", lw + 1.5, 3, 0.95),
                (color, lw, 4, 1.0),
            )
        for col, width, z, a in layers:
            (line,) = ax.plot(
                xs, ys, color=col, lw=width, linestyle=ls,
                solid_capstyle="round", solid_joinstyle="round",
                zorder=z, alpha=a,
                clip_path=clip_path, clip_on=True,
            )
            if clip_path is not None:
                line.set_clip_path(clip_path)
                line.set_clip_on(True)

    r2_lim = inset * inset + 1e-10
    drawn = 0
    for lo, la in lon_segs:
        # Mollweide in normalized oval coords, then pixels.
        x_m, y_m = mollweide_xy(lo, la)
        Xn = x_m / (2.0 * np.sqrt(2.0))
        Yn = y_m / np.sqrt(2.0)
        px = cx + Xn * Rx
        py = cy - Yn * Ry
        # Mild inset so thick strokes stay inside published oval (no rim-snap).
        px = cx + (px - cx) * inset
        py = cy + (py - cy) * inset
        Xn2 = (px - cx) / Rx
        Yn2 = (cy - py) / Ry
        inside = (Xn2 * Xn2 + Yn2 * Yn2) <= r2_lim
        if color_oval is not None:
            ocx, ocy, oRx, oRy, o_r = color_oval
            inside = inside & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)

        # Pixel-space split on NaN / large jumps / outside oval.
        cur_x, cur_y = [], []
        for x, y, ok in zip(px, py, inside):
            if (not ok) or (not np.isfinite(x)) or (not np.isfinite(y)):
                if len(cur_x) >= 2:
                    stroke_xy(cur_x, cur_y)
                    drawn += 1
                cur_x, cur_y = [], []
                continue
            if cur_x and np.hypot(x - cur_x[-1], y - cur_y[-1]) > jump_px:
                if len(cur_x) >= 2:
                    stroke_xy(cur_x, cur_y)
                    drawn += 1
                cur_x, cur_y = [float(x)], [float(y)]
            else:
                cur_x.append(float(x))
                cur_y.append(float(y))
        if len(cur_x) >= 2:
            # Close only if this single segment is already nearly a loop
            # (non-wrapping ring such as MAX 30°).
            if (
                drawn == 0
                and len(lon_segs) == 1
                and np.hypot(cur_x[0] - cur_x[-1], cur_y[0] - cur_y[-1]) < jump_px
            ):
                cur_x.append(cur_x[0])
                cur_y.append(cur_y[0])
            stroke_xy(cur_x, cur_y)
            drawn += 1



def _ring_px_arrays(mu, cx, cy, Rx, Ry, inset, n_pts=721):
    """Small-circle n·D=μ → pixel (px, py) arrays, azimuth-ordered, with inset."""
    pts = small_circle_about_axis(D_HAT, float(mu), n=n_pts)
    phi = np.arctan2(pts[:, 2], pts[:, 1])
    order = np.argsort(phi)
    pts = pts[order]
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
    return lon, px, py


def fill_highlight_band(
    ax, theta_lo, theta_hi, theta_mid, cx, cy, Rx, Ry, color,
    inset=1.0, clip_path=None, color_oval=None, n_pts=721, jump_px=80.0,
):
    """
    One highlight band matching three-panel Panel B:
      continuous transparent polygon strips between ±BAND_HALF edges
      (broken only at Mollweide antimeridian / large pixel jumps — NOT a dense
      multi-stroke wash) + ONE solid mid-stroke at θ_mid.
    """
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    lon_a, px_a, py_a = _ring_px_arrays(mu_a, cx, cy, Rx, Ry, inset, n_pts=n_pts)
    lon_b, px_b, py_b = _ring_px_arrays(mu_b, cx, cy, Rx, Ry, inset, n_pts=n_pts)

    r2_lim = inset * inset + 1e-10

    def _inside(px, py):
        Xn2 = (px - cx) / Rx
        Yn2 = (cy - py) / Ry
        ok = (Xn2 * Xn2 + Yn2 * Yn2) <= r2_lim
        if color_oval is not None:
            ocx, ocy, oRx, oRy, o_r = color_oval
            ok = ok & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)
        return ok

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
        elif np.hypot(px_a[i + 1] - px_a[i], py_a[i + 1] - py_a[i]) > jump_px:
            bad[i] = True
        elif np.hypot(px_b[i + 1] - px_b[i], py_b[i + 1] - py_b[i]) > jump_px:
            bad[i] = True
        elif (not _inside(px_a[i], py_a[i])) or (not _inside(px_b[i], py_b[i])):
            bad[i] = True
        elif (not _inside(px_a[i + 1], py_a[i + 1])) or (not _inside(px_b[i + 1], py_b[i + 1])):
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
            xs = np.concatenate([px_a[i0:i1 + 1], px_b[i0:i1 + 1][::-1]])
            ys = np.concatenate([py_a[i0:i1 + 1], py_b[i0:i1 + 1][::-1]])
            polys.append(np.column_stack([xs, ys]))
        i = j + 1

    if polys:
        # Soft transparent highlight (clean pre-halo look).
        coll = PolyCollection(
            polys, facecolors=color, edgecolors="none",
            alpha=0.30, zorder=2, linewidths=0,
        )
        ax.add_collection(coll)
        if clip_path is not None:
            coll.set_clip_path(clip_path)
            coll.set_clip_on(True)

    # ONE solid mid-stroke only (thin white underlay lw+1.5). No boundary-edge
    # strokes — dual edges with halos multiply visual lines into "extra loops".
    plot_closed_ring(
        ax, mu_mid, cx, cy, Rx, Ry, color, lw=2.8,
        inset=inset, clip_path=clip_path, color_oval=color_oval,
    )



# --- Opacity gradient (same rule as make_coe_rings_migkas_pin_axis_matched_grad*.py):
# alpha *= (1+cos psi)/2, psi = ring azimuth from the θ-ray point (ring point most
# along sign(sin θ)·Ŷ, Ŷ = Galactic +z projected ⊥ D), rolled by ROLL (Y' = cY − sZ).
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba


def _ring_basis_roll(D, roll_deg):
    a = D / np.linalg.norm(D)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    y = tmp - np.dot(tmp, a) * a
    y /= np.linalg.norm(y)
    z = np.cross(a, y)
    z /= np.linalg.norm(z)
    r = np.deg2rad(float(roll_deg))
    c, s_ = np.cos(r), np.sin(r)
    return c * y - s_ * z, s_ * y + c * z


def _grad_ring_px(mu, theta_mid, cx, cy, Rx, Ry, inset, n):
    a = D_HAT / np.linalg.norm(D_HAT)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    rho = np.sqrt(1.0 - mu * mu)
    y, z = _ring_basis_roll(a, ROLL)
    phi0 = 0.5 * np.pi if np.sin(np.deg2rad(theta_mid)) >= 0 else 1.5 * np.pi
    phi = phi0 + np.linspace(-np.pi, np.pi, n)
    pts = (mu * a)[None, :] + rho * (np.outer(np.sin(phi), y) + np.outer(np.cos(phi), z))
    ell = np.array([n_to_lb(p)[0] for p in pts]); b = np.array([n_to_lb(p)[1] for p in pts])
    lon = np.deg2rad(lon_disp_from_ell(ell)); lat = np.deg2rad(b)
    x_m, y_m = mollweide_xy(lon, lat)
    px = cx + x_m / (2.0 * np.sqrt(2.0)) * Rx
    py = cy - y_m / np.sqrt(2.0) * Ry
    px = cx + (px - cx) * inset; py = cy + (py - cy) * inset
    return lon, px, py, 0.5 * (1.0 + np.cos(phi - phi0))


def _rgba(color, alphas):
    r, g, b_, _ = to_rgba(color)
    al = np.clip(np.asarray(alphas, dtype=float), 0.0, 1.0)
    out = np.empty((al.size, 4)); out[:, 0], out[:, 1], out[:, 2], out[:, 3] = r, g, b_, al
    return out


def _runs(seg_sc, valid, n_levels=64):
    lev = np.round(np.clip(seg_sc, 0.0, 1.0) * n_levels) / n_levels
    ok = (lev > 0) & valid
    runs, i, n = [], 0, lev.size
    while i < n:
        if not ok[i]:
            i += 1; continue
        j = i + 1
        while j < n and ok[j] and lev[j] == lev[i]:
            j += 1
        runs.append((i, j, float(lev[i]))); i = j
    return runs


def fill_highlight_band_grad(
    ax, theta_lo, theta_hi, theta_mid, cx, cy, Rx, Ry, color,
    inset=1.0, clip_path=None, color_oval=None, n_pts=721, jump_px=80.0,
):
    """Same band (alpha 0.30) + mid-stroke (white lw+1.5 @0.95 under colour lw 2.8)
    as fill_highlight_band, but with the (1+cos psi)/2 opacity gradient."""
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    r2_lim = inset * inset + 1e-10

    def _inside(px, py):
        ok = ((px - cx) / Rx) ** 2 + ((cy - py) / Ry) ** 2 <= r2_lim
        if color_oval is not None:
            ocx, ocy, oRx, oRy, o_r = color_oval
            ok = ok & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)
        return ok

    def _clip(coll):
        if clip_path is not None:
            coll.set_clip_path(clip_path); coll.set_clip_on(True)
        ax.add_collection(coll)

    gap = 0.9 * np.pi
    lon_a, pxa, pya, sc = _grad_ring_px(mu_a, theta_mid, cx, cy, Rx, Ry, inset, n_pts)
    lon_b, pxb, pyb, _ = _grad_ring_px(mu_b, theta_mid, cx, cy, Rx, Ry, inset, n_pts)
    qa = np.stack([lon_a[:-1], lon_a[1:], lon_b[1:], lon_b[:-1]])
    ins_a, ins_b = _inside(pxa, pya), _inside(pxb, pyb)
    valid = ((qa.max(0) - qa.min(0)) <= gap)
    valid &= np.hypot(np.diff(pxa), np.diff(pya)) <= jump_px
    valid &= np.hypot(np.diff(pxb), np.diff(pyb)) <= jump_px
    valid &= ins_a[:-1] & ins_a[1:] & ins_b[:-1] & ins_b[1:]
    runs = _runs(0.5 * (sc[:-1] + sc[1:]), valid)
    if runs:
        polys = [np.column_stack([np.concatenate([pxa[i0:i1 + 1], pxb[i0:i1 + 1][::-1]]),
                                  np.concatenate([pya[i0:i1 + 1], pyb[i0:i1 + 1][::-1]])])
                 for i0, i1, _ in runs]
        _clip(PolyCollection(polys, facecolors=_rgba(color, [0.30 * r[2] for r in runs]),
                             edgecolors="none", linewidths=0, zorder=2))
    lon_m, pxm, pym, scm = _grad_ring_px(mu_mid, theta_mid, cx, cy, Rx, Ry, inset, 4001)
    ins = _inside(pxm, pym)
    valid = (np.abs(np.diff(lon_m)) <= np.pi) & (np.hypot(np.diff(pxm), np.diff(pym)) <= jump_px)
    valid &= ins[:-1] & ins[1:]
    runs = _runs(0.5 * (scm[:-1] + scm[1:]), valid)
    if runs:
        segs = [np.column_stack([pxm[i0:i1 + 1], pym[i0:i1 + 1]]) for i0, i1, _ in runs]
        _clip(LineCollection(segs, colors=_rgba("white", [0.95 * r[2] for r in runs]),
                             linewidths=2.8 + 1.5, zorder=3, capstyle="butt", joinstyle="round"))
        _clip(LineCollection(segs, colors=_rgba(color, [1.0 * r[2] for r in runs]),
                             linewidths=2.8, zorder=4, capstyle="butt", joinstyle="round"))


def draw_overlay(rgb, reg, seed, source):

    h, w, _ = rgb.shape
    cx, cy, Rx, Ry = reg["cx"], reg["cy"], reg["Rx"], reg["Ry"]
    # Published colored Mollweide (seed detect) — used to keep overlays inside
    # the Migkas map even when registration ellipse is slightly larger.
    scx = float(seed["cx"])
    scy = float(seed["cy"])
    sRx = float(seed["Rx"])
    sRy = float(seed["Ry"])

    dpi = 120
    _extra_h = 1.35 if SHOW_LEGEND else 0.35  # tighten when legend omitted
    fig, ax = plt.subplots(figsize=(w / dpi, (h / dpi) + _extra_h), dpi=dpi)
    ax.imshow(rgb, origin="upper", extent=(0, w, h, 0), interpolation="nearest")
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    # clip_path = registered oval (cx,cy,Rx,Ry). Registration itself is unchanged.
    # Because the refined ellipse can sit a bit outside the published colored
    # map, we ALSO inset strokes and cull/inset markers against the seed oval.
    clip_oval = Ellipse(
        (cx, cy), width=2.0 * Rx, height=2.0 * Ry,
        facecolor="none", edgecolor="none", transform=ax.transData,
    )
    ax.add_patch(clip_oval)

    # Stroke inset so thick white+color underlays do not AA-bleed past the rim.
    STROKE_INSET = 0.97
    # Cull points outside this fraction of the published colored seed oval.
    COLOR_R = 0.978
    color_oval = (scx, scy, sRx, sRy, COLOR_R)

    def stroke_curve(ell, b, color, lw, ls="-"):
        # White underlay then color — stands out on magma without looking oversized.
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=STROKE_INSET, clip_path=clip_oval,
            color_oval=color_oval,
            color="white", lw=lw + 1.6, ls=ls, alpha=0.95,
            solid_capstyle="round", zorder=3,
        )
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=STROKE_INSET, clip_path=clip_oval,
            color_oval=color_oval,
            color=color, lw=lw, ls=ls, alpha=1.0,
            solid_capstyle="round", zorder=4,
        )

    def mark_point(ell, b, marker, color, ms):
        """Place marker; pull limb symbols inward so halo stays in colored oval."""
        px, py = as_xy(*lb_to_px(ell, b, cx, cy, Rx, Ry))
        # Approx marker+halo radius in pixels (ms is points).
        halo_ms = ms + 1.8
        ms_px = 0.5 * halo_ms * (dpi / 72.0)
        # Normalized margin inside the published colored oval.
        margin = ms_px / min(sRx, sRy) + 0.012
        r_lim = max(0.55, COLOR_R - margin)

        Xn = (px - scx) / sRx
        Yn = (py - scy) / sRy
        r_col = float(np.hypot(Xn, Yn))

        if r_col > 1e-9:
            # Scale toward colored-oval center so the whole symbol fits.
            if r_col > r_lim:
                scale = r_lim / r_col
                px = scx + (px - scx) * scale
                py = scy + (py - scy) * scale
                Xn = (px - scx) / sRx
                Yn = (py - scy) / sRy
                r_col = float(np.hypot(Xn, Yn))
            # Mild extra inset for everyone (matches stroke inset spirit).
            mild = min(STROKE_INSET, 0.995)
            px = scx + (px - scx) * mild
            py = scy + (py - scy) * mild
            Xn = (px - scx) / sRx
            Yn = (py - scy) / sRy
            r_col = float(np.hypot(Xn, Yn))

        # Skip if still too close to the rim for this marker size.
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

    # CoE axis at geometric inset=1.0. Do NOT overlay Dec=0 — Migkas already
    # printed the celestial equator on Fig. 8; a second green dashed stroke fights it.
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

    # CoE↔anti-CoE axis = great circle through ±D (tilted +30°; not Galactic equator)
    _ax_ell, _ax_b = great_circle_through_poles(D_HAT, n=721)
    stroke_curve_geom(_ax_ell, _ax_b, C_AXIS, lw=1.15)

    # Exactly four highlight bands (match three-panel Panel B): each ±9°
    # (= 18° end-to-end = 5% of 360°) + one solid mid-stroke. No companion.
    for (lo, hi), mid, color in zip(MAX_BANDS, THETA_MAX, COLORS_MAX):
        (fill_highlight_band_grad if GRAD else fill_highlight_band)(
            ax, lo, hi, mid, cx, cy, Rx, Ry, color,
            inset=STROKE_INSET, clip_path=clip_oval, color_oval=color_oval,
        )
    for (lo, hi), mid, color in zip(MIN_BANDS, THETA_MIN, COLORS_MIN):
        (fill_highlight_band_grad if GRAD else fill_highlight_band)(
            ax, lo, hi, mid, cx, cy, Rx, Ry, color,
            inset=STROKE_INSET, clip_path=clip_oval, color_oval=color_oval,
        )

    # Markers (limb ones get stronger inset inside mark_point); per-loop colors
    mark_point(ANTI_COE_LB[0], ANTI_COE_LB[1], "o", C_ANTI, ms=5.5)
    mark_point(COE_LB[0], COE_LB[1], "*", C_COE, ms=7.5)
    for (ell, b), color in zip(MIGKAS_MAXES, COLORS_MAX):
        mark_point(ell, b, "^", color, ms=6.5)
    for (ell, b), color in zip(MIGKAS_MINS, COLORS_MIN):
        mark_point(ell, b, "s", color, ms=6.0)

    handles = [
        Patch(
            facecolor=C_MAX_30, alpha=0.55, edgecolor=C_MAX_30,
            label=(
                rf"MAX $\theta\!=\!{THETA_MAX[0]:g}^\circ$ band "
                rf"[{MAX_BANDS[0][0]:g},{MAX_BANDS[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$=5% of $360^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MAX_185, alpha=0.55, edgecolor=C_MAX_185,
            label=(
                rf"MAX $\theta\!=\!{THETA_MAX[1]:g}^\circ$ band "
                rf"[{MAX_BANDS[1][0]:g},{MAX_BANDS[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_115, alpha=0.55, edgecolor=C_MIN_115,
            label=(
                rf"MIN $\theta\!=\!{THETA_MIN[0]:g}^\circ$ band "
                rf"[{MIN_BANDS[0][0]:g},{MIN_BANDS[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_265, alpha=0.55, edgecolor=C_MIN_265,
            label=(
                rf"MIN $\theta\!=\!{THETA_MIN[1]:g}^\circ$ band "
                rf"[{MIN_BANDS[1][0]:g},{MIN_BANDS[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$)"
            ),
        ),
        Line2D([0], [0], color=C_AXIS, lw=2.0, label=r"CoE axis tilted $+15^\circ$"),
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
        Line2D([0], [0], marker="o", color=C_ANTI, markeredgecolor="k",
               linestyle="None", ms=6,
               label=r"anti-CoE $+x$ $(0^\circ,+15^\circ)$"),
        Line2D([0], [0], marker="*", color=C_COE, markeredgecolor="k",
               linestyle="None", ms=9,
               label=r"CoE $(180^\circ,-15^\circ)$"),
    ]
    if SHOW_LEGEND:
        leg = ax.legend(
            handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.012),
            ncol=4, fontsize=7.5, frameon=True, fancybox=False,
            edgecolor="#222222", framealpha=1.0, facecolor="white",
            labelcolor="#111111", borderpad=0.6, handlelength=2.0,
        )
        leg.get_frame().set_linewidth(1.2)
    fig.suptitle(
        r"Migkas Fig. 8 $\theta=75^\circ$ $A/A_{\rm all}$"
        r" + four CoE highlight bands, each $18^\circ$ end-to-end ($\pm9^\circ$ = 5% of $360^\circ$)"
        "\n"
        r"(two MAX + two MIN; cyan/green/magenta/blue; anti-CoE at $(0^\circ,+15^\circ)$; "
        r"CoE axis tilted $+15^\circ$; $\ell$ increases LEFT"
        + (r"; gradient: darkest at $\theta$ ray, transparent opposite" if GRAD else "")
        + (rf"; ${ROLL:+g}^\circ$ CCW about CoE axis" if abs(ROLL) > 1e-9 else "")
        + ")",
        fontsize=10.5, y=0.985,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=(0.13 if SHOW_LEGEND else 0.02))
    fig.savefig(PNG_PATH, dpi=dpi, facecolor="white")
    plt.close(fig)
    print(f"wrote {PNG_PATH}")


def main():
    source = choose_source()
    print(f"source: {source}")
    rgb = np.asarray(Image.open(source).convert("RGB"))
    seed = detect_oval_bbox(rgb)
    print("oval seed:", {k: seed[k] for k in ("cx", "cy", "Rx", "Ry", "x0", "x1", "y0", "y1")})
    # Restore previous visually-good oval (Chris preferred). No medial-axis Rx stretch.
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
    print("restored reg:", {k: reg[k] for k in ("cx", "cy", "Rx", "Ry", "ceq_median_px")})
    checks = validate_landmarks(rgb, reg)
    print(
        "validation OK — max0 RGB", checks["migkas_max0_rgb_mean"],
        "min1 RGB", checks["migkas_min1_rgb_mean"],
    )
    global PNG_PATH, JSON_PATH, GRAD, ROLL
    for name, g, r in VERSIONS:
        PNG_PATH = SM_DIR / f"fig8_overlay_{name}.png"
        JSON_PATH = SM_DIR / f"fig8_overlay_{name}.json"
        GRAD, ROLL = g, r
        draw_overlay(rgb, reg, seed, source)
        img = Image.open(PNG_PATH).convert("RGB")
        img.resize((1700, round(img.height * 1700 / img.width)), Image.LANCZOS).save(
            SM_DIR / f"fig8_overlay_{name}_web.jpg", quality=88)

    meta = {
        "png": str(PNG_PATH),
        "script": str(Path(__file__).resolve()),
        "source_image": str(source),
        "convention": (
            "Galactic Mollweide; lon_disp=-wrap(ell); ell increases LEFT "
            "(Migkas Fig.8)"
        ),
        "axis_tilt_b_deg": AXIS_TILT_B_DEG,
        "anti_CoE_lb_deg": [float(ANTI_COE_LB[0]), float(ANTI_COE_LB[1])],
        "CoE_lb_deg": [float(COE_LB[0]), float(COE_LB[1])],
        "D_hat_galactic_xyz": D_HAT.tolist(),
        "axis_tilt_reason": (
            "Chris 2026-09-22 — markers fixed, rings rotate up 30°"
        ),
        "migkas_max_lb_deg": [list(p) for p in MIGKAS_MAXES],
        "migkas_min_lb_deg": [list(p) for p in MIGKAS_MINS],
        "theta_max_from_plus_x_deg": list(THETA_MAX),
        "theta_min_from_plus_x_deg": list(THETA_MIN),
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "band_width_note": "5% of 360° = 18° end-to-end = ±9° about each MAX/MIN center",
        "max_bands_deg": [list(b) for b in MAX_BANDS],
        "min_bands_deg": [list(b) for b in MIN_BANDS],
        "band_intervals": {
            "MAX_30": [float(MAX_BANDS[0][0]), float(MAX_BANDS[0][1])],
            "MAX_185": [float(MAX_BANDS[1][0]), float(MAX_BANDS[1][1])],
            "MIN_115": [float(MIN_BANDS[0][0]), float(MIN_BANDS[0][1])],
            "MIN_265": [float(MIN_BANDS[1][0]), float(MIN_BANDS[1][1])],
        },
        "colors": {
            "MAX_30": C_MAX_30,
            "MAX_185": C_MAX_185,
            "MIN_115": C_MIN_115,
            "MIN_265": C_MIN_265,
            "palette": "magma-popping cyan/magenta/lime/sky (no reds/oranges)",
        },
        "mu_max_bands": [list(mu_pair_for_band(lo, hi)) for lo, hi in MAX_BANDS],
        "mu_min_bands": [list(mu_pair_for_band(lo, hi)) for lo, hi in MIN_BANDS],
        "scan_stop_source": (
            "best_fit_2ell_multiplier.npz A_after extrema (two max, two min); "
            f"each drawn as one transparent {BAND_FULL:.0f}° end-to-end highlight "
            f"(±{BAND_HALF:.0f}° = 5% of 360°) + one solid mid-stroke; "
            "matches three-panel Panel B; no companion"
        ),
        "matches_three_panel": str(OUT_DIR / "make_coe_rings_migkas_pin_axis.py"),
        "oval_seed_bbox": {
            "x0": seed["x0"], "x1": seed["x1"], "y0": seed["y0"], "y1": seed["y1"],
            "cx": seed["cx"], "cy": seed["cy"], "Rx": seed["Rx"], "Ry": seed["Ry"],
        },
        "registration": {
            "cx": reg["cx"], "cy": reg["cy"], "Rx": reg["Rx"], "Ry": reg["Ry"],
            "ceq_median_px": reg["ceq_median_px"],
            "registration_source": reg.get("registration_source"),
            "ceq_residual_before_registration": dict(GOOD_REG),
            "rejected_medial_axis_reg": {
                "cx": 865.6, "cy": 766.4, "Rx": 706.19, "Ry": 393.70,
                "note": "medial-axis / achromatic DT-ridge retune; blew up Rx; reverted",
            },
            "dofs": ["cx", "cy", "Rx", "Ry"],
            "drawing": {
                "dec0_drawn": False,
                "dec0_note": "not drawn; use Migkas published celestial equator",
                "coe_axis_inset": 1.0,
                "rings_inset": 0.97,
                "note": "CoE axis geometric inset=1; rings mild inset; no Dec=0 overlay",
            },
            "method": (
                "oval seed detect + hard-set to pre-retune visually-good (cx,cy,Rx,Ry); "
                "medial-axis Rx stretch rejected; Dec=0 NOT overlaid (published CEQ kept); "
                "CoE-axis inset=1.0; rings inset~0.97"
            ),
        },
        "landmark_checks": checks,
        "drawing": (
            "Exactly 4 features. Each = ONE smooth transparent PolyCollection "
            "band (θ±9°, 18° end-to-end) + ONE solid mid-stroke at center θ. "
            "Fill alpha≈0.30. Mid-stroke: thin white underlay (lw+1.5) then color "
            "(lw=2.8). NO black outer rim; NO thick white halo stacks; NO dual "
            "boundary-edge strokes. Colors: cyan/magenta/lime/sky for magma. "
            "No 283° companion. Dec=0 NOT overlaid."
        ),
        "legend_summary": (
            "four highlight bands, each 18° end-to-end (±9° = 5% of 360%); "
            "MAX cyan+magenta, MIN lime+sky; intervals [21,39]/[176,194]/[106,124]/[256,274]"
        ),
        "notes": [
            "Exactly 4 features matching three-panel Panel B geometry/width",
            "Each = ONE smooth transparent PolyCollection ribbon (±9°) + ONE solid mid-stroke only",
            "Clean draw: fill alpha=0.30; mid-stroke thin white underlay (lw+1.5) + color lw=2.8; no black rim; no boundary-edge halos",
            "No multi-stroke scribble wash; no 283° companion; no Dec=0 overlay",
            "MAX bands [21,39] & [176,194]; MIN bands [106,124] & [256,274]",
            "Colors: MAX_30=#00e5ff, MAX_185=#ff2d95, MIN_115=#b8f200, MIN_265=#4fc3f7",
            "band_half_width_deg=9, band_full_width_deg=18 (5% of 360°)",
            "Exact n·x̂=μ small-circles; seam-split polygon fill; mid via plot_closed_ring",
            "No ellipse-fit, no interior chords, no decorative oval-rim loops",
            "No blocking text box on the map; legend below",
            "Colorbar preserved from source (overlays stay on oval)",
            "Overlays clip_path = registered Ellipse(cx,cy,Rx,Ry); rings stroke inset 0.97",
            "Dec=0 NOT drawn (Migkas Fig.8 already has celestial equator); CoE axis = great circle through ±D (tilted +30°) at inset=1.0",
            "axis_tilt_b_deg=30; anti-CoE (0°,+30°), CoE (180°,−30°); Migkas markers fixed at true sky positions",
            "Registration hard-set to pre medial-axis params (863,763.25,693.9927,392.542425); DT-ridge fit rejected",
            "Limb markers inset toward colored-oval center so halo stays inside",
            "Points culled at COLOR_R=0.978 of seed oval; limb markers inset to fit halo",
        ],
    }
    pass  # per-version JSON not written (PNG only)
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
