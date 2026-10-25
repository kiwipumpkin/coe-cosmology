#!/usr/bin/env python3
"""
Authoritative spherical small-circle + Mollweide geometry.

Equation
--------
Center n0 = n(l0, b0), angular radius r ∈ [0, 180°]:
    n(l,b) · n0 = cos(r)

Stable parameterization (orthonormal frame about n0)
----------------------------------------------------
Pick any tmp not parallel to n0; then
    e1 = normalize(tmp × n0)   # right-handed
    e2 = n0 × e1
    n(φ) = cos(r) n0 + sin(r) (cos φ e1 + sin φ e2),  φ ∈ [0, 2π)

(Equivalently μ = cos r, ρ = sin r; for r > 90°, μ < 0.)

Mollweide (standard; Snyder / matplotlib)
-----------------------------------------
Solve  2θ + sin(2θ) = π sin(φ)  for auxiliary θ ∈ [-π/2, π/2]
(φ = latitude). Then
    x = (2√2 / π) λ cos θ
    y = √2 sin θ
with λ = longitude in (-π, π] relative to the map center.

Galactic display used here (Migkas): λ = lon_disp = -wrap(l)
so Galactic longitude increases LEFT.

Antimeridian / seam
-------------------
A single spherical small circle is topologically one closed curve. On one
Mollweide chart cut at λ=±π it projects to one or more polyline segments.
If the circle crosses the cut, segments terminate on BOTH left and right
boundaries (duplicated edge points); they must NOT be joined by an interior
chord or by a fake oval-rim bridge. They are also not a single closed
Euclidean ellipse in the map interior.

theta > 180° and μ = cos(θ)
---------------------------
Spherical angular distance is only in [0, 180°]. Using μ = cos(θ) for
θ ∈ (180°, 360°) yields the SAME cone as μ = cos(360°−θ) = cos(r) with
r = 360°−θ ∈ (0°, 180°), or equivalently a circle of radius (θ−180°) about −n0.
"""

from __future__ import annotations

import numpy as np


def wrap_lon_deg(lon):
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def lon_disp_from_ell(ell_deg):
    """Migkas / Galactic Mollweide: ℓ increases LEFT."""
    return -wrap_lon_deg(ell_deg)


def lb_to_n(l_deg, b_deg):
    l = np.deg2rad(np.asarray(l_deg, dtype=float))
    b = np.deg2rad(np.asarray(b_deg, dtype=float))
    return np.stack(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)],
        axis=-1,
    )


def n_to_lb(n):
    n = np.asarray(n, dtype=float)
    if n.ndim == 1:
        n = n / np.linalg.norm(n)
        b = float(np.rad2deg(np.arcsin(np.clip(n[2], -1.0, 1.0))))
        l = float(np.rad2deg(np.arctan2(n[1], n[0])) % 360.0)
        return l, b
    nn = n / np.linalg.norm(n, axis=-1, keepdims=True)
    b = np.rad2deg(np.arcsin(np.clip(nn[..., 2], -1.0, 1.0)))
    l = np.rad2deg(np.arctan2(nn[..., 1], nn[..., 0])) % 360.0
    return l, b


def orthonormal_frame(n0):
    """Right-handed ONB (e1, e2, n0) with e1 ⟂ n0, e2 = n0 × e1."""
    a = np.asarray(n0, dtype=float)
    a = a / np.linalg.norm(a)
    tmp = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(tmp, a)  # note: tmp × a so (e1,e2,a) is RH with e2 = a × e1
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    return e1, e2, a


def small_circle_points(l0_deg, b0_deg, r_deg, n=1441):
    """
    Exact locus n·n0 = cos(r) for r in [0, 180].
    Returns unit vectors shape (n, 3), closed (first==last conceptually via full 2π).
    """
    r = float(r_deg)
    if not (0.0 <= r <= 180.0):
        raise ValueError(f"angular radius must be in [0,180], got {r}")
    if r < 1e-12 or r > 180.0 - 1e-12:
        raise ValueError("degenerate pole/antipode; not a circle")
    n0 = lb_to_n(l0_deg, b0_deg).reshape(3)
    e1, e2, a = orthonormal_frame(n0)
    mu = np.cos(np.deg2rad(r))
    rho = np.sin(np.deg2rad(r))
    phi = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return (mu * a)[None, :] + rho * (
        np.outer(np.cos(phi), e1) + np.outer(np.sin(phi), e2)
    )


def small_circle_about_axis(axis, mu, n=1441):
    """Unit vectors on n·axis = mu with |mu|<1 (same ONB construction)."""
    a = np.asarray(axis, dtype=float)
    a = a / np.linalg.norm(a)
    mu = float(np.clip(mu, -1.0 + 1e-12, 1.0 - 1e-12))
    e1, e2, a = orthonormal_frame(a)
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    t = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return (mu * a)[None, :] + rho * (
        np.outer(np.cos(t), e1) + np.outer(np.sin(t), e2)
    )


def mollweide_xy(lon_rad, lat_rad, n_bisect=64):
    """
    Standard Mollweide. Auxiliary θ from 2θ + sin(2θ) = π sin(lat),
    solved by bisection (stable, residual ~1e-15). Matches matplotlib's
    formula; matplotlib's near-pole Taylor is slightly less accurate.
    """
    lon = np.asarray(lon_rad, dtype=float)
    lat = np.asarray(lat_rad, dtype=float)
    target = np.pi * np.sin(np.clip(lat, -np.pi / 2, np.pi / 2))
    lo = np.full_like(lat, -np.pi / 2, dtype=float)
    hi = np.full_like(lat, np.pi / 2, dtype=float)
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        g = 2.0 * mid + np.sin(2.0 * mid)
        lo = np.where(g < target, mid, lo)
        hi = np.where(g >= target, mid, hi)
    theta = 0.5 * (lo + hi)
    # Exact endpoints
    theta = np.where(lat >= np.pi / 2 - 1e-15, np.pi / 2, theta)
    theta = np.where(lat <= -np.pi / 2 + 1e-15, -np.pi / 2, theta)
    x = (2.0 * np.sqrt(2.0) / np.pi) * lon * np.cos(theta)
    y = np.sqrt(2.0) * np.sin(theta)
    return x, y


def split_antimeridian_lonlat(lon_rad, lat_rad, gap=np.pi):
    """
    Split a closed-ish polyline where display longitude jumps across the cut.
    Returns list of (lon, lat) segments; each segment is continuous in λ.
    """
    lon = np.asarray(lon_rad, dtype=float)
    lat = np.asarray(lat_rad, dtype=float)
    if lon.size < 2:
        return []
    # ensure we catch wrap between last and first for closed curves by
    # working on the open sample; caller may pass already-open 2π samples.
    d = np.diff(lon)
    cuts = np.where(np.abs(d) > gap)[0] + 1
    idx = np.concatenate([[0], cuts, [len(lon)]])
    segs = []
    for i0, i1 in zip(idx[:-1], idx[1:]):
        if i1 - i0 >= 2:
            segs.append((lon[i0:i1].copy(), lat[i0:i1].copy()))
    # The input circle is cyclic. If the artificial sampling join (last→first)
    # is NOT itself a seam crossing, merge the last and first linear pieces.
    # This removes a false third segment / tiny plotting gap.
    if len(segs) > 1 and abs(lon[0] - lon[-1]) <= gap:
        first_lo, first_la = segs[0]
        last_lo, last_la = segs[-1]
        merged = (
            np.concatenate([last_lo, first_lo]),
            np.concatenate([last_la, first_la]),
        )
        segs = [merged] + segs[1:-1]
    return segs


def snap_segment_to_seam(lon, lat, seam=np.pi):
    """
    If a segment endpoint is near ±seam (within a few degrees of the cut
    after a wrap split), pin that endpoint longitude exactly to ±seam so the
    stroke visibly meets the Mollweide boundary. Does not invent new points
    in the interior.
    """
    lon = lon.copy()
    lat = lat.copy()
    # After wrap-split, endpoints that crossed are typically near ±π already.
    if abs(abs(lon[0]) - seam) < np.deg2rad(8):
        lon[0] = np.sign(lon[0]) * seam if lon[0] != 0 else lon[0]
    if abs(abs(lon[-1]) - seam) < np.deg2rad(8):
        lon[-1] = np.sign(lon[-1]) * seam if lon[-1] != 0 else lon[-1]
    return lon, lat


def project_small_circle_mollweide(l0, b0, r_deg, n=2881, galactic_left=True):
    """
    Exact small circle → list of Mollweide (x,y) polyline segments.
    Seam-crossing circles become multiple segments ending on both edges.
    """
    pts = small_circle_points(l0, b0, r_deg, n=n)
    ell, b = n_to_lb(pts)
    if galactic_left:
        lon = np.deg2rad(lon_disp_from_ell(ell))
    else:
        lon = np.deg2rad(wrap_lon_deg(ell))
    lat = np.deg2rad(b)
    segs_lb = split_antimeridian_lonlat(lon, lat)
    out = []
    for lo, la in segs_lb:
        lo, la = snap_segment_to_seam(lo, la)
        x, y = mollweide_xy(lo, la)
        out.append((x, y, lo, la))
    return out


def theta_scan_to_spherical(theta_deg):
    """
    Map a 0..360 scan angle (along a great circle from +x) to the unique
    spherical small-circle radius in [0,180] implied by μ=cos(θ).
    Returns (r_from_plus_x, r_from_minus_x, mu).
    """
    th = float(theta_deg) % 360.0
    mu = float(np.cos(np.deg2rad(th)))
    r = float(np.rad2deg(np.arccos(np.clip(mu, -1.0, 1.0))))
    r_anti = float(np.rad2deg(np.arccos(np.clip(-mu, -1.0, 1.0))))
    return r, r_anti, mu
