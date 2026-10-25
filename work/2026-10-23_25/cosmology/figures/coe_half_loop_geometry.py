#!/usr/bin/env python3
"""
Shared θ-ray-side CoE ring geometry (Chris 2026-09-23 whiteout — CORRECTED).

Half-loop rule (θ-ray side — NOT outer limb):
  Panel A foreshortened_ellipse_xy uses
    ex = mu + 0.045*rho*cos(phi);  ey = rho*sin(phi)
  θ-rays go to (cos θ, sin θ). Keep the half on that same side:
    keep ey * sin(θ_mid) >= 0
  i.e. θ∈(0,180°) → TOP (y≥0); θ∈(180,360°) → BOTTOM (y≤0).

3D ring parameterization matching that phi:
  n(phi) = mu * D + rho * (sin(phi) * Y_hat + cos(phi) * Z_hat)
  so (n · Y_hat) = rho * sin(phi)
  Keep: (n · Y_hat) * sin(θ_mid) >= 0

DISCARD the old outer-limb rule sign(μ)·cos(φ)≥0.

Axis pose (Panels B/C): anti-CoE (+D) at (0°, +b_tilt), CoE at (180°, -b_tilt).
Panel A: model frame, CoE axis horizontal (no Galactic yet).
"""

from __future__ import annotations

import numpy as np

BAND_HALF = 9.0
BAND_FULL = 2.0 * BAND_HALF
STRIP_THETA_MAX = [30.0, 185.0]
STRIP_THETA_MIN = [115.0, 265.0]

MIGKAS_MAXES = [(17.0, -9.0), (185.0, 0.0)]
MIGKAS_MINS = [(115.0, 0.0), (274.0, -22.0)]

HALF_LOOP_RULE = (
    "θ-ray side: keep (n·Ŷ)*sin(θ_mid)>=0 "
    "(Panel A: ey*sin(θ_mid)>=0; TOP if sinθ>0, BOTTOM if sinθ<0)"
)


def wrap_lon_deg(lon):
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def lon_disp_from_ell(ell_deg):
    return -wrap_lon_deg(ell_deg)


def lb_to_n(l_deg, b_deg):
    l = np.deg2rad(float(l_deg))
    b = np.deg2rad(float(b_deg))
    return np.array(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], dtype=float
    )


def n_to_lb(n):
    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    b = float(np.rad2deg(np.arcsin(np.clip(n[2], -1.0, 1.0))))
    l = float(np.rad2deg(np.arctan2(n[1], n[0])) % 360.0)
    return l, b


def mu_of_theta(theta_deg):
    return float(np.cos(np.deg2rad(theta_deg)))


def theta_ray_sign(theta_mid_deg, eps=1e-9):
    """
    sign(sin θ_mid): +1 → keep TOP (ey≥0 / n·Ŷ≥0);
    −1 → keep BOTTOM (ey≤0 / n·Ŷ≤0).
    """
    s = float(np.sin(np.deg2rad(float(theta_mid_deg))))
    if s > eps:
        return 1.0
    if s < -eps:
        return -1.0
    return 1.0


# Back-compat alias (unused by θ-ray rule; kept so old imports don't explode)
def mu_sign(mu, eps=1e-9):
    if mu > eps:
        return 1.0
    if mu < -eps:
        return -1.0
    return 1.0


def ring_basis(D, roll_deg=0.0):
    """
    Orthonormal Y_hat, Z_hat ⊥ D matching Panel A:
      n(phi) = mu*D + rho*(sin(phi)*Y_hat + cos(phi)*Z_hat)
    In Panel A model frame, +Y is the θ-ray plane "up" (direction of θ-rays).
    roll_deg: CCW about +D when looking along +D (right-hand).
    """
    a = np.asarray(D, dtype=float)
    a = a / np.linalg.norm(a)
    # Prefer Galactic +z projected into plane ⊥ D as unrolled Y ("up")
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    y0 = tmp - np.dot(tmp, a) * a
    y0 /= np.linalg.norm(y0)
    # Right-handed: X=D, Y=y0, Z = X × Y
    z0 = np.cross(a, y0)
    z0 /= np.linalg.norm(z0)
    r = np.deg2rad(float(roll_deg))
    c, s = np.cos(r), np.sin(r)
    # CCW about +D looking along +D: Y' = c Y - s Z, Z' = s Y + c Z
    y = c * y0 - s * z0
    z = s * y0 + c * z0
    return y, z


def _phi_keep_interval(theta_mid_deg, n, endpoint=True):
    """Open 180° phi arc on the θ-ray side: sin(phi)*sin(θ_mid) >= 0."""
    sgn = theta_ray_sign(theta_mid_deg)
    if sgn >= 0:
        # TOP: ey = rho*sin(phi) >= 0 → phi ∈ [0, π]
        return np.linspace(0.0, np.pi, n, endpoint=endpoint)
    # BOTTOM: ey <= 0 → phi ∈ [π, 2π]
    return np.linspace(np.pi, 2.0 * np.pi, n, endpoint=endpoint)


def small_circle_theta_ray_half(D, mu, theta_mid_deg, roll_deg=0.0, n=721, endpoint=True):
    """
    Unit vectors on the θ-ray-side half of the small circle n·D = mu.
    Returns points ordered along the open 180° arc (not a closed loop).
    Keep: (n · Ŷ) * sin(θ_mid) >= 0.
    """
    a = np.asarray(D, dtype=float)
    a = a / np.linalg.norm(a)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    y_hat, z_hat = ring_basis(a, roll_deg=roll_deg)
    phi = _phi_keep_interval(theta_mid_deg, n, endpoint=endpoint)
    pts = (
        (mu * a)[None, :]
        + rho * (np.outer(np.sin(phi), y_hat) + np.outer(np.cos(phi), z_hat))
    )
    return pts, phi


# Alias: old name now means θ-ray half (NOT outer limb)
def small_circle_outer_half(D, mu, roll_deg=0.0, n=721, endpoint=True, theta_mid_deg=None):
    """
    Back-compat wrapper. Prefer small_circle_theta_ray_half.
    If theta_mid_deg is omitted, infer from mu via θ = arccos(mu) in [0,180]
    (WRONG for θ>180 bands — always pass theta_mid_deg explicitly).
    """
    if theta_mid_deg is None:
        theta_mid_deg = float(np.rad2deg(np.arccos(np.clip(float(mu), -1.0, 1.0))))
    return small_circle_theta_ray_half(
        D, mu, theta_mid_deg, roll_deg=roll_deg, n=n, endpoint=endpoint
    )


def small_circle_full(D, mu, roll_deg=0.0, n=721):
    """Full small circle with the SAME (Y,Z) basis (for reference / tests)."""
    a = np.asarray(D, dtype=float)
    a = a / np.linalg.norm(a)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    y_hat, z_hat = ring_basis(a, roll_deg=roll_deg)
    phi = np.linspace(0.0, 2.0 * np.pi, n)
    pts = (
        (mu * a)[None, :]
        + rho * (np.outer(np.sin(phi), y_hat) + np.outer(np.cos(phi), z_hat))
    )
    return pts, phi


def half_keep_mask_on_full(pts, D, theta_mid_deg, roll_deg=0.0):
    """Boolean keep mask: (n · Ŷ) * sin(θ_mid) >= 0."""
    a = np.asarray(D, dtype=float)
    a = a / np.linalg.norm(a)
    y_hat, _ = ring_basis(a, roll_deg=roll_deg)
    dots = np.asarray(pts, dtype=float) @ y_hat
    return theta_ray_sign(theta_mid_deg) * dots >= -1e-12


def foreshortened_half_ellipse_xy(mu, theta_mid_deg, roll_deg=0.0, n=256):
    """
    Panel A θ-ray-side half-arc of the foreshortened ellipse at n·x = mu.
    Returns (ex, ey) for the open 180° arc only.
    Keep: ey * sin(θ_mid) >= 0.
    roll_deg≠0 rotates the keep window in phi (rarely used; this deliverable is roll=0).
    """
    mu = float(mu)
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    r = np.deg2rad(float(roll_deg))
    # Sample psi in keep interval, then phi = psi + roll
    psi = _phi_keep_interval(theta_mid_deg, n, endpoint=True)
    phi = psi + r
    ex = mu + 0.045 * rho * np.cos(phi)
    ey = rho * np.sin(phi)
    return ex, ey, rho


def foreshortened_half_band_polygon(mu_lo, mu_hi, theta_mid_deg, roll_deg=0.0, n=128):
    """
    Closed polygon covering the θ-ray-side half-annulus between two mu edges
    in the Panel A foreshortened drawing. Used as the highlight fill.
    """
    ex_a, ey_a, _ = foreshortened_half_ellipse_xy(
        mu_lo, theta_mid_deg, roll_deg=roll_deg, n=n
    )
    ex_b, ey_b, _ = foreshortened_half_ellipse_xy(
        mu_hi, theta_mid_deg, roll_deg=roll_deg, n=n
    )
    xs = np.concatenate([ex_a, ex_b[::-1]])
    ys = np.concatenate([ey_a, ey_b[::-1]])
    return np.column_stack([xs, ys])


def great_circle_through_poles(D, n=361):
    a = np.asarray(D, dtype=float)
    a = a / np.linalg.norm(a)
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


def celestial_equator_lb(n=721):
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


def axis_pose(tilt_b_deg=15.0):
    anti_l, anti_b = 0.0, float(tilt_b_deg)
    coe_l, coe_b = 180.0, -float(tilt_b_deg)
    n_anti = lb_to_n(anti_l, anti_b)
    n_coe = lb_to_n(coe_l, coe_b)
    return {
        "anti_lb": (anti_l, anti_b),
        "coe_lb": (coe_l, coe_b),
        "n_anti": n_anti,
        "n_coe": n_coe,
        "D": n_anti,
        "tilt_b_deg": float(tilt_b_deg),
    }


def band_edges(centers, half=BAND_HALF):
    return [(c - half, c + half) for c in centers]


def mu_pair_for_band(lo, hi):
    return float(np.cos(np.deg2rad(lo))), float(np.cos(np.deg2rad(hi)))
