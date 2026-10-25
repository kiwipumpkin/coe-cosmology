#!/usr/bin/env python3
"""
Full-loop opacity-gradient helpers (Chris 2026-09-25).

Instead of erasing the opposite half of each CoE ring, the FULL loop is drawn
with an opacity gradient:
  alpha_scale(psi) = (1 + cos psi) / 2
where psi is the ring azimuth measured from the θ-ray point, i.e. the point of
the ring most in the direction sign(sin θ_mid)·Ŷ (Ŷ from ring_basis, rolled for
Panels B/C; Panel A always roll 0).

With the ring parameterization of coe_half_loop_geometry
  n(phi) = mu*D + rho*(sin(phi)*Ŷ + cos(phi)*Ẑ)
the θ-ray point is phi0 = +π/2 (sinθ>0) or 3π/2 (sinθ<0), and
  cos(psi) = cos(phi - phi0) = sign(sin θ_mid) * sin(phi)
so alpha_scale = (1 + sign(sin θ_mid)*sin(phi)) / 2 — 1 at the θ-ray point,
0 at the diametrically opposite point (middle of the formerly erased half).
The existing half-loop keep rule (n·Ŷ)·sinθ ≥ 0 is exactly alpha_scale ≥ 1/2.

This module only ADDS functions; coe_half_loop_geometry is unchanged.
"""

from __future__ import annotations

import numpy as np

from coe_half_loop_geometry import ring_basis, theta_ray_sign

GRADIENT_RULE = (
    "full loops; alpha *= (1+cos psi)/2, psi = ring azimuth from the θ-ray point "
    "(ring point most along sign(sinθ_mid)·Ŷ, Ŷ rolled for B/C; A unrolled); "
    "equivalently alpha_scale = (1 + sign(sinθ_mid)·sin(phi))/2"
)


def theta_ray_phi0(theta_mid_deg):
    """phi of the θ-ray point on the ring (psi = 0)."""
    return 0.5 * np.pi if theta_ray_sign(theta_mid_deg) > 0 else 1.5 * np.pi


def gradient_alpha_scale(phi, theta_mid_deg):
    """(1 + cos psi)/2 with psi = phi - phi0; in [0, 1]."""
    psi = np.asarray(phi, dtype=float) - theta_ray_phi0(theta_mid_deg)
    return 0.5 * (1.0 + np.cos(psi))


def full_phi_from_theta_ray(theta_mid_deg, n):
    """Closed phi samples starting AND ending at the antipode of the θ-ray point
    (psi = -180°..+180°), so the only wrap-around point is the invisible one."""
    phi0 = theta_ray_phi0(theta_mid_deg)
    return phi0 + np.linspace(-np.pi, np.pi, n)


def small_circle_full_gradient(D, mu, theta_mid_deg, roll_deg=0.0, n=721):
    """Full closed small circle n·D = mu with per-point alpha_scale.
    Returns (pts[n,3], phi[n], scale[n])."""
    a = np.asarray(D, dtype=float)
    a = a / np.linalg.norm(a)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    y_hat, z_hat = ring_basis(a, roll_deg=roll_deg)
    phi = full_phi_from_theta_ray(theta_mid_deg, n)
    pts = (mu * a)[None, :] + rho * (
        np.outer(np.sin(phi), y_hat) + np.outer(np.cos(phi), z_hat)
    )
    return pts, phi, gradient_alpha_scale(phi, theta_mid_deg)


def foreshortened_full_ellipse_gradient(mu, theta_mid_deg, n=361):
    """Panel A full foreshortened ellipse (roll 0 always) with alpha_scale.
    ex = mu + 0.045*rho*cos(phi), ey = rho*sin(phi)."""
    mu = float(mu)
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    phi = full_phi_from_theta_ray(theta_mid_deg, n)
    ex = mu + 0.045 * rho * np.cos(phi)
    ey = rho * np.sin(phi)
    return ex, ey, gradient_alpha_scale(phi, theta_mid_deg)


def segment_scale(scale):
    """Per-segment alpha scale = mean of the two endpoint scales."""
    s = np.asarray(scale, dtype=float)
    return 0.5 * (s[:-1] + s[1:])


def rgba_array(color, alphas):
    from matplotlib.colors import to_rgba

    r, g, b, _ = to_rgba(color)
    al = np.clip(np.asarray(alphas, dtype=float), 0.0, 1.0)
    out = np.empty((al.size, 4))
    out[:, 0], out[:, 1], out[:, 2], out[:, 3] = r, g, b, al
    return out


def quantized_runs(seg_scale, valid=None, n_levels=64):
    """Group consecutive segments into runs sharing one quantized alpha level.

    Drawing each run as ONE polyline / strip polygon (instead of one element per
    segment) removes the beading at segment joints and most anti-aliasing seams,
    while the gradient stays visually smooth (step <= 1/n_levels of base alpha).
    Runs never include invalid segments (e.g. Mollweide antimeridian crossings)
    and fully transparent (level 0) segments are skipped.
    Returns list of (i0, i1, level) meaning segments i0..i1-1, i.e. points i0..i1.
    """
    s = np.asarray(seg_scale, dtype=float)
    lev = np.round(np.clip(s, 0.0, 1.0) * n_levels) / n_levels
    ok = lev > 0 if valid is None else (lev > 0) & np.asarray(valid, dtype=bool)
    runs = []
    i, n = 0, s.size
    while i < n:
        if not ok[i]:
            i += 1
            continue
        j = i + 1
        while j < n and ok[j] and lev[j] == lev[i]:
            j += 1
        runs.append((i, j, float(lev[i])))
        i = j
    return runs
