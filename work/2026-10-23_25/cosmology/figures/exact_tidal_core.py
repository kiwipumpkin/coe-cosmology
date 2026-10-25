#!/usr/bin/env python3
"""Exact (no eps expansion) line-of-sight tidal term for the CoE model, sector-averaged like the analytic sheets.
Geometry (same as r2_vs_Dr_density_families.py): CoE at origin, observer at x_o = D xhat, source at x_s = x_o + d n,
d = r (vision shell), eps = r/D, mu = n.xhat = cos b cos b0 cos(l-l0) + sin b sin b0.
Radial peculiar velocity u(R) Rhat. Line-of-sight term: dH/H = [u(x_s).n - u(x_o).n]/(d H).
Write u(R) = u(D) f(R/D) (f = shape, f=1/rho^2 for a point mass). With K = -u(D)/(H D):
    dH/H = K T(mu, eps),   T = (1/eps) [ mu - f(rho_s) (mu+eps)/rho_s ],  rho_s = |x_s|/D = sqrt(1+2 eps mu+eps^2).
Point mass: T = [mu - (mu+eps) rho_s^-3]/eps = (3mu^2-1) + eps(9/2 mu - 15/2 mu^3) + O(eps^2)  (the series form).
Model: A/A_all = N [1 + K T]; sector average (90-deg longitude sectors, |b|>20, cos b weight) is linear, so
    <A/A_all>_sector(psi) = N [1 + K tau_0] + N K sum_m tau_m cos(m psi),  psi = l - l0,
tau_m = W_m x (Fourier cos coefficient of the latitude-averaged T), W_m = sin(m pi/4)/(m pi/4) (exact 90-deg boxcar).
Series equivalents: tau_0 = q0, tau_1 = eps o1, tau_2 = q2, tau_3 = eps o3.
"""
import numpy as np
d2r = np.pi / 180
BCUT = 20.0

def lat_nodes(nb):
    # Gauss-Legendre in b over [20,90] deg, weight cos b, both hemispheres
    x, w = np.polynomial.legendre.leggauss(nb)
    lo, hi = BCUT * d2r, np.pi / 2
    b = 0.5 * (hi - lo) * x + 0.5 * (hi + lo); wb = 0.5 * (hi - lo) * w * np.cos(b)
    b = np.concatenate([b, -b]); wb = np.concatenate([wb, wb]); return b, wb / wb.sum()

_cache = {}
def grid(nb, npsi):
    key = (nb, npsi)
    if key not in _cache:
        b, wb = lat_nodes(nb); psi = np.arange(npsi) * 2 * np.pi / npsi
        _cache[key] = (b, wb, psi)
    return _cache[key]

def T_point(mu, eps):
    rho2 = 1 + 2 * eps * mu + eps * eps
    return (mu - (mu + eps) * rho2 ** -1.5) / eps

def T_series(mu, eps):
    return (3 * mu * mu - 1) + eps * (4.5 * mu - 7.5 * mu ** 3)

def T_shape(mu, eps, fshape):
    rho = np.sqrt(1 + 2 * eps * mu + eps * eps)
    return (mu - fshape(rho) * (mu + eps) / rho) / eps

MMAX = 30
W = np.array([1.0] + [np.sin(m * np.pi / 4) / (m * np.pi / 4) for m in range(1, MMAX + 1)])

def taus(b0, eps, Tfun=T_point, nb=300, npsi=720, extra=None):
    """Return tau_0..tau_MMAX (signed cos coefficients after latitude + 90-deg sector averaging)."""
    b, wb, psi = grid(nb, npsi)
    B0 = b0 * d2r
    mu = np.cos(b)[:, None] * np.cos(B0) * np.cos(psi)[None, :] + np.sin(b)[:, None] * np.sin(B0)
    Tv = Tfun(mu, eps) if extra is None else Tfun(mu, eps, extra)
    prof = wb @ Tv                                   # latitude average vs psi
    F = np.fft.rfft(prof) / npsi
    t = np.empty(MMAX + 1); t[0] = F[0].real; t[1:] = 2 * F[1:MMAX + 1].real
    return t * W

def sector_profile_direct(b0, eps, Tfun=T_point, nb=300, step=0.25):
    """Brute-force check: explicit 90-deg sector averages at the 72 centres, |b|>20, cos b weight."""
    b, wb = lat_nodes(nb); B0 = b0 * d2r
    L = np.arange(-45 + step / 2, 45, step) * d2r
    out = []
    for c in np.arange(0, 360, 5.0) * d2r:
        mu = np.cos(b)[:, None] * np.cos(B0) * np.cos(c + L)[None, :] + np.sin(b)[:, None] * np.sin(B0)
        out.append((wb @ Tfun(mu, eps)).mean())
    return np.array(out)

# ---------- extended density: reuse the numeric code's profile + ST ideal-drag time integration ----------
import sys
sys.path.insert(0, "/workspace/cosmology/figures")
import r2_vs_Dr_density_families as NM
RGRID = np.linspace(0.004, 2.0, 500)
_shape_cache = {}
def shape_exp(x, fam="d_pure_exp", fu=0.0):
    """Pure exponential core rho = rho0 exp(-R/L) (family d of the numeric code), x = D/L.
    Returns (fshape, alpha, kappa_base): f(rho)=u(rho D)/u(D); alpha = K/kappa = -u(D)/(H D kappa) for mass_scale 1."""
    key = (fam, round(float(x), 12), fu)
    if key not in _shape_cache:
        p = [x] if fam == "d_pure_exp" else [x, fu]
        rg, u, kb = NM.u_table(fam, p, RGRID)
        uD = float(np.interp(1.0, rg, u))
        fr = u / uD
        rr = np.concatenate([[0.0], rg]); ff = np.concatenate([[0.0], fr])   # u -> 0 at the centre (extended mass)
        f = lambda rho, rr=rr, ff=ff: np.interp(rho, rr, ff)
        _shape_cache[key] = (f, -uD / kb, kb)
    return _shape_cache[key]

def shape_point():
    """Point mass through the same time integration (check: f = rho^-2, alpha ~ 2/3)."""
    kb = 1.0
    rg, u, _ = NM.u_table(None, None, RGRID, menc=lambda r: np.ones_like(r) / NM.G)
    uD = float(np.interp(1.0, rg, u))
    return rg, u / uD, -uD / kb
