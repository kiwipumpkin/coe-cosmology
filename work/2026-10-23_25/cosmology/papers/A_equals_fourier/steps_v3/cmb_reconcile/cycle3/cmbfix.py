"""Exact-LTB CMB for models with a bang-time bump (decaying mode). Same null-geodesic equations as Chris's rays.py (ltb.rhs,
Alnes-Amarzguioui), but each ray stops on an event (local density > 30 x the recombination-density guess, or local age < 1e-7 t0)
instead of a fixed t_end, so rays through young regions never step past the local bang. Emission at fixed local density rho_rec
(the recombination condition), rho_rec set so the z monopole is 1089; dipole removed with analyze.boost_multipoles.
Null-condition error is reported per ray."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
from ltb_common import rays, analyze, LT
ZREC = 1089.0
def tB(r, p): return p[6] * np.exp(-(r / p[7]) ** 2) if p[7] > 0 else 0.0
def trace(p, t0, D, psi, rho_stop):
    y0, J = rays.start(p, t0, D, psi)
    ev1 = lambda x, y: np.log(LT.density(y[0], t0 * np.exp(-x), p)) - np.log(rho_stop); ev1.terminal = True
    ev2 = lambda x, y: (t0 * np.exp(-x) - tB(y[0], p)) / t0 - 1e-7; ev2.terminal = True
    s = solve_ivp(lambda x, y: rays.ltb.rhs(x, y, p, J, t0), (0.0, np.log(1e7)), y0, method="DOP853", rtol=1e-12, atol=1e-13,
                  dense_output=True, events=(ev1, ev2))
    s.J = J; return s
def z_at_age(sol, p, t0, tau):
    f = lambda x: np.log(max(t0 * np.exp(-x) - tB(sol.sol(x)[0], p), 1e-300)) - np.log(tau)
    xs = np.linspace(sol.t[0], sol.t[-1], 400); fv = np.array([f(x) for x in xs])
    i = np.where(np.diff(np.sign(fv)) != 0)[0]
    if len(i) == 0: return np.nan
    x = brentq(f, xs[i[0]], xs[i[0] + 1], xtol=1e-14); return -sol.sol(x)[2] - 1.0
def multip(mu, w, z):
    T = 1 / (1 + z); T = T / (0.5 * np.sum(w * T)); c = rays.legendre_coeffs(mu, w, T - 1, 8)
    beta, cb = analyze.boost_multipoles(mu, w, T, lmax=8); return c, beta, cb
def cmb_age(p, t0, D, N=24, psis=None, weights=None, sols=None):
    """(a) emission at fixed local age tau_rec = t - t_B(r) (fixed proper time since the local bang), monopole z = 1089."""
    if psis is None: mu, w = leggauss(N); psis = np.arccos(mu)
    else: mu, w = np.cos(psis), weights
    if sols is None: sols = [trace(p, t0, D, ps, 1e3 * p[8] * (1 + ZREC) ** 3) for ps in psis]
    zs = lambda lt: np.array([z_at_age(s, p, t0, np.exp(lt)) for s in sols])
    def g(lt):
        z = zs(lt); ok = np.isfinite(z)
        if not ok.any(): return np.nan
        m = np.mean(z[ok]) if w is None else np.sum(w[ok] * z[ok]) / np.sum(w[ok]); return np.log1p(m) - np.log1p(ZREC)
    grid = np.linspace(np.log(t0 * 2e-5), np.log(t0 * 1e-3), 25); gv = np.array([g(x) for x in grid])
    k = np.where(np.isfinite(gv[:-1]) & np.isfinite(gv[1:]) & (gv[:-1] * gv[1:] <= 0))[0]
    lt = brentq(g, grid[k[0]], grid[k[0] + 1], xtol=1e-13); z = zs(lt)
    out = dict(psi=psis, mu=mu, z=z, tau_rec_over_t0=float(np.exp(lt) / t0), nerr=np.array([rays.null_err(p, t0, s) for s in sols]),
               nsteps=[len(s.t) for s in sols], stop=[int(np.argmax([len(e) for e in s.t_events])) if s.status == 1 else -1 for s in sols], sols=sols)
    if w is not None and np.all(np.isfinite(z)): out["c"], out["beta"], out["cb"] = multip(mu, w, z)
    return out
def cmb(p, t0, D, N=24, psis=None, weights=None):
    if psis is None: mu, w = leggauss(N); psis = np.arccos(mu)
    else: mu, w = np.cos(psis), weights
    rho_guess = p[8] * (1 + ZREC) ** 3
    sols = [trace(p, t0, D, ps, 30 * rho_guess) for ps in psis]
    zs = lambda lr: np.array([rays.z_at_rho(s, p, t0, np.exp(lr))[0] for s in sols])
    def g(lr):
        z = zs(lr); ok = np.isfinite(z)
        if not ok.any(): return np.nan
        m = np.mean(z[ok]) if w is None else np.sum(w[ok] * z[ok]) / np.sum(w[ok])
        return np.log1p(m) - np.log1p(ZREC)
    lr = brentq(g, np.log(rho_guess / 20), np.log(rho_guess * 20), xtol=1e-13); z = zs(lr)
    ne = np.array([rays.null_err(p, t0, s) for s in sols])
    out = dict(psi=psis, mu=mu, z=z, rho_rec=float(np.exp(lr)), nerr=ne, nsteps=[len(s.t) for s in sols])
    if w is not None and np.all(np.isfinite(z)):
        T = 1 / (1 + z); T = T / (0.5 * np.sum(w * T))
        out["c"] = rays.legendre_coeffs(mu, w, T - 1, 8)
        out["beta"], out["cb"] = analyze.boost_multipoles(mu, w, T, lmax=8)
    return out
