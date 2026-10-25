"""Shared pieces for dq_vs_model_match.py: the 'ideal drag' physical model copied
verbatim (arithmetic identical) from make_model_param_heatmaps.py, plus harmonic tools.
Importing make_model_param_heatmaps would re-run all its grids, so the model is copied."""
import sys
import numpy as np
sys.path.insert(0, "/workspace/cosmology/density_local_match")
from run_density_local_match import params_twoscale_local_match, m_enc_at_r, H, G, D  # noqa: E402

CSV = "/workspace/cosmology/migkas/fig6_digitized.csv"
FU, ETA, NT = 0.94, 1e-4, 800
T = 2.0 / (5.0 * H)
dat = np.genfromtxt(CSV, delimiter=",", names=True)
ell_m, A_m = dat["l_deg"].astype(float), dat["A_over_Aall"].astype(float)
A_lo, A_hi = dat["A_err_low"].astype(float), dat["A_err_high"].astype(float)
best = np.load("/workspace/cosmology/map_compare/amp_mean_param_noise/ideal_fit_M_d_D_L_migkas.npz")
MS0, EPS0, X0, L0DEG = float(best["mass_scale"]), float(best["eps"]), float(best["x"]), float(best["phi0_deg"])

_ucache = {}
def u_base(x, nrad=50):
    key = (round(float(x), 6), nrad)
    if key in _ucache:
        return _ucache[key]
    p = params_twoscale_local_match(float(x), FU)
    r_grid = np.linspace(0.2 * D, 2.0 * D, nrad)
    chi = r_grid / T ** 0.4
    t = np.linspace(ETA * T, T, NT)
    dt = t[1] - t[0]
    u = np.zeros(nrad)
    for ti in t[1:]:
        Hi = 2.0 / (5.0 * ti)
        ri = chi * ti ** 0.4
        gi = G * m_enc_at_r(ri, p) / np.maximum(ri, 1e-30) ** 2
        u = u + dt * (-Hi * u - gi)
    _ucache[key] = (r_grid, u, float(p["kappa"]))
    return _ucache[key]

def kbase(x): return u_base(x)[2]

def boxcar_dl90(y):
    n = len(y); half = int(round(45.0 / (360.0 / n)))
    k = np.ones(2 * half + 1) / (2 * half + 1)
    return np.convolve(np.pad(y, (half, half), mode="wrap"), k, mode="valid")

def model_fine(kappa, eps, x, l0, b0, n_phi=181, nrad=50, smooth=True):
    """Returns (l_deg grid, A on grid). kappa -> mass_scale = kappa/kappa_base(x)."""
    r_grid, ub, kb = u_base(x, nrad)
    ms = kappa / kb
    u = lambda r: ms * np.interp(r, r_grid, ub)
    l_deg = np.linspace(0.0, 360.0, n_phi, endpoint=False)
    mu = np.cos(np.radians(b0)) * np.cos(np.radians(l_deg - l0))
    s = np.sqrt(np.clip(1 - mu ** 2, 0, None))
    n_hat = np.stack([mu, s, np.zeros_like(mu)], axis=-1)
    d = eps * D
    r_vec = np.array([D, 0, 0]) + d * n_hat
    r = np.linalg.norm(r_vec, axis=-1)
    u_r = u(r)[:, None] * r_vec / r[:, None]
    u_D = u(D) * np.array([1.0, 0, 0])
    Hobs = H + np.sum((u_r - u_D) * n_hat, axis=-1) / d
    AF = Hobs / H; AF = AF / AF.mean()
    if smooth:
        AF = boxcar_dl90(AF); AF = AF / AF.mean()
    return l_deg, AF

def model_at_data(kappa, eps, x, l0, b0, **kw):
    l, A = model_fine(kappa, eps, x, l0, b0, **kw)
    return np.interp(ell_m % 360.0, l, A, period=360.0)

# ---------------- harmonic tools ----------------
def design(ell, mmax=2):
    L = np.radians(ell); cols = [np.ones_like(L)]
    for m in range(1, mmax + 1):
        cols += [np.cos(m * L), np.sin(m * L)]
    return np.stack(cols, axis=1)

def lsq(ell, A, ms=(1, 2)):
    L = np.radians(ell); cols = [np.ones_like(L)]
    for m in ms:
        cols += [np.cos(m * L), np.sin(m * L)]
    X = np.stack(cols, 1)
    beta, *_ = np.linalg.lstsq(X, A, rcond=None)
    pred = X @ beta
    return beta, X, pred

def amp_phase(c, s, m):
    return float(np.hypot(c, s)), float((np.degrees(np.arctan2(s, c)) / m) % (360.0 / m))

def r2(A, pred, Aref=None):
    Aref = A if Aref is None else Aref
    return float(1 - ((A - pred) ** 2).sum() / ((Aref - Aref.mean()) ** 2).sum())

def harm71(A):
    """Same LSQ as the data fit (a0 + dipole + quadrupole at the 71 data longitudes)."""
    b, _, _ = lsq(ell_m, A)
    a1, l1 = amp_phase(b[1], b[2], 1); a2, l2 = amp_phase(b[3], b[4], 2)
    return dict(a0=float(b[0]), c1=float(b[1]), s1=float(b[2]), c2=float(b[3]), s2=float(b[4]),
                a1=a1, l1=l1, a2=a2, l2=l2)

def harm_fine(l_deg, A, mmax=8):
    """Exact Fourier on the uniform full-circle model grid."""
    n = len(A); F = np.fft.rfft(A) / n
    out = dict(a0=float(F[0].real))
    for m in range(1, mmax + 1):
        c, s = 2 * F[m].real, -2 * F[m].imag
        out[f"c{m}"], out[f"s{m}"] = float(c), float(s)
        out[f"a{m}"], out[f"l{m}"] = amp_phase(c, s, m)
    P = np.abs(F[1:]) ** 2
    out["frac_power_m1"] = float(P[0] / P.sum()); out["frac_power_m2"] = float(P[1] / P.sum())
    out["frac_power_m_ge3"] = float(P[2:].sum() / P.sum())
    return out

def model_harm(kappa, eps, x, l0, b0, **kw):
    l, A = model_fine(kappa, eps, x, l0, b0, **kw)
    h = harm_fine(l, A); A71 = np.interp(ell_m, l, A, period=360.0)
    h71 = harm71(A71)
    return dict(fine=h, lsq71=h71, R2_direct=r2(A_m, A71),
                R2_corr2=float(np.corrcoef(A_m, A71)[0, 1] ** 2))

def nelder_mead(f, x0, step, tol=1e-10, maxit=4000):
    n = len(x0); pts = [np.array(x0, float)]
    for i in range(n):
        p = np.array(x0, float); p[i] += step[i]; pts.append(p)
    vals = [f(p) for p in pts]
    for it in range(maxit):
        o = np.argsort(vals); pts = [pts[i] for i in o]; vals = [vals[i] for i in o]
        if abs(vals[-1] - vals[0]) < tol and it > 50: break
        c = np.mean(pts[:-1], 0); xr = c + (c - pts[-1]); fr = f(xr)
        if fr < vals[0]:
            xe = c + 2 * (c - pts[-1]); fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = c + 0.5 * (pts[-1] - c); fc = f(xc)
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:
                pts = [pts[0]] + [pts[0] + 0.5 * (p - pts[0]) for p in pts[1:]]
                vals = [vals[0]] + [f(p) for p in pts[1:]]
    i = int(np.argmin(vals)); return pts[i], vals[i]
