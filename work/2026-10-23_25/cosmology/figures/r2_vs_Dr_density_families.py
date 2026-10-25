#!/usr/bin/env python3
"""R^2 vs D/r (d = r, so d/D = r/D) for the Sedov-Taylor + Gauss 'ideal drag' CoE model,
with full re-optimisation of kappa, D/L (or density scale), l0, b0 (|b0|<=30 deg) and the
density-shape parameters, for four density families.

Physics / conventions are IDENTICAL to make_model_param_heatmaps.py::model_A:
  * ST background R ∝ t^{2/5}, H=2/(5t), T=2/(5H); units H=G=D=1.
  * du/dt + H u = -G M(<r)/r^2 integrated (explicit Euler, 800 steps, t_i=1e-4 T) along
    homologous paths r_i = chi t^{2/5}; static physical-radius mass profile (as in the code).
  * H_obs(n) = H + [(u(r)-u(D)).n]/d ; A = H_obs/H mean-renormalised, 90-deg boxcar (181-pt grid),
    renormalised, interpolated to the 71 Migkas longitudes; mu = cos b0 cos(l-l0).
  * R^2 PRIMARY = 'direct' (1-SS_res/SS_tot, no rescale) -> reproduces the 0.776 best fit.
    corr2 (Pearson^2) also reported.
  * Local-density convention (same as code): the density SHAPE is normalised so that
    rho(D) = 0.3 rho_crit at mass_scale = 1; kappa = G M(<D)/(D^3 H^2) is then FREE
    (kappa = mass_scale * kappa_base), so the implied local density is
    rho(D)/rho_crit = 0.3 * kappa / kappa_base. This implied value is reported.
Changes vs the heat-map code (numerically checked, see JSON 'sanity'):
  * M(<r) from cumulative quadrature of rho on a fine grid (any family), not analytic.
  * u(r) table on r in [0.05, 2.0] D (79 pts) instead of [0.2, 2.0] D (50 pts), because
    d/D up to 0.833 needs r down to 0.167 D.
  * Euler recursion written as an exact linear sum (same recursion, different summation order).
Families (rho/rho(D) = (1-f_u) core(R)/core(D) + f_u):
  a) exp core  + floor : core = exp(-R/L),              params D/L, f_u in [0,0.99]
  b) power law + floor : core = (1+R/R_c)^-n,           params n in [0,10], R_c/D, f_u
  c) Gauss core + floor: core = exp(-R^2/(2 s^2)),      params D/s, f_u
  d) pure exponential  : core = exp(-R/L), f_u = 0,     params D/L
"""
import json, sys, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize

OUT = "/workspace/cosmology/figures"
CSV = "/workspace/cosmology/migkas/fig6_digitized.csv"
H = G = D = 1.0
OM = 0.3
RHOC = 3 * H ** 2 / (8 * np.pi * G)
ETA, NT = 1e-4, 800
T = 2.0 / (5.0 * H)
dat = np.genfromtxt(CSV, delimiter=",", names=True)
ell_m, A_m = dat["l_deg"].astype(float), dat["A_over_Aall"].astype(float)
SST = float(((A_m - A_m.mean()) ** 2).sum())
R_PHYS_GPC = 14.3
GLY_PER_GPC = 3.26156

# ---- time integration weights (exact rewrite of the Euler loop) ----
_t = np.linspace(ETA * T, T, NT); _dt = _t[1] - _t[0]
_tk = _t[1:]; _Hk = 2.0 / (5.0 * _tk)
_fac = 1.0 - _dt * _Hk                      # u <- u*fac_k - dt*g_k
_W = np.append(np.cumprod(_fac[::-1])[::-1][1:], 1.0)   # prod_{j>k} fac_j
_W = -_dt * _W

RG_LO, RG_HI, NRAD = 0.05, 2.0, 79
MF = np.linspace(0.0, 2.05, 8201)         # fine grid for M(<r)

FAM = {
    "a_exp_floor":   dict(names=["D_over_L", "f_u"], lo=[0.05, 0.0], hi=[40.0, 0.99], log=[True, False]),
    "b_powerlaw_floor": dict(names=["n", "Rc_over_D", "f_u"], lo=[0.0, 0.005, 0.0], hi=[10.0, 20.0, 0.99], log=[False, True, False]),
    "c_gauss_floor": dict(names=["D_over_sigma", "f_u"], lo=[0.05, 0.0], hi=[20.0, 0.99], log=[True, False]),
    "d_pure_exp":    dict(names=["D_over_L"], lo=[0.05], hi=[40.0], log=[True]),
}

def core(fam, p, r):
    if fam in ("a_exp_floor", "d_pure_exp"):
        return np.exp(-p[0] * r)
    if fam == "b_powerlaw_floor":
        return (1.0 + r / p[1]) ** (-p[0])
    if fam == "c_gauss_floor":
        return np.exp(-0.5 * (p[0] * r) ** 2)
    raise ValueError(fam)

def fu_of(fam, p):
    return 0.0 if fam == "d_pure_exp" else p[-1]

def rho_shape(fam, p, r):
    fu = fu_of(fam, p)
    return OM * RHOC * ((1 - fu) * core(fam, p, r) / core(fam, p, D) + fu)

def menc_table(fam, p):
    rho = rho_shape(fam, p, MF)
    f = 4 * np.pi * MF ** 2 * rho
    M = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(MF))])
    return M

def u_table(fam, p, r_grid=None, menc=None):
    r_grid = np.linspace(RG_LO, RG_HI, NRAD) if r_grid is None else r_grid
    chi = r_grid / T ** 0.4
    ri = _tk[:, None] ** 0.4 * chi[None, :]
    if menc is None:
        M = menc_table(fam, p); m = np.interp(ri, MF, M)
        kb = float(np.interp(D, MF, M)) * G / (D ** 3 * H ** 2)
    else:
        m, kb = menc(ri), None
    g = G * m / np.maximum(ri, 1e-30) ** 2
    u = _W @ g
    return r_grid, u, kb

def boxcar_dl90(y):
    n = len(y); half = int(round(45.0 / (360.0 / n)))
    k = np.ones(2 * half + 1) / (2 * half + 1)
    return np.convolve(np.pad(y, (half, half), mode="wrap"), k, mode="valid")

LG = np.linspace(0.0, 360.0, 181, endpoint=False)
def model_fine(ms, eps, l0, b0, tab):
    r_grid, ub = tab[0], tab[1]
    mu = np.cos(np.radians(b0)) * np.cos(np.radians(LG - l0))
    s = np.sqrt(np.clip(1 - mu ** 2, 0, None))
    n_hat = np.stack([mu, s, np.zeros_like(mu)], axis=-1)
    d = eps * D
    r_vec = np.array([D, 0, 0]) + d * n_hat
    r = np.linalg.norm(r_vec, axis=-1)
    ur = ms * np.interp(r, r_grid, ub); uD = ms * np.interp(D, r_grid, ub)
    Hobs = H + (ur * np.sum(r_vec * n_hat, -1) / r - uD * mu) / d
    AF = Hobs / H; AF = AF / AF.mean(); AF = boxcar_dl90(AF); AF = AF / AF.mean()
    return AF

def at_data(AF):
    return np.interp(ell_m % 360.0, LG, AF, period=360.0)

def r2_direct(m): return float(1 - ((A_m - m) ** 2).sum() / SST)
def r2_corr(m): return float(np.corrcoef(m, A_m)[0, 1] ** 2)

# ---- parameter vector: [log kappa, l0, b0, density...] ----
def unpack(fam, v):
    F = FAM[fam]
    dens = [np.exp(v[3 + i]) if F["log"][i] else v[3 + i] for i in range(len(F["names"]))]
    return np.exp(v[0]), v[1], v[2], np.array(dens)

def bounds(fam):
    F = FAM[fam]
    b = [(np.log(0.005), np.log(20.0)), (-180.0, 180.0), (-30.0, 30.0)]
    for lo, hi, lg in zip(F["lo"], F["hi"], F["log"]):
        b.append((np.log(lo), np.log(hi)) if lg else (lo, hi))
    return b

def evaluate(fam, eps, v, strict=False):
    kap, l0, b0, dens = unpack(fam, v)
    tab = u_table(fam, dens)
    if strict: kap = tab[2]
    AF = model_fine(kap / tab[2], eps, l0, b0, tab)
    return AF, kap, tab[2]

def objective(v, fam, eps, strict=False):
    b = bounds(fam)
    pen = 0.0
    vv = np.array(v, float)
    for i, (lo, hi) in enumerate(b):
        if i == 1: continue
        if vv[i] < lo: pen += (lo - vv[i]) * 10; vv[i] = lo
        if vv[i] > hi: pen += (vv[i] - hi) * 10; vv[i] = hi
    try:
        AF, _, _ = evaluate(fam, eps, vv, strict)
        m = at_data(AF)
        if not np.all(np.isfinite(m)): return 10.0
        return -r2_direct(m) + pen
    except Exception:
        return 10.0

def harmonics(AF):
    n = len(AF); F = np.fft.rfft(AF) / n; out = {}
    P = np.abs(F[1:]) ** 2
    for m in (1, 2, 3):
        c, s = 2 * F[m].real, -2 * F[m].imag
        out[f"a{m}"] = float(np.hypot(c, s)); out[f"l{m}_deg"] = float((np.degrees(np.arctan2(s, c)) / m) % (360.0 / m))
    out["frac_power_m1"] = float(P[0] / P.sum()); out["frac_power_m2"] = float(P[1] / P.sum()); out["frac_power_m3"] = float(P[2] / P.sum())
    # same LSQ as the data (a0 + m=1,2,3 at the 71 data longitudes)
    A71 = at_data(AF); L = np.radians(ell_m)
    X = np.stack([np.ones_like(L)] + [f(m * L) for m in (1, 2, 3) for f in (np.cos, np.sin)], 1)
    b, *_ = np.linalg.lstsq(X, A71, rcond=None)
    for m in (1, 2, 3):
        c, s = b[2 * m - 1], b[2 * m]
        out[f"lsq71_a{m}"] = float(np.hypot(c, s)); out[f"lsq71_l{m}_deg"] = float((np.degrees(np.arctan2(s, c)) / m) % (360.0 / m))
    return out

def describe(fam, eps, v, strict=False):
    AF, kap, kb = evaluate(fam, eps, v, strict)
    kap_, l0, b0, dens = unpack(fam, v)
    m = at_data(AF)
    ms = kap / kb
    fu = fu_of(fam, dens)
    rho0 = rho_shape(fam, dens, np.array([0.0]))[0] * ms / RHOC
    rhoD = OM * ms
    M = menc_table(fam, dens)
    rho_mean = ms * np.interp(D, MF, M) / (4 / 3 * np.pi) / RHOC
    Dr = 1.0 / eps
    out = dict(family=fam, D_over_r=Dr, d_over_D=eps, R2_direct=r2_direct(m), R2_corr2=r2_corr(m),
               kappa=float(kap), l0_deg=float(l0 % 360), b0_deg=float(b0),
               density={k: float(x) for k, x in zip(FAM[fam]["names"], dens)},
               kappa_base=float(kb), mass_scale=float(ms),
               rho_D_over_rhocrit=float(rhoD), rho_center_over_rhocrit=float(rho0),
               rho_center_over_rho_D=float(rho0 / rhoD), rho_mean_encl_over_rhocrit=float(rho_mean),
               f_u=float(fu), D_Gpc=float(Dr * R_PHYS_GPC), D_Gly=float(Dr * R_PHYS_GPC * GLY_PER_GPC),
               harmonics=harmonics(AF), x=[float(t) for t in v])
    return out

def fit_one(args):
    fam, Dr, seed, strict = args
    eps = 1.0 / Dr
    b = bounds(fam)
    rng = np.random.default_rng(seed)
    sols = []
    t0 = time.time()
    for trial in range(2):
        de = differential_evolution(objective, b, args=(fam, eps, strict), seed=int(rng.integers(1 << 30)),
                                    popsize=18, maxiter=250, tol=1e-7, polish=False, init="sobol",
                                    mutation=(0.5, 1.0), recombination=0.8)
        sols.append(de.x)
        # also keep a few diverse population members
        pop = de.population if hasattr(de, "population") else None
    # random multi-starts with Nelder-Mead
    for k in range(10):
        x0 = np.array([rng.uniform(lo, hi) for lo, hi in b])
        sols.append(x0)
    polished = []
    for x0 in sols:
        r = minimize(objective, x0, args=(fam, eps, strict), method="Nelder-Mead",
                     options=dict(xatol=1e-6, fatol=1e-9, maxiter=6000, maxfev=12000, adaptive=True))
        r = minimize(objective, r.x, args=(fam, eps, strict), method="Nelder-Mead",
                     options=dict(xatol=1e-7, fatol=1e-10, maxiter=6000, maxfev=12000, adaptive=True))
        polished.append((float(-r.fun), r.x))
    polished.sort(key=lambda t: -t[0])
    best = describe(fam, eps, polished[0][1], strict)
    near = [describe(fam, eps, x, strict) for (R, x) in polished if R > polished[0][0] - 0.01]
    keys = ["kappa", "l0_deg", "b0_deg", "rho_center_over_rho_D", "f_u"] + [f"density.{k}" for k in FAM[fam]["names"]]
    def get(d, k): return d["density"][k.split(".")[1]] if k.startswith("density.") else d[k]
    spread = {k: [float(min(get(d, k) for d in near)), float(max(get(d, k) for d in near))] for k in keys}
    best["n_local_optima_within_0.01"] = len(near)
    best["param_range_within_dR2_0.01"] = spread
    best["all_starts_R2"] = [p[0] for p in polished]
    best["strict_local_density"] = bool(strict)
    best["runtime_s"] = time.time() - t0
    print(f"{fam:18s} D/r={Dr:5.2f} strict={strict} R2={best['R2_direct']:.4f} ({best['runtime_s']:.0f}s)", flush=True)
    return best

def sanity():
    """Reproduce the stored 0.776 best fit (original analytic M, original grid, then new machinery)."""
    sys.path.insert(0, "/workspace/cosmology/density_local_match")
    from run_density_local_match import params_twoscale_local_match, m_enc_at_r
    best = np.load("/workspace/cosmology/map_compare/amp_mean_param_noise/ideal_fit_M_d_D_L_migkas.npz")
    MS0, EPS0, X0, L00 = float(best["mass_scale"]), float(best["eps"]), float(best["x"]), float(best["phi0_deg"])
    p = params_twoscale_local_match(X0, 0.94)
    # (i) original grid & analytic M, original Euler loop
    rg = np.linspace(0.2, 2.0, 50); chi = rg / T ** 0.4; u = np.zeros(50)
    for ti in _t[1:]:
        Hi = 2.0 / (5.0 * ti); ri = chi * ti ** 0.4
        u = u + _dt * (-Hi * u - G * m_enc_at_r(ri, p) / np.maximum(ri, 1e-30) ** 2)
    AF = model_fine(MS0, EPS0, L00, 0.0, (rg, u)); R_orig = r2_direct(at_data(AF))
    # (ii) same grid, analytic M, linear-sum integration
    rg2, u2, _ = u_table(None, None, rg, menc=lambda r: m_enc_at_r(r, p))
    R_lin = r2_direct(at_data(model_fine(MS0, EPS0, L00, 0.0, (rg2, u2))))
    # (iii) new machinery: numeric M, extended grid, kappa-parameterisation
    kap = MS0 * p["kappa"]
    v = [np.log(kap), L00, 0.0, np.log(X0), 0.94]
    AF3, _, kb = evaluate("a_exp_floor", EPS0, v)
    R_new = r2_direct(at_data(AF3))
    return dict(stored_params=dict(mass_scale=MS0, d_over_D=EPS0, D_over_L=X0, l0=L00, kappa=kap),
                R2_direct_original_code=R_orig, R2_direct_linear_sum=R_lin, R2_direct_new_machinery=R_new,
                R2_corr2_new_machinery=r2_corr(at_data(AF3)), kappa_base_analytic=p["kappa"], kappa_base_numeric=kb,
                heatmap_json_direct_best=0.7755772297783201, heatmap_json_best_fit_ref_corr2=0.7778974654481913)

DRS = [1.2, 1.32, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]

if __name__ == "__main__":
    s = sanity(); print(json.dumps(s, indent=1), flush=True)
    if len(sys.argv) > 1 and sys.argv[1] == "sanity": sys.exit()
    jobs = [(f, Dr, 1000 + i * 17 + j, False) for i, f in enumerate(FAM) for j, Dr in enumerate(DRS)]
    jobs += [(f, Dr, 5000 + i * 17 + j, True) for i, f in enumerate(FAM) for j, Dr in enumerate(DRS)]
    jobs.sort(key=lambda j: -len(FAM[j[0]]["names"]))
    with Pool(8) as pool:
        res = pool.map(fit_one, jobs, chunksize=1)
    json.dump(dict(doc=__doc__, sanity=s, results=res), open(f"{OUT}/R2_vs_Dr_density_families_raw.json", "w"), indent=1)
    print("done")
