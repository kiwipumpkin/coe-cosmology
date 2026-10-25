#!/usr/bin/env python3
"""Term (harmonic-amplitude) matching of the Sedov-Taylor CoE model to Migkas+2020 Fig.6 A(l) at D/r = 2, 3, 4.
Sources on the vision shell: d = r, eps = r/D = 1/(D/r).
A. ANALYTIC point-mass tidal form (latitude-averaged, 90-deg sectors, |b|>20 cut, cos b weight; as solve_v2.py),
   FIXED eps. Unknowns N, K, b0 (l0 = data m=2 peak 14.96 deg). Amplitude eqs:
     m0: A0 = N(1+K q0) ; m1: |A1| = NK eps|o1| ; m2: |A2| = NK|q2| ; m3: |A3| = NK eps|o3|
   Weighted LSQ (weights 1/sigma_m^2, sigma_m = amplitude errors from the 71-pt LSQ Fourier fit with per-point
   errors (lo+hi)/2, correlations between overlapping sectors ignored). With P=NK: m0 fixes N exactly for any
   (P,b0); (P,b0) from LSQ on m1..m3 (P linear -> closed form at each b0, 1-D scan + refine in b0).
   Transforms: phi1 = data dipole peak - model natural dipole peak; phi2 = l2 - l0 = 0 by construction.
B. NUMERIC ideal-drag model (r2_vs_Dr_density_families.py, family a_exp_floor, ring geometry, kappa free),
   curve on 181-pt grid, FFT m=0..3; m=1 and m=2 components phase-shifted by phi1, phi2 (from the curve's own
   natural peaks), m>=3 untouched. Params (log kappa, l0, b0, log D/L, f_u, phi1, phi2) fitted to maximise R^2
   of the transformed curve vs the 71 points. |b0|<=30 and D/L<=300 as in the stored fits.
   PRIMARY ('fixl0'): l0 fixed at the data m=2 peak 14.96 deg, same as the analytic solve, so the octupole
   sits at its natural phase relative to that axis; phi1, phi2 free. REFERENCE ('freel0'): l0 also free; since
   phi1, phi2 already set the m=1, m=2 phases, a free l0 only moves the octupole, i.e. it is in effect a third
   free phase (reported for completeness, not as the term-matching answer).
   Warm starts: stored best fits (phi from data-model phase difference) + DE + Nelder-Mead polish."""
import json, sys, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize, minimize_scalar
sys.path.insert(0, "/workspace/cosmology/figures")
import r2_vs_Dr_density_families as M
M.FAM["a_exp_floor"]["hi"][0] = 300.0
M.RG_LO, M.NRAD = 0.02, 80
OUT = "/workspace/cosmology/figures/termmatch_Dr234"
DRS = [2.0, 3.0, 4.0]
R_GPC, GLY = 14.3, 3.26156
d2r = np.pi / 180
wrap = lambda x, P: (x + P / 2) % P - P / 2

dat = np.genfromtxt(M.CSV, delimiter=",", names=True)
l_d, A_d, elo, ehi = dat["l_deg"], dat["A_over_Aall"], dat["A_err_low"], dat["A_err_high"]
SST = np.sum((A_d - A_d.mean()) ** 2)
r2 = lambda p: float(1 - np.sum((A_d - p) ** 2) / SST)

# ---------------- data harmonics ----------------
coef = np.array([0.9965, 0.0296, 0.0434, 0.0636, 0.0366, -0.0299, -0.0114])
X = np.column_stack([np.ones_like(l_d)] + [f(m * l_d * d2r) for m in (1, 2, 3) for f in (np.cos, np.sin)])
fourier = lambda l: np.column_stack([np.ones_like(l)] + [f(m * l * d2r) for m in (1, 2, 3) for f in (np.cos, np.sin)]) @ coef
R2_fourier = r2(X @ coef)
amp = {0: coef[0]}; ph = {}
for m in (1, 2, 3):
    c, s = coef[2 * m - 1], coef[2 * m]; amp[m] = float(np.hypot(c, s)); ph[m] = float((np.degrees(np.arctan2(s, c)) / m) % (360 / m))
sig_i = 0.5 * (elo + ehi)
Wt = np.diag(1 / sig_i ** 2); cov = np.linalg.inv(X.T @ Wt @ X)
sig = {0: float(np.sqrt(cov[0, 0]))}
for m in (1, 2, 3):
    i, j = 2 * m - 1, 2 * m; c, s = coef[i], coef[j]; g = np.array([c, s]) / np.hypot(c, s)
    sig[m] = float(np.sqrt(g @ cov[np.ix_([i, j], [i, j])] @ g))
l0 = ph[2]

# ---------------- analytic coefficients (identical to solve_v2.py) ----------------
bc, hp = 20 * d2r, np.pi / 2
Fi = lambda f: f(hp) - f(bc); nrm = Fi(np.sin)
c1 = Fi(lambda b: b / 2 + np.sin(2 * b) / 4) / nrm
c2 = Fi(lambda b: np.sin(b) - np.sin(b) ** 3 / 3) / nrm
c3 = Fi(lambda b: 3 * b / 8 + np.sin(2 * b) / 4 + np.sin(4 * b) / 32) / nrm
s2 = Fi(lambda b: np.sin(b) ** 3 / 3) / nrm
cs2 = Fi(lambda b: b / 8 - np.sin(4 * b) / 32) / nrm
W = {m: np.sin(m * np.pi / 4) / (m * np.pi / 4) for m in (1, 2, 3)}
def coefs(b0):
    a, s0 = np.cos(b0 * d2r), np.sin(b0 * d2r)
    return dict(q0=1.5 * a * a * c2 + 3 * s0 * s0 * s2 - 1, q2=1.5 * a * a * c2 * W[2],
                o1=a * W[1] * (4.5 * c1 - 45 / 8 * a * a * c3 - 22.5 * s0 * s0 * cs2), o3=-15 / 8 * a ** 3 * c3 * W[3])

def analytic(Dr):
    eps = 1 / Dr
    w = {m: 1 / sig[m] ** 2 for m in (1, 2, 3)}
    def solveP(b0):
        t = coefs(b0); cm = {1: eps * abs(t["o1"]), 2: abs(t["q2"]), 3: eps * abs(t["o3"])}
        P = sum(w[m] * amp[m] * cm[m] for m in cm) / sum(w[m] * cm[m] ** 2 for m in cm)
        chi2 = sum(w[m] * (amp[m] - P * cm[m]) ** 2 for m in cm)
        return P, chi2, cm, t
    grid = np.arange(0, 89.5, 0.05); ch = np.array([solveP(b)[1] for b in grid])
    # local minima
    mins = []
    for i in range(len(grid)):
        if (i == 0 or ch[i] <= ch[i - 1]) and (i == len(grid) - 1 or ch[i] <= ch[i + 1]):
            lo_, hi_ = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
            if lo_ == hi_: b = lo_
            else: b = minimize_scalar(lambda b: solveP(b)[1], bounds=(lo_, hi_), method="bounded", options=dict(xatol=1e-9)).x
            mins.append((solveP(b)[1], b))
    mins.sort()
    sols = []
    for chi2, b0 in mins:
        P, chi2, cm, t = solveP(b0)
        N = amp[0] - P * t["q0"]; K = P / N
        Ms = {0: N * (1 + K * t["q0"]), 1: P * eps * t["o1"], 2: P * t["q2"], 3: P * eps * t["o3"]}
        nat = {m: (l0 + (0 if Ms[m] > 0 else 180 / m)) % (360 / m) for m in (1, 2, 3)}
        phi = {m: wrap(ph[m] - nat[m], 360 / m) for m in (1, 2, 3)}
        def curve(l, p1, p2, p3=0.0):
            L = l * d2r
            return (Ms[0] + abs(Ms[1]) * np.cos(L - (nat[1] + p1) * d2r) + abs(Ms[2]) * np.cos(2 * (L - (nat[2] + p2) * d2r))
                    + abs(Ms[3]) * np.cos(3 * (L - (nat[3] + p3) * d2r)))
        res = {m: float(abs(Ms[m]) - amp[m]) for m in (0, 1, 2, 3)}
        sols.append(dict(D_over_r=Dr, eps=eps, b0=float(b0), N=float(N), K=float(K), NK=float(P), chi2=float(chi2),
                         q0=t["q0"], q2=t["q2"], o1=t["o1"], o3=t["o3"], c_m={str(k): float(v) for k, v in cm.items()},
                         model_amp={str(m): float(abs(Ms[m])) for m in Ms}, model_signed={str(m): float(Ms[m]) for m in Ms},
                         resid={str(m): res[m] for m in res}, resid_over_sigma={str(m): res[m] / sig[m] for m in res},
                         nat_peak={str(m): float(nat[m]) for m in nat}, phi1=float(phi[1]), phi2=float(phi[2]),
                         delta3=float(phi[3]), delta3_3l=float(3 * phi[3]),
                         R2_transformed=r2(curve(l_d, phi[1], phi[2])), R2_notransform=r2(curve(l_d, 0, 0)),
                         R2_all3phases=r2(curve(l_d, phi[1], phi[2], phi[3])),
                         MoverMcrit_alpha1=float(2 * K), D_Gpc=Dr * R_GPC, D_Gly=Dr * R_GPC * GLY,
                         phiC_model=float(wrap(90 - nat[1], 360)), phiC_data=float(wrap(90 - ph[1], 360)),
                         phi2C=float(wrap(90 - 2 * (nat[2] + phi[2]), 360))))
        sols[-1]["_curve"] = curve
    return sols

# ---------------- numeric ----------------
def num_harm(AF):
    n = len(AF); F = np.fft.rfft(AF) / n; out = {0: float(F[0].real)}
    for m in (1, 2, 3):
        c, s = 2 * F[m].real, -2 * F[m].imag; out[m] = (float(np.hypot(c, s)), float((np.degrees(np.arctan2(s, c)) / m) % (360 / m)))
    return out

def transform(AF, p1, p2):
    h = num_harm(AF); y = AF.copy(); L = M.LG * d2r
    for m, p in ((1, p1), (2, p2)):
        a, lm = h[m]
        y = y - a * np.cos(m * (L - lm * d2r)) + a * np.cos(m * (L - (lm + p) * d2r))
    return y

BND = M.bounds("a_exp_floor") + [(-180.0, 180.0), (-90.0, 90.0)]
def num_obj(v, eps, use_tf=True, fixl0=None):
    if fixl0 is not None: v = np.r_[v[0], fixl0, v[1:]]
    pen = 0.0; vv = np.array(v, float)
    for i, (lo, hi) in enumerate(BND):
        if i in (1, 5, 6): continue
        if vv[i] < lo: pen += (lo - vv[i]) * 10; vv[i] = lo
        if vv[i] > hi: pen += (vv[i] - hi) * 10; vv[i] = hi
    try:
        AF, _, _ = M.evaluate("a_exp_floor", eps, vv[:5])
        y = transform(AF, vv[5], vv[6]) if use_tf else AF
        m = M.at_data(y)
        if not np.all(np.isfinite(m)): return 10.0
        return -r2(m) + pen
    except Exception:
        return 10.0

def numeric(job):
    Dr, mode = job
    t0 = time.time(); eps = 1 / Dr; fl = l0 if mode == "fixl0" else None
    idx = [0, 2, 3, 4, 5, 6] if fl is not None else list(range(7))
    bnd = [BND[i] for i in idx]
    stored = [r for r in json.load(open("/workspace/cosmology/figures/R2_vs_Dr_1to10.json"))["results"]
              if r["family"] == "a_exp_floor" and not r["strict_local_density"] and abs(r["D_over_r"] - Dr) < 1e-9][0]
    xs = np.array(stored["x"])
    if fl is not None: xs[1] = fl
    AF0, _, _ = M.evaluate("a_exp_floor", eps, xs); h0 = num_harm(AF0)
    R2_stored_check = r2(M.at_data(AF0))
    p1s = wrap(ph[1] - h0[1][1], 360); p2s = wrap(ph[2] - h0[2][1], 180)
    starts = [np.r_[xs, p1s, p2s][idx]]
    rng = np.random.default_rng(int(Dr * 100))
    de = differential_evolution(num_obj, bnd, args=(eps, True, fl), seed=int(Dr * 7), popsize=14, maxiter=150, tol=1e-8,
                                polish=False, init="sobol", x0=starts[0])
    starts.append(de.x)
    for k in range(4):
        x = np.array([rng.uniform(lo, hi) for lo, hi in bnd]); starts.append(x)
    pol = []
    for x0 in starts:
        r = minimize(num_obj, x0, args=(eps, True, fl), method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-10, maxiter=8000, maxfev=16000, adaptive=True))
        r = minimize(num_obj, r.x, args=(eps, True, fl), method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-11, maxiter=8000, maxfev=16000, adaptive=True))
        pol.append((-r.fun, r.x))
    pol.sort(key=lambda t: -t[0])
    full = lambda x: np.r_[x[0], fl, x[1:]] if fl is not None else np.array(x)
    def summ(R, x):
        x = full(x).copy(); x[1] %= 360; x[5] = wrap(x[5], 360); x[6] = wrap(x[6], 180)
        AF, kap, kb = M.evaluate("a_exp_floor", eps, x[:5]); h = num_harm(AF)
        return dict(R2_transformed=float(R), R2_notransform=r2(M.at_data(AF)), kappa=float(kap), MoverMcrit=float(2 * kap),
                    l0=float(x[1]), b0=float(x[2]), D_over_L=float(np.exp(x[3])), f_u=float(x[4]), phi1=float(x[5]), phi2=float(x[6]),
                    harm_natural={str(m): h[m] for m in h}, delta3=float(wrap(ph[3] - h[3][1], 120)))
    branches = {}
    for R, xx in pol:   # dipole-on-axis (natural m=1 peak within 90 deg of l0) vs flipped
        sm = summ(R, xx); key = "on_axis" if abs(wrap(sm["harm_natural"]["1"][1] - sm["l0"], 360)) < 90 else "flipped"
        if key not in branches: branches[key] = sm
    x = full(pol[0][1])
    x[1] = x[1] % 360; x[5] = wrap(x[5], 360); x[6] = wrap(x[6], 180)
    AF, kap, kb = M.evaluate("a_exp_floor", eps, x[:5]); h = num_harm(AF); y = transform(AF, x[5], x[6])
    ht = num_harm(y)
    kap_, l0n, b0n, dens = M.unpack("a_exp_floor", x[:5])
    out = dict(D_over_r=Dr, mode=mode, eps=eps, x=[float(t) for t in x], kappa=float(kap), MoverMcrit=float(2 * kap), l0=float(l0n), b0=float(b0n),
               D_over_L=float(dens[0]), f_u=float(dens[1]), L_Gpc=float(Dr * R_GPC / dens[0]), phi1=float(x[5]), phi2=float(x[6]),
               R2_transformed=float(pol[0][0]), R2_notransform=r2(M.at_data(AF)),
               R2_stored=stored["R2_direct"], R2_stored_recomputed=R2_stored_check,
               stored_x=stored["x"], stored_harm={str(m): h0[m] for m in h0},
               start_phi1=float(p1s), start_phi2=float(p2s),
               harm_natural={str(m): h[m] for m in h}, harm_transformed={str(m): ht[m] for m in ht},
               delta3=float(wrap(ph[3] - h[3][1], 120)), all_starts_R2=[float(p[0]) for p in pol], dipole_branches=branches,
               b0_at_bound=bool(abs(abs(b0n) - 30) < 1e-3), DL_at_bound=bool(dens[0] > 299.9),
               kappa_base=float(kb), runtime_s=time.time() - t0,
               curve_transformed=list(map(float, y)), curve_natural=list(map(float, AF)))
    print(f"numeric {mode} D/r={Dr}: R2_tf={out['R2_transformed']:.4f} R2_nat={out['R2_notransform']:.4f} ({out['runtime_s']:.0f}s)", flush=True)
    return out

if __name__ == "__main__":
    print("data amps", amp, "phases", ph, "sigma", sig, "R2 fourier", R2_fourier)
    with Pool(6) as pool:
        nres_all = pool.map(numeric, [(Dr, md) for md in ("fixl0", "freel0") for Dr in DRS])
    nres = [r for r in nres_all if r["mode"] == "fixl0"]; nfree = [r for r in nres_all if r["mode"] == "freel0"]
    ares = {Dr: analytic(Dr) for Dr in DRS}
    res = dict(doc=__doc__, data=dict(amp={str(k): v for k, v in amp.items()}, peak={str(k): v for k, v in ph.items()},
               sigma={str(k): v for k, v in sig.items()}, l0=l0, R2_fourier=R2_fourier, coef=list(coef)),
               lat=dict(c1=c1, c2=c2, c3=c3, s2=s2, cs2=cs2), W={str(k): v for k, v in W.items()},
               analytic={str(Dr): [{k: v for k, v in s.items() if k != "_curve"} for s in ares[Dr]] for Dr in DRS},
               numeric={str(r["D_over_r"]): r for r in nres}, numeric_freel0={str(r["D_over_r"]): r for r in nfree}, LG=list(map(float, M.LG)))
    json.dump(res, open(OUT + "_results.json", "w"), indent=1, default=float)
    for Dr in DRS:
        for s in ares[Dr]:
            print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in s.items() if k != "_curve"})
        for n in [r for r in nres_all if r["D_over_r"] == Dr]:
            print({k: v for k, v in n.items() if not k.startswith("curve")})
