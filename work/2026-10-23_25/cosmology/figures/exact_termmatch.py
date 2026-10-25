#!/usr/bin/env python3
"""Exact (no eps expansion) term matching of the CoE model to Migkas+2020 Fig.6 A(l).
Point mass (exact 1/R^2 field) and extended pure-exponential density (numeric code's profile + ST drag integration).
Same sector averaging, data harmonics, sigmas and transforms as fixedDr_octmiss.py / termmatch_Dr234.py.
Writes figures/exact_termmatch_results.json (plots made by exact_termmatch_plot.py)."""
import sys, json, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import brentq, minimize_scalar, least_squares
sys.path.insert(0, "/workspace/cosmology/figures")
import exact_tidal_core as E
import termmatch_Dr234 as T

OUT = "/workspace/cosmology/figures/exact_termmatch"
d2r = np.pi / 180; wrap = T.wrap
amp = {m: float(T.amp[m]) for m in range(4)}; ph = T.ph; l0 = T.l0
OM = json.load(open("/workspace/cosmology/figures/fixedDr_octmiss_results.json"))
SIG = {k: {int(m): v for m, v in OM["sigmas"][k].items()} for k in ("corr", "naive", "neff16")}
R12, R31, R32 = amp[1] / amp[2], amp[3] / amp[1], amp[3] / amp[2]
R_GPC = 14.3
_, _, ALPHA_PM = E.shape_point()

def tau_pm(b0, eps, **kw): return E.taus(b0, eps, E.T_point, **kw)
def tau_ser(b0, eps): return E.taus(b0, eps, E.T_series)
def tau_ext(b0, eps, x, **kw): return E.taus(b0, eps, E.T_shape, extra=E.shape_exp(x)[0], **kw)

# ---------------- model curve / packaging ----------------
def package(label, tau, N, K, eps, b0, extra=None):
    P = N * K
    Ms = np.r_[N * (1 + K * tau[0]), P * tau[1:]]
    nat = {m: (l0 + (0 if Ms[m] > 0 else 180 / m)) % (360 / m) for m in (1, 2, 3)}
    phi = {m: float(wrap(ph[m] - nat[m], 360 / m)) for m in (1, 2, 3)}
    def curve(l, p3=False):
        L = l * d2r
        y = Ms[0] + sum(abs(Ms[m]) * np.cos(m * (L - (nat[m] + (phi[m] if (m < 3 or p3) else 0)) * d2r)) for m in (1, 2, 3))
        y = y + sum(Ms[m] * np.cos(m * (L - l0 * d2r)) for m in range(4, len(Ms)))
        return y
    d = dict(label=label, D_over_r=1 / eps, eps=eps, b0=float(b0), N=float(N), K=float(K), NK=float(P),
             tau={m: float(tau[m]) for m in range(4)}, tau_hi_max=float(np.abs(tau[4:]).max()),
             model_amp={m: float(abs(Ms[m])) for m in range(4)}, model_signed={m: float(Ms[m]) for m in range(4)},
             pull={k: {m: float((abs(Ms[m]) - amp[m]) / S[m]) for m in range(4)} for k, S in SIG.items()},
             chi2={k: float(sum(((abs(Ms[m]) - amp[m]) / S[m]) ** 2 for m in (1, 2, 3))) for k, S in SIG.items()},
             nat_peak=nat, phi1=phi[1], phi2=phi[2], delta3=phi[3],
             R2=T.r2(curve(T.l_d)), R2_oct_free=T.r2(curve(T.l_d, True)),
             D_Gpc=R_GPC / eps, D_Gly=R_GPC / eps * 3.26156, MoverMcrit_alpha1=2 * K,
             phiC=float(wrap(90 - (nat[1] + phi[1]), 360)), phiC_model=float(wrap(90 - nat[1], 360)),
             phi2C=float(wrap(90 - 2 * (nat[2] + phi[2]), 360)))
    if extra: d.update(extra)
    d["_curve"] = curve
    return d

def solve_amp(tau):
    """Given tau, exact m0,m2: P = A2/|tau2|, N = A0 - P tau0."""
    P = amp[2] / abs(tau[2]); N = amp[0] - P * tau[0]; return N, P / N

# ---------------- LSQ at fixed eps (m0 exact via N; weights 1/sigma^2 on m1..m3) ----------------
def lsq_fixed(taufun, eps, S, bgrid=np.arange(0, 89.51, 0.25), label="LSQ", extra=None):
    w = {m: 1 / S[m] ** 2 for m in (1, 2, 3)}
    def solve(b0):
        t = taufun(b0, eps); c = {m: abs(t[m]) for m in (1, 2, 3)}
        P = sum(w[m] * amp[m] * c[m] for m in c) / sum(w[m] * c[m] ** 2 for m in c)
        return P, sum(w[m] * (amp[m] - P * c[m]) ** 2 for m in c), t
    ch = np.array([solve(b)[1] for b in bgrid]); mins = []
    for i in range(len(bgrid)):
        if (i == 0 or ch[i] <= ch[i - 1]) and (i == len(bgrid) - 1 or ch[i] <= ch[i + 1]):
            lo_, hi_ = bgrid[max(i - 1, 0)], bgrid[min(i + 1, len(bgrid) - 1)]
            b = minimize_scalar(lambda b: solve(b)[1], bounds=(lo_, hi_), method="bounded", options=dict(xatol=1e-6)).x
            mins.append((solve(b)[1], b))
    mins.sort(); out = []
    for c2, b0 in mins[:3]:
        P, c2, t = solve(b0); N = amp[0] - P * t[0]
        ex = dict(chi2_fit=float(c2), branch="on-axis" if t[1] > 0 else "flipped"); ex.update(extra or {})
        out.append(package(label, t, N, P / N, eps, b0, ex))
    return out

def exact012(taufun, eps, bgrid=np.arange(0, 89.51, 0.25), extra=None):
    """Match m0,m1,m2 exactly (|tau1|/|tau2| = A1/A2 -> b0), m3 is a prediction."""
    g = lambda b: np.log(abs(taufun(b, eps)[1]) / abs(taufun(b, eps)[2])) - np.log(R12)
    gv = np.array([g(b) for b in bgrid]); out = []
    for i in range(len(bgrid) - 1):
        if np.isfinite(gv[i]) and np.isfinite(gv[i + 1]) and gv[i] * gv[i + 1] < 0:
            b0 = brentq(g, bgrid[i], bgrid[i + 1], xtol=1e-9); t = taufun(b0, eps); N, K = solve_amp(t)
            ex = dict(branch="on-axis" if t[1] > 0 else "flipped"); ex.update(extra or {})
            out.append(package("exact m0-m2", t, N, K, eps, b0, ex))
    return out

# ---------------- free-eps exact solve (point mass): scan (b0, eps) ----------------
def scan_row_pm(eps):
    return [(tau_pm(b, eps)) for b in B_SCAN]
B_SCAN = np.arange(0, 89.01, 0.5)
EPS_SCAN = np.r_[np.linspace(0.12, 0.98, 87), np.linspace(1.02, 1.6, 30)]

def ratio_resid(t):
    return np.array([np.log(abs(t[3]) / abs(t[1])) - np.log(R31), np.log(abs(t[1]) / abs(t[2])) - np.log(R12)])

def find_roots_2d(F1, F2, xs, ys, fun):
    """F1,F2 on grid [ix,iy]; candidate cells where both change sign; refine with least_squares on fun(x,y)."""
    cands = []
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            c1 = F1[i:i + 2, j:j + 2]; c2 = F2[i:i + 2, j:j + 2]
            if not (np.all(np.isfinite(c1)) and np.all(np.isfinite(c2))): continue
            if c1.min() < 0 < c1.max() and c2.min() < 0 < c2.max():
                cands.append((0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[j] + ys[j + 1]), (xs[i], xs[i + 1]), (ys[j], ys[j + 1])))
    roots = []
    for x0, y0, bx, by in cands:
        try:
            r = least_squares(lambda v: fun(v[0], v[1]), [x0, y0], bounds=([bx[0] - 2 * (bx[1] - bx[0]), by[0] - 2 * (by[1] - by[0])],
                                                                         [bx[1] + 2 * (bx[1] - bx[0]), by[1] + 2 * (by[1] - by[0])]), xtol=1e-12, ftol=1e-14)
        except Exception:
            continue
        if np.max(np.abs(r.fun)) < 1e-7 and not any(abs(r.x[0] - q[0]) < 1e-3 and abs(r.x[1] - q[1]) < 1e-4 for q in roots):
            roots.append((float(r.x[0]), float(r.x[1])))
    return roots

# ---------------- extended ----------------
X_SCAN = np.exp(np.linspace(np.log(0.3), np.log(200), 41))
def ext_grid(args):
    eps, = args
    F1 = np.full((len(B_SCAN), len(X_SCAN)), np.nan); F2 = F1.copy()
    for j, x in enumerate(X_SCAN):
        for i, b in enumerate(B_SCAN[::2]):
            pass
    # coarse b0 grid (1 deg) for speed
    bg = np.arange(0, 89.01, 1.0)
    F1 = np.full((len(bg), len(X_SCAN)), np.nan); F2 = F1.copy()
    for j, x in enumerate(X_SCAN):
        for i, b in enumerate(bg):
            t = tau_ext(b, eps, x); F1[i, j], F2[i, j] = ratio_resid(t)
    fun = lambda b, lx: ratio_resid(tau_ext(b, eps, np.exp(lx)))
    roots = find_roots_2d(F1, F2, bg, np.log(X_SCAN), fun)
    return eps, roots, F1, F2, bg

def ext_lsq(eps, S):
    """No exact root: minimise chi2 over (b0, x) with P closed form (m0 exact via N)."""
    w = {m: 1 / S[m] ** 2 for m in (1, 2, 3)}
    LO, HI = np.log(X_SCAN[0]), np.log(X_SCAN[-1])
    def chi(v):
        pen = 0.0
        if v[1] < LO: pen += 1e3 * (LO - v[1]); v = [v[0], LO]
        if v[1] > HI: pen += 1e3 * (v[1] - HI); v = [v[0], HI]
        if v[0] < 0: pen += 1e3 * (-v[0]); v = [0.0, v[1]]
        t = tau_ext(v[0], eps, np.exp(v[1])); c = {m: abs(t[m]) for m in (1, 2, 3)}
        P = sum(w[m] * amp[m] * c[m] for m in c) / sum(w[m] * c[m] ** 2 for m in c)
        return sum(w[m] * (amp[m] - P * c[m]) ** 2 for m in c) + pen
    best = None
    for b in np.arange(0, 89, 3.0):
        for lx in np.log(X_SCAN[::4]):
            v = chi([b, lx])
            if best is None or v < best[0]: best = (v, b, lx)
    from scipy.optimize import minimize
    r = minimize(chi, [best[1], best[2]], method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-10))
    x = np.array([max(r.x[0], 0.0), np.clip(r.x[1], LO, HI)])
    return x, r.fun

def cont_row(Dr):
    eps = 1 / Dr
    out = dict(Dr=Dr)
    for b0 in (0.0, 19.1, 40.0):
        te, ts = tau_pm(b0, eps), tau_ser(b0, eps)
        out[f"ratio_b{b0}"] = [float(te[m] / ts[m]) for m in range(4)]
        out[f"tau_exact_b{b0}"] = [float(v) for v in te[:4]]; out[f"tau_series_b{b0}"] = [float(v) for v in ts[:4]]
    L = lsq_fixed(tau_pm, eps, SIG["naive"], bgrid=np.arange(0, 89.51, 0.5))[0]
    Lc = lsq_fixed(tau_pm, eps, SIG["corr"], bgrid=np.arange(0, 89.51, 0.5))[0]
    ex = exact012(tau_pm, eps, bgrid=np.arange(0, 89.51, 0.5))
    out["m3_lsq_naive"] = L["model_amp"][3]; out["b0_lsq_naive"] = L["b0"]; out["R2_lsq_naive"] = L["R2"]
    out["m3_lsq_corr"] = Lc["model_amp"][3]; out["chi2_corr_profile"] = Lc["chi2"]["corr"]
    out["exact012"] = [(s["b0"], s["model_amp"][3], s["branch"]) for s in ex]
    return out

def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items() if k != "_curve"}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, np.ndarray): return clean(o.tolist())
    if isinstance(o, (np.floating, np.integer)): return float(o)
    return o

if __name__ == "__main__":
    t0 = time.time()
    res = dict(doc=__doc__, alpha_point_mass_numeric=ALPHA_PM, data=dict(amp=amp, peak=ph, l0=l0, R12=R12, R31=R31, R32=R32), sigmas=SIG)
    # 0. limit check
    lim = []
    for b0 in (0.0, 19.1, 40.0):
        for eps in (0.3, 0.1, 0.03, 0.01):
            te, ts = tau_pm(b0, eps), tau_ser(b0, eps)
            lim.append(dict(b0=b0, eps=eps, ratio=[float(te[m] / ts[m]) for m in range(4)]))
    mu = np.linspace(-1, 1, 2001)
    res["limit"] = dict(tau_ratios=lim, pointwise_maxdiff_over_eps2={str(e): float(np.max(np.abs(E.T_point(mu, e) - E.T_series(mu, e))) / e ** 2) for e in (0.1, 0.01, 0.001)})
    # check sector average: FFT route vs brute-force 72 sectors
    prof = E.sector_profile_direct(30, 0.6, step=0.05); ps = np.arange(0, 360, 5.) * d2r; tt = E.taus(30, 0.6, npsi=2880)
    res["sector_check_maxdiff"] = float(np.max(np.abs(prof - sum(tt[m] * np.cos(m * ps) for m in range(len(tt))))))
    print("limit done", time.time() - t0, flush=True)

    with Pool(6) as pool:
        # 1. free eps point-mass scan
        rows = pool.map(scan_row_pm, list(EPS_SCAN))
        F1 = np.full((len(B_SCAN), len(EPS_SCAN)), np.nan); F2 = F1.copy()
        for j, row in enumerate(rows):
            for i, t in enumerate(row): F1[i, j], F2[i, j] = ratio_resid(t)
        roots = find_roots_2d(F1, F2, B_SCAN, EPS_SCAN, lambda b, e: ratio_resid(tau_pm(b, e)))
        free_pm = []
        for b0, eps in sorted(roots, key=lambda r: r[1]):
            t = tau_pm(b0, eps); N, K = solve_amp(t)
            # convergence check with finer grid
            tf = tau_pm(b0, eps, nb=800, npsi=2880)
            free_pm.append(package("free eps exact point mass", t, N, K, eps, b0, dict(conv=float(np.max(np.abs(tf[:4] - t[:4]))))))
        res["free_pm"] = free_pm; res["free_scan"] = dict(b0=B_SCAN, eps=EPS_SCAN, F1=F1, F2=F2)
        print("free pm roots", [(s["b0"], s["D_over_r"], s["K"]) for s in free_pm], time.time() - t0, flush=True)
        # 2. fixed D/r point mass
        DRS = [1.3, 2.0, 3.0, 4.0]
        res["fixed_pm"] = {}
        for Dr in DRS:
            eps = 1 / Dr
            res["fixed_pm"][Dr] = dict(lsq_naive=lsq_fixed(tau_pm, eps, SIG["naive"]), lsq_corr=lsq_fixed(tau_pm, eps, SIG["corr"]),
                                       exact012=exact012(tau_pm, eps),
                                       series_lsq_naive=lsq_fixed(lambda b, e: tau_ser(b, e), eps, SIG["naive"], label="series LSQ"))
        print("fixed pm done", time.time() - t0, flush=True)
        # 3. extended: exact 4x4 at fixed D/r (unknowns N,K,b0,D/L)
        EXT_DRS = [1.3, 2.0, 3.0, 4.0, 1.1, 1.5, 1.75, 2.5, 3.5, 5.0, 6.0]
        ext = pool.map(ext_grid, [(1 / Dr,) for Dr in EXT_DRS])
        res["fixed_ext"] = {}
        for (eps, roots, F1e, F2e, bg), Dr in zip(ext, EXT_DRS):
            sols = []
            for b0, lx in roots:
                x = float(np.exp(lx)); t = tau_ext(b0, eps, x); N, K = solve_amp(t); f, al, kb = E.shape_exp(x)
                sols.append(package("extended exact 4x4", t, N, K, eps, b0, dict(D_over_L=x, L_Gpc=R_GPC / eps / x, alpha=al,
                                    kappa=K / al, MoverMcrit_ST=2 * K / al, branch="on-axis" if t[1] > 0 else "flipped")))
            lsq = None
            if not sols:
                v, c2 = ext_lsq(eps, SIG["naive"]); x = float(np.exp(v[1])); t = tau_ext(v[0], eps, x)
                w = {m: 1 / SIG["naive"][m] ** 2 for m in (1, 2, 3)}; c = {m: abs(t[m]) for m in (1, 2, 3)}
                P = sum(w[m] * amp[m] * c[m] for m in c) / sum(w[m] * c[m] ** 2 for m in c); N = amp[0] - P * t[0]; f, al, kb = E.shape_exp(x)
                lsq = package("extended LSQ", t, N, P / N, eps, v[0], dict(D_over_L=x, L_Gpc=R_GPC / eps / x, alpha=al, kappa=P / N / al, MoverMcrit_ST=2 * P / N / al))
            res["fixed_ext"][Dr] = dict(exact=sols, lsq=lsq, grid=dict(b0=bg, x=X_SCAN, F1=F1e, F2=F2e))
            print("ext", Dr, [(round(s["b0"], 2), round(s["D_over_L"], 3), round(s["R2"], 4)) for s in sols], "lsq" if lsq else "", time.time() - t0, flush=True)
        # 4. continuous comparisons (point mass exact vs series)
        DRC = list(np.round(np.r_[np.linspace(1.05, 2.0, 20), np.linspace(2.1, 6.0, 40)], 4))
        res["cont"] = pool.map(cont_row, DRC)
        print("cont done", time.time() - t0, flush=True)
    json.dump(clean(res), open(OUT + "_results.json", "w"))
    print("saved", time.time() - t0)
