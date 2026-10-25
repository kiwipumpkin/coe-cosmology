#!/usr/bin/env python3
"""Continuous D/r fit with the exact model: profile chi^2(D/r) (sector-correlated and naive errors) and R^2(D/r),
minimised jointly over b0, K, N (and D/L for the extended density), transforms included. 1 < D/r <= 6.
chi^2: amplitude objective of the term-matching sheets (m0 exact through N; m1..m3 weighted by 1/sigma^2).
R^2: transformed curve (m1, m2 at the data phases; m3 and m>=4 at natural phases) fitted to the 71 points;
     for fixed (b0, D/L, eps) the curve is linear in (N, P=NK), solved in closed form.
Writes figures/exact_profile_results.json."""
import sys, json, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize_scalar, minimize, brentq
sys.path.insert(0, "/workspace/cosmology/figures")
from exact_termmatch import tau_pm, tau_ext, amp, ph, l0, SIG, package, clean, E, R_GPC, T, d2r, X_SCAN

DR = np.round(np.r_[np.arange(1.01, 1.5001, 0.01), np.arange(1.55, 3.001, 0.05), np.arange(3.1, 6.0001, 0.1)], 4)
BPM = np.arange(0, 89.76, 0.25)
BEX = np.arange(0, 89.1, 1.5)
XEX = np.exp(np.linspace(np.log(0.3), np.log(200), 25))
Ld = T.l_d * d2r; yd = T.A_d; SST = np.sum((yd - yd.mean()) ** 2)
MM = np.arange(len(E.W))

def chi2(t, S):
    w = np.array([1 / S[m] ** 2 for m in (1, 2, 3)]); c = np.abs(t[1:4]); a = np.array([amp[1], amp[2], amp[3]])
    P = np.sum(w * a * c) / np.sum(w * c * c); return float(np.sum(w * (a - P * c) ** 2)), P

def gcurve(t, l=Ld):
    g = t[0] * np.ones_like(l)
    for m in (1, 2, 3):
        pk = ph[m] if m < 3 else (l0 + (0 if t[3] > 0 else 60)) % 120
        g = g + abs(t[m]) * np.cos(m * (l - pk * d2r))
    for m in range(4, len(t)):
        g = g + t[m] * np.cos(m * (l - l0 * d2r))
    return g

def r2fit(t):
    g = gcurve(t); X = np.column_stack([np.ones_like(g), g]); coef, *_ = np.linalg.lstsq(X, yd, rcond=None)
    return float(1 - np.sum((yd - X @ coef) ** 2) / SST), coef   # coef = (N, P)

def objs(t):
    c, _ = chi2(t, SIG["corr"]); n, _ = chi2(t, SIG["naive"]); r, _ = r2fit(t); return np.array([c, n, -r])

LOX, HIX = np.log(XEX[0]), np.log(XEX[-1])
def tfun(kind, eps, v):
    b0 = min(max(v[0], 0.0), 89.9)
    if kind == "pm": return tau_pm(b0, eps)
    return tau_ext(b0, eps, float(np.exp(min(max(v[1], LOX), HIX))))

def profile_one(args):
    kind, Dr = args; eps = 1 / Dr
    if kind == "pm":
        vals = np.array([objs(tau_pm(b, eps)) for b in BPM])
        out = {}
        for k, name in enumerate(("chi2_corr", "chi2_naive", "negR2")):
            i = int(np.argmin(vals[:, k])); lo, hi = BPM[max(i - 1, 0)], BPM[min(i + 1, len(BPM) - 1)]
            r = minimize_scalar(lambda b: objs(tau_pm(b, eps))[k], bounds=(lo, hi), method="bounded", options=dict(xatol=1e-6))
            out[name] = (float(r.fun), [float(r.x)])
        return kind, Dr, out
    vals = np.full((len(BEX), len(XEX), 3), np.nan)
    for j, x in enumerate(XEX):
        for i, b in enumerate(BEX):
            vals[i, j] = objs(tau_ext(b, eps, x))
    out = {}
    for k, name in enumerate(("chi2_corr", "chi2_naive", "negR2")):
        i, j = np.unravel_index(np.nanargmin(vals[:, :, k]), vals.shape[:2])
        f = lambda v: objs(tfun("ext", eps, v))[k] + 1e3 * (max(0, -v[0]) + max(0, v[1] - HIX) + max(0, LOX - v[1]))
        r = minimize(f, [BEX[i], np.log(XEX[j])], method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-10, initial_simplex=[[BEX[i], np.log(XEX[j])], [BEX[i] + 1.5, np.log(XEX[j])], [BEX[i], np.log(XEX[j]) + 0.3]]))
        v = [min(max(r.x[0], 0.0), 89.9), float(min(max(r.x[1], LOX), HIX))]
        best = min((float(vals[i, j, k]), [float(BEX[i]), float(np.log(XEX[j]))]), (float(objs(tfun("ext", eps, v))[k]), v))
        out[name] = best
    return kind, Dr, out

if __name__ == "__main__":
    t0 = time.time()
    jobs = [("ext", Dr) for Dr in DR] + [("pm", Dr) for Dr in DR]
    with Pool(6) as pool:
        res = pool.map(profile_one, jobs, chunksize=2)
    prof = {"pm": {}, "ext": {}}
    for kind, Dr, out in res: prof[kind][float(Dr)] = out
    print("grid done", time.time() - t0, flush=True)
    R = dict(doc=__doc__, DR=list(map(float, DR)), profile={k: {str(d): v for d, v in prof[k].items()} for k in prof})
    # best fits and intervals
    summ = {}
    for kind in ("pm", "ext"):
        for name in ("chi2_corr", "chi2_naive", "negR2"):
            y = np.array([prof[kind][float(d)][name][0] for d in DR])
            if kind == "ext":   # point mass is the D/L->inf member of the family
                ypm = np.array([prof["pm"][float(d)][name][0] for d in DR]); y = np.minimum(y, ypm)
            i = int(np.argmin(y)); ymin = float(y[i])
            # refine best D/r (pm only: 1-D inner problem)
            best_Dr = float(DR[i])
            if kind == "pm" and 0 < i < len(DR) - 1:
                k = ("chi2_corr", "chi2_naive", "negR2").index(name)
                def inner(Dr):
                    eps = 1 / Dr; b = prof["pm"][float(DR[i])][name][1][0]
                    r = minimize_scalar(lambda bb: objs(tau_pm(bb, eps))[k], bounds=(max(b - 3, 0), min(b + 3, 89.9)), method="bounded", options=dict(xatol=1e-7))
                    return r.fun
                r = minimize_scalar(inner, bounds=(DR[i - 1], DR[i + 1]), method="bounded", options=dict(xatol=1e-6))
                best_Dr, ymin = float(r.x), float(min(r.fun, ymin))
            iv = {}
            for dth in (1.0, 4.0):
                if name == "negR2": continue
                inside = y <= ymin + dth
                # contiguous region around best
                lo_i = i
                while lo_i > 0 and inside[lo_i - 1]: lo_i -= 1
                hi_i = i
                while hi_i < len(DR) - 1 and inside[hi_i + 1]: hi_i += 1
                lo = "boundary(1)" if lo_i == 0 else float(np.interp(ymin + dth, [y[lo_i], y[lo_i - 1]], [DR[lo_i], DR[lo_i - 1]]))
                hi = "beyond 6" if hi_i == len(DR) - 1 else float(np.interp(ymin + dth, [y[hi_i], y[hi_i + 1]], [DR[hi_i], DR[hi_i + 1]]))
                iv[str(dth)] = [lo, hi, bool(np.all(inside))]
            summ[f"{kind}_{name}"] = dict(best_Dr=best_Dr, min=ymin, intervals=iv, grid_min_Dr=float(DR[i]),
                                          zero_region=[float(DR[0]), float(DR[np.max(np.where(y < 1e-6)[0])])] if np.any(y < 1e-6) else None)
            print(kind, name, summ[f"{kind}_{name}"], flush=True)
    R["summary"] = summ
    # best-fit solutions packaged (pm: chi2_corr best, chi2_naive best, R2 best)
    best = {}
    for name in ("chi2_corr", "chi2_naive", "negR2"):
        s = summ[f"pm_{name}"]; Dr = s["best_Dr"]; eps = 1 / Dr; k = ("chi2_corr", "chi2_naive", "negR2").index(name)
        b0g = prof["pm"][min(prof["pm"].keys(), key=lambda d: abs(d - Dr))][name][1][0]
        r = minimize_scalar(lambda bb: objs(tau_pm(bb, eps))[k], bounds=(max(b0g - 3, 0), min(b0g + 3, 89.9)), method="bounded", options=dict(xatol=1e-8))
        t = tau_pm(r.x, eps)
        if name == "negR2":
            R2v, (N, P) = r2fit(t)
        else:
            _, P = chi2(t, SIG["corr" if name == "chi2_corr" else "naive"]); N = amp[0] - P * t[0]
        s_ = package(f"best pm ({name})", t, N, P / N, eps, float(r.x), dict(alpha=E.shape_point()[2], branch="on-axis" if t[1] > 0 else "flipped"))
        s_["R2_fitNP"] = r2fit(t)[0]
        best[name] = clean(s_)
    # extended R2 best
    ye = np.array([prof["ext"][float(d)]["negR2"][0] for d in DR]); i = int(np.argmin(ye)); Dr = float(DR[i]); v = prof["ext"][Dr]["negR2"][1]
    t = tfun("ext", 1 / Dr, v); R2v, (N, P) = r2fit(t); x = float(np.exp(v[1])); f, al, kb = E.shape_exp(x)
    best["ext_negR2"] = clean(package("best ext (R2)", t, N, P / N, 1 / Dr, v[0], dict(D_over_L=x, L_Gpc=R_GPC * Dr / x, alpha=al, MoverMcrit_ST=2 * P / N / al)))
    best["ext_negR2"]["R2_fitNP"] = R2v
    R["best"] = best
    json.dump(clean(R), open("/workspace/cosmology/figures/exact_profile_results.json", "w"))
    print("saved", time.time() - t0)
