import sys, json, numpy as np, time
sys.path.insert(0, "."); from param_scan_core import *
from scipy.optimize import minimize
rng = np.random.default_rng(1)
M0 = CUR["M"]
BOUNDS = dict(DR=(1.05, 10.0), LD=(0.02, 3.0), M=(0.1*M0, 1000*M0), amp=(-0.15, 0.15), phi=(0, 2*np.pi), l0=(0, 360))
LOG = {"DR": False, "LD": True, "M": True, "amp": False, "phi": False, "l0": False}
PER = {"phi": 2*np.pi, "l0": 360.0}
def to_u(n, v):
    lo, hi = BOUNDS[n]
    if LOG[n]: lo, hi, v = np.log(lo), np.log(hi), np.log(v)
    return (v - lo) / (hi - lo)
def from_u(n, u):
    lo, hi = BOUNDS[n]
    if n in PER: u = u % 1.0
    else: u = min(max(u, 0.0), 1.0)
    if LOG[n]: return float(np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo))))
    return float(lo + u * (hi - lo))
def fit(target, free, fixed, starts=12, seeds=()):
    f = TARGETS[target]
    def obj(u):
        p = dict(fixed); pen = 0
        for n, ui in zip(free, u):
            if n not in PER: pen += 10 * (max(0, -ui) + max(0, ui - 1)) ** 2
            p[n] = from_u(n, ui)
        v = f(p); return (1 - v if np.isfinite(v) else 1e3) + pen
    if not free:
        return f(fixed), dict(fixed)
    # coarse random screen
    U = rng.random((400, len(free)))
    for s in seeds: U = np.vstack([U, [to_u(n, s[n]) for n in free]])
    vals = np.array([obj(u) for u in U]); order = np.argsort(vals)
    cand = [U[i] for i in order[:starts]] + [U[-j-1] for j in range(len(seeds))]
    best = (np.inf, None)
    for u0 in cand:
        r = minimize(obj, u0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-10, maxiter=4000, maxfev=8000))
        r = minimize(obj, r.x, method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-12, maxiter=4000))
        if r.fun < best[0]: best = (r.fun, r.x)
    p = dict(fixed); [p.__setitem__(n, from_u(n, ui)) for n, ui in zip(free, best[1])]
    return f(p), p
BASE = ["DR", "LD", "M", "amp", "phi"]
out = {"current": {t: TARGETS[t](CUR) for t in TARGETS}, "global": {}}
t0 = time.time()
for target in ["mig", "four", "sheet"]:
    for with_l0 in ([False, True] if target != "sheet" else [False]):
        free = BASE + (["l0"] if with_l0 else [])
        fixed = {} if with_l0 else {"l0": 0.0}
        v, p = fit(target, free, fixed, starts=25, seeds=[CUR])
        other = {t: TARGETS[t](p) for t in TARGETS}
        edges = [n for n in free if n not in PER and (abs(to_u(n, p[n])) < 1e-3 or abs(to_u(n, p[n]) - 1) < 1e-3)]
        key = target + ("_l0" if with_l0 else "")
        out["global"][key] = dict(R2=v, params=p, R2_all=other, at_edge=edges)
        print(key, round(v, 4), {n: (round(p[n], 4) if n != "M" else "%.3e" % p[n]) for n in p}, other, edges, round(time.time() - t0))
json.dump(out, open("param_scan_results.json", "w"), indent=1)
