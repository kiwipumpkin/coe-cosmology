#!/usr/bin/env python3
"""Numeric (ideal-drag ST + Gauss, exp core + floor, kappa free, non-strict) fits at D/r = 0.9, 1.01, 1.3,
which are not in R2_vs_Dr_1to10.json. Same model/objective as r2_vs_Dr_1to10.py (D/L upper bound 300),
but the u(r) table is extended down to r = 0.001 D (600 pts) because for D/r < ~1.05 the sightlines
toward the CoE pass within |d - D| of the centre (0.0099 D at D/r = 1.01, 0.111 D at 0.9).
Longitude grid kept at 181 pts (model definition). Optimiser: 2 DE runs (popsize 15, maxiter 200)
+ warm starts from all a_exp_floor/d_pure_exp solutions in R2_vs_Dr_1to10.json at the 2 nearest D/r
+ 6 random starts, each Nelder-Mead polished twice."""
import json, time, numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
import r2_vs_Dr_density_families as M
M.RG_LO, M.NRAD = 0.001, 600
for f in ("a_exp_floor", "d_pure_exp"): M.FAM[f]["hi"][0] = 300.0
FAM = "a_exp_floor"
old = json.load(open(f"{M.OUT}/R2_vs_Dr_1to10.json"))["results"]

def warm(Dr):
    Ds = sorted(set(r["D_over_r"] for r in old), key=lambda d: abs(np.log(d / Dr)))[:2]
    S = []
    for r in old:
        if r["D_over_r"] in Ds and r["family"] in ("a_exp_floor", "d_pure_exp"):
            x = list(r["x"]); x = x + [0.0] if len(x) == 4 else x
            S.append(np.array(x))
    return S

def fit(args):
    Dr, seed = args; eps = 1 / Dr; b = M.bounds(FAM); rng = np.random.default_rng(seed); t0 = time.time()
    S = []
    for k in range(2):
        de = differential_evolution(M.objective, b, args=(FAM, eps, False), seed=seed + k, popsize=15, maxiter=200,
                                    tol=1e-7, polish=False, init="sobol")
        S.append(de.x)
    S += warm(Dr) + [np.array([rng.uniform(lo, hi) for lo, hi in b]) for _ in range(6)]
    P = []
    for x0 in S:
        r = minimize(M.objective, x0, args=(FAM, eps, False), method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9, maxfev=4000, adaptive=True))
        r = minimize(M.objective, r.x, args=(FAM, eps, False), method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-10, maxfev=2000, adaptive=True))
        P.append((-r.fun, r.x))
    P.sort(key=lambda t: -t[0])
    d = M.describe(FAM, eps, P[0][1], False)
    d["L_Gpc"] = d["D_Gpc"] / d["density"]["D_over_L"]; d["strict_local_density"] = False
    d["all_starts_R2"] = [float(p[0]) for p in P]; d["runtime_s"] = time.time() - t0; d["seed"] = seed
    d["grid"] = dict(RG_LO=M.RG_LO, NRAD=M.NRAD, NLG=len(M.LG))
    print(f"D/r={Dr} seed={seed} R2={d['R2_direct']:.4f} {d['runtime_s']:.0f}s", flush=True)
    return d

if __name__ == "__main__":
    jobs = [(Dr, s) for Dr in (0.9, 1.01, 1.3) for s in (11, 29)]
    with Pool(6) as p: res = p.map(fit, jobs, chunksize=1)
    best = {}
    for d in res:
        if d["D_over_r"] not in best or d["R2_direct"] > best[d["D_over_r"]]["R2_direct"]: best[d["D_over_r"]] = d
    json.dump(dict(doc=__doc__, results=list(best.values()), all_runs=res), open(f"{M.OUT}/R2_numeric_Dr_0p9_1p01_1p3.json", "w"), indent=1)
    print("done")
