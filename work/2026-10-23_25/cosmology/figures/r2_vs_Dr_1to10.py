#!/usr/bin/env python3
"""Extension of r2_vs_Dr_density_families.py to 1 < D/r <= 10 (same model, same direct R^2).
Families: exp core + floor (kappa free, code convention), exp core + floor with rho(D)=0.3 rho_crit
exactly (kappa = kappa_base), pure exponential (kappa free).
Changes for this run: u(r) table on [0.02, 2.0] D (80 pts) because D/r=1.05 needs r down to 0.048 D;
D/L upper bound widened 40 -> 300 (the earlier fits hit 40 at D/r 1.5-2).
Reduced optimiser (time budget): 1 differential-evolution run (popsize 12, maxiter 120) +
warm starts from R2_vs_Dr_density_families_raw.json (nearest D/r, all exp-type solutions) +
4 random starts, each Nelder-Mead polished twice."""
import json, time, numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
import r2_vs_Dr_density_families as M
M.RG_LO, M.NRAD = 0.02, 80
for f in ("a_exp_floor", "d_pure_exp"): M.FAM[f]["hi"][0] = 300.0
OUT = M.OUT
DRS = [1.05, 1.1, 1.2, 1.32, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10]
old = json.load(open(f"{OUT}/R2_vs_Dr_density_families_raw.json"))["results"]

def warm(fam, Dr):
    near = min(set(r["D_over_r"] for r in old), key=lambda d: abs(np.log(d / Dr)))
    S = []
    for r in old:
        if r["D_over_r"] != near or r["family"] not in ("a_exp_floor", "d_pure_exp"): continue
        x = list(r["x"])
        if fam == "a_exp_floor" and len(x) == 4: x = x + [0.0]
        if fam == "d_pure_exp": x = x[:4]
        S.append(np.array(x))
    return S

def fit(args):
    fam, Dr, strict, seed = args
    eps = 1 / Dr; b = M.bounds(fam); rng = np.random.default_rng(seed); t0 = time.time()
    de = differential_evolution(M.objective, b, args=(fam, eps, strict), seed=seed, popsize=12, maxiter=120,
                                tol=1e-6, polish=False, init="sobol")
    S = [de.x] + warm(fam, Dr) + [np.array([rng.uniform(lo, hi) for lo, hi in b]) for _ in range(4)]
    P = []
    for x0 in S:
        r = minimize(M.objective, x0, args=(fam, eps, strict), method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9, maxfev=4000, adaptive=True))
        r = minimize(M.objective, r.x, args=(fam, eps, strict), method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-10, maxfev=2000, adaptive=True))
        P.append((-r.fun, r.x))
    P.sort(key=lambda t: -t[0])
    d = M.describe(fam, eps, P[0][1], strict)
    d["L_Gpc"] = d["D_Gpc"] / d["density"]["D_over_L"]
    d["strict_local_density"] = strict; d["all_starts_R2"] = [float(p[0]) for p in P]; d["runtime_s"] = time.time() - t0
    print(f"{fam} strict={strict} D/r={Dr} R2={d['R2_direct']:.4f} {time.time()-t0:.0f}s", flush=True)
    return d

if __name__ == "__main__":
    jobs = [("a_exp_floor", Dr, s, 100 + i) for i, Dr in enumerate(DRS) for s in (False, True)] + \
           [("d_pure_exp", Dr, False, 300 + i) for i, Dr in enumerate(DRS)]
    with Pool(8) as p: res = p.map(fit, jobs, chunksize=1)
    json.dump(dict(doc=__doc__, results=res), open(f"{OUT}/R2_vs_Dr_1to10.json", "w"), indent=1)
    print("done")
