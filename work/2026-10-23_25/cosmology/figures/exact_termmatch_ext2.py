#!/usr/bin/env python3
"""Second pass for the extended-density part: bounded LSQ (D/L in [0.3,200]; D/L->inf is the point mass) at the fixed D/r,
and a trace of the exact 4x4 extended solution family near D/r ~ 1.0-1.2. Merges into exact_termmatch_results.json."""
import sys, json, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, "/workspace/cosmology/figures")
from exact_termmatch import *

if __name__ == "__main__":
    t0 = time.time()
    R = json.load(open(OUT + "_results.json"))
    FAM_DRS = [1.01, 1.03, 1.05, 1.07, 1.09, 1.11, 1.12, 1.13, 1.14, 1.145, 1.15, 1.16, 1.18, 1.2, 1.25]
    with Pool(6) as pool:
        ext = pool.map(ext_grid, [(1 / Dr,) for Dr in FAM_DRS])
    fam = []
    for (eps, roots, F1e, F2e, bg), Dr in zip(ext, FAM_DRS):
        for b0, lx in roots:
            x = float(np.exp(lx)); t = tau_ext(b0, eps, x); N, K = solve_amp(t); f, al, kb = E.shape_exp(x)
            s = package("extended exact 4x4", t, N, K, eps, b0, dict(D_over_L=x, L_Gpc=R_GPC / eps / x, alpha=al, kappa=K / al,
                        MoverMcrit_ST=2 * K / al, branch="on-axis" if t[1] > 0 else "flipped"))
            fam.append(clean(s))
        print(Dr, [(round(b, 2), round(float(np.exp(l)), 2)) for b, l in roots], time.time() - t0, flush=True)
    R["ext_family"] = fam
    for Dr in list(R["fixed_ext"].keys()):
        v = R["fixed_ext"][Dr]
        if v["exact"]: continue
        eps = 1 / float(Dr)
        xv, c2 = ext_lsq(eps, SIG["naive"]); x = float(np.exp(xv[1])); t = tau_ext(xv[0], eps, x)
        w = {m: 1 / SIG["naive"][m] ** 2 for m in (1, 2, 3)}; c = {m: abs(t[m]) for m in (1, 2, 3)}
        P = sum(w[m] * amp[m] * c[m] for m in c) / sum(w[m] * c[m] ** 2 for m in c); N = amp[0] - P * t[0]; f, al, kb = E.shape_exp(x)
        v["lsq"] = clean(package("extended LSQ", t, N, P / N, eps, xv[0], dict(D_over_L=x, L_Gpc=R_GPC / eps / x, alpha=al, kappa=P / N / al,
                                  MoverMcrit_ST=2 * P / N / al, chi2_fit=float(c2), at_upper_DL_bound=bool(x > 199.0))))
        s = v["lsq"]; print("lsq", Dr, round(s["b0"], 2), round(x, 2), s["model_amp"]["3"], s["R2"], s["chi2"]["naive"], flush=True)
    json.dump(R, open(OUT + "_results.json", "w"))
    print("saved", time.time() - t0)
