"""Stage 2: per Mc, Nelder-Mead over (log D, log(D/R - 1)) from the best valid grid points; l0 (and Amp, phi) fitted inside.
Valid = rays finite and 1 + 2E > 0 on every matter shell. Objective: R2 at z = 0.1 without sine ('ns') and with sine ('s').
Also mode 'along': same but only solutions with H_r > H_perp at Earth (peaks along the CoE axis) are allowed."""
import json, sys, numpy as np, multiprocessing as mp, collections
from scipy.optimize import minimize
import rm_core as C
def ev(Mc, x):
    D = 10**x[0]; R = D / (1 + 10**x[1])
    try: return C.evaluate(Mc, D, R)
    except Exception as e: return dict(fail=str(e))
def good(o, along): return ("R2_ns" in o) and o["twoE_min"] > -1 and (o["along"] or not along)
def run(a):
    Mc, starts, key, along = a; best = None
    for x0 in starts:
        cache = {}
        def f(x):
            o = ev(Mc, x); cache[tuple(x)] = o
            return -o[key] if good(o, along) else 5.0
        r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=2e-3, fatol=1e-4, maxiter=70, initial_simplex=[x0, x0 + [0.15, 0], x0 + [0, 0.25]]))
        o = ev(Mc, r.x)
        if good(o, along) and (best is None or o[key] > best[key]): best = o
    return dict(Mc=Mc, key=key, along=along, best=best)
if __name__ == "__main__":
    G = json.load(open("grid_mc.json")); by = collections.defaultdict(list)
    for o in G:
        if good(o, False): by[o["Mc_kg"]].append(o)
    tasks = []
    for Mc, L in by.items():
        for key in ("R2_ns", "R2_s"):
            L2 = sorted(L, key=lambda o: -o[key])[:3]
            tasks.append((Mc, [np.array([np.log10(o["D"]), np.log10(o["D"] / o["R"] - 1)]) for o in L2], key, False))
    with mp.Pool(6) as P: res = P.map(run, tasks, chunksize=1)
    json.dump(res, open("refine_mc.json", "w"), default=float)
    for r in sorted(res, key=lambda r: (r["key"], r["Mc"])):
        b = r["best"]
        print("%s Mc %.2e  R2ns %.3f R2s %.3f D %.2f R %.2f l0ns %.1f l0s %.1f along %s 2E %.2f age %.2f Hp %.1f Hr %.1f m2 %.3f" % (r["key"], r["Mc"], b["R2_ns"], b["R2_s"], b["D"], b["R"], b["fits"]["0.1"]["l0_ns"], b["fits"]["0.1"]["l0_s"], b["along"], b["twoE"], b["age_Gyr"], b["H_perp"], b["H_r"], b["fits"]["0.1"]["m2"]))
