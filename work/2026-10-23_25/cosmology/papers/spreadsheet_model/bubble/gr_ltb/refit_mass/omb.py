"""Omega_s fixed at 0.049 (ordinary matter, no dark matter); refit D/R, Mc/Mcrit, l0 (+ Amp, phi). Same definitions/grid as wgrid.py."""
import json, numpy as np, multiprocessing as mp
from scipy.optimize import minimize
import wgrid as WG
OS = 0.049; DRS = WG.DRS; MU = np.geomspace(0.005, 3.0, 16)
def one(a): return WG.ev(a[0], a[1], OS, a[2])
def refine(a):
    case, o0 = a; x0 = np.array([np.log(o0["DR"] - 1), np.log(o0["mu"])])
    def f(x):
        o = WG.ev(case, 1 + np.exp(x[0]), OS, np.exp(x[1])); return -o["R2_s"] if WG.valid(o) else 5.0
    r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=0.02, fatol=5e-4, maxiter=60, initial_simplex=[x0, x0 + [0.3, 0], x0 + [0, 0.3]]))
    return WG.ev(case, 1 + np.exp(r.x[0]), OS, np.exp(r.x[1]))
if __name__ == "__main__":
    T = [(c, d, m) for c in WG.CASES for d in DRS for m in MU]
    with mp.Pool(7) as P:
        res = P.map(one, T, chunksize=3)
        starts = [(c, max([o for o in res if o["case"] == c and WG.valid(o)], key=lambda o: -o["R2_s"] * -1)) for c in WG.CASES]
        starts = [(c, max([o for o in res if o["case"] == c and WG.valid(o)], key=lambda o: o["R2_s"])) for c in WG.CASES]
        ref = P.map(refine, starts)
    json.dump(dict(grid=res, refined=ref), open("omb.json", "w"), default=float)
    for c in WG.CASES:
        V = sorted([o for o in res if o["case"] == c and WG.valid(o)], key=lambda o: -o["R2_s"])
        print(c, "valid", len(V), "of", sum(o["case"] == c for o in res))
        for o in V[:5]: print("  grid DR %.2f mu %.4f R2s %.3f R2ns %.3f age %.2f 2E %.3f" % (o["DR"], o["mu"], o["R2_s"], o["R2_ns"], o["age_Gyr"], o["twoE"]))
        o = [r for r in ref if r["case"] == c][0]; print("  REF DR %.3f mu %.4f R2s %.3f R2ns %.3f" % (o["DR"], o["mu"], o["R2_s"], o["R2_ns"]))
