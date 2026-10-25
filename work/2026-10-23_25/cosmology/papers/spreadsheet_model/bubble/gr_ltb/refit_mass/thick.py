"""Thick shell: Omega_s = 0.049 (ordinary matter), shell D-h .. D+h with half-thickness h >= 14.3 Gpc (decoupled from sampling),
empty interior + core, vacuum outside. A sampled at z = 0.075 (headline) / 0.3 with the linear model A, A_all = mean over directions.
Free: D, h, Mc/Mcrit(D), l0 (+ Amp, phi). LTB GR-A, H = 70, t_B = 0."""
import json, numpy as np, multiprocessing as mp
from scipy.optimize import minimize
import rm_core as C
OS = 0.049; ZS = {"median": 0.075, "zmax": 0.3}
HS = [14.3, 20.0, 35.0, 70.0]; DH = [1.05, 1.2, 1.5, 2.0, 3.0, 6.0]; MU = np.geomspace(0.002, 2.0, 12)
def ev(case, h, dh, mu):
    z = ZS[case]; C.ZL = (z,); D = dh * h; Mc = mu * 4 * np.pi / 3 * C.RHOC * D**3 * C.KG
    try: o = C.evaluate(Mc, D, h, Os=OS)
    except Exception as e: o = dict(fail=str(e))
    o.update(case=case, h=h, dh=dh, mu=mu, z=z, D=D, Mc_kg=Mc); return o
def valid(o): return "R2_s" in o and o["twoE_min"] > -1 and o["h"] >= 14.3 and (o["D"] - o["h"]) > 5 * 0.005 * o["D"]
def one(a): return ev(*a)
def refine(a):
    case, o0 = a; x0 = np.array([np.log(o0["h"] / 14.3 + 1e-3), np.log(o0["dh"] - 1), np.log(o0["mu"])])
    def f(x):
        o = ev(case, 14.3 * np.exp(x[0]), 1 + np.exp(x[1]), np.exp(x[2])); return -o["R2_s"] if valid(o) else 5.0
    r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=0.03, fatol=5e-4, maxiter=50, initial_simplex=[x0, x0 + [0.3, 0, 0], x0 + [0, 0.3, 0], x0 + [0, 0, 0.3]]))
    return ev(case, 14.3 * np.exp(r.x[0]), 1 + np.exp(r.x[1]), np.exp(r.x[2]))
if __name__ == "__main__":
    T = [(c, h, d, m) for c in ZS for h in HS for d in DH for m in MU]
    with mp.Pool(7) as P:
        res = P.map(one, T, chunksize=3)
        starts = [(c, max([o for o in res if o["case"] == c and valid(o)], key=lambda o: o["R2_s"])) for c in ZS]
        ref = P.map(refine, starts)
    json.dump(dict(grid=res, refined=ref), open("thick.json", "w"), default=float)
    for c in ZS:
        V = sorted([o for o in res if o["case"] == c and valid(o)], key=lambda o: -o["R2_s"])
        print(c, "valid", len(V), "of", sum(o["case"] == c for o in res), "along&valid", sum(o["along"] for o in V))
        for o in V[:5]: print("  grid h %.1f D %.1f mu %.4f R2s %.3f R2ns %.3f along %s age %.2f 2E %.2f" % (o["h"], o["D"], o["mu"], o["R2_s"], o["R2_ns"], o["along"], o["age_Gyr"], o["twoE"]))
        o = [r for r in ref if r["case"] == c][0]
        if "R2_s" in o: print("  REF h %.2f D %.2f mu %.4f R2s %.3f R2ns %.3f along %s valid %s" % (o["h"], o["D"], o["mu"], o["R2_s"], o["R2_ns"], o["along"], valid(o)))
