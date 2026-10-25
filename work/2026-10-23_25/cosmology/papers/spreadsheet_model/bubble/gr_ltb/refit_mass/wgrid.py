"""Widened coarse grid (D/R down to 1.05, up to 20) + local Nelder-Mead refine around the optimum. Same definitions as rgrid.py."""
import json, numpy as np, multiprocessing as mp
from scipy.optimize import minimize
import rm_core as C
CASES = {"median": (0.316, 0.075), "zmax": (1.19, 0.3)}
DRS = [1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 6.55, 12.0, 20.0]
OSS = [0.005, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]; MU = np.geomspace(0.005, 3.0, 11)
def ev(case, DR, Os, mu):
    R, z = CASES[case]; C.ZL = (z,); D = DR * R; Mc = mu * 4 * np.pi / 3 * C.RHOC * D**3 * C.KG
    try: o = C.evaluate(Mc, D, R, Os=Os)
    except Exception as e: o = dict(fail=str(e))
    o.update(case=case, DR=DR, Os=Os, mu=mu, z=z, D=D, R=R, Mc_kg=Mc); return o
def valid(o): return "R2_s" in o and o["twoE_min"] > -1 and (o["D"] - o["R"]) > 5 * 0.005 * o["D"]
def one(a): return ev(*a)
def refine(a):
    case, o0 = a
    x0 = np.array([np.log(o0["DR"] - 1), np.log(o0["Os"]), np.log(o0["mu"])])
    def f(x):
        o = ev(case, 1 + np.exp(x[0]), np.exp(x[1]), np.exp(x[2])); return -o["R2_s"] if valid(o) else 5.0
    r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=0.02, fatol=5e-4, maxiter=80, initial_simplex=[x0, x0 + [0.3, 0, 0], x0 + [0, 0.5, 0], x0 + [0, 0, 0.3]]))
    return ev(case, 1 + np.exp(r.x[0]), np.exp(r.x[1]), np.exp(r.x[2]))
if __name__ == "__main__":
    T = [(c, d, o, m) for c in CASES for d in DRS for o in OSS for m in MU]
    with mp.Pool(7) as P:
        res = P.map(one, T, chunksize=4)
        starts = []
        for c in CASES:
            V = sorted([o for o in res if o["case"] == c and valid(o)], key=lambda o: -o["R2_s"])
            # distinct starts: best overall + best at other Omega_s values
            seen = set()
            for o in V:
                if o["Os"] in seen: continue
                seen.add(o["Os"]); starts.append((c, o))
                if len(seen) == 3: break
        ref = P.map(refine, starts)
    json.dump(dict(grid=res, refined=ref), open("wgrid.json", "w"), default=float)
    for c in CASES:
        V = [o for o in res if o["case"] == c and valid(o)]
        print(c, "valid", len(V), "of", sum(o["case"] == c for o in res), " invalid reasons:",
              sum(1 for o in res if o["case"] == c and "R2_s" in o and o["twoE_min"] <= -1), "1+2E<=0;", sum(1 for o in res if o["case"] == c and "fail" in o), "no solution")
        for o in sorted(V, key=lambda o: -o["R2_s"])[:8]:
            f = o["fits"][str(o["z"])]
            print("  grid DR %.2f Os %.3f mu %.4f | R2s %.3f R2ns %.3f l0 %.1f | along %s 2E %.3f age %.2f" % (o["DR"], o["Os"], o["mu"], o["R2_s"], o["R2_ns"], f["l0_s"], o["along"], o["twoE"], o["age_Gyr"]))
        for o in [r for r in ref if r["case"] == c]:
            if "R2_s" not in o: print("  refined fail", o.get("fail")); continue
            f = o["fits"][str(o["z"])]
            print("  REF DR %.3f D %.3f Os %.4f mu %.4f Mc %.2e | R2s %.3f R2ns %.3f l0 %.1f amp %.3f phi %.2f | along %s 2E %.3f age %.2f Hp %.1f Hr %.1f m2 %.3f corr %.2f valid %s" % (
                o["DR"], o["D"], o["Os"], o["mu"], o["Mc_kg"], o["R2_s"], o["R2_ns"], f["l0_s"], f["amp"], f["phi"], o["along"], o["twoE"], o["age_Gyr"], o["H_perp"], o["H_r"], f["m2"], f["corr"], valid(o)))
