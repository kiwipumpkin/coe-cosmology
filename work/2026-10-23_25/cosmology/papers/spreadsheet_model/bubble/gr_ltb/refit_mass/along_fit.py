"""Along-axis (H_r > H_perp) valid solutions: ray-evaluate the flipsearch sample points that are along & valid, then NM-refine the best."""
import json, numpy as np, multiprocessing as mp
from scipy.optimize import minimize
import rm_core as C
from refine_mc import ev, good
def one(r):
    try: o = C.evaluate(r["Mc"], r["D"], r["R"])
    except Exception as e: o = dict(fail=str(e))
    return o
def refine(a):
    o0, key = a; x0 = np.array([np.log10(o0["Mc_kg"]), np.log10(o0["D"]), np.log10(o0["D"] / o0["R"] - 1)])
    def f(x):
        o = ev(10**x[0], x[1:]); return -o[key] if good(o, True) else 5.0
    r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=2e-3, fatol=1e-4, maxiter=120))
    return ev(10**r.x[0], r.x[1:])
if __name__ == "__main__":
    S = [r for r in json.load(open("flipsearch.json"))["sample"] if r["ok"] and r["Hr"] > r["Hp"]]
    with mp.Pool(6) as P:
        R = P.map(one, S, chunksize=4)
        V = [o for o in R if good(o, True)]
        starts = [(o, "R2_ns") for o in sorted(V, key=lambda o: -o["R2_ns"])[:3]] + [(o, "R2_s") for o in sorted(V, key=lambda o: -o["R2_s"])[:3]]
        F = P.map(refine, starts)
    json.dump(dict(sample=R, refined=F), open("along_fit.json", "w"), default=float)
    print("along&valid with rays:", len(V))
    for o in F:
        if "R2_ns" in o: print("Mc %.3e D %.3f R %.3f R2ns %.3f R2s %.3f l0 %.1f along %s 2E %.3f 2Emin %.3f age %.2f Hp %.1f Hr %.1f m2 %.3f" % (o["Mc_kg"], o["D"], o["R"], o["R2_ns"], o["R2_s"], o["fits"]["0.1"]["l0_ns"], o["along"], o["twoE"], o["twoE_min"], o["age_Gyr"], o["H_perp"], o["H_r"], o["fits"]["0.1"]["m2"]))
