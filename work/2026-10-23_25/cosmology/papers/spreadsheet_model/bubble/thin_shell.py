"""Chris's sketch: constant-density shell hugging the vision sphere, R_in = D - kin R, R_out = D + kout R, optional compact core M_c (Lc = 0.005 D).
Built on ptshell (same integrator); rin/rout are tied to D/R."""
import json, sys, time, numpy as np
import ptshell as S, prof_core as P
P.LB.update(Os=(0.003, 3.0), mcc=(1e-5, 3.0)); P.LOGV.update(Os=True, mcc=True)
def geom(p, kin, kout):
    q = dict(p); q["rin"] = 1 - kin / p["DR"]; q["rout"] = 1 + kout / p["DR"]; q.setdefault("Lc", 0.005); return q
def fit(free, fixed, kin, kout, norm="mean", seeds=(), starts=4, nscreen=160):
    def full(u): p = dict(fixed); p.update(P.unpack(free, u)); return geom(p, kin, kout)
    def obj(u):
        pen = sum(10*(max(0, -x) + max(0, x-1))**2 for n, x in zip(free, u) if n != "l0")
        v = S.score(full(u), "mig", norm, inside=False)[0]
        return (1 - v if np.isfinite(v) else 50.0) + pen
    rng = np.random.default_rng(7)
    U = [rng.random(len(free)) for _ in range(nscreen)] + [P.pack(free, dict(fixed, **s)) for s in seeds]
    vals = np.array([obj(u) for u in U]); order = np.argsort(vals); best = (np.inf, None)
    from scipy.optimize import minimize
    for i in order[:starts]:
        if vals[i] >= 50: continue
        r = minimize(obj, U[i], method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-7, maxiter=700))
        if r.fun < best[0]: best = (r.fun, r.x)
    if best[1] is None: return None
    p = full(best[1]); v, a, ph, rr = S.score(p, "mig", norm, inside=False); p.update(amp=a, phi=ph, R2=v)
    d = S.describe(p, norm); d.update(kin=kin, kout=kout, norm=norm); return d
JOBS = {
 "a|free|hug": (["DR", "Os", "l0"], {"mcc": 0.0}, 1.0, 1.0, "mean"),
 "a|0.05|hug": (["DR", "l0"], {"mcc": 0.0, "Os": 0.05}, 1.0, 1.0, "mean"),
 "b|free|hug": (["DR", "Os", "mcc", "l0"], {}, 1.0, 1.0, "mean"),
 "b|0.05|hug": (["DR", "mcc", "l0"], {"Os": 0.05}, 1.0, 1.0, "mean"),
 "a|free|marg": (["DR", "Os", "l0"], {"mcc": 0.0}, 1.1, 1.1, "mean"),
 "a|0.05|marg": (["DR", "l0"], {"mcc": 0.0, "Os": 0.05}, 1.1, 1.1, "mean"),
 "b|free|marg": (["DR", "Os", "mcc", "l0"], {}, 1.1, 1.1, "mean"),
 "b|0.05|marg": (["DR", "mcc", "l0"], {"Os": 0.05}, 1.1, 1.1, "mean"),
 "a|free|hug|p2": (["DR", "Os", "l0"], {"mcc": 0.0}, 1.0, 1.0, "mean2"),
 "b|0.05|hug|p2": (["DR", "mcc", "l0"], {"Os": 0.05}, 1.0, 1.0, "mean2"),
}
if __name__ == "__main__":
    key = sys.argv[1]; free, fixed, kin, kout, norm = JOBS[key]; t0 = time.time()
    seeds = []
    for DR in (1.5, 2.9, 5.0):
        for l0 in (13.0, 103.0, 193.0, 283.0):
            seeds.append(dict(DR=DR, Os=0.3, mcc=0.002, l0=l0)); seeds.append(dict(DR=DR, Os=0.05, mcc=0.02, l0=l0))
    d = fit(free, fixed, kin, kout, norm, seeds=seeds)
    json.dump(d, open("thin_%s.json" % key.replace("|", "_"), "w"), indent=1)
    print(key, round(time.time() - t0), {k: round(v, 4) for k, v in d.items() if isinstance(v, float)} if d else None, flush=True)
