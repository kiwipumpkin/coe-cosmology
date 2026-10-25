import json, sys, time
import prof_core as P
P.LB["DR"] = (1.05, 400.0); P.LOGV["DR"] = True
key = sys.argv[1]; Ob = float(sys.argv[2]); cmax = None if sys.argv[3] == "none" else float(sys.argv[3])
t0 = time.time()
seeds = [dict(DR=dr, LD=ld, mc=mc, l0=13.0) for dr in (3.0, 10.0, 50.0, 200.0) for ld in (0.01, 0.1) for mc in (0.01, 0.05, 0.3)]
p = P.fit("cored", "mig", ["DR", "LD", "mc", "l0"], {"Ob": Ob}, starts=5, nscreen=200, seeds=seeds, norm="mean", cmax=cmax)
d = None if p is None else P.describe("cored", p)
json.dump(d, open("cored_%s.json" % key, "w"), indent=1)
print(key, round(time.time() - t0), None if d is None else {k: round(v, 4) if isinstance(v, float) else v for k, v in d.items()}, flush=True)
