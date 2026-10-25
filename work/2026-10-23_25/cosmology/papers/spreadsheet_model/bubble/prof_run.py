import json, time
from prof_core import *
t0 = time.time(); out = {}
try: out = json.load(open("prof_results.json"))
except Exception: pass
PAR = dict(pl=["LD", "n"], gauss=["LD"], cored=["LD"], shell=["rs", "w"])
SEED = dict(pl=[dict(DR=2.8, LD=0.03, n=4.0, mc=0.04, l0=13), dict(DR=2.8, LD=0.1, n=5.5, mc=0.04, l0=13)],
            gauss=[dict(DR=2.8, LD=0.1, mc=0.04, l0=13), dict(DR=2.8, LD=0.2, mc=0.04, l0=13)],
            cored=[dict(DR=2.8, LD=0.03, mc=0.04, l0=13), dict(DR=2.8, LD=0.01, mc=0.04, l0=13)],
            shell=[dict(DR=2.8, rs=0.2, w=0.05, mc=0.04, l0=13), dict(DR=2.8, rs=0.4, w=0.05, mc=0.04, l0=13), dict(DR=2.8, rs=0.1, w=0.1, mc=0.04, l0=13)])
for cmax in (0.02, 0.005):
    for kind in ("pl", "gauss", "cored", "shell"):
        for l0 in (False, True):
            key = "%s|c=%s|%s" % (kind, cmax, "mig_l0" if l0 else "mig")
            if key in out: continue
            free = ["DR", "mc"] + PAR[kind] + (["l0"] if l0 else [])
            fixed = {"Ob": 0.3}; fixed.update({} if l0 else {"l0": 0.0})
            p = fit(kind, "mig", free, fixed, starts=5, nscreen=200, seeds=SEED[kind], cmax=cmax)
            out[key] = None if p is None else describe(kind, p)
            print(key, None if p is None else {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[key].items()}, round(time.time() - t0), flush=True)
            json.dump(out, open("prof_results.json", "w"), indent=1)
print("done")
