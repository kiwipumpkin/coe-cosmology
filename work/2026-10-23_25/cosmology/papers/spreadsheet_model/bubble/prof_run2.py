import json, time
from prof_core import *
t0 = time.time(); out = {}
try: out = json.load(open("prof_results.json"))
except Exception: pass
PAR = dict(pl=["LD", "n"], gauss=["LD"], cored=["LD"], shell=["rs", "w"], st53=["rst"], st43=["rst"], exp=["LD"])
SEED = dict(pl=[dict(DR=2.8, LD=0.02, n=6.0, mc=0.04, l0=13), dict(DR=2.8, LD=0.1, n=7.5, mc=0.04, l0=13)],
            gauss=[dict(DR=2.8, LD=0.1, mc=0.04, l0=13), dict(DR=2.8, LD=0.2, mc=0.04, l0=13)],
            cored=[dict(DR=2.8, LD=0.01, mc=0.04, l0=13), dict(DR=2.8, LD=0.005 * 2, mc=0.03, l0=13)],
            shell=[dict(DR=2.8, rs=0.2, w=0.05, mc=0.04, l0=13), dict(DR=2.8, rs=0.4, w=0.05, mc=0.04, l0=13)],
            st53=[dict(DR=2.8, rst=0.3, mc=0.04, l0=13), dict(DR=2.8, rst=0.55, mc=0.04, l0=13), dict(DR=2.8, rst=3.0, mc=0.04, l0=13)],
            st43=[dict(DR=2.8, rst=0.3, mc=0.04, l0=13), dict(DR=2.8, rst=0.55, mc=0.04, l0=13), dict(DR=2.8, rst=3.0, mc=0.04, l0=13)],
            exp=[dict(DR=2.83, LD=0.0877, mc=0.039, l0=13.1), dict(DR=2.83, LD=0.0877, mc=0.02, l0=13.1)])
jobs = [(k, 0.02, True, "mean", 0.3) for k in ("cored", "st53", "st43", "shell")]
jobs += [("exp", 0.02, True, "mean", 0.05), ("exp", 0.02, False, "mean", 0.05), ("exp", 0.02, True, "mean2", 0.3)]
for kind, cmax, l0, norm, ob in jobs:
    key = "%s|c=%s|%s|%s" % (kind, cmax, "mig_l0" if l0 else "mig", norm) + ("" if ob == 0.3 else "|Ob=%g" % ob)
    if key in out: continue
    free = ["DR", "mc"] + PAR[kind] + (["l0"] if l0 else [])
    fixed = {"Ob": ob}; fixed.update({} if l0 else {"l0": 0.0})
    import prof_core
    p = fit(kind, "mig", free, fixed, starts=3, nscreen=120, seeds=SEED[kind], cmax=cmax, norm=norm)
    out[key] = None if p is None else describe(kind, p, norm)
    print(key, None if p is None else {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[key].items()}, round(time.time() - t0), flush=True)
    json.dump(out, open("prof_results.json", "w"), indent=1)
print("done")
