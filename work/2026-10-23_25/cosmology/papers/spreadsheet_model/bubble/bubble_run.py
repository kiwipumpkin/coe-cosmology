import json, time, numpy as np
from bubble_core import *
t0 = time.time(); out = {}
try: out = json.load(open("bubble_results.json"))
except Exception: pass
seeds = [dict(DR=2.9, LD=0.07, mc=0.0273, l0=13.2, Ob=0.3), dict(DR=2.7, LD=0.072, mc=0.025, l0=13.0, Ob=0.3), dict(DR=2.0, LD=0.12, mc=0.1, l0=13.0, Ob=0.3), dict(DR=1.6, LD=0.1, mc=0.3, l0=13.0, Ob=0.3)]
cases = []
for cmax in (0.02, 0.005):
    for mode in ("a", "b"):
        for tg, l0 in (("mig", False), ("mig", True), ("four", True)):
            cases.append(("mean", cmax, mode, tg, l0, "0.3"))
        cases.append(("mean", cmax, mode, "mig", True, "free"))
cases += [("kR", 0.02, "a", "mig", False, "0.3"), ("kR", 0.02, "a", "mig", True, "0.3"), ("mean", None, "a", "mig", True, "0.3"), ("mean", None, "a", "mig", True, "0")]
for norm, cmax, mode, tg, l0, ob in cases:
    key = "%s|c=%s|%s|Ob=%s|%s%s" % (norm, cmax, mode, ob, tg, "_l0" if l0 else "")
    if key in out: continue
    free = ["DR", "LD", "mc"] + (["l0"] if l0 else []) + (["Ob"] if ob == "free" else [])
    fixed = {} if l0 else {"l0": 0.0}
    if ob != "free": fixed["Ob"] = float(ob)
    kw = dict(norm=norm, cmax=cmax, need_bound=(mode == "b"), need_age=(mode == "b"))
    r = fit(tg, free, fixed, starts=6, nscreen=250, seeds=seeds, **kw)
    if r is None: out[key] = None; print(key, "no feasible point", flush=True); continue
    p, rr = r; d = describe(p, norm); out[key] = d
    print(key, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in d.items()}, round(time.time() - t0), flush=True)
    json.dump(out, open("bubble_results.json", "w"), indent=1)
# previous best (Omega_b = 0) with Omega_b = 0.3 added, same M/Mcrit today
p = dict(DR=2.904374576719836, LD=0.07027609488588729, mc=0.02726, l0=13.166)
for ob in (0.0, 0.3):
    for norm in ("mean", "kR"):
        q = dict(p, Ob=ob); v, a, ph, rr = score(q, "mig", norm); q.update(amp=a, phi=ph, R2=v)
        out["prev|Ob=%g|%s" % (ob, norm)] = describe(q, norm)
# background radius check at the main best fit
json.dump(out, open("bubble_results.json", "w"), indent=1)
print("done", round(time.time() - t0))
