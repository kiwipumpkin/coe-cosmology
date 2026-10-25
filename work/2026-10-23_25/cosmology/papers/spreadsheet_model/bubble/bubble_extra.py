import json, time
from bubble_core import *
out = json.load(open("bubble_results.json")); t0 = time.time()
seeds = [dict(DR=2.9, LD=0.07, mc=0.0273, l0=13.2, Ob=0.3), dict(DR=2.83, LD=0.0877, mc=0.039, l0=13.1, Ob=0.3), dict(DR=1.5, LD=0.5, mc=1.0, l0=12.5, Ob=0.3)]
jobs = [("kR|c=0.02|a|Ob=0.3|mig", "kR", 0.02, False, 0.3), ("kR|c=0.02|a|Ob=0.3|mig_l0", "kR", 0.02, True, 0.3),
        ("mean|c=None|a|Ob=0.3|mig_l0", "mean", None, True, 0.3), ("mean|c=None|a|Ob=0|mig_l0", "mean", None, True, 0.0)]
for key, norm, cmax, l0, ob in jobs:
    if key in out: continue
    free = ["DR", "LD", "mc"] + (["l0"] if l0 else []); fixed = {"Ob": ob}; fixed.update({} if l0 else {"l0": 0.0})
    r = fit("mig", free, fixed, starts=4, nscreen=150, seeds=seeds, norm=norm, cmax=cmax)
    out[key] = None if r is None else describe(r[0], norm)
    print(key, None if r is None else {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[key].items()}, round(time.time() - t0), flush=True)
    json.dump(out, open("bubble_results.json", "w"), indent=1)
p = dict(DR=2.904374576719836, LD=0.07027609488588729, mc=0.02726, l0=13.166)
for ob in (0.0, 0.3):
    for norm in ("mean", "kR"):
        q = dict(p, Ob=ob); v, a, ph, rr = score(q, "mig", norm); q.update(amp=a, phi=ph, R2=v)
        out["prev|Ob=%g|%s" % (ob, norm)] = describe(q, norm); print("prev", ob, norm, out["prev|Ob=%g|%s" % (ob, norm)]["R2_mig"], flush=True)
json.dump(out, open("bubble_results.json", "w"), indent=1)
print("done")
