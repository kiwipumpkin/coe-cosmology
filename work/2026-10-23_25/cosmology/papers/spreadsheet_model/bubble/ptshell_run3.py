import json, time
import ptshell as S
t0 = time.time(); out = json.load(open("ptshell_results.json"))
seeds = [dict(DR=2.74, mcc=0.03, Os=0.3, rin=0.2, rout=3.0, Lc=0.06, l0=13.1), dict(DR=3.2, mcc=0.003, Os=0.3, rin=0.05, rout=3.0, Lc=0.006, l0=13.2),
         dict(DR=3.0, mcc=0.01, Os=0.05, rin=0.1, rout=3.0, Lc=0.02, l0=13.1), dict(DR=3.2, mcc=0.002, Os=0.05, rin=0.05, rout=3.0, Lc=0.005, l0=13.2)]
jobs = [("inside|Os=0.3|mig_l0|c=0.02", ["DR", "mcc", "rin", "rout", "Lc", "l0"], {"Os": 0.3}),
        ("inside|Os=0.05|mig_l0|c=0.02", ["DR", "mcc", "rin", "rout", "Lc", "l0"], {"Os": 0.05}),
        ("inside|Os=0.05|mig|c=0.02", ["DR", "mcc", "rin", "rout", "Lc"], {"Os": 0.05, "l0": 0.0})]
for key, free, fixed in jobs:
    if key in out: continue
    p = S.fit(free, fixed, inside=True, starts=4, nscreen=200, seeds=seeds, cmax=0.02)
    out[key] = None if p is None else S.describe(p)
    print(key, None if p is None else {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[key].items()}, round(time.time() - t0), flush=True)
    json.dump(out, open("ptshell_results.json", "w"), indent=1)
print("done")
