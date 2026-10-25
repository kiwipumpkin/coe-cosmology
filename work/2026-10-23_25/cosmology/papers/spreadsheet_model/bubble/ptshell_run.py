import json, time
import ptshell as S
t0 = time.time(); out = {}
try: out = json.load(open("ptshell_results.json"))
except Exception: pass
seeds = [dict(DR=2.83, mcc=0.03, Os=0.3, rin=0.3, rout=2.5, Lc=0.085, l0=13.1), dict(DR=2.83, mcc=0.05, Os=0.3, rin=0.5, rout=2.0, Lc=0.06, l0=13.1),
         dict(DR=2.5, mcc=0.02, Os=0.5, rin=0.2, rout=3.0, Lc=0.1, l0=13.1), dict(DR=2.0, mcc=0.1, Os=0.3, rin=0.4, rout=1.2, Lc=0.08, l0=13.1)]
jobs = [("inside|Os=free|mig_l0", ["DR", "mcc", "Os", "rin", "rout", "Lc", "l0"], {}, True),
        ("inside|Os=0.3|mig_l0", ["DR", "mcc", "rin", "rout", "Lc", "l0"], {"Os": 0.3}, True),
        ("inside|Os=0.05|mig_l0", ["DR", "mcc", "rin", "rout", "Lc", "l0"], {"Os": 0.05}, True),
        ("rout_free|Os=free|mig_l0", ["DR", "mcc", "Os", "rin", "rout", "Lc", "l0"], {}, False),
        ("inside|Os=free|mig", ["DR", "mcc", "Os", "rin", "rout", "Lc"], {"l0": 0.0}, True),
        ("inside|Os=free|mig_l0|AF2", ["DR", "mcc", "Os", "rin", "rout", "Lc", "l0"], {}, True)]
for key, free, fixed, inside in jobs:
    if key in out: continue
    norm = "mean2" if key.endswith("AF2") else "mean"
    p = S.fit(free, fixed, inside=inside, norm=norm, starts=4, nscreen=250, seeds=seeds)
    out[key] = None if p is None else S.describe(p, norm)
    print(key, None if p is None else {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[key].items()}, round(time.time() - t0), flush=True)
    json.dump(out, open("ptshell_results.json", "w"), indent=1)
print("done")
