import json
from bubble_core import *
out = json.load(open("bubble_results.json"))
seeds = [dict(DR=2.86, LD=0.085, mc=0.039, l0=13.1, Ob=0.5), dict(DR=2.83, LD=0.1, mc=0.05, l0=13.1, Ob=0.7)]
key = "mean|c=0.005|b|Ob=free|mig_l0"
r = fit("mig", ["DR", "LD", "mc", "l0", "Ob"], {}, starts=4, nscreen=150, seeds=seeds, norm="mean", cmax=0.005, need_bound=True, need_age=True)
out[key] = None if r is None else describe(r[0], "mean"); print(key, out[key], flush=True)
json.dump(out, open("bubble_results.json", "w"), indent=1)
