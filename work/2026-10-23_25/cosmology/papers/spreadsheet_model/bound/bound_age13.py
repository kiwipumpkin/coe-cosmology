import json, numpy as np
import bound_core as B
from bound_core import *
TH0 = 13.967499547347456; HTMIN = 13.0 / TH0
_orig = B.score
def score_age(p, target, norm="kR", need_bound=False):
    v, a, ph, rr = _orig(p, target, norm, need_bound)
    if np.isfinite(v) and rr["vE"] / rr["rE"] < HTMIN: return -np.inf, None, None, rr
    return v, a, ph, rr
B.score = score_age
res = json.load(open("bound_results.json"))
p = B.fit("mig", ["DR", "LD", "mc", "l0"], {"ob": 0.0}, "mean", True, starts=8, nscreen=400, seeds=[dict(DR=1.65, LD=0.44, mc=0.45, l0=12.5), dict(DR=2, LD=1.0, mc=0.3, l0=12.5), dict(DR=3, LD=2.5, mc=0.05, l0=12.5)])
if p is None: print("none"); res["fits"]["age13|mean|ob=0|mig_l0"] = None
else:
    rr = run(p["DR"], p["LD"], p["mc"], 0.0); Ht = rr["vE"]/rr["rE"]
    p.update(t_turn=t_turn(p["DR"], p["LD"], p["mc"], 0.0), Ht=Ht, age_Gyr=Ht*TH0, M_over_Mcrit_today=2*p["mc"]/MCRIT_FAC/Ht**2)
    res["fits"]["age13|mean|ob=0|mig_l0"] = {a: float(b) for a, b in p.items()}; print({a: round(float(b), 4) for a, b in p.items()})
json.dump(res, open("bound_results.json", "w"), indent=1)
