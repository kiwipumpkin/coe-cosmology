import json, numpy as np
import bound_core as B
from bound_core import *
TH0 = 13.967499547347456; HTMIN = 12.0 / TH0
_orig = B.score
def score_age(p, target, norm="kR", need_bound=False):
    v, a, ph, rr = _orig(p, target, norm, need_bound)
    if np.isfinite(v) and rr["vE"] / rr["rE"] < HTMIN: return -np.inf, None, None, rr
    return v, a, ph, rr
B.score = score_age
res = json.load(open("bound_results.json"))
for tg, l0, ob in [("mig", False, 0.0), ("mig", True, 0.0), ("four", True, 0.0), ("mig", False, 0.19), ("mig", True, 0.19), ("four", True, 0.19)]:
    free = ["DR", "LD", "mc"] + (["l0"] if l0 else []); fixed = {"ob": ob}; fixed.update({} if l0 else {"l0": 0.0})
    seeds = [dict(DR=1.5, LD=0.5, mc=0.35, l0=12.5), dict(DR=1.4, LD=0.4, mc=0.2, l0=12.5), dict(DR=1.5, LD=0.6, mc=0.45, l0=12.5)]
    p = B.fit(tg, free, fixed, "mean", True, starts=8, nscreen=300, seeds=seeds)
    rr = run(p["DR"], p["LD"], p["mc"], ob); Ht = rr["vE"] / rr["rE"]
    af = rr["af"] / rr["af"].mean(); afs = shifted(af, p["l0"]) if p["l0"] else af; A = A_from(afs, TH1 - p["l0"], p["amp"], p["phi"])
    p.update(R2_mig=float(r2(MA, A[IMIG])), R2_four=float(r2(AFfour(TH1), A)), t_turn=t_turn(p["DR"], p["LD"], p["mc"], ob), Ht=Ht, age_Gyr=Ht * TH0,
             M_over_Mcrit_today=2 * p["mc"] / MCRIT_FAC / Ht**2, ob_today=2 * ob / MCRIT_FAC / Ht**2)
    k = "age12|mean|ob=%g|%s%s" % (ob, tg, "_l0" if l0 else ""); res["fits"][k] = {a: float(b) for a, b in p.items()}
    print(k, {a: round(float(b), 4) for a, b in p.items()}, flush=True)
json.dump(res, open("bound_results.json", "w"), indent=1)
