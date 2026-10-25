import json, numpy as np
from bound_core import *
res = json.load(open("bound_results.json")); TH0 = 13.967499547347456
for norm, tg, l0 in [("mean","mig",False),("mean","mig",True),("mean","four",False),("mean","four",True),("kR","mig",False)]:
    ob = 0.19
    for it in range(2):
        free = ["DR","LD","mc"] + (["l0"] if l0 else []); fixed = {"ob": ob}; fixed.update({} if l0 else {"l0": 0.0})
        src = res["fits"]["bound|%s|ob=0.3|%s%s" % (norm, tg, "_l0" if l0 else "")] if norm == "mean" else dict(DR=5.5, LD=2.6, mc=0.03, l0=0)
        p = fit(tg, free, fixed, norm, True, starts=6, nscreen=200, seeds=[src])
        rr = run(p["DR"], p["LD"], p["mc"], ob); Ht = rr["vE"]/rr["rE"]; obnow = 2*ob/MCRIT_FAC/Ht**2
        ob = ob * 0.3 / obnow
    af = rr["af"] / (rr["af"].mean() if norm == "mean" else 1.0)
    afs = shifted(af, p["l0"]) if p["l0"] else af; A = A_from(afs, TH1 - p["l0"], p["amp"], p["phi"])
    p.update(R2_mig=float(r2(MA, A[IMIG])), R2_four=float(r2(AFfour(TH1), A)), t_turn=t_turn(p["DR"], p["LD"], p["mc"], p["ob"]), Ht=Ht, age_Gyr=Ht*TH0,
             M_over_Mcrit_today=2*p["mc"]/MCRIT_FAC/Ht**2, ob_today=obnow)
    p["t_turn_Gyr"] = p["t_turn"]*p["age_Gyr"]
    k = "bound|%s|obnow=0.3|%s%s" % (norm, tg, "_l0" if l0 else ""); res["fits"][k] = {a: float(b) for a, b in p.items()}
    print(k, {a: round(float(b), 4) for a, b in p.items()}, flush=True)
json.dump(res, open("bound_results.json", "w"), indent=1)
