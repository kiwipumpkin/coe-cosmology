import json, numpy as np
from bound_core import *
res = json.load(open("bound_results.json"))
# extra: literal normalisation with Omega_b = 0.3
if "bound|kR|ob=0.3|mig" not in res["fits"]:
    p = fit("mig", ["DR", "LD", "mc"], {"l0": 0.0, "ob": 0.3}, "kR", True, starts=8, nscreen=400,
            seeds=[dict(DR=5.5, LD=2.6, mc=0.03), dict(DR=2.7, LD=0.07, mc=0.02)])
    if p is not None:
        p["t_turn"] = t_turn(p["DR"], p["LD"], p["mc"], 0.3)
        A = A_from(run(p["DR"], p["LD"], p["mc"], 0.3)["af"], TH1, p["amp"], p["phi"])
        p["R2_mig"] = float(r2(MA, A[IMIG])); p["R2_four"] = float(r2(AFfour(TH1), A))
    res["fits"]["bound|kR|ob=0.3|mig"] = None if p is None else {k: float(v) for k, v in p.items()}
TH0 = 1 / (70 / 3.0857e19) / 3.156e16          # 1/H0 in Gyr
G = GSI; Rm = RKM * 1e3
for k, p in res["fits"].items():
    if p is None: continue
    rr = run(p["DR"], p["LD"], p["mc"], p.get("ob", 0.0))
    Ht = rr["vE"] / rr["rE"]; age = Ht * TH0
    Dm = p["DR"] * Rm; t0s = age * 3.156e16
    mu = p["mc"] / MCRIT_FAC
    p.update(Ht=Ht, age_Gyr=age, t_turn_Gyr=p["t_turn"] * age if np.isfinite(p["t_turn"]) else float("inf"),
             M_kg=mu * Dm**3 / (G * t0s**2), M_over_Mcrit_today=2 * mu / Ht**2, Mb_over_Mcrit_today=2 * p.get("ob", 0) / MCRIT_FAC / Ht**2,
             Dgpc=p["DR"] * RKM / GPC, Lgpc=p["LD"] * p["DR"] * RKM / GPC, young=age < 12)
    print(f'{k:28s} R2m {p["R2_mig"]:.3f} R2f {p["R2_four"]:.3f} DR {p["DR"]:.2f} LD {p["LD"]:.3f} mc {p["mc"]:.3g} Mnow {p["M_over_Mcrit_today"]:.3g} ob {p.get("ob",0):.2f} l0 {p.get("l0",0):.1f} amp {p["amp"]:.3f} phi {p["phi"]:.2f} v2 {p["vesc2"]:.3f} Ht {Ht:.3f} age {age:.2f} tturn {p["t_turn"]:.2f} ({p["t_turn_Gyr"]:.1f} Gyr) M {p["M_kg"]:.3g} meanAF {p["mean_AF"] if "mean_AF" in p else 0:.3f}')
# profiles: add age
for key, rows in res["prof"].items():
    for p in rows:
        if p.get("R2") is None: continue
        rr = run(p["DR"], p["LD"], p["mc"], p.get("ob", 0.0)); p["Ht"] = rr["vE"] / rr["rE"]; p["age_Gyr"] = p["Ht"] * TH0
json.dump(res, open("bound_results.json", "w"), indent=1)
print("1/H0 Gyr", TH0)
