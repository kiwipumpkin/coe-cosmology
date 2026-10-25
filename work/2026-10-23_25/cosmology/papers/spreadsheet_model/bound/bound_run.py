import json, time, numpy as np
from bound_core import *
t0 = time.time(); out = {"fits": {}, "prof": {}}
def finish(p, norm):
    if p is None: return None
    for tg in ("mig", "four"):
        q = dict(p); v, a, ph, rr = score(q, tg, norm)   # R2 with sine refitted per target? keep fitted sine:
    af = run(p["DR"], p["LD"], p["mc"], p.get("ob", 0.0))["af"]
    if norm == "mean": af = af / af.mean()
    l0 = p.get("l0", 0.0); afs = shifted(af, l0) if l0 else af; th = TH1 - l0
    A = A_from(afs, th, p["amp"], p["phi"])
    p["R2_mig"] = float(r2(MA, A[IMIG])); p["R2_four"] = float(r2(AFfour(TH1), A))
    p["t_turn"] = t_turn(p["DR"], p["LD"], p["mc"], p.get("ob", 0.0))
    p["mean_AF"] = float(run(p["DR"], p["LD"], p["mc"], p.get("ob", 0.0))["af"].mean())
    return {k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in p.items()}
cases = []
for norm in ("mean", "kR"):
    for ob in ("0", "0.3", "free"):
        for tg in ("mig", "four"):
            for l0 in (False, True):
                if norm == "kR" and (ob != "0"): continue
                cases.append((norm, ob, tg, l0))
for bound in (False, True):
    for norm, ob, tg, l0 in cases:
        if not bound and ob != "0": continue
        free = ["DR", "LD", "mc"] + (["l0"] if l0 else []) + (["ob"] if ob == "free" else [])
        fixed = {} if l0 else {"l0": 0.0}
        if ob == "0.3": fixed["ob"] = 0.3
        if ob == "0": fixed["ob"] = 0.0
        seeds = [dict(DR=1.5, LD=1.0, mc=1.0, l0=13.0, ob=0.3), dict(DR=2.7, LD=0.072, mc=0.022, l0=13.0, ob=0.0), dict(DR=1.5, LD=3.0, mc=2.0, l0=13.0, ob=0.3)]
        p = fit(tg, free, fixed, norm, bound, starts=8, nscreen=400, seeds=seeds)
        key = "%s|%s|ob=%s|%s%s" % ("bound" if bound else "free", norm, ob, tg, "_l0" if l0 else "")
        out["fits"][key] = finish(p, norm); print(key, out["fits"][key], round(time.time() - t0), flush=True)
json.dump(out, open("bound_results.json", "w"), indent=1)
# profiles of best R2 vs M/Mcrit, constraint (b) only, with (v/vesc)^2 recorded
MCG = np.geomspace(0.01, 10, 19)
for norm, tg, l0, ob in [("mean", "mig", False, 0.0), ("mean", "mig", True, 0.0), ("mean", "four", False, 0.0), ("mean", "four", True, 0.0), ("kR", "mig", False, 0.0), ("mean", "mig", False, 0.3), ("mean", "mig", True, 0.3)]:
    rows = []; prev = None
    for mc in MCG:
        free = ["DR", "LD"] + (["l0"] if l0 else []); fixed = dict(mc=float(mc), ob=ob); fixed.update({} if l0 else {"l0": 0.0})
        sd = [dict(DR=1.5, LD=1.0, l0=13.0), dict(DR=2.7, LD=0.072, l0=13.0), dict(DR=1.5, LD=3.0, l0=13.0)] + ([prev] if prev else [])
        p = fit(tg, free, fixed, norm, False, starts=5, nscreen=150, seeds=sd)
        if p is None: rows.append(dict(mc=float(mc), R2=None)); continue
        prev = p; rows.append({k: float(v) for k, v in p.items() if isinstance(v, (int, float, np.floating))})
    key = "%s|%s%s|ob=%g" % (norm, tg, "_l0" if l0 else "", ob); out["prof"][key] = rows
    print(key, [(round(r["mc"], 3), None if r.get("R2") is None else round(r["R2"], 3), None if r.get("vesc2") is None else round(r["vesc2"], 3)) for r in rows], round(time.time() - t0), flush=True)
    json.dump(out, open("bound_results.json", "w"), indent=1)
