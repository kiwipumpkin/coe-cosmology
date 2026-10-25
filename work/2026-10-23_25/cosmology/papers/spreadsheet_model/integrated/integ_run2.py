import sys, json, time, numpy as np
sys.path.insert(0, "."); import integ_scan as S
from integ_scan import *
t0 = time.time()
exec(open("integ_scan.py").read().split('if __name__ == "__main__":')[1].split("def vrmin_of")[0].replace("\n    ", "\n"))
res = json.load(open("integ_results.json"))
for k in list(res): 
    if isinstance(res[k], dict) and "DR" in res[k]: res[k]["vr_min"] = vrmin_of(res[k])
# shot variant (final positions on his geometry) and comoving variants
cfg_s = dict(CFG, shoot=6)
p = fit("mig", ["DR", "LD", "mc"], {"l0": 0.0}, starts=12, nscreen=300, seeds=[res["mig"]], cfg=cfg_s)
p["R2_four"] = score(p["DR"], p["LD"], p["mc"]/MCRIT_FAC, "four", cfg=cfg_s, amp=p["amp"], phi=p["phi"])[0]; p["R2_mig"] = p["R2"]
res["mig_shoot"] = p; print("shoot", p, round(time.time()-t0), flush=True)
for t in (0.01, 0.001):
    cfg = dict(como=True, taui=t, nsteps=20000)
    p = fit("mig", ["DR", "LD", "mc"], {"l0": 0.0}, starts=8, nscreen=200, seeds=[dict(res["mig"], mc=res["mig"]["mc"]*t)], cfg=cfg)
    p["R2_four"] = score(p["DR"], p["LD"], p["mc"]/MCRIT_FAC, "four", cfg=cfg, amp=p["amp"], phi=p["phi"])[0]; p["R2_mig"] = p["R2"]
    res["como_%g" % t] = p; print("como", t, p, round(time.time()-t0), flush=True)
json.dump({k: ({kk: float(vv) for kk, vv in v.items()} if isinstance(v, dict) else v) for k, v in res.items()}, open("integ_results.json", "w"), indent=1)
