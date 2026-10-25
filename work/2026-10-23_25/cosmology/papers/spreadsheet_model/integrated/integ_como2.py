import sys, json, numpy as np
sys.path.insert(0, "."); import integ_scan as S
from integ_scan import *
S.LB["mc"] = (1e-7, 10.0)
res = json.load(open("integ_results.json"))
c1 = res["como_0.01"]
for t, ns in ((0.001, 20000), (0.01, 20000)):
    cfg = dict(como=True, taui=t, nsteps=ns)
    seed = dict(c1, mc=c1["mc"] * t / 0.01)
    p = fit("mig", ["DR", "LD", "mc"], {"l0": 0.0}, starts=2, nscreen=30, seeds=[seed], cfg=cfg)
    p["R2_four"] = score(p["DR"], p["LD"], p["mc"]/MCRIT_FAC, "four", cfg=cfg, amp=p["amp"], phi=p["phi"])[0]; p["R2_mig"] = p["R2"]
    # step check
    a = af_grid(p["DR"], p["LD"], p["mc"]/MCRIT_FAC, cfg); b = af_grid(p["DR"], p["LD"], p["mc"]/MCRIT_FAC, dict(cfg, nsteps=4*ns))
    p["step_err"] = float(np.nanmax(np.abs(a - b)))
    res["como_%g" % t] = {k: float(v) for k, v in p.items()}; print(t, res["como_%g" % t], flush=True)
json.dump(res, open("integ_results.json", "w"), indent=1)
