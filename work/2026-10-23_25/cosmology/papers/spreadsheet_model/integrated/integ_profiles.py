import sys, json, time, numpy as np
sys.path.insert(0, "."); from integ_scan import *
res = json.load(open("integ_results.json")); gb = res["mig"]; hisM = res["his_geom_Mfit"]["mc"]
CUR = dict(DR=5.0, LD=0.18, mc=hisM, amp=-0.06, phi=0.35, l0=0.0)
GR = dict(DR=np.linspace(1.05, 10, 30), LD=np.geomspace(0.02, 3, 30), mc=np.geomspace(1e-3, 10, 30),
          amp=np.linspace(-0.15, 0.15, 31), phi=np.linspace(0, 2*np.pi, 31), l0=np.arange(0, 361, 10.0))
prof = {}; t0 = time.time()
for n, g in GR.items():
    fx, pr, pars = [], [], []; prev = gb
    for v in g:
        c = dict(CUR); c[n] = v
        fx.append(score(c["DR"], c["LD"], c["mc"]/MCRIT_FAC, "mig", c["l0"], amp=c["amp"], phi=c["phi"])[0])
        if n in ("DR", "LD", "mc", "l0"):
            free = [m for m in ("DR", "LD", "mc") if m != n] + ([] if n == "l0" else [])
            fixed = {n: v}; fixed.setdefault("l0", 0.0)
            p = fit("mig", free, fixed, starts=3, nscreen=60, seeds=[prev, gb])
        else:
            kw = {"amp": v} if n == "amp" else {"phi": v}
            p = fit("mig", ["DR", "LD", "mc"], {"l0": 0.0}, starts=3, nscreen=60, seeds=[prev, gb], **kw)
        pr.append(p["R2"]); pars.append({k: float(x) for k, x in p.items()}); prev = p
    prof[n] = dict(grid=list(map(float, g)), fixed=list(map(float, fx)), profile=list(map(float, pr)), pars=pars)
    print(n, "prof max %.4f at %.4g; fixed max %.4f at %.4g" % (max(pr), g[int(np.argmax(pr))], max(fx), g[int(np.argmax(fx))]), round(time.time()-t0), flush=True)
json.dump(prof, open("integ_profiles.json", "w"))
# heat map D/R x L/D, M refitted (1D) with amp/phi analytic
DRg = np.linspace(1.1, 10, 28); LDg = np.geomspace(0.02, 3, 28); Z = np.zeros((len(LDg), len(DRg))); MC = np.zeros_like(Z)
for j, ld in enumerate(LDg):
    for i, dr in enumerate(DRg):
        lm = np.linspace(np.log(1e-3), np.log(10), 25)
        v = [score(dr, ld, np.exp(x)/MCRIT_FAC, "mig")[0] for x in lm]; k = int(np.argmax(v))
        lo, hi = lm[max(k-1, 0)], lm[min(k+1, len(lm)-1)]
        r = minimize_scalar(lambda x: -score(dr, ld, np.exp(x)/MCRIT_FAC, "mig")[0], bounds=(lo, hi), method="bounded", options=dict(xatol=1e-4))
        Z[j, i] = max(-r.fun, v[k]); MC[j, i] = np.exp(r.x)
print("map max", Z.max(), round(time.time()-t0))
np.savez("integ_map.npz", DR=DRg, LD=LDg, Z=Z, MC=MC)
