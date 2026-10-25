import sys, json, numpy as np
sys.path.insert(0, "."); from integ_scan import *
sys.path.insert(0, ".."); import param_scan_core as PS          # impulse model, R = 14.3 Gpc
res = json.load(open("integ_results.json")); gb = res["mig"]; out = {}
mu = gb["mc"]/MCRIT_FAC
def r2mig(cfg, p=gb, refit=False):
    if refit:
        q = fit("mig", ["DR", "LD", "mc"], {"l0": 0.0}, starts=4, nscreen=80, seeds=[p], cfg=cfg); return q
    return score(p["DR"], p["LD"], p["mc"]/MCRIT_FAC, "mig", cfg=cfg)[0]
# t_i sensitivity (fixed profile) and step check
base = af_grid(gb["DR"], gb["LD"], mu)
ti = {}
for t in (0.05, 0.01, 0.001):
    cfg = dict(como=False, taui=t, nsteps=1600)
    a = af_grid(gb["DR"], gb["LD"], mu, cfg)
    q = r2mig(cfg, refit=True)
    ti[str(t)] = dict(max_rel_dAF=float(np.max(np.abs(a - base) / base)), R2_same_params=float(r2mig(cfg)),
                      refit=dict(DR=q["DR"], LD=q["LD"], mc=q["mc"], amp=float(q["amp"]), phi=float(q["phi"]), R2=float(q["R2"])))
    print("ti", t, ti[str(t)], flush=True)
out["taui"] = ti
# shooting: land final positions on his geometry
cfg_s = dict(CFG); cfg_s["shoot"] = 6
a_s = af_grid(gb["DR"], gb["LD"], mu, cfg_s)
_, X0 = AF(TH1, gb["DR"], gb["LD"], mu, **CFG); _, Xs = AF(TH1, gb["DR"], gb["LD"], mu, **cfg_s)
q0 = geometry(TH1, gb["DR"])[0]
qs = fit("mig", ["DR", "LD", "mc"], {"l0": 0.0}, starts=4, nscreen=80, seeds=[gb], cfg=cfg_s)
out["shoot"] = dict(max_rel_dAF=float(np.max(np.abs(a_s - base)/base)), R2_same=float(r2mig(cfg_s)),
                    pos_err_unshot=float(np.abs(X0 - q0).max()), pos_err_shot=float(np.abs(Xs - q0).max()),
                    refit=dict(DR=qs["DR"], LD=qs["LD"], mc=qs["mc"], amp=float(qs["amp"]), phi=float(qs["phi"]), R2=float(qs["R2"])))
print("shoot", out["shoot"], flush=True)
# harmonic content, 1-degree grid
def harm(y):
    c = np.fft.rfft(y) / len(y) * 2; return [float(abs(c[m])) for m in (1, 2, 3)] + [float(np.sqrt(np.sum(abs(c[4:])**2)))]
A_int = A_from(base, TH1, gb["amp"], gb["phi"])
ib = json.load(open("../param_scan_results.json"))["global"]["mig"]["params"]
A_imp = PS.BH(TH1, ib["DR"], ib["LD"], ib["M"], ib["amp"], ib["phi"])
AFF = AFfour(TH1)
out["harmonics"] = dict(A_F=harm(AFF), integrated=harm(A_int), impulse=harm(A_imp),
                        AF_only_integrated=harm(base),
                        shot=harm(A_from(af_grid(res["mig_shoot"]["DR"], res["mig_shoot"]["LD"], res["mig_shoot"]["mc"]/MCRIT_FAC, dict(CFG, shoot=6)), TH1, res["mig_shoot"]["amp"], res["mig_shoot"]["phi"])),
                        comoving=harm(A_from(af_grid(res["como_0.01"]["DR"], res["como_0.01"]["LD"], res["como_0.01"]["mc"]/MCRIT_FAC, dict(como=True, taui=0.01, nsteps=20000)), TH1, res["como_0.01"]["amp"], res["como_0.01"]["phi"])), AF_only_impulse=harm(PS.BH(TH1, ib["DR"], ib["LD"], ib["M"], 0.0, 0.0)))
out["shape"] = dict(R2_int_vs_imp=float(r2(A_imp, A_int)), rms_int_minus_imp=float(np.sqrt(np.mean((A_int-A_imp)**2))))
# impulse AF shape at the integrated geometry, best Gamma: can impulse reproduce integrated AF?
from scipy.optimize import minimize_scalar as ms
f = lambda lm: np.sum((PS.BH(TH1, gb["DR"], gb["LD"], np.exp(lm), 0, 0) - base)**2)
r = ms(f, bounds=(np.log(1e66), np.log(1e72)), method="bounded")
out["shape"]["impulse_AF_at_int_geom_best_R2"] = float(r2(base, PS.BH(TH1, gb["DR"], gb["LD"], np.exp(r.x), 0, 0)))
print(out["harmonics"], out["shape"])
json.dump(out, open("integ_checks.json", "w"), indent=1)
