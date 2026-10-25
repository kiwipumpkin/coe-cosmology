from dscan_core import *
import multiprocessing as mp, time, sys
T = json.load(open(BUB + "/universe_thin.json")); O = json.load(open(BUB + "/universe_newt_omb.json"))["Rs14.3|sine"]
rng = np.random.default_rng(3)
def unpack(kind, u):
    if kind == "free": lm, Os, l0 = u; return 10**lm, Os, 1.0, 1.0, l0
    lm, ki, ko, l0 = u; return 10**lm, 0.049, ki, ko, l0
def obj(u, kind, D):
    mcc, Os, ki, ko, l0 = unpack(kind, u)
    if kind == "free" and not (0.001 < Os < 20): return 5.0
    if kind == "omb" and not (1.0 <= ki <= 4 and 1.0 <= ko <= 4): return 5.0
    if not (-7 < np.log10(mcc) < 1.5) or not (-20 < l0 < 50): return 5.0
    return 1 - newt_eval(D, mcc, Os, ki, ko, l0)
def fit_newt(args):
    kind, DR = args; D = DR * RV; t = time.time()
    if kind == "free": seeds = [[lm, Os, l0] for lm in np.arange(-5, 1.01, 0.5) for Os in (0.05, 0.3, 1.0, 3.0) for l0 in (10, 14)]
    else: seeds = [[lm, ki, ko, l0] for lm in np.arange(-5, 1.01, 0.5) for ki in (1.0, 1.5) for ko in (1.0, 1.5, 2.5) for l0 in (10, 14)]
    vals = sorted([(obj(s, kind, D), s) for s in seeds], key=lambda x: x[0])
    best = None
    for v, s in vals[:4]:
        r = minimize(obj, s, args=(kind, D), method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-6, maxiter=600))
        if best is None or r.fun < best.fun: best = r
    mcc, Os, ki, ko, l0 = unpack(kind, best.x); R2, ex = newt_eval(D, mcc, Os, ki, ko, l0, True)
    m = masses(D, mcc, Os, ki, ko)
    out = dict(fit=kind, DR=DR, D=D, R2=R2, mcc=mcc, Os=Os, kin=ki, kout=ko, l0=l0 % 360, **m, **{k: float(v) for k, v in ex.items()},
               VO_c=70 * ex["rE"] * D * MPC_PER_GPC / c, t0_Gyr=ex["Ht"] * 977.8 / 70, sec=time.time() - t)
    print(kind, "D/R %.1f R2 %.4f mcc %.4g Os %.3g kin %.2f kout %.2f share %.4f (%.0fs)" % (DR, R2, mcc, Os, ki, ko, m["share"], out["sec"]), flush=True)
    return out
def fit_rel(DR):
    import relfree as RF
    D = DR * RV; t = time.time()
    def unp(u): return 10**u[0] * 4 * np.pi / 3 * 9.2042e-27 * (D * 3.0857e25)**3, 1 + np.exp(u[1]), 1 + np.exp(u[2])
    def ob(u):
        Mc, ki, ko = unp(u); o = RF.ev(Mc, D, ki, ko); return -o["R2_s"] if "R2_s" in o else 5.0
    starts = [[np.log10(m), np.log(0.04), np.log(0.17)] for m in (0.02, 0.08, 0.3)]
    sv = sorted([(ob(s), s) for s in starts]); best = None
    for v, s in sv[:2]:
        r = minimize(ob, s, method="Nelder-Mead", options=dict(xatol=0.01, fatol=3e-4, maxiter=90))
        if best is None or r.fun < best.fun: best = r
    Mc, ki, ko = unp(best.x); o = RF.ev(Mc, D, ki, ko)
    rs = 0.049 * 9.2042e-27; Rin, Rout = D - ki * RV, D + ko * RV; Msh = 4 * np.pi / 3 * rs * ((Rout * GPCm)**3 - (Rin * GPCm)**3)
    out = dict(fit="rel", DR=DR, D=D, R2=o.get("R2_s"), R2geo=o.get("R2_ns"), mcc=o.get("mcc"), Os=0.049, kin=ki, kout=ko, l0=o.get("l0_s"), amp=o.get("amp"), phi=o.get("phi"),
               Mc=Mc, Mshell=Msh, share=Mc / (Mc + Msh), R_in=Rin, R_out=Rout, VO_c=70 * D * MPC_PER_GPC / c, t0_Gyr=o.get("age_Gyr"), H_r=o.get("H_r"), H_perp=o.get("H_perp"),
               fail=o.get("fail"), sec=time.time() - t)
    print("rel D/R %.1f R2 %s mcc %s share %.4f (%.0fs)" % (DR, out["R2"], out["mcc"], out["share"], out["sec"]), flush=True)
    return out
if __name__ == "__main__":
    which = sys.argv[1]
    if which == "newt":
        jobs = [("free", d) for d in (T["free"]["DR"], 5, 10, 20, 50, 100, 200)] + [("omb", d) for d in (O["D"] / RV, 5, 10, 20, 50, 100, 200)]
        with mp.Pool(7) as P: res = P.map(fit_newt, jobs, chunksize=1)
    else:
        R = json.load(open(BUB + "/universe_rel_mcfree.json"))
        with mp.Pool(7) as P: res = P.map(fit_rel, [R["D"] / RV, 5, 10, 20, 50, 100, 200], chunksize=1)
    json.dump(res, open("/workspace/cosmology/papers/A_equals_fourier/steps_v3/dscan/dscan_%s.json" % which, "w"), indent=1, default=float)
