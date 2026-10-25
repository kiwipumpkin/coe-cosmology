"""Chris-geometry variants, Migkas consistency, joint SN+Migkas grid."""
import numpy as np, json, sys, warnings, multiprocessing as mp
from scipy.optimize import minimize_scalar
warnings.filterwarnings("ignore")
import sncore as s
ZMIN = float(sys.argv[1]) if len(sys.argv) > 1 else 0.023
LDc, MCc = 2.919 / 41.53, 0.02395
COE = s.DIRS["CoE"]
def table(D, LD, mc, Om):
    p, t0, ro = s.geom_from(D, LD, mc, Om); return s.dl_table(p, t0, ro), t0 * s.ltb.GYR
def variant(args):
    name, Om, fix = args; S = s.load(ZMIN)
    if Om == "free":
        best = {}
        def f(om):
            T, t0 = table(41.53, LDc, MCc, om); DLn = s.dl_nodes(T, S["z"])
            l, b, c = s.best_direction(S, T, DLn, ngrid=400, fix=COE if fix else None)
            best[om] = (c, l, b, t0); return c
        minimize_scalar(f, bounds=(0.0, 0.9), method="bounded", options=dict(xatol=0.01))
        om = min(best, key=lambda k: best[k][0]); c, l, b, t0 = best[om]
    else:
        om = Om; T, t0 = table(41.53, LDc, MCc, om); DLn = s.dl_nodes(T, S["z"])
        l, b, c = s.best_direction(S, T, DLn, ngrid=1200, fix=COE if fix else None)
    T, _ = table(41.53, LDc, MCc, om)
    mig = s.migkas_chi2(T, l, b)
    return dict(name=name, Om=om, fixdir=fix, chi2=c, l=l, b=b, t0_Gyr=t0, mig_chi2=mig[0], mig_r=mig[1])
def gridpt(args):
    LD, mc = args; S = s.load(ZMIN)
    T, t0 = table(41.53, LD, mc, 0.0); DLn = s.dl_nodes(T, S["z"])
    ls = np.arange(0, 360, 5.0); bs = np.arange(-60, 61, 10.0)
    sn = np.array([[s.chi2_marg(S, s.mu_ltb(T, DLn, S["n"], l, b))[0] for b in bs] for l in ls])
    mg = np.array([[s.migkas_chi2(T, l, b)[0] for b in bs] for l in ls])
    i = np.unravel_index(np.argmin(sn), sn.shape); j = np.unravel_index(np.argmin(sn + mg), sn.shape)
    return dict(LD=LD, mc=mc, t0=t0, sn_best=float(sn[i]), sn_dir=[ls[i[0]], bs[i[1]]], mig_at_snbest=float(mg[i]),
                joint=float((sn + mg)[j]), joint_sn=float(sn[j]), joint_mig=float(mg[j]), joint_dir=[ls[j[0]], bs[j[1]]],
                mig_min=float(mg.min()))
if __name__ == "__main__":
    V = [("Chris Om0 CoE", 0.0, True), ("Chris Om0 free", 0.0, False), ("Chris Om0.05 CoE", 0.05, True), ("Chris Om0.05 free", 0.05, False), ("Chris Om0.3 CoE", 0.3, True),
         ("Chris Om0.3 free", 0.3, False), ("Chris Omfree CoE", "free", True), ("Chris Omfree free", "free", False)]
    G = [(ld, mc) for ld in [0.03, 0.07, 0.15, 0.3, 0.6, 1.2] for mc in [0.001, 0.003, 0.01, 0.024, 0.06, 0.15]]
    with mp.Pool(8) as pool:
        v = pool.map(variant, V); print(json.dumps(v), flush=True)
        g = pool.map(gridpt, G)
    json.dump(dict(variants=v, grid=g), open("post_z%s.json" % ZMIN, "w"), indent=1)
