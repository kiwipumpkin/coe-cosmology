import numpy as np, json, warnings, multiprocessing as mp
from scipy.optimize import minimize
warnings.filterwarnings("ignore")
import sncore as s
def run(Om):
    S = s.load(0.023); best = [np.inf, None]
    def f(x):
        LD, mc = 10 ** x[0], 10 ** x[1]
        if not (0.01 < LD < 10 and 1e-4 < mc < 1): return 1e6
        try:
            p, t0, ro = s.geom_from(41.53, LD, mc, Om); T = s.dl_table(p, t0, ro)
        except Exception: return 1e6
        c = s.chi2_marg(S, s.mu_ltb(T, s.dl_nodes(T, S["z"]), S["n"], x[2], x[3]))[0]; m = s.migkas_chi2(T, x[2], x[3])
        if c + m[0] < best[0]: best[:] = [c + m[0], dict(Om=Om, LD=LD, mc=mc, l=x[2] % 360, b=x[3], sn=c, mig=m[0], mig_r=m[1], t0=t0 * s.ltb.GYR)]
        return c + m[0]
    minimize(f, [np.log10(0.6), np.log10(0.024), 190, -20], method="Nelder-Mead",
             options=dict(maxfev=140, initial_simplex=None, xatol=0.01, fatol=0.05))
    return best[1]
if __name__ == "__main__":
    with mp.Pool(2) as p: r = p.map(run, [0.0, 0.05])
    json.dump(r, open("joint.json", "w"), indent=1); print(r)
