"""Free LTB fits to Pantheon+ (Lambda = 0), several starts in parallel."""
import numpy as np, json, sys, warnings, multiprocessing as mp
from scipy.optimize import minimize
warnings.filterwarnings("ignore")
import sncore as s
ZMIN = float(sys.argv[1]) if len(sys.argv) > 1 else 0.023
def evaluate(S, D, LD, mc, Om, fixdir=None):
    p, t0, ro = s.geom_from(D, LD, mc, Om)
    T = s.dl_table(p, t0, ro); DLn = s.dl_nodes(T, S["z"])
    l, b, c = s.best_direction(S, T, DLn, ngrid=400, fix=fixdir)
    return c, l, b, t0
def run(start):
    S = s.load(ZMIN); tag, x0, mode = start
    def unpack(x):
        D = 10 ** x[2] if mode == "Dfree" else 41.53
        Om = float(np.clip(x[2], 0, 0.9)) if mode == "Omfree" else 0.0
        return D, 10 ** x[0], 10 ** x[1], Om
    best = [np.inf]
    def f(x):
        D, LD, mc, Om = unpack(x)
        if not (0.01 < LD < 30 and 1e-5 < mc < 3 and 1 < D < 300): return 1e6
        try: c, l, b, t0 = evaluate(S, D, LD, mc, Om)
        except Exception: return 1e6
        if c < best[0]: best[:] = [c, dict(D=D, LD=LD, mc=mc, Om=Om, l=l, b=b, t0_Gyr=t0 * s.ltb.GYR)]
        return c
    r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=0.01, fatol=0.02, maxfev=90))
    out = dict(tag=tag, mode=mode, chi2=best[0], **best[1], nfev=r.nfev)
    print(json.dumps(out), flush=True); return out
if __name__ == "__main__":
    starts = [("A", [np.log10(0.07), np.log10(0.003)], "pure"), ("B", [np.log10(0.5), np.log10(0.03)], "pure"),
              ("C", [np.log10(1.0), np.log10(0.024)], "pure"), ("D", [np.log10(3.0), np.log10(0.05)], "pure"),
              ("E", [np.log10(0.2), np.log10(0.01)], "pure"), ("F", [np.log10(1.0), np.log10(0.024), 0.05], "Omfree"),
              ("G", [np.log10(0.07), np.log10(0.003), np.log10(41.53)], "Dfree"), ("H", [np.log10(10.0), np.log10(0.3)], "pure")]
    with mp.Pool(8) as pool: res = pool.map(run, starts)
    json.dump(res, open("fit_z%s.json" % ZMIN, "w"), indent=1)
