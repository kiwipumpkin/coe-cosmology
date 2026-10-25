import numpy as np, json, itertools, sys, warnings, multiprocessing as mp
warnings.filterwarnings("ignore")
import sncore as s
ZMIN = float(sys.argv[1]) if len(sys.argv) > 1 else 0.023
S = None
def work(args):
    global S
    if S is None: S = s.load(ZMIN)
    D, LD, mc, Om = args
    try:
        p, t0, ro = s.geom_from(D, LD, mc, Om)
        T = s.dl_table(p, t0, ro); DLn = s.dl_nodes(T, S["z"])
        l, b, c = s.best_direction(S, T, DLn, ngrid=600)
        return dict(D=D, LD=LD, mc=mc, Om=Om, l=l, b=b, chi2=c, t0_Gyr=t0 * s.ltb.GYR)
    except Exception as e:
        return dict(D=D, LD=LD, mc=mc, Om=Om, err=str(e)[:80])
if __name__ == "__main__":
    grid = list(itertools.product([10.0, 20.0, 41.53, 80.0], [0.07, 0.2, 0.5, 1.0], [0.003, 0.024, 0.1, 0.3], [0.0, 0.3]))
    with mp.Pool(8) as pool:
        res = []
        for r in pool.imap_unordered(work, grid):
            res.append(r); print(json.dumps(r), flush=True)
    json.dump(res, open("grid_z%s.json" % ZMIN, "w"), indent=1)
