"""Stage 1: for each Mc, grid over D (log) and D/R; every point evaluated with rays; l0 fitted analytically."""
import json, sys, numpy as np, multiprocessing as mp, time
import rm_core as C
MCS = 2.472258e54 * 10 ** (np.arange(0, 14) / 2.0)
DS = np.geomspace(0.5, 400, 24); DRS = np.array([1.1, 1.3, 1.6, 2.0, 3.0, 4.5, 6.55, 10.0, 16.0])
def job(a):
    Mc, D, DR = a
    try: o = C.evaluate(Mc, D, D / DR)
    except Exception as e: o = dict(Mc_kg=Mc, D=D, R=D / DR, fail="exc %s" % e)
    o["DR"] = DR; return o
if __name__ == "__main__":
    t = time.time(); tasks = [(m, d, r) for m in MCS for d in DS for r in DRS]
    with mp.Pool(6) as P: res = P.map(job, tasks, chunksize=4)
    json.dump(res, open("grid_mc.json", "w"), default=float)
    print("done", len(res), round(time.time() - t))
