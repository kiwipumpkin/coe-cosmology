"""Fast version of relfree.py: same ev(); per-evaluation time cap 15 s (slow far-D ray solves skipped), D 15-400 Gpc log-uniform,
Mc/Mcrit(D) 1e-5..1e3 log-uniform, kin,kout>=1 independent. Screen -> Nelder-Mead from the best 4 (R2 with sine). Incremental output."""
import json, signal, numpy as np, multiprocessing as mp
from scipy.optimize import minimize
import relfree as R
class TO(Exception): pass
def _h(*a): raise TO()
def evs(u):
    signal.signal(signal.SIGALRM, _h); signal.alarm(15)
    try: o = R.ev(*R.unpack(u))
    except TO: o = dict(fail="timeout")
    finally: signal.alarm(0)
    o["u"] = list(u); return o
def nm(u0):
    f = lambda u: -evs(u).get("R2_s", -5.0)
    r = minimize(f, np.array(u0), method="Nelder-Mead", options=dict(xatol=0.02, fatol=5e-4, maxiter=70))
    return evs(r.x)
if __name__ == "__main__":
    rng = np.random.default_rng(11); n = 560
    U = np.c_[rng.uniform(np.log10(15), np.log10(400), n), rng.uniform(-5, 3, n), rng.uniform(-4, 3, n), rng.uniform(-4, 3, n)]
    S = []
    with mp.Pool(7) as P:
        for i, o in enumerate(P.imap_unordered(evs, U, chunksize=2)):
            S.append(o)
            if i % 80 == 0: print("screen", i, flush=True)
        json.dump(dict(screen=S), open("relfree2.json", "w"), default=float)
        ok = sorted([o for o in S if "R2_s" in o], key=lambda o: -o["R2_s"])
        for o in ok[:6]: print(" scr D %.1f mcc %.3g Rin %.1f Rout %.1f R2s %.3f R2ns %.3f" % (o["D"], o["mcc"], o["R_in"], o["R_out"], o["R2_s"], o["R2_ns"]), flush=True)
        okn = sorted([o for o in S if "R2_ns" in o], key=lambda o: -o["R2_ns"])
        for o in okn[:3]: print(" scr-ns D %.1f mcc %.3g Rin %.1f Rout %.1f R2s %.3f R2ns %.3f" % (o["D"], o["mcc"], o["R_in"], o["R_out"], o["R2_s"], o["R2_ns"]), flush=True)
        F = P.map(nm, [o["u"] for o in ok[:4]] + [o["u"] for o in okn[:3]], chunksize=1)
    json.dump(dict(screen=S, refined=F), open("relfree2.json", "w"), default=float)
    for o in F:
        if "R2_s" in o: print(" NM D %.2f mcc %.4g Mc %.3e Rin %.2f Rout %.2f R2s %.4f R2ns %.4f l0 %.1f 2E %.2f age %.2f" % (o["D"], o["mcc"], o["Mc_kg"], o["R_in"], o["R_out"], o["R2_s"], o["R2_ns"], o["l0_s"], o["twoE"], o["age_Gyr"]), flush=True)
