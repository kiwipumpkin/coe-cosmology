"""Ballpark: Omega_s free + Mc (as Mc/Mcrit(D)) + D, D/R = 6.55 kept, l0 (and Amp, phi) refit; ray-traced LTB GR-A, H = 70, t_B = 0."""
import json, numpy as np, multiprocessing as mp
import rm_core as C
OSS = [0.02, 0.05, 0.2, 1.012]; MU = np.geomspace(0.005, 1.0, 10); DS = [93.72, 500.0]
def one(a):
    Os, mu, D = a; Mcrit = 4 * np.pi / 3 * C.RHOC * D**3 * C.KG
    try: o = C.evaluate(mu * Mcrit, D, D / 6.5536, Os=Os)
    except Exception as e: o = dict(fail=str(e))
    o.update(Os=Os, mu=mu, D=D); return o
if __name__ == "__main__":
    T = [(o, m, d) for o in OSS for m in MU for d in DS]
    with mp.Pool(7) as P: R = P.map(one, T, chunksize=2)
    json.dump(R, open("os_grid.json", "w"), default=float)
    V = [o for o in R if "R2_s" in o and o["twoE_min"] > -1]
    for Os in OSS:
        for D in DS:
            L = [o for o in V if o["Os"] == Os and o["D"] == D]
            if not L: print("Os %.3f D %.0f: no valid" % (Os, D)); continue
            b = max(L, key=lambda o: o["R2_s"])
            print("Os %.3f D %6.1f best: mu %.3f Mc %.2e R2s %.3f R2ns %.3f l0s %.1f along %s 2E %.1f age %.2f Hp %.1f Hr %.1f m2 %.3f | valid n=%d, along n=%d" % (
                Os, D, b["mu"], b["Mc_kg"], b["R2_s"], b["R2_ns"], b["fits"]["0.1"]["l0_s"], b["along"], b["twoE"], b["age_Gyr"], b["H_perp"], b["H_r"], b["fits"]["0.1"]["m2"], len(L), sum(o["along"] for o in L)))
