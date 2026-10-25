"""No-ray search: is there any (Mc, D, R) with Omega_s = 1.012, H = 70, single bang, where H_r > H_perp at Earth (peaks along the
CoE axis) AND 1 + 2E > 0 on every matter shell?  Dense log-random sample + summary; also threshold Mc(flip) at the Newtonian D, R."""
import json, numpy as np, multiprocessing as mp
import rm_core as C, ltb_fill as LT
def one(a):
    lMc, lD, lq = a; Mc = 10**lMc; D = 10**lD; R = D / (1 + 10**lq)
    p = C.params(Mc, D, R); t0 = C.solve_t0(p, D)
    if t0 is None: return None
    g = LT.geom(D, t0, p); v = C.validity(p, t0)
    return dict(Mc=Mc, D=D, R=R, age=t0 * C.GYR, Hp=g[3] / g[0] * C.K, Hr=g[4] / g[1] * C.K, twoE=2 * g[5], twoE_min=v["twoE_min"], ok=v["twoE_min"] > -1)
if __name__ == "__main__":
    rng = np.random.default_rng(3); n = 6000
    T = np.c_[rng.uniform(50, 62, n), rng.uniform(-0.5, 3.0, n), rng.uniform(-1.5, 1.5, n)]
    with mp.Pool(6) as P: res = [r for r in P.map(one, T, chunksize=20) if r]
    ok = [r for r in res if r["ok"]]; al = [r for r in res if r["Hr"] > r["Hp"]]; both = [r for r in ok if r["Hr"] > r["Hp"]]
    print("solved", len(res), "valid", len(ok), "along", len(al), "along&valid", len(both))
    if al: print("along: max twoE_min", max(r["twoE_min"] for r in al), min(al, key=lambda r: -r["twoE_min"]))
    if ok: print("valid: max (Hr-Hp)", max(ok, key=lambda r: r["Hr"] - r["Hp"]))
    # threshold at Newtonian D, R
    thr = []
    for lMc in np.linspace(54, 58, 81):
        r = one((lMc, np.log10(93.72014521956349), np.log10(93.72014521956349 / 14.30048287260589 - 1)))
        if r: thr.append(r)
    json.dump(dict(sample=res, thr=thr), open("flipsearch.json", "w"), default=float)
    for r in thr[::4]: print("%.3e  Hp %.2f Hr %.2f 2E %.2f 2Emin %.2f age %.2f" % (r["Mc"], r["Hp"], r["Hr"], r["twoE"], r["twoE_min"], r["age"]))
