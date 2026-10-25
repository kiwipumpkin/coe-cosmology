"""Scan core mass Mc (1e54..1e61 kg), Earth distance D and shell density Os in the comoving-dust LTB (uniform bang, label = areal radius today,
local H_loc = 70). Shell edges hug the vision sphere as in the Newtonian fit: Rin = D - R, Rout = D + R (R = 14.30 Gpc; local observables don't depend on R).
Records the local shear s = H_r/H_perp - 1, age, boundness, 1+2E>0 on all shells."""
import json, numpy as np, multiprocessing as mp
import ltb_local as LL
R = 14.30048
def one(a):
    Mc, D, Os = a
    m = LL.Model(Mc, D, Os, max(D - R, 0.5 * D), D + R)
    try: r = m.solve()
    except Exception as e: r = None
    if r is None: return dict(Mc=Mc, D=D, Os=Os, ok=False)
    try: ok, mn = m.shells_ok(r["t0"], 25)
    except Exception: ok, mn = False, np.nan
    Mcrit = 3 * LL.H70**2 / (8 * np.pi) * 4 * np.pi / 3 * D**3 * LL.KG
    return dict(Mc=Mc, D=D, Os=Os, ok=True, shells_ok=ok, min_1p2E=mn, s=r["H_r"] / r["H_perp"] - 1, Mc_over_Mcrit=Mc / Mcrit,
                MD_over_Mcrit=m.M(D) * LL.KG / Mcrit, **r)
if __name__ == "__main__":
    grid = [(Mc, D, Os) for Mc in np.geomspace(1e54, 1e61, 15) for D in np.geomspace(20, 3000, 12) for Os in (1.012, 0.05)]
    with mp.Pool(8) as P: res = P.map(one, grid, chunksize=4)
    json.dump(res, open("scan_mc.json", "w"), default=float)
    print(sum(r["ok"] for r in res), "solved of", len(res))
