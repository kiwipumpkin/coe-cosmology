"""Coarse free refit in comoving-dust LTB: for each (Os, D) find Mc giving the Migkas-optimal local shear, then R^2 (A~H^2, l0 refit)."""
import json, numpy as np, warnings
warnings.filterwarnings("ignore")
from scipy.optimize import brentq
import ltb_local as LL, quad_fit as Q
R = 14.30048
def model(Mc, D, Os): return LL.Model(Mc, D, Os, max(D - R, 0.5 * D), D + R)
def s_of(lMc, D, Os):
    r = model(10**lMc, D, Os).solve(); return np.nan if r is None else r["H_r"] / r["H_perp"] - 1
rows = []
L = np.arange(0, 180, 2.0)
for Os in (0.05, 0.3, 1.012):
    for D in (30, 60, 94, 150, 300, 800, 2000):
        Mcrit = 3 * LL.H70**2 / (8*np.pi) * 4*np.pi/3 * D**3 * LL.KG
        for target in (0.0813, -0.069):
            lg = np.linspace(np.log10(Mcrit) - 6, np.log10(Mcrit) + 0.5, 27); sv = np.array([s_of(x, D, Os) for x in lg])
            idx = [i for i in range(len(lg)-1) if np.isfinite(sv[i]) and np.isfinite(sv[i+1]) and (sv[i]-target)*(sv[i+1]-target) < 0]
            if not idx:
                rows.append(dict(Os=Os, D=D, target_s=target, found=False, s_range=[float(np.nanmin(sv)), float(np.nanmax(sv))])); continue
            lM = brentq(lambda x: s_of(x, D, Os) - target, lg[idx[0]], lg[idx[0]+1], xtol=1e-3)
            m = model(10**lM, D, Os); r = m.solve(); ok, mn = m.shells_ok(r["t0"], 25); s = r["H_r"]/r["H_perp"] - 1
            fits = [(Q.fit_s_l0(s, l0)[1], Q.fit_s_l0(s, l0)[0], l0) for l0 in L]; b = max(fits)
            rows.append(dict(Os=Os, D=D, target_s=target, found=True, Mc_kg=10**lM, Mc_over_Mcrit=10**lM/Mcrit, s=s, R2_sine=b[0], R2_geom=b[1], l0=b[2],
                             age_Gyr=r["age_Gyr"], bound=bool(r["bound"]), twoE=r["twoE"], shells_ok=ok, min_1p2E=mn, H_r=r["H_r"], H_perp=r["H_perp"]))
        print(rows[-2:], flush=True)
json.dump(rows, open("refit.json", "w"), indent=1, default=float)
