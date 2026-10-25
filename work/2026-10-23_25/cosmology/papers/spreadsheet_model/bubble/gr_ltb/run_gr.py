"""GR (LTB) version of the thin shell + core best fit (universe_thin.json 'free'). Two constructions (see setup.py):
 GR-A: today's areal mass distribution = Newtonian exactly, uniform bang, t0 from local H = 70.
 GR-B: earlier convention (ltb_cmb case B): Newtonian profile shape on the Lagrangian label, E = r^2/(2 t0^2), uniform bang.
Nearby A(l) as in papers/cmb_redshift (analyze.nearby: D_L(z, psi) from ray bundles, A propto D_L^-2 at fixed z, b = 0),
compared with the Migkas points with the sine factor off and refitted (best_sine, |amp| <= 0.15), geometry and l0 fixed."""
import sys, os, json, time, warnings, numpy as np
warnings.filterwarnings("ignore")
import setup as S, ltb_thin as LT
sys.path.insert(0, os.path.join(S.HERE, "..")); sys.path.insert(0, os.path.join(S.HERE, "../../integrated"))
from bubble_core import MA, IMIG, TH1, r2, A_from, best_sine
from scipy.optimize import brentq
analyze = S.analyze
l0 = S.N["l0"]; analyze.L_COE = (l0 + 180.0) % 360       # CoE at l0 + 180 (l0 = anti-CoE, as in ltb_cmb)
K = 70 / S.H0
ZL = (0.05, 0.1, 0.2, 0.3)

def info(p, t0, ro, name):
    R, Rr, Rrr, Rt, Rtr, E, E1, M, M1 = LT.geom(ro, t0, p)
    Hp, Hr = Rt / R, Rtr / Rr
    d = dict(name=name, t0_Gpc=t0, age_Gyr=t0 * S.GYR, r_label=ro, R_obs=R, M_obs_kg=M * S.KG_PER_GPC, rs_enc=2 * M, twoM_over_R=2 * M / R,
             twoE=2 * E, bound=bool(E < 0), Rdot=Rt, H_perp=Hp * K, H_r=Hr * K, H_loc=(Hr + 2 * Hp) / 3 * K,
             shear=(Hr - Hp) / ((Hr + 2 * Hp) / 3), rho_obs_over_crit=LT.density(ro, t0, p) / S.RHOC)
    if E < 0:
        d["R_max"] = M / (-E); d["t_turn_Gyr"] = np.pi * M / (-2 * E) ** 1.5 * S.GYR; d["Rmax_over_rs"] = d["R_max"] / (2 * M)
    # matter shells: sign/limit of 2E, shell crossing (R' > 0) for 1e-4 t0 < t <= 3 t0 (past and future), edges today
    a, b = p[9], p[10]
    rr = np.linspace(a - p[11], b + p[11], 241)
    twoE = np.array([2 * LT.geom(r, t0, p)[5] for r in rr])
    d["twoE_min_matter"] = float(twoE.min()); d["twoE_max_matter"] = float(twoE.max()); d["any_1p2E_neg"] = bool((1 + twoE <= 0).any())
    minRr, where = np.inf, None
    for r in rr[::4]:
        for t in t0 * np.r_[np.geomspace(1e-4, 1, 80), np.linspace(1.02, 3, 40)]:
            g = LT.geom(r, t, p)
            if np.isfinite(g[1]) and g[1] < minRr: minRr, where = g[1], (float(r), float(t / t0))
    d["min_Rprime"] = float(minRr); d["min_Rprime_at_r_t/t0"] = where
    d["R_inner_edge_today"] = float(LT.geom(a, t0, p)[0]); d["R_outer_edge_today"] = float(LT.geom(b, t0, p)[0])
    Hr_rows = []
    for r in np.linspace(a + 2 * p[11], b - 2 * p[11], 7):
        Hp_, Hr_ = LT.local_rates(r, t0, p); Hr_rows.append((float(LT.geom(r, t0, p)[0]), float(Hp_ * K), float(Hr_ * K)))
    d["H_perp_H_r_across_shell"] = Hr_rows
    return d

def fits(ell, A):
    """A on 1-deg grid in Galactic longitude, mean-normalised; th = l - l0 as in the Newtonian fit."""
    A = A / A.mean(); x = A[IMIG]; th = (TH1 - l0)[IMIG]
    amp, ph = best_sine(x, th, MA)
    return dict(R2_nosine=float(r2(MA, x)), R2_sine=float(r2(MA, A_from(x, th, amp, ph))), amp=float(amp), phi=float(ph),
                four=analyze.fourier(ell, A).tolist(), corr=float(np.corrcoef(x, MA)[0, 1]), ptp=float(A.max() - A.min()))

if __name__ == "__main__":
    out = {}
    for name, fn in (("GR-A", S.model_areal), ("GR-B", S.model_lagr)):
        t1 = time.time(); p, t0, ro = fn()
        d = info(p, t0, ro, name)
        nb_ = analyze.nearby(p, t0, ro, zlist=ZL, N=24)
        A = 1.0 / nb_["A"]                               # A propto D_L^-2 (H0^2) at fixed z
        d["fits"] = {str(z): fits(nb_["ell"], A[:, k]) for k, z in enumerate(ZL)}
        d["p"] = p.tolist(); d["secs"] = time.time() - t1
        np.save("A_%s.npy" % name, A)
        out[name] = d
        print(json.dumps({k: v for k, v in d.items() if k != "p"}, indent=1, default=float), flush=True)
    json.dump(out, open("gr_results.json", "w"), indent=1, default=float)
