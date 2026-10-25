"""Exact-LTB CMB for the full-mass Omega_s=0.049 fits (rays traced to t0*5e-5 only, which is enough for z=1089; analyze.cmb's
default 1e-6 t0 end hits the core singularity on the rays that pass it). Same emission rules: (a) fixed proper time with z monopole = 1089,
dipole removed by an axial boost (analyze.boost_multipoles). Writes A_full_<name>.json."""
import sys, json, time, numpy as np
from scipy.optimize import brentq
from ltb_common import *
F = fitparams(); k, cav, N = sys.argv[1], sys.argv[2] == "cav", int(sys.argv[3]); f = dict(F[k]); SC = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0; f["Mc"] *= SC; f["mcc"] *= SC
p = make_p(f["Mc"], f["D"], f["Os"], f["Rin"], 5000.0, cavity=cav); t0 = setup(p, f["D"])
mu, w = leggauss(N); T = time.time(); sols = []
for m in mu:
    s = rays.trace(p, t0, f["D"], np.arccos(m), t0 * 5e-5); sols.append(s); print("psi %.1f status %d %.0fs" % (np.degrees(np.arccos(m)), s.status, time.time() - T), flush=True)
zbar = lambda tr: 0.5 * np.sum(w * np.array([rays.z_at_t(s, t0, tr) for s in sols]))
tr = np.exp(brentq(lambda lt: np.log1p(zbar(np.exp(lt))) - np.log1p(1089.0), np.log(t0 * 6e-5), np.log(t0 * 0.01)))
z = np.array([rays.z_at_t(s, t0, tr) for s in sols]); Ta = 1 / (1 + z); Ta /= 0.5 * np.sum(w * Ta)
ca = rays.legendre_coeffs(mu, w, Ta - 1, 6); beta, cb = analyze.boost_multipoles(mu, w, Ta)
# emission at the fixed local density reached at tr on the far (anti-CoE) side  -> fixed-rho variant
out = dict(name=f"{k}|{'cavity' if cav else 'core-only'}", mcc=f["mcc"], D=f["D"], N=N, t_rec_over_t0=tr / t0, c_a=ca.tolist(), cb_a=cb.tolist(), beta=beta,
           z=z.tolist(), psi_deg=np.degrees(np.arccos(mu)).tolist(), status=[int(s.status) for s in sols], sec=time.time() - T)
print(out["name"], "a2..a4 boosted", cb[2:5], flush=True)
json.dump(out, open(f"A_full_{k}_{'cav' if cav else 'core'}_x{SC:g}.json", "w"), indent=1, default=float)
