"""B: Chris's geometry + a tilt (matter bulk velocity relative to the CMB frame) carrying the dipole/sine part.
The tilt enters A(l) as the multiplicative dipole factor [1 - Amp sin(psi + phi)] (same form as the headline sine term; a bulk flow v
changes the apparent rate by v cos(angle)/(c z) -> a dipole in l). Refit with Mc/Mcrit fixed on a ladder: D, kin, kout, l0, Amp, phi free
(Omega_s = 0.049 rule). Also: the same with NO cavity (uniform dust from 0.01 D outward + core), Mc free.
Uses Chris's Newtonian integrator unchanged (via newt_omb helpers). Writes B_hybrid_fit.json."""
import sys, os, json, numpy as np
from scipy.optimize import minimize
BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"; sys.path.insert(0, BUB); os.chdir(BUB)
import newt_omb as N
OUTF = "/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle2/B_hybrid_fit.json"
def geo(D, mcc, Os, rin, rout): return dict(DR=D / N.RV, Os=Os, mcc=mcc, rin=rin, rout=rout, Lc=0.005, mc=0.0, Ob=0.0)
def score(D, mcc, kin, kout, cavity=True):
    if cavity:
        if kin * N.RV >= 0.99 * D: return -5, None
        rin, rout = 1 - kin * N.RV / D, 1 + kout * N.RV / D
    else: rin, rout = 0.03, 1 + kout * N.RV / D
    r = N.af_of(geo(D, mcc, 0.049, rin, rout))
    if r is None: return -5, None
    v, l0, am, ph = N.fit_l0(r[0], True); return v, dict(D=D, mcc=mcc, kin=kin if cavity else None, kout=kout, l0=l0, amp=am, phi=ph, R2=v, cavity=cavity)
def fit(mcc, cavity=True, seed=1, n=80):
    rng = np.random.default_rng(seed)
    if mcc is None:   # mcc free
        U = np.c_[rng.uniform(1.25, 3.5, n), rng.uniform(-4, 2, n), rng.uniform(-4, 2, n), rng.uniform(-5, 0, n)]
        obj = lambda u: -score(10**u[0], 10**u[3], 1 + np.exp(u[1]), 1 + np.exp(u[2]), cavity)[0]
    else:
        U = np.c_[rng.uniform(1.25, 3.5, n), rng.uniform(-4, 2, n), rng.uniform(-4, 2, n)]
        obj = lambda u: -score(10**u[0], mcc, 1 + np.exp(u[1]), 1 + np.exp(u[2]), cavity)[0]
    vals = np.array([obj(u) for u in U]); best = (-9, None)
    for i in np.argsort(vals)[:3]:
        r = minimize(obj, U[i], method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-5, maxiter=150))
        if -r.fun > best[0]: best = (-r.fun, r.x)
    u = best[1]; m = mcc if mcc is not None else 10**u[3]
    return score(10**u[0], m, 1 + np.exp(u[1]), 1 + np.exp(u[2]), cavity)[1]
if __name__ == "__main__":
    out = {"ladder_cavity": [], "no_cavity": None}
    for mcc in (1.7e-3, 1e-3, 5e-4, 2e-4, 1e-4, 2e-5, 2e-6):
        d = fit(mcc); out["ladder_cavity"].append(d); print(d, flush=True)
    out["no_cavity"] = fit(None, cavity=False); print("no cavity", out["no_cavity"], flush=True)
    out["no_cavity_ladder"] = []
    for mcc in (1.7e-3, 5e-4, 1e-4, 2e-6):
        d = fit(mcc, cavity=False); out["no_cavity_ladder"].append(d); print("nocav", d, flush=True)
    json.dump(out, open(OUTF, "w"), indent=1, default=float)
