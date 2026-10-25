"""Candidate 5: kinematic centre (Mc -> 0). Refit the Migkas A(l) profile with Mc/Mcrit fixed tiny.
Uses Chris's Newtonian integrator (ptshell/prof_core via newt_omb helpers) unchanged; A = |dV|/sep, mean-normalised, linear vs Migkas A/A_all,
with the same geometry + sine (Amp, phi) + l0 as the headline fits.
Cases:
  core-only kinematic: no core, no cavity -> pure V = k r -> A = const, only the sine term survives (= best single dipole).
  integrator, mcc fixed (2e-6 and 1e-9), shell rules as in the fits (Os=0.049 with free kin,kout; Os free grid with kin=kout=1).
    NOTE: these keep the shell's empty cavity (R_in = D - kin R), which the user excludes from the CMB model; any R2 above the
    dipole-only value comes from that cavity.
Writes c5_kinematic_fit.json."""
import sys, json, numpy as np
from scipy.optimize import minimize
BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"; sys.path.insert(0, BUB)
import os; os.chdir(BUB)
import newt_omb as N
from bubble_core import MA, IMIG, TH1, r2, A_from, shifted
ML = TH1[IMIG]

def geo(D, mcc, Os, kin, kout):
    return dict(DR=D / N.RV, Os=Os, mcc=mcc, rin=1 - kin * N.RV / D, rout=1 + kout * N.RV / D, Lc=0.005, mc=0.0, Ob=0.0)

def score(D, mcc, Os, kin, kout):
    if kin * N.RV >= 0.99 * D: return -5, None
    r = N.af_of(geo(D, mcc, Os, kin, kout))
    if r is None: return -5, None
    v, l0, am, ph = N.fit_l0(r[0], True); return v, dict(D=D, mcc=mcc, Os=Os, kin=kin, kout=kout, l0=l0, amp=am, phi=ph, R2=v)

def dipole_only():
    # best A = C0 + C1 cos l + S1 sin l on the 71 Migkas points (what a pure kinematic centre + sine can do)
    X = np.c_[np.ones_like(ML), np.cos(np.radians(ML)), np.sin(np.radians(ML))]
    c, *_ = np.linalg.lstsq(X, MA, rcond=None); out = {"dipole": r2(MA, X @ c)}
    for L in (2, 3):
        cols = [np.ones_like(ML)] + [f(np.radians(m * ML)) for m in range(1, L + 1) for f in (np.cos, np.sin)]
        X = np.c_[tuple(cols)]; c, *_ = np.linalg.lstsq(X, MA, rcond=None); out[f"harm<=%d" % L] = r2(MA, X @ c)
    X = np.c_[np.ones_like(ML), np.cos(np.radians(ML)), np.sin(np.radians(ML)), np.cos(np.radians(2 * ML)), np.sin(np.radians(2 * ML))]
    return out

def fit_omb(mcc, seed=0, n=120):
    rng = np.random.default_rng(seed); best = (-9, None)
    U = np.c_[rng.uniform(1.25, 3.5, n), rng.uniform(-4, 2.0, n), rng.uniform(-4, 2.0, n)]
    obj = lambda u: -score(10**u[0], mcc, 0.049, 1 + np.exp(u[1]), 1 + np.exp(u[2]))[0]
    vals = np.array([obj(u) for u in U])
    for i in np.argsort(vals)[:4]:
        r = minimize(obj, U[i], method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-5, maxiter=150))
        if -r.fun > best[0]: best = (-r.fun, r.x)
    u = best[1]; return score(10**u[0], mcc, 0.049, 1 + np.exp(u[1]), 1 + np.exp(u[2]))[1]

def fit_free(mcc):
    best = None
    for Os in (0.049, 0.3, 1.0, 3.0):
        for D in (20, 33, 49, 72, 94, 143, 286, 715):
            v, d = score(D, mcc, Os, 1.0, 1.0)
            if d and (best is None or v > best["R2"]): best = d
    return best

if __name__ == "__main__":
    out = {"dipole_only_and_harmonic_refs": dipole_only()}
    print(out, flush=True)
    for mcc in (2e-6, 1e-9):
        out[f"omb_rule|mcc={mcc:g}"] = fit_omb(mcc); print(mcc, out[f"omb_rule|mcc={mcc:g}"], flush=True)
        out[f"free_rule_grid|mcc={mcc:g}"] = fit_free(mcc); print(mcc, out[f"free_rule_grid|mcc={mcc:g}"], flush=True)
    # reference: headline omb parameters with mcc swept down (how R2 decays with core mass)
    sweep = []
    for mcc in (1.7e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 2e-6, 1e-9):
        v, d = score(49.1746, mcc, 0.049, 1.0528, 1.0815); sweep.append(d); print("sweep", mcc, v, flush=True)
    out["omb_headline_geometry_mcc_sweep"] = sweep
    json.dump(out, open("/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle1/c5_kinematic_fit.json", "w"), indent=1, default=float)
