"""Exact check of the formula roots of d(ell)/d(input) = 0: evaluate the exact
no-boxcar-with-tilt extremum derivative at each formula root and locate the nearby
exact zero (Brent root-finding in the same variable). No fits, no sweeps."""
import json, sys, numpy as np
from scipy.optimize import brentq
from nobox_tilt_deriv_check import exact_vals, BASE, D, L, d, M
SYM = {"D": D, "L": L, "d": d, "M": M}
def dl_exact(name, X, v, rel=1e-3):
    h = rel*v[X]; vp = dict(v); vm = dict(v); vp[X] += h; vm[X] -= h
    return (exact_vals(vp)[name] - exact_vals(vm)[name])/(2*h)
def check(name, X, Y, y0, over=None):
    """d ell_name/dX (exact) at Y = y0, and the exact zero in Y nearest y0."""
    v = dict(BASE); v.update(over or {}); X, Y = SYM[X], SYM[Y]
    def g(y):
        u = dict(v); u[Y] = y
        try: return dl_exact(name, X, u)
        except AssertionError: return np.nan
    g0 = g(y0); out = {"at_root": g0}
    ks = [k for k in range(-25, 26) if 1 + 0.02*k > 0]
    ys = np.array([y0*(1 + 0.02*k) for k in ks]); gs = np.array([g(y) for y in ys])
    out["n_valid"] = int(np.isfinite(gs).sum()); out["y_valid"] = [float(ys[np.isfinite(gs)].min()), float(ys[np.isfinite(gs)].max())] if np.isfinite(gs).any() else None
    cands = [i for i in range(len(ys) - 1) if np.isfinite(gs[i]) and np.isfinite(gs[i+1]) and np.sign(gs[i]) != np.sign(gs[i+1])]
    root = None
    if cands:
        i = min(cands, key=lambda i: abs(ys[i] - y0))
        root = brentq(g, ys[i], ys[i+1], xtol=1e-7*y0)
        out["g_at_exact_root_neighbours"] = [float(gs[i]), float(gs[i+1])]
    out["exact_root"] = root
    return out
CASES = json.load(open("nobox_tilt_stat_roots.json"))["checks"]
if __name__ == "__main__":
    res = {}
    for key, c in CASES.items():
        over = {SYM[k]: val for k, val in c.get("over", {}).items()}
        r = check(c["name"], c["X"], c["Y"], c["y0"], over)
        res[key] = r; print(key, c, r, flush=True)
    json.dump(res, open("nobox_tilt_stat_exact.json", "w"), indent=1, default=float)
