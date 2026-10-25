"""Follow-ups: (a) where the exact second minimum disappears for the cases whose formula
root lies outside it; (b) wider local search for lmin2 dD-root in M."""
import json, numpy as np
from scipy.optimize import brentq
from nobox_tilt_deriv_check import exact_vals, BASE, D, L, d, M
from nobox_tilt_stat_exact import dl_exact, SYM
def exists(v):
    try: exact_vals(v); return True
    except AssertionError: return False
def edge(Y, a, b, over=None):
    """bisection for the boundary of existence in Y between a (exists) and b (not)."""
    v = dict(BASE); v.update(over or {})
    for _ in range(40):
        m = 0.5*(a + b); u = dict(v); u[Y] = m
        if exists(u): a = m
        else: b = m
    return a
out = {}
out["Mc_d0.1"] = edge(M, 0.045, 0.0393)
out["Mc_dline"] = edge(M, 0.03, 0.025, {d: 4/3*np.tan(0.35)})
out["Dc"] = edge(D, 2.39, 2.5)
out["Lc"] = edge(L, 0.4172, 0.45)
for k, Y, over in (("Mc_d0.1", M, {}), ("Mc_dline", M, {d: 4/3*np.tan(0.35)}), ("Dc", D, {}), ("Lc", L, {})):
    for eps in (1e-2, 3e-2):
        v = dict(BASE); v.update(over); v[Y] = out[k]*(1 + (eps if Y in (M,) else -eps))
        try: out[k + "_dldd_near_%g" % eps] = dl_exact("lmin2", d, v, rel=1e-4)
        except AssertionError: out[k + "_dldd_near_%g" % eps] = None
        out[k + "_lmin2_near_%g" % eps] = exact_vals(v)["lmin2"]
# (b) lmin2 dD root in M: wider bracket
def g(m):
    v = dict(BASE); v[M] = m
    try: return dl_exact("lmin2", D, v)
    except AssertionError: return np.nan
ms = np.geomspace(0.04, 0.95, 40); gs = np.array([g(m) for m in ms])
out["lmin2_D_M_scan_signs"] = [(float(m), float(x)) for m, x in zip(ms, gs)]
c = [i for i in range(len(ms) - 1) if np.isfinite(gs[i]) and np.isfinite(gs[i+1]) and np.sign(gs[i]) != np.sign(gs[i+1])]
out["lmin2_D_M_exact_roots"] = [brentq(g, ms[i], ms[i+1], xtol=1e-8) for i in c]
print(json.dumps(out, indent=1, default=float))
json.dump(out, open("nobox_tilt_stat_exact2.json", "w"), indent=1, default=float)
