"""Full-domain LTB regularity check (rm_core.validity only scans the shell interval): 1+2E > 0 (metric signature), R' > 0
(no shell crossing) at t0 and at earlier epochs, for r in [1e-3, 60] Gpc."""
import numpy as np
from ltb_common import LT
RS = np.r_[np.geomspace(1e-3, 1, 30), np.linspace(1.05, 60, 400)]
def full(p, t0, fracs=(1.0, 0.5, 0.2, 0.05, 0.01, 1e-3)):
    tB = lambda r: p[6] * np.exp(-(r / p[7]) ** 2) if p[7] > 0 else 0.0
    e = np.array([1 + 2 * LT.geom(r, t0, p)[5] for r in RS]); bad = RS[e <= 0]
    sc = []
    for f in fracs:
        t = t0 * f; rr = [r for r in RS if t > tB(r) * 1.0001]
        if not rr: continue
        rp = np.array([LT.geom(r, t, p)[1] for r in rr]); sc.append(int(np.sum(~(rp > 0))))
    return dict(ok=bool(len(bad) == 0 and sum(sc) == 0), min_1p2E=float(np.nanmin(e)), r_bad_max=float(bad.max()) if len(bad) else 0.0,
                n_shell_cross=int(sum(sc)))
