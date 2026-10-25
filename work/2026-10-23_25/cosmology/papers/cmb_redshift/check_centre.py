"""Centred-observer check: in model B, dipole and quadrupole must vanish as r_obs -> 0 (dipole ~ r_obs)."""
import numpy as np, json, models as m, analyze as an, ltb, rays, warnings
warnings.filterwarnings("ignore")
from numpy.polynomial.legendre import leggauss
p, t0, ro = m.model_B()
mu, wts = leggauss(16); out = {}
for r in (1e-3, 1e-2, 1e-1):
    tr = t0 * 1e-4
    z = np.array([rays.z_at_t(rays.trace(p, t0, r, np.arccos(x), tr * 0.99), t0, tr) for x in mu])
    T = 1 / (1 + z); T /= 0.5 * np.sum(wts * T)
    c = rays.legendre_coeffs(mu, wts, T - 1, 3)
    out[str(r)] = dict(c1=c[1], c2=c[2], c3=c[3])
print(json.dumps(out, indent=1)); json.dump(out, open("check_centre.json", "w"), indent=1)
