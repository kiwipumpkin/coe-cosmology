import numpy as np, json, ltb, rays
from scipy.optimize import brentq
H0 = ltb.H0_70; Om = 0.3
p = np.zeros(8); p[0] = 0; p[1] = Om * H0**2 / 2; p[3] = 0; p[4] = (1 - Om) * H0**2 / 2
# age: R(r,t0) = r
f = lambda t0: ltb.geom(1.0, t0, p)[0] - 1.0
t0 = brentq(f, 0.5, 10, xtol=1e-15); p[5] = t0
Hp, Hr = ltb.local_rates(3.0, t0, p)
print("t0 [Gyr]", t0 * ltb.GYR, "H/H0", Hp / H0, Hr / H0)
robs = 3.0; mu, wts = np.polynomial.legendre.leggauss(16); res = {}
trec = t0 * 2e-5
zs = []; errs = []; zth = []
for m in mu:
    s = rays.trace(p, t0, robs, np.arccos(m), trec)
    zs.append(rays.z_at_t(s, t0, trec)); errs.append(rays.null_err(p, t0, s))
    # FRW relation along ray at mid time
    tm = t0 * 0.3; zth.append(rays.z_at_t(s, t0, tm))
zs = np.array(zs); a_rec = ltb.geom(1.0, trec, p)[0]; ztrue = 1 / a_rec - 1
zm = 1 / ltb.geom(1.0, t0 * 0.3, p)[0] - 1
c = rays.legendre_coeffs(mu, wts, (1 + zs.mean()) / (1 + zs) - 1)
out = dict(t0_Gyr=t0 * ltb.GYR, z_rec_true=ztrue, z_rec_rays_max_rel_err=float(np.max(abs(zs - ztrue)) / (1 + ztrue)),
           z_mid_max_rel_err=float(np.max(abs(np.array(zth) - zm)) / (1 + zm)), null_err=float(max(errs)),
           multipoles=c[:4].tolist())
# D_A check at z~0.1..1 against FRW open formula
zq, DA = rays.area_distance(p, t0, robs, 1.0, t0 * np.array([0.9, 0.7, 0.5]))
from scipy.integrate import quad
def DA_frw(z):
    Ok = 1 - Om; chi = quad(lambda x: 1 / np.sqrt(Om * (1 + x)**3 + Ok * (1 + x)**2), 0, z, epsabs=1e-13)[0]
    return np.sinh(np.sqrt(Ok) * chi) / np.sqrt(Ok) / H0 / (1 + z)
out["DA_rel_err"] = [float(DA[i] / DA_frw(zq[i]) - 1) for i in range(3)]
print(json.dumps(out, indent=1)); json.dump(out, open("check_frw.json", "w"), indent=1)
