import numpy as np, json, sys, warnings
warnings.filterwarnings("ignore")
import sncore as s
from scipy.optimize import minimize_scalar
out = {}
for zmin in (0.023, 0.01):
    S = s.load(zmin); key = "z%s" % zmin; o = out[key] = {}
    Om = minimize_scalar(lambda x: s.chi2_marg(S, s.mu_frw(S["z"], x, 1 - x))[0], bounds=(0.01, 1), method="bounded").x
    refs = dict(LCDM=s.mu_frw(S["z"], Om, 1 - Om), Milne=s.mu_frw(S["z"], 0, 0))
    for name, mu in refs.items():
        o["dip_" + name] = s.residual_dipole(S, mu)
        for lo, hi in ((0, 0.1), (0.1, 0.3), (0.3, 10)):
            m = (S["z"] > lo) & (S["z"] <= hi); T = s.sub(S, m)
            o["dip_%s_%g_%g" % (name, lo, hi)] = s.residual_dipole(T, mu[m]); o["dip_%s_%g_%g" % (name, lo, hi)]["N"] = T["N"]
    for dn, (l, b) in s.DIRS.items():
        o["colin_" + dn] = s.colin_fit(S, l, b)
        o["hemi_" + dn] = [list(map(float, x)) for x in s.hemi_q0(S, l, b)]
    dl = o["dip_LCDM"]; o["colin_bestdip"] = s.colin_fit(S, dl["l"], dl["b"]); o["colin_bestdip"]["dir"] = [dl["l"], dl["b"]]
    o["Om_LCDM"] = Om
    print(key, json.dumps(o, default=float)[:3000], flush=True)
json.dump(out, open("direction.json", "w"), indent=1, default=float)
