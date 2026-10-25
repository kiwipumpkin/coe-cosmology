import numpy as np, json, warnings
from scipy.optimize import minimize_scalar
warnings.filterwarnings("ignore")
import sncore as s
out = {}
l, A, e = s.migkas_data(); w = 1 / e**2; k = np.sum(w * A) / w.sum()
out["mig_iso_chi2"] = float(np.sum(w * (A - k)**2)); out["mig_N"] = len(l)
for zmin in [0.023, 0.01]:
    S = s.load(zmin); o = {"N": S["N"]}
    r = minimize_scalar(lambda om: s.chi2_marg(S, s.mu_frw(S["z"], om, 1 - om))[0], bounds=(0.05, 0.6), method="bounded")
    o["LCDM"] = [float(r.fun), float(r.x)]
    for nm, om in [("Milne", 0.0), ("FRWb", 0.05), ("EdS", 1.0)]: o[nm] = float(s.chi2_marg(S, s.mu_frw(S["z"], om, 0.0))[0])
    for nm, om in [("ChrisCoE", 0.0), ("ChrisCoE_b", 0.05), ("ChrisCoE_03", 0.3)]:
        p, t0, ro = s.geom_from(41.53, 2.919 / 41.53, 0.02395, om); T = s.dl_table(p, t0, ro)
        o[nm] = float(s.chi2_marg(S, s.mu_ltb(T, s.dl_nodes(T, S["z"]), S["n"], *s.DIRS["CoE"]))[0])
    fit = min(json.load(open("fit_z%s.json" % zmin)), key=lambda r: r["chi2"])
    p, t0, ro = s.geom_from(fit["D"], fit["LD"], fit["mc"], fit["Om"]); T = s.dl_table(p, t0, ro)
    mu = s.mu_ltb(T, s.dl_nodes(T, S["z"]), S["n"], fit["l"], fit["b"])
    o["best"] = fit; o["best_dip"] = {k: v for k, v in s.residual_dipole(S, mu).items() if k != "vec"}
    o["best_mig"] = s.migkas_chi2(T, fit["l"], fit["b"])
    # quadrupole-to-mean of best model at z=0.1
    out["z%s" % zmin] = o
json.dump(out, open("extra.json", "w"), indent=1, default=float); print(json.dumps(out, default=float, indent=0))
