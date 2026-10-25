"""Robustness check only: same objective as thin_nosine, wider search (600 random screens, 16 starts)."""
import json, numpy as np, thin_nosine as N, thin_shell as T
free, fixed, kin, kout, norm = T.JOBS["b|free|hug"]
old = json.load(open("thin_b_free_hug.json"))
d, runs = N.fit(free, fixed, kin, kout, norm, seeds=[dict(DR=old["DR"], Os=old["Os"], mcc=old["mcc"], l0=old["l0"])], starts=16, nscreen=600)
json.dump(dict(best={k: d[k] for k in ("DR", "Os", "mcc", "l0", "R2", "D_Gpc")}, starts=runs), open("thin_nosine_check.json", "w"), indent=1, default=float)
print("BEST", {k: round(d[k], 5) for k in ("DR", "Os", "mcc", "l0", "R2", "D_Gpc")})
