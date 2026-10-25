import numpy as np
from nobox_check import exact_ext, formula, K0
def cd(a, b): return abs((a-b+180) % 360-180)
SW = [("x", 1.01, 15, {}), ("alpha", 0.01, 8, {}), ("kappa", 0.001, 2, {}), ("eps", 0.005, 0.999, {"alpha": 1.3})]
for key, lo, hi, over in SW:
    rows = []
    for v in np.linspace(lo, hi, 60):
        kw = dict(x=5.9, eps=0.1, alpha=1.0, kappa=K0); kw.update(over); kw[key] = v
        ex, hbar, _ = exact_ext(**kw); f = formula(**kw)
        mx = [e for e in ex if e[0] == "max"]; mn = [e for e in ex if e[0] == "min"]
        ok = len(mx) == 2 and len(mn) == 2
        if not ok:
            rows.append((v, len(mx), len(mn), np.nan, np.nan)); continue
        dl = 0; dA = 0
        for (t, l, A) in ex:
            j = int(np.argmin([cd(l, lf) for lf in f["l"]]))
            dl = max(dl, cd(l, f["l"][j])); dA = max(dA, abs(A-f["A"][j]))
        rows.append((v, 2, 2, dl, dA))
    print(f"== {key}")
    for r in rows: print("  %.4g  nmax=%d nmin=%d  dl=%.3f  dA=%.4f" % r)
