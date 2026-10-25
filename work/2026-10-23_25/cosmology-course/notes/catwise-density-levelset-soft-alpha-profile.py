"""Soft core+tail α profile at ℓ₀=15.2° (β=1): refit {D, mc} per α under δ lock."""
import sys, os, json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
code = Path('/workspace/cosmology-course/notes/catwise-density-levelset-migkas.py').read_text()
code = code.replace('if __name__ == "__main__":\n    main()', 'pass')
ns = {}; exec(compile(code, 'levelset', 'exec'), ns)
rows = []
for alpha in [0.31, 0.5, 0.7, 1.0, 1.5, 2.0]:
    sh = dict(alpha=alpha, beta=1.0)
    obj = lambda u: ns['chi2_of']('soft', 10**u[0], 10**u[1], 15.2, ns['D_LOCK'], sh)
    best = None
    for D0, mc0 in [(39.0, 4.4), (48.0, 4.5), (60.0, 4.6)]:
        r = minimize(obj, [np.log10(D0), np.log10(mc0)], method='Nelder-Mead', options=dict(maxiter=200, xatol=1e-4, fatol=1e-4))
        if best is None or r.fun < best.fun: best = r
    D, mc = 10**best.x[0], 10**best.x[1]
    full, rr, LD, L = ns['model_bins']('soft', D, mc, 15.2, ns['D_LOCK'], sh)
    fp = ns['_pack_fit']('soft', D, mc, 15.2, L, LD, rr, full, ns['D_LOCK'], sh, 3, 'fixed15.2')  # nu=3: D, mc, alpha
    c = ns['density_card']('soft', 58.2, ns['D_LOCK'], sh)
    row = dict(alpha=alpha, R2=fp['R2_unweighted'], chi2=fp['chi2'], sigma_nu3=fp['sigma_nested_one_sided'], D=D, L=L, mc=mc,
               delta_check_D58=c['delta_check'], cliff_D58=c['cliff_index'], smooth_at_centre=alpha >= 2.0 or alpha == 1.0 and False)
    row['central_slope'] = 'cusp (rho\' -> -inf)' if alpha < 1 else ('finite kink (rho\'(0)<0)' if alpha < 2 else 'smooth (rho\'(0)=0)')
    rows.append(row); print(row, flush=True)
Path('/workspace/cosmology-course/notes/catwise-density-levelset-soft-alpha-profile.json').write_text(json.dumps(rows, indent=1))
