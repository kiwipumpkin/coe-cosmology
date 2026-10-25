import setup as S, ltb_thin as LT, numpy as np, json
rays = S.rays
out = {}
for name, fn in (("GR-A", S.model_areal), ("GR-B", S.model_lagr)):
    p, t0, ro = fn()
    rr = np.linspace(p[9] - p[11], p[10] + p[11], 400)
    ts = t0 * np.r_[np.geomspace(1e-5, 1, 200), np.linspace(1.0, 5, 400)]
    first = None
    for t in ts:
        m = min(LT.geom(r, t, p)[1] for r in rr)
        if m <= 0: first = t / t0; break
    past_min = min(LT.geom(r, t, p)[1] for r in rr[::2] for t in t0 * np.geomspace(1e-5, 1, 120))
    ne = []
    for psi in (0.1, 1.0, 2.0, 3.0):
        sol = rays.trace(p, t0, ro, psi, 0.6 * t0); ne.append(rays.null_err(p, t0, sol))
    out[name] = dict(first_shell_crossing_t_over_t0=first, first_crossing_after_today_Gyr=None if first is None else (first - 1) * t0 * S.GYR,
                     min_Rprime_past=past_min, null_err_max=max(ne))
    print(name, out[name])
json.dump(out, open("checks.json", "w"), indent=1, default=float)
