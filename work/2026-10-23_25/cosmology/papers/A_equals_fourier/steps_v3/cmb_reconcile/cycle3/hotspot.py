"""Part 2b: resolve the CoE hot spot in the regular two-mode models. N = 96 Gauss-Legendre rays for multipoles (vs 24), plus 96 rays
uniform in psi over 0.1-60 deg around the CoE; emission at fixed local age (tau_rec from the 96-ray monopole). Also the single-bang
relativistic headline (irregular: 1+2E < 0 inside r ~ 18 Gpc) with per-ray null-condition errors, to show which directions break."""
import json, numpy as np
from numpy.polynomial.legendre import leggauss
import mix3, cmbfix
D = mix3.D; out = {}
for name, (ms, Ab, w) in {"reg_ms0.7_w12": (0.7, 0.99999, 12.0), "reg_ms1.0_w13.5": (1.0, 0.99999, 13.5), "reg_ms1.0_w14_Ab0.999": (1.0, 0.999, 14.0),
                          "headline_single_bang": (1.0, 0.0, 1.0)}.items():
    p = mix3.build(ms, Ab, w); t0 = p[5]
    c = cmbfix.cmb_age(p, t0, D, N=96)
    tau = c["tau_rec_over_t0"] * t0
    ps = np.radians(np.linspace(0.1, 60, 96)); sols = [cmbfix.trace(p, t0, D, x, 1e3 * p[8] * (1 + cmbfix.ZREC) ** 3) for x in ps]
    zd = np.array([cmbfix.z_at_age(s, p, t0, tau) for s in sols]); ned = np.array([cmbfix.rays.null_err(p, t0, s) for s in sols])
    Tm = 0.5 * np.sum(leggauss(96)[1] / (1 + c["z"]))
    e = dict(ms=ms, Ab_frac=Ab, w=w, t0=t0, cb96=np.asarray(c.get("cb", [np.nan] * 9)).tolist(), beta96=c.get("beta"), nerr96_max=float(c["nerr"].max()),
             psi96=np.degrees(c["psi"]).tolist(), z96=c["z"].tolist(), nerr96=c["nerr"].tolist(),
             psi_dense=np.degrees(ps).tolist(), z_dense=zd.tolist(), dT_dense=((1 / (1 + zd)) / Tm - 1).tolist(), nerr_dense=ned.tolist())
    out[name] = e
    print(name, "a2-4 (96)", np.round(e["cb96"][2:5], 4), "beta", e["beta96"], "nerr96 %.1e" % e["nerr96_max"], "dense nerr max %.1e" % ned.max(),
          "dT/T at psi 0.1/5/20/40/60:", np.round(np.interp([0.1, 5, 20, 40, 60], np.degrees(ps), e["dT_dense"]), 3), flush=True)
json.dump(out, open("hotspot.json", "w"), indent=1, default=float)
