#!/usr/bin/env python3
"""Term-by-term equality between the m=0..3 harmonic fit of the 71 digitised Migkas+20 Fig.6
points and the paper's analytic tidal equation (main.tex eq:tidal):
  A/A_all = N[1 + ak(3mu^2-1) + ak(d/D)(9/2 mu - 15/2 mu^3)],  mu = cos b0 cos(l-l0),  ak = alpha*kappa
in-plane, 90-deg sector window W_m = sin(m pi/4)/(m pi/4):
  m0: N[1+ak(3c^2/2-1)]   m1: N ak eps c(9/2-45c^2/8) W1   m2: N ak (3c^2/2) W2   m3: -N ak eps (15/8) c^3 W3
Also compares b0 = 0/15/33.5 deg in the full numerical ideal-drag model.
Outputs dq_model_equality.{png,pdf,json}."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dq_model_lib import ell_m, A_m, A_lo, A_hi, MS0, EPS0, X0, L0DEG, kbase, model_at_data, model_fine, harm_fine, lsq, amp_phase, r2, nelder_mead

OUT = "/workspace/cosmology/figures"
R = {}
# ---- data m=0..3 fit
b, X, pred = lsq(ell_m, A_m, (1, 2, 3))
c = {"c0": b[0], "c1": b[1], "s1": b[2], "c2": b[3], "s2": b[4], "c3": b[5], "s3": b[6]}
R["data_fit_m0_3"] = {"coef": {k: float(v) for k, v in c.items()}, "R2": r2(A_m, pred),
                      "amp_phase": {f"m{m}": amp_phase(c[f"c{m}"], c[f"s{m}"], m) for m in (1, 2, 3)}}
def proj(l0):
    P = {0: c["c0"]}; Q = {}
    for m in (1, 2, 3):
        t = np.radians(m * l0)
        P[m] = c[f"c{m}"] * np.cos(t) + c[f"s{m}"] * np.sin(t)
        Q[m] = -c[f"c{m}"] * np.sin(t) + c[f"s{m}"] * np.cos(t)
    return P, Q

W_cont = {m: np.sin(m * np.pi / 4) / (m * np.pi / 4) for m in (1, 2, 3)}
W_disc = {1: 0.8927632947856, 2: 0.6119081150141242, 3: 0.2615587166275956}  # 47-pt boxcar used in the numeric code

def ratio_solve(P, W):
    r13 = P[1] / P[3]            # = -[W1 (9/2 - 45y/8)] / [(15/8) y W3]
    y = (4.5 * W[1]) / (45 / 8 * W[1] - r13 * 15 / 8 * W[3])
    cc = np.sqrt(y)
    r32 = P[3] / P[2]            # = -(15/8) c eps W3 / ((3/2) W2)
    eps = -r32 * 1.5 * W[2] / (15 / 8 * cc * W[3])
    Nak = P[2] / (1.5 * y * W[2])
    # m0: N + Nak*(3y/2 - 1) = P0
    N = P[0] - Nak * (1.5 * y - 1)
    return {"c2": float(y), "b0_deg": float(np.degrees(np.arccos(cc))), "d_over_D": float(eps),
            "alpha_kappa": float(Nak / N), "N": float(N), "N_alpha_kappa": float(Nak)}

def coeffs(N, ak, eps, cc, W):
    return {0: N * (1 + ak * (1.5 * cc ** 2 - 1)), 1: N * ak * eps * cc * (4.5 - 45 * cc ** 2 / 8) * W[1],
            2: N * ak * 1.5 * cc ** 2 * W[2], 3: -N * ak * eps * 15 / 8 * cc ** 3 * W[3]}
def curve(ell, N, ak, eps, cc, l0, W):
    k = coeffs(N, ak, eps, cc, W); L = np.radians(ell - l0)
    return k[0] + k[1] * np.cos(L) + k[2] * np.cos(2 * L) + k[3] * np.cos(3 * L)

P15, Q15 = proj(15.0)
R["projections_l0_15"] = {"parallel": {f"m{m}": float(P15[m]) for m in P15},
                          "orthogonal_unmodelled": {f"m{m}": float(Q15[m]) for m in Q15}}
sol = ratio_solve(P15, W_cont); sol_d = ratio_solve(P15, W_disc)
R["ratio_solve_l0_15_Wcont"] = sol; R["ratio_solve_l0_15_Wdisc47pt"] = sol_d
cc = np.sqrt(sol["c2"])
mod15 = curve(ell_m, sol["N"], sol["alpha_kappa"], sol["d_over_D"], cc, 15.0, W_cont)
R["ratio_solve_l0_15_Wcont"]["R2_vs_71pts"] = r2(A_m, mod15)
R["ratio_solve_l0_15_Wcont"]["model_coeffs"] = {f"m{m}": float(v) for m, v in coeffs(sol["N"], sol["alpha_kappa"], sol["d_over_D"], cc, W_cont).items()}
# ratio solve as a function of l0 (sensitivity)
R["ratio_solve_vs_l0"] = {}
for l0 in (5.0, 10.0, 14.4, 15.0, 20.0, 25.0):
    P, Q = proj(l0); s = ratio_solve(P, W_cont)
    s["R2"] = r2(A_m, curve(ell_m, s["N"], s["alpha_kappa"], s["d_over_D"], np.sqrt(max(s["c2"], 0)), l0, W_cont)) if 0 < s["c2"] <= 1 else None
    R["ratio_solve_vs_l0"][str(l0)] = s
# free LSQ of the analytic model to the 71 points
def neg(p):
    N, ak, eps, cb, l0 = p
    if not (0 < cb <= 1): return 1e3
    return float(((A_m - curve(ell_m, N, ak, eps, cb, l0, W_cont)) ** 2).sum())
pf, _ = nelder_mead(neg, [sol["N"], sol["alpha_kappa"], sol["d_over_D"], cc, 15.0], [0.01, 0.02, 0.1, 0.05, 3.0], tol=1e-14, maxit=5000)
Nf, akf, ef, cf, lf = pf
R["free_LSQ_analytic"] = {"N": Nf, "alpha_kappa": akf, "d_over_D": ef, "c2": cf ** 2, "b0_deg": float(np.degrees(np.arccos(min(cf, 1)))),
                          "l0_deg": lf, "R2": r2(A_m, curve(ell_m, Nf, akf, ef, cf, lf, W_cont))}
R["free_LSQ_analytic"] = {k: float(v) for k, v in R["free_LSQ_analytic"].items()}

# ---- numerical ideal-drag model: b0 = 0 / 15 / 33.5
K0 = MS0 * kbase(X0)
def num_summary(k, e, x, l0, b0):
    A71 = model_at_data(k, e, x, l0, b0); l, Af = model_fine(k, e, x, l0, b0); h = harm_fine(l, Af)
    pr = {}
    for m in (1, 2, 3):   # signed coefficient along the model axis: A_m cos m(l-l0)
        t = np.radians(m * l0); pr[f"m{m}_along_axis"] = float(h[f"c{m}"] * np.cos(t) + h[f"s{m}"] * np.sin(t))
    return {"kappa": float(k), "d_over_D": float(e), "D_over_L": float(x), "l0": float(l0), "b0": float(b0),
            "R2_direct": r2(A_m, A71), "R2_corr2": float(np.corrcoef(A_m, A71)[0, 1] ** 2), **pr,
            "frac_m_ge3": h["frac_power_m_ge3"]}
num = {"b0_only_changed": [num_summary(K0, EPS0, X0, L0DEG, b0) for b0 in (0.0, 15.0, 33.5)], "reoptimized": []}
for b0 in (0.0, 15.0, 33.5):
    def f(p):
        k, e, x, l0 = p
        if not (0.01 < k < 3 and 0.02 <= e <= 0.8 and 3 <= x <= 14): return 10.0
        return -r2(A_m, model_at_data(k, e, x, l0, b0))
    best = None
    for start in ([K0, EPS0, X0, L0DEG], [0.3, 0.5, 5.0, 15.0], [0.6, 0.7, 10.0, 15.0]):
        p, v = nelder_mead(f, start, [0.05, 0.05, 0.5, 3.0], tol=1e-10, maxit=1500)
        if best is None or v < best[1]: best = (p, v)
    num["reoptimized"].append(num_summary(*best[0], b0))
R["numeric_ideal_drag_b0_check"] = num
R["notes"] = ["W_cont = continuous 90-deg window (paper); the numeric code's boxcar is 47 samples of 360/181 deg = 93.5 deg, W_disc.",
              "Projections: parallel = c_m cos(m l0) + s_m sin(m l0); orthogonal part cannot be produced by an axisymmetric model.",
              "Numeric model: d/D restricted to <= 0.8 (u table range); D/L allowed 3-14 in the re-optimisation.",
              "71 points from overlapping 90-deg sectors are not independent."]
json.dump(R, open(f"{OUT}/dq_model_equality.json", "w"), indent=1)

# ---- figure
fig = plt.figure(figsize=(15, 12.5)); fig.patch.set_facecolor("white")
def tx(y, s, size=16, x=0.03, **kw): fig.text(x, y, s, fontsize=size, va="top", **kw)
s = sol; N, ak, e = s["N"], s["alpha_kappa"], s["d_over_D"]; y = s["c2"]; W = W_cont
mc = R["ratio_solve_l0_15_Wcont"]["model_coeffs"]
tx(0.985, r"Data $m$=0–3 fit $=$ analytic tidal model (main.tex eq:tidal), term by term", 20, weight="bold")
tx(0.945, r"Data (71 pts):  $A = %.4f %+.4f\cos\ell %+.4f\sin\ell %+.4f\cos2\ell %+.4f\sin2\ell %+.4f\cos3\ell %+.4f\sin3\ell$,  $R^2=%.3f$"
   % (c["c0"], c["c1"], c["s1"], c["c2"], c["s2"], c["c3"], c["s3"], R["data_fit_m0_3"]["R2"]), 14.5)
tx(0.905, r"Model:  $A/A_{\rm all} = N\,[\,1+\alpha\kappa(3\mu^2-1)+\alpha\kappa\,\frac{d}{D}(\frac{9}{2}\mu-\frac{15}{2}\mu^3)\,]$,"
          r"   $\mu=\cos b_0\cos(\ell-\ell_0)$,  $c\equiv\cos b_0$,  $W_m=\frac{\sin(m\pi/4)}{m\pi/4}$", 14.5)
tx(0.862, r"Projected on $\ell_0=15^\circ$ (data coefficient of $\cos m(\ell-\ell_0)$)  =  model expression  =  value at solution:", 14.5, color="#333")
lines = [
 (r"$m=0$:   $%.4f$" % P15[0], r"$N[1+\alpha\kappa(\frac{3}{2}c^2-1)]$", r"$%.3f\,[1+%.3f(\frac{3}{2}\cdot%.3f-1)] = %.4f$" % (N, ak, y, mc["m0"])),
 (r"$m=1$:   $%+.4f$" % P15[1], r"$N\alpha\kappa\frac{d}{D}\,c\,(\frac{9}{2}-\frac{45}{8}c^2)\,W_1$",
  r"$%.3f\cdot%.3f\cdot%.3f\cdot%.3f\,(4.5-5.625\cdot%.3f)\cdot%.3f = %+.4f$" % (N, ak, e, np.sqrt(y), y, W[1], mc["m1"])),
 (r"$m=2$:   $%+.4f$" % P15[2], r"$N\alpha\kappa\,\frac{3}{2}c^2\,W_2$", r"$%.3f\cdot%.3f\cdot1.5\cdot%.3f\cdot%.3f = %+.4f$" % (N, ak, y, W[2], mc["m2"])),
 (r"$m=3$:   $%+.4f$" % P15[3], r"$-N\alpha\kappa\frac{d}{D}\,\frac{15}{8}c^3\,W_3$", r"$-%.3f\cdot%.3f\cdot%.3f\cdot1.875\cdot%.3f\cdot%.3f = %+.4f$" % (N, ak, e, np.sqrt(y) ** 3, W[3], mc["m3"])),
]
yy = 0.822
for a_, b_, c_ in lines:
    tx(yy, a_, 15.5); tx(yy, "=  " + b_, 15.5, x=0.17); tx(yy, "=  " + c_, 15.5, x=0.47); yy -= 0.042
tx(yy - 0.002, (r"Solve:  $m_1/m_3 \Rightarrow c^2=%.3f$ ($|b_0|=%.1f^\circ$);   $m_3/m_2 \Rightarrow d/D=%.3f$;   $m_2 \Rightarrow N\alpha\kappa=%.4f$;   $m_0 \Rightarrow N=%.4f$, $\alpha\kappa=%.4f$"
     % (y, s["b0_deg"], e, s["N_alpha_kappa"], N, ak)), 14.5, color="#8b0000")
Qs = R["projections_l0_15"]["orthogonal_unmodelled"]
tx(yy - 0.042, (r"Not produced by an axisymmetric model (data $\sin m(\ell-15^\circ)$ parts):  $m_1$: $%+.4f$,  $m_2$: $%+.4f$,  $m_3$: $%+.4f$"
                % (Qs["m1"], Qs["m2"], Qs["m3"])), 13.5, color="#444")
tx(yy - 0.074, (r"With the numeric code's 93.5° window ($W_{1,2,3}$=0.893, 0.612, 0.262): $c^2=%.3f$ ($|b_0|=%.1f^\circ$), $d/D=%.3f$, $\alpha\kappa=%.4f$."
                % (sol_d["c2"], sol_d["b0_deg"], sol_d["d_over_D"], sol_d["alpha_kappa"])), 13.5, color="#444")
nr = num["reoptimized"]
tx(yy - 0.112, (r"Full numeric ideal-drag model, re-optimised ($\kappa$, d/D, D/L, $\ell_0$) at fixed $b_0$:  $R^2$ = %.3f ($b_0$=0°), %.3f (15°), %.3f (33.5°, D/L at 14 bound)"
                % (nr[0]["R2_direct"], nr[1]["R2_direct"], nr[2]["R2_direct"])), 13.5, color="#00008b")
tx(yy - 0.146, ("→ b0≈33° does not help the numeric model (its m1 is already >0 along ℓ0 at b0=0, from higher orders); "
                "the truncated analytic form needs cos²b0<0.8 to get m1>0."), 12, color="#00008b")
ax = fig.add_axes([0.06, 0.06, 0.62, 0.40])
ax.errorbar(ell_m, A_m, yerr=[A_lo, A_hi], fmt="o", ms=4, color="k", ecolor="0.82", elinewidth=0.8, label="Migkas+20 Fig. 6 (digitised)")
lg = np.linspace(0, 360, 721); Lg = np.radians(lg)
fit_curve = c["c0"] + sum(c[f"c{m}"] * np.cos(m * Lg) + c[f"s{m}"] * np.sin(m * Lg) for m in (1, 2, 3))
ax.plot(lg, fit_curve, color="C0", lw=2.4, label=rf"data $m$=0–3 fit, $R^2$={R['data_fit_m0_3']['R2']:.3f}")
ax.plot(lg, curve(lg, N, ak, e, np.sqrt(y), 15.0, W), color="C3", lw=2.4,
        label=rf"analytic tidal model, solved ($\ell_0$=15°, $|b_0|$={s['b0_deg']:.1f}°, $d/D$={e:.2f}, $\alpha\kappa$={ak:.3f}): $R^2$={R['ratio_solve_l0_15_Wcont']['R2_vs_71pts']:.3f}")
F = R["free_LSQ_analytic"]
ax.plot(lg, curve(lg, F["N"], F["alpha_kappa"], F["d_over_D"], np.sqrt(F["c2"]), F["l0_deg"], W), color="C2", lw=1.5, ls="--",
        label=rf"analytic model, free LSQ ($\ell_0$={F['l0_deg']:.1f}°, $|b_0|$={F['b0_deg']:.1f}°, $d/D$={F['d_over_D']:.2f}, $\alpha\kappa$={F['alpha_kappa']:.3f}): $R^2$={F['R2']:.3f}")
ax.set_xlim(0, 360); ax.set_xticks(range(0, 361, 45)); ax.set_xlabel(r"Galactic longitude $\ell$ [deg]", fontsize=13)
ax.set_ylabel(r"$A/A_{\rm all}$", fontsize=13); ax.legend(fontsize=9.5, loc="lower left", framealpha=0.92)
ax2 = fig.add_axes([0.74, 0.06, 0.24, 0.40])
ms_ = [1, 2, 3]; xx = np.arange(3); w = 0.36
ax2.bar(xx - w / 2, [P15[m] for m in ms_], w, color="C0", label=r"data, along $\ell_0$=15°")
ax2.bar(xx + w / 2, [mc[f"m{m}"] for m in ms_], w, color="C3", label="analytic model")
ax2.bar(xx - w / 2, [Qs[f"m{m}"] for m in ms_], w * 0.45, color="none", edgecolor="k", hatch="//", label=r"data $\perp$ (unmodelled)")
ax2.axhline(0, color="k", lw=0.7); ax2.set_xticks(xx); ax2.set_xticklabels([r"$m=1$", r"$m=2$", r"$m=3$"], fontsize=13)
ax2.set_ylabel("coefficient", fontsize=12); ax2.legend(fontsize=9, loc="upper right")
fig.text(0.01, 0.008, "71 points come from overlapping 90° sectors (not independent). The analytic form is a small-d/D expansion evaluated at d/D≈0.8.", fontsize=9, color="#444")
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/dq_model_equality.{ext}", dpi=300)
print(json.dumps({k: R[k] for k in R if k not in ("numeric_ideal_drag_b0_check", "notes")}, indent=1))
for k in num:
    for r_ in num[k]: print(k, {a: round(v, 4) for a, v in r_.items()})
