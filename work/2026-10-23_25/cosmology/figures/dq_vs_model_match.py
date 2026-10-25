#!/usr/bin/env python3
"""Dipole+quadrupole ("dq") fit of the 71 digitised Migkas+2020 Fig.6 A/A_all points,
harmonic decomposition of the physical ideal-drag Sedov-Taylor + Gauss model, and
coefficient-by-coefficient matching of the model to the data dq fit.
Outputs: dq_vs_model_match.{png,pdf,json} in /workspace/cosmology/figures."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dq_model_lib import (u_base, ell_m, A_m, A_lo, A_hi, MS0, EPS0, X0, L0DEG, kbase, model_fine,
                          model_at_data, harm71, harm_fine, model_harm, lsq, amp_phase, r2, nelder_mead)

OUT = "/workspace/cosmology/figures"
K0 = MS0 * kbase(X0)
res = {"notes": ["Model = ideal-drag physical model copied from make_model_param_heatmaps.py (identical arithmetic).",
                 "'lsq71' = same LSQ (a0 + cos/sin l + cos/sin 2l) at the 71 data longitudes as the data fit; "
                 "'fine' = exact Fourier series of the model on its 181-point full-circle grid.",
                 "Phases: l1 = direction of dipole maximum; l2 = quadrupole maximum (mod 180).",
                 "71 points come from overlapping 90-deg sectors stepped by 5 deg: not independent."]}

# ------------------------------------------------------------------ 1. data fit
N = len(A_m)
fits = {}
for name, ms in (("dq", (1, 2)), ("quad_only", (2,)), ("dip_only", (1,)), ("dip_quad_oct", (1, 2, 3))):
    b, X, pred = lsq(ell_m, A_m, ms)
    fits[name] = {"coef": b.tolist(), "R2": r2(A_m, pred), "RMSE": float(np.sqrt(((A_m - pred) ** 2).mean()))}
b, X, pred = lsq(ell_m, A_m, (1, 2))
resid = A_m - pred
sig2 = (resid ** 2).sum() / (N - 5)
cov = sig2 * np.linalg.inv(X.T @ X)
# autocorrelation of residuals (5-deg steps) -> effective N
rc = resid - resid.mean()
rho = np.array([np.sum(rc[:-k] * rc[k:]) / np.sum(rc * rc) for k in range(1, 30)])
kz = int(np.argmax(rho <= 0)) if np.any(rho <= 0) else len(rho)
tau = 1 + 2 * rho[:kz].sum()
infl_ac = np.sqrt(tau)
infl_overlap = np.sqrt(90.0 / 5.0)   # one 90-deg sector spans 18 steps

def derived(beta):
    a1, l1 = amp_phase(beta[1], beta[2], 1); a2, l2 = amp_phase(beta[3], beta[4], 2)
    return np.array([beta[0], a1, l1, a2, l2])
d0 = derived(b)
J = np.zeros((5, 5)); h = 1e-7
for i in range(5):
    bp = b.copy(); bp[i] += h; J[:, i] = (derived(bp) - d0) / h
sd = np.sqrt(np.diag(J @ cov @ J.T))
names = ["a0", "a1", "l1_deg", "a2", "l2_deg"]
data_dq = {n: float(v) for n, v in zip(names, d0)}
data_dq.update({"c1": float(b[1]), "s1": float(b[2]), "c2": float(b[3]), "s2": float(b[4])})
res["task1_data_fit"] = {
    "values": data_dq, "R2": fits["dq"]["R2"], "RMSE": fits["dq"]["RMSE"],
    "sigma_formal_from_scatter": dict(zip(names, sd.tolist())),
    "residual_autocorr_lag1..": rho[:kz + 1].tolist(), "tau_int": float(tau),
    "N_eff_from_residual_autocorr": float(N / tau),
    "sigma_inflated_autocorr": dict(zip(names, (sd * infl_ac).tolist())),
    "sigma_inflated_sector_overlap_x4.24": dict(zip(names, (sd * infl_overlap).tolist())),
    "other_fits": fits,
    "median_quoted_Migkas_errbar": float(np.median(0.5 * (A_lo + A_hi))),
}
dl12 = (data_dq["l1_deg"] - data_dq["l2_deg"])
res["task1_data_fit"]["l1_minus_l2_deg"] = dl12
T1 = res["task1_data_fit"]

# ------------------------------------------------------------------ 2. model decompositions
def pack(kappa, eps, x, l0, b0, label):
    m = model_harm(kappa, eps, x, l0, b0)
    return {"label": label, "params": {"kappa": kappa, "d_over_D": eps, "D_over_L": x, "l0": l0, "b0": b0}, **m}
poses = [pack(K0, EPS0, X0, L0DEG, 0.0, "best fit"),
         pack(0.15, EPS0, X0, 0.0, 15.0, "paper pose (0,+15), kappa=0.15"),
         pack(K0, EPS0, X0, 0.0, 15.0, "paper pose (0,+15), kappa=0.458")]
res["task2_model_harmonics"] = poses

# ------------------------------------------------------------------ 3. scalings
def row(kappa, eps, x, l0, b0, nrad=50):
    l, A = model_fine(kappa, eps, x, l0, b0, nrad=nrad)
    hf = harm_fine(l, A)
    return {"a1": hf["a1"], "l1": hf["l1"], "a2": hf["a2"], "l2": hf["l2"], "a3": hf["a3"],
            "a1_over_a2": hf["a1"] / hf["a2"], "frac_m_ge3": hf["frac_power_m_ge3"]}
scan = {}
scan["d_over_D"] = [{"v": e, **row(K0, e, X0, L0DEG, 0.0), "nrad400": row(K0, e, X0, L0DEG, 0.0, nrad=400)}
                    for e in (0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.759, 0.8)]
scan["kappa"] = [{"v": k, **row(k, EPS0, X0, L0DEG, 0.0)} for k in (0.05, 0.1, 0.15, 0.3, 0.458, 0.6, 0.8, 1.0)]
scan["D_over_L"] = [{"v": x, **row(K0, EPS0, x, L0DEG, 0.0)} for x in (3.0, 4.0, 5.0, 6.0, 7.62, 9.0, 10.0)]
scan["l0"] = [{"v": l0, **row(K0, EPS0, X0, l0, 0.0)} for l0 in (0.0, 14.38, 45.0, 90.0, 135.0, 194.38, 270.0)]
scan["b0"] = [{"v": b0, **row(K0, EPS0, X0, L0DEG, b0)} for b0 in (0.0, 15.0, 30.0, 45.0, 60.0)]
# unsmoothed pre-boxcar amplitudes at best fit (boxcar transfer check)
l, Araw = model_fine(K0, EPS0, X0, L0DEG, 0.0, smooth=False)
hr = harm_fine(l, Araw); hs = poses[0]["fine"]
scan["boxcar_transfer_best_fit"] = {f"m{m}": hs[f"a{m}"] / hr[f"a{m}"] for m in (1, 2, 3, 4, 5)}
# small-d power laws (nrad=400 to avoid linear-interp artefacts)
es = np.array([0.05, 0.1, 0.2]); a1s = [row(K0, e, X0, L0DEG, 0, 400)["a1"] for e in es]
a2s = [row(K0, e, X0, L0DEG, 0, 400)["a2"] for e in es]
scan["loglog_slope_small_d(0.05-0.2,nrad400)"] = {"a1": float(np.polyfit(np.log(es), np.log(a1s), 1)[0]),
                                                  "a2": float(np.polyfit(np.log(es), np.log(a2s), 1)[0])}
res["task3_scalings"] = scan

# ------------------------------------------------------------------ 4. matching
TA2, TL2, TA1, TL1 = data_dq["a2"], data_dq["l2_deg"], data_dq["a1"], data_dq["l1_deg"]
TA1proj = TA1 * np.cos(np.radians(TL1 - TL2))
res["task4_targets"] = {"a2": TA2, "l2": TL2, "a1": TA1, "l1": TL1,
                        "a1_projected_on_axis_l2": float(TA1proj)}

def solve_k_l0(eps, x, b0, k=0.4, l0=TL2):
    """match lsq71 a2 and l2 by adjusting kappa and l0 (others fixed)."""
    for _ in range(60):
        h = harm71(model_at_data(k, eps, x, l0, b0))
        dl = ((TL2 - h["l2"] + 90) % 180) - 90
        k_new = k * TA2 / h["a2"]; l0 += dl
        if abs(k_new - k) < 1e-10 and abs(dl) < 1e-9: k = k_new; break
        k = k_new
    h = harm71(model_at_data(k, eps, x, l0, b0))
    return k, l0, h

def summarize(k, eps, x, l0, b0, extra=None):
    A71 = model_at_data(k, eps, x, l0, b0); h = harm71(A71)
    l, A = model_fine(k, eps, x, l0, b0); hf = harm_fine(l, A)
    d = {"kappa": k, "d_over_D": eps, "D_over_L": x, "l0": l0, "b0": b0,
         "a0": h["a0"], "a1": h["a1"], "l1": h["l1"], "a2": h["a2"], "l2": h["l2"],
         "frac_m_ge3": hf["frac_power_m_ge3"], "a3_fine": hf["a3"],
         "R2_direct": r2(A_m, A71), "R2_corr2": float(np.corrcoef(A_m, A71)[0, 1] ** 2)}
    if extra: d.update(extra)
    return d

# (A) match a2, l2 only: family over (d/D, D/L, b0); kappa & l0 solved
famA = []
for b0 in (0.0, 15.0, 30.0):
    for x in (3.0, 5.0, 7.62, 10.0):
        for eps in (0.1, 0.3, 0.5, 0.759):
            k, l0, _ = solve_k_l0(eps, x, b0)
            famA.append(summarize(k, eps, x, l0, b0))
# member of family A with max R2 (optimize eps, x, b0; kappa,l0 solved)
def negR2A(p):
    eps, x, b0 = p
    if not (0.02 <= eps <= 0.8 and 3.0 <= x <= 10.0 and -60 <= b0 <= 60): return 10.0
    k, l0, _ = solve_k_l0(eps, x, b0)
    return -r2(A_m, model_at_data(k, eps, x, l0, b0))
pA, _ = nelder_mead(negR2A, [EPS0, X0, 0.0], [0.05, 0.5, 5.0], tol=1e-9, maxit=400)
kA, l0A, _ = solve_k_l0(*pA)
bestA = summarize(kA, *pA[:2], l0A, pA[2])
# the kappa=0.15 question: which D/L (at best d/D, b0=0) or which b0 gives a2 with kappa=0.15?
res["task4A_match_a2_l2"] = {"family_samples": famA, "max_R2_member": bestA}

# (B) match a2, l2 AND a1 amplitude: solve eps by bisection (a1 rises with eps)
def solve_eps_for_a1(x, b0, target):
    def f(e):
        k, l0, h = solve_k_l0(e, x, b0); return h["a1"] - target, k, l0
    lo, hi = 0.02, 0.8
    flo, fhi = f(lo)[0], f(hi)[0]
    if flo * fhi > 0: return None, (flo + target, fhi + target)
    for _ in range(50):
        mid = 0.5 * (lo + hi); fm = f(mid)[0]
        if flo * fm <= 0: hi = mid
        else: lo, flo = mid, fm
    e = 0.5 * (lo + hi); _, k, l0 = f(e)
    return (k, e, l0), None
famB, famBp = [], []
for b0 in (0.0, 15.0):
    for x in (3.0, 4.0, 5.0, 6.0, 7.62, 9.0, 10.0):
        for tgt, fam in ((TA1, famB), (TA1proj, famBp)):
            sol, rng = solve_eps_for_a1(x, b0, tgt)
            if sol is None:
                fam.append({"D_over_L": x, "b0": b0, "solution": None, "a1_range_eps0.02_0.8": list(rng)})
            else:
                k, e, l0 = sol; fam.append(summarize(k, e, x, l0, b0))
res["task4B_match_a2_l2_a1amp"] = famB
res["task4Bp_match_a2_l2_a1projected"] = famBp

# (C) unconstrained LSQ in harmonic space: minimize coefficient mismatch of (c1,s1,c2,s2)
tv = np.array([data_dq["c1"], data_dq["s1"], data_dq["c2"], data_dq["s2"]])
def coefmis(p):
    k, eps, x, l0, b0 = p
    if not (0.01 < k < 2 and 0.02 <= eps <= 0.8 and 3 <= x <= 10 and -60 <= b0 <= 60): return 10.0
    h = harm71(model_at_data(k, eps, x, l0, b0))
    return float(((np.array([h["c1"], h["s1"], h["c2"], h["s2"]]) - tv) ** 2).sum())
pC, fC = nelder_mead(coefmis, [K0, EPS0, X0, L0DEG, 0.0], [0.05, 0.05, 0.5, 5, 5], tol=1e-14, maxit=3000)
bestC = summarize(*pC, extra={"coef_rms_mismatch": float(np.sqrt(fC / 4)),
                              "min_possible_dipole_mismatch_given_phase_lock": float(TA1 * abs(np.sin(np.radians(TL1 - TL2))))})
res["task4C_min_coef_mismatch"] = bestC

def flag(fam, need_a1=None):
    for f in fam:
        if f.get("kappa") is None: f["converged"] = False; continue
        ok = abs(f["a2"] - TA2) < 2e-4 and abs(((f["l2"] - TL2 + 90) % 180) - 90) < 0.05 and abs(((f["l0"] - TL2 + 90) % 180) - 90) < 1
        if need_a1 is not None: ok = ok and abs(f["a1"] - need_a1) < 2e-4
        f["converged"] = bool(ok)
flag(famA); flag(famB, TA1); flag(famBp, TA1proj)
for nm, fam in (("A", famA), ("B", famB), ("Bp", famBp)):
    ok = [f["R2_direct"] for f in fam if f["converged"]]
    res[f"task4{nm}_R2_range_converged"] = [min(ok), max(ok), len(ok), len(fam)]

# kappa = 0.15 with a2,l2 matched: which D/L is needed? (bisection in D/L, b0=0)
def dl_for_kappa(kt, eps, b0=0.0):
    g = lambda x: solve_k_l0(eps, x, b0)[0] - kt
    lo, hi = 7.0, 16.0
    glo, ghi = g(lo), g(hi)
    if glo * ghi > 0: return None
    for _ in range(40):
        mid = 0.5 * (lo + hi); gm = g(mid)
        if glo * gm <= 0: hi = mid
        else: lo, glo = mid, gm
    x = 0.5 * (lo + hi); k, l0, _ = solve_k_l0(eps, x, b0)
    return summarize(k, eps, x, l0, b0)
res["task4A_kappa0.15_needs_DL"] = [dl_for_kappa(0.15, e) for e in (0.3, 0.5, 0.759)]

# chosen matched model for the figure: a2,l2,|a1| matched at the best-fit D/L=7.62, b0=0
chosen = [f for f in famB if f.get("converged") and f["b0"] == 0.0 and abs(f["D_over_L"] - 7.62) < 1e-9][0]
res["figure_matched_model"] = chosen
# l1 test: rotate the matched model so its axis sits on the data dipole (l0 = l1_data)
res["task4_l1_test_axis_at_l1data"] = summarize(chosen["kappa"], chosen["d_over_D"], chosen["D_over_L"], TL1, chosen["b0"])
# data octupole (for the m>=3 comparison)
bo = fits["dip_quad_oct"]["coef"]; a3d, l3d = amp_phase(bo[5], bo[6], 3)
res["task1_data_fit"]["octupole_fit"] = {"a3": a3d, "l3_deg": l3d, "R2": fits["dip_quad_oct"]["R2"],
                                         "note": "axisymmetric model predicts l3 = l0 + 60 (mod 120) when its m=3 coefficient is negative"}

# ---- analytic small-d/D expansion (Taylor of H_obs in eps = d/D), f(r) = u(r)/r, units H=D=1, per unit kappa
#  dH(mu) = f + D f' mu^2 + eps [ (3/2) D f' mu + (1/2)(D^2 f'' - D f') mu^3 ] + O(eps^2),  mu = cos b0 cos(l - l0)
#  => Q2 = T2 kappa cos^2 b0 g1/2 / N,  D1 = T1 kappa eps [ (3/2) cos b0 g1 + (3/8) cos^3 b0 (g2 - g1) ] / N
#     N = 1 + kappa (f + g1 cos^2 b0 / 2);  T_m = discrete boxcar transfer (T1=0.8928, T2=0.6119)
def gfun(x):
    rg, ub, kb = u_base(x); u = ub / kb
    m = (rg > 0.6) & (rg < 1.4); p = np.polyfit(rg[m] - 1, u[m] / rg[m], 4)
    return float(p[-1]), float(p[-2]), float(2 * p[-3])
Tm = scan["boxcar_transfer_best_fit"]; T1_, T2_ = Tm["m1"], Tm["m2"]
def analytic(k, eps, x, b0):
    f0, g1, g2 = gfun(x); cb = np.cos(np.radians(b0)); Nn = 1 + k * (f0 + g1 * cb ** 2 / 2)
    Q2 = T2_ * k * cb ** 2 * g1 / 2 / Nn
    D1 = T1_ * k * eps * (1.5 * cb * g1 + 0.375 * cb ** 3 * (g2 - g1)) / Nn
    return {"f": f0, "g1_Dfprime": g1, "g2_D2fpp": g2, "N": Nn, "Q2_leading": Q2, "D1_leading_signed": D1}
anal = {"formula": "Q2 = T2*kappa*cos^2(b0)*g1/2/N ; D1 = T1*kappa*(d/D)*[1.5 cos b0 g1 + 0.375 cos^3 b0 (g2-g1)]/N ; "
                   "g1 = D f'(D)/kappa, g2 = D^2 f''(D)/kappa, f = u_pec/r ; D1<0 means dipole max at l0+180",
        "T1": T1_, "T2": T2_,
        "g_of_DL": {str(x): gfun(x) for x in (3.0, 4.0, 5.0, 6.0, 7.62, 9.0, 10.0, 12.0, 14.0)},
        "check_small_d": [{"d_over_D": e, **analytic(K0, e, X0, 0.0), "numeric": row(K0, e, X0, L0DEG, 0.0, 400)}
                          for e in (0.05, 0.1, 0.2)],
        "at_matched_solution": {**analytic(chosen["kappa"], chosen["d_over_D"], chosen["D_over_L"], chosen["b0"]),
                                "numeric_a1": chosen["a1"], "numeric_a2": chosen["a2"]}}
res["analytic_expansion"] = anal
json.dump(res, open(f"{OUT}/dq_vs_model_match.json", "w"), indent=1)

# ------------------------------------------------------------------ 5. figure
lg = np.linspace(0, 360, 721)
Xg = np.stack([np.ones_like(lg), np.cos(np.radians(lg)), np.sin(np.radians(lg)),
               np.cos(2 * np.radians(lg)), np.sin(2 * np.radians(lg))], 1)
fig = plt.figure(figsize=(12.5, 4.8))
ax = fig.add_axes([0.06, 0.13, 0.52, 0.78])
ax.errorbar(ell_m, A_m, yerr=[A_lo, A_hi], fmt="o", ms=4, color="k", ecolor="0.8", elinewidth=0.8,
            label="Migkas+20 Fig. 6 (digitised, 90° sectors)")
ax.plot(lg, Xg @ b, color="C0", lw=2.2, label=f"data dipole+quadrupole fit, $R^2$={T1['R2']:.3f}")
lf, Af = model_fine(chosen["kappa"], chosen["d_over_D"], chosen["D_over_L"], chosen["l0"], chosen["b0"])
ax.plot(np.append(lf, 360), np.append(Af, Af[0]), color="C3", lw=2.2,
        label=(rf"matched physical model ($a_2,\ell_2,|a_1|$): $\kappa$={chosen['kappa']:.3f}, "
               rf"$d/D$={chosen['d_over_D']:.3f}, $D/L$={chosen['D_over_L']:.2f}, "
               rf"$\ell_0$={chosen['l0']:.1f}°, $b_0$={chosen['b0']:.0f}°; $R^2$={chosen['R2_direct']:.3f}"))
lb, Ab = model_fine(K0, EPS0, X0, L0DEG, 0.0)
ax.plot(np.append(lb, 360), np.append(Ab, Ab[0]), color="C1", lw=1.3, ls="--",
        label=f"physical best fit (R²-optimal), $R^2$={poses[0]['R2_direct']:.3f}")
ax.set_xlim(0, 360); ax.set_xticks(range(0, 361, 45))
ax.set_xlabel(r"Galactic longitude $\ell$ [deg]"); ax.set_ylabel(r"$A/A_{\rm all}$")
ax.legend(fontsize=7.3, loc="lower left", framealpha=0.9)
ax.set_title("Migkas A(ℓ): data dq fit vs matched physical model", fontsize=10)
ax2 = fig.add_axes([0.66, 0.13, 0.32, 0.78])
keys = ["c1", "s1", "c2", "s2", "a1", "a2"]
dv = [data_dq[k] for k in keys]
hm = harm71(model_at_data(chosen["kappa"], chosen["d_over_D"], chosen["D_over_L"], chosen["l0"], chosen["b0"]))
hb = poses[0]["lsq71"]
mv = [hm[k] for k in keys]; bv = [hb[k] for k in keys]
err = [np.nan, np.nan, np.nan, np.nan, T1["sigma_inflated_autocorr"]["a1"], T1["sigma_inflated_autocorr"]["a2"]]
xx = np.arange(len(keys)); w = 0.27
ax2.bar(xx - w, dv, w, color="C0", label="data dq fit", yerr=err, capsize=3)
ax2.bar(xx, mv, w, color="C3", label="matched model")
ax2.bar(xx + w, bv, w, color="C1", alpha=0.6, label="best-fit model")
ax2.axhline(0, color="k", lw=0.6)
ax2.set_xticks(xx); ax2.set_xticklabels([r"$c_1$", r"$s_1$", r"$c_2$", r"$s_2$", r"$a_1$", r"$a_2$"])
ax2.set_ylabel(r"coefficient of $A/A_{\rm all}$")
ax2.set_title(rf"phases  data: $\ell_1$={data_dq['l1_deg']:.0f}°, $\ell_2$={data_dq['l2_deg']:.1f}°" "\n"
              rf"matched: $\ell_1$={hm['l1']:.0f}°, $\ell_2$={hm['l2']:.1f}°  (model $\ell_1$ locked to $\ell_0$)",
              fontsize=8.5)
ax2.legend(fontsize=7.5, loc="upper left")
fig.text(0.01, 0.01, "71 points from overlapping 90° sectors are not independent; error bars on $a_1,a_2$ "
         "inflated by residual autocorrelation.", fontsize=7, color="#444")
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/dq_vs_model_match.{ext}", dpi=300)
print("done")

json.dump(res, open(f"{OUT}/dq_vs_model_match.json", "w"), indent=1)
print("equality figure: see dq_model_equality.py")
