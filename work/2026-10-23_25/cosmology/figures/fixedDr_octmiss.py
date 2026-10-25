#!/usr/bin/env python3
"""Fixed-D/r 'let the octupole miss' analysis for the analytic Sedov-Taylor CoE model vs Migkas+2020 Fig.6 A(l).
Reuses data, Fourier coefficients, latitude/sector coefficients and the weighted-LSQ solver of termmatch_Dr234.py.
Outputs figures/fixedDr_octmiss_results.json and figures/fixedDr_octmiss.{png,pdf}.
Harmonic-amplitude sigmas: the quoted Fourier coefficients are the OLS fit to the 71 points. Three error models:
  naive : 71 points independent, per-point sigma=(lo+hi)/2, WLS covariance (the weights used in termmatch_Dr234.py)
  neff16: naive x sqrt(71/16) (rule of thumb: ~16 independent sectors)
  corr  : sandwich covariance of the OLS estimator, B C B^T, with C_ij = s_i s_j max(0, 1-|l_i-l_j|/90)
          (overlap fraction of two 90-deg sectors; uniform clusters). PRIMARY.
"""
import sys, json
import numpy as np
from scipy.optimize import brentq, minimize_scalar
sys.path.insert(0, "/workspace/cosmology/figures")
import termmatch_Dr234 as T
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/workspace/cosmology/figures/fixedDr_octmiss"
d2r = np.pi / 180; wrap = T.wrap
amp, ph, l0 = T.amp, T.ph, T.l0
X, y, s_i, l_d = T.X, T.A_d, T.sig_i, T.l_d
R_GPC = 14.3

# ---------- sigmas ----------
def amp_sig(cov):
    out = {0: float(np.sqrt(cov[0, 0]))}
    for m in (1, 2, 3):
        i, j = 2 * m - 1, 2 * m; g = np.array([T.coef[i], T.coef[j]]) / amp[m]
        out[m] = float(np.sqrt(g @ cov[np.ix_([i, j], [i, j])] @ g))
    return out
sig_naive = dict(T.sig)
sig_n16 = {m: v * np.sqrt(71 / 16) for m, v in sig_naive.items()}
dl = np.abs((l_d[:, None] - l_d[None, :] + 180) % 360 - 180)
C = np.outer(s_i, s_i) * np.clip(1 - dl / 90, 0, None)
B = np.linalg.solve(X.T @ X, X.T)
cov_corr = B @ C @ B.T
sig_corr = amp_sig(cov_corr)
SIGS = dict(corr=sig_corr, neff16=sig_n16, naive=sig_naive)
resid_rms = float(np.std(y - X @ T.coef))

coefs = T.coefs
def cm_of(b0, eps):
    t = coefs(b0); return t, {1: eps * abs(t["o1"]), 2: abs(t["q2"]), 3: eps * abs(t["o3"])}

# ---------- model curve with Chris's transforms ----------
def build(N, K, b0, eps):
    t = coefs(b0); P = N * K
    Ms = {0: N * (1 + K * t["q0"]), 1: P * eps * t["o1"], 2: P * t["q2"], 3: P * eps * t["o3"]}
    nat = {m: (l0 + (0 if Ms[m] > 0 else 180 / m)) % (360 / m) for m in (1, 2, 3)}
    phi = {m: wrap(ph[m] - nat[m], 360 / m) for m in (1, 2, 3)}
    def curve(l, p3=False):
        L = l * d2r
        return (Ms[0] + abs(Ms[1]) * np.cos(L - (nat[1] + phi[1]) * d2r) + abs(Ms[2]) * np.cos(2 * (L - (nat[2] + phi[2]) * d2r))
                + abs(Ms[3]) * np.cos(3 * (L - (nat[3] + (phi[3] if p3 else 0)) * d2r)))
    return Ms, nat, phi, curve, t

def package(label, N, K, b0, eps, extra=None):
    Ms, nat, phi, curve, t = build(N, K, b0, eps)
    d = dict(label=label, D_over_r=1 / eps, eps=eps, N=N, K=K, NK=N * K, b0=b0, a=float(np.cos(b0 * d2r)),
             q0=t["q0"], q2=t["q2"], o1=t["o1"], o3=t["o3"],
             model_amp={m: float(abs(Ms[m])) for m in Ms}, model_signed={m: float(Ms[m]) for m in Ms},
             resid={m: float(abs(Ms[m]) - amp[m]) for m in Ms},
             pull={k: {m: float((abs(Ms[m]) - amp[m]) / S[m]) for m in Ms} for k, S in SIGS.items()},
             chi2={k: float(sum(((abs(Ms[m]) - amp[m]) / S[m]) ** 2 for m in (1, 2, 3))) for k, S in SIGS.items()},
             phi1=float(phi[1]), phi2=float(phi[2]), delta3=float(phi[3]), nat_peak={m: float(nat[m]) for m in nat},
             R2=T.r2(curve(l_d)), R2_oct_phase_free=T.r2(curve(l_d, True)),
             D_Gpc=R_GPC / eps, D_Gly=R_GPC / eps * 3.26156, MoverMcrit_alpha1=2 * K,
             phiC_model=float(wrap(90 - nat[1], 360)), phiC=float(wrap(90 - (nat[1] + phi[1]), 360)),
             phi2C=float(wrap(90 - 2 * (nat[2] + phi[2]), 360)))
    if extra: d.update(extra)
    d["_curve"] = curve
    return d

# ---------- 1. free eps: exact 4x4 ----------
ratio31 = amp[3] / amp[1]
f31 = lambda b: abs(coefs(b)["o3"]) / abs(coefs(b)["o1"]) - ratio31
bflip = np.degrees(np.arccos(np.sqrt(3.1958 / 3.9948)))  # o1 = 0 (approx; refine)
bflip = brentq(lambda b: coefs(b)["o1"], 20, 35)
free = []
for lo, hi, lab in ((0.0, bflip - 1e-6, "A (on-axis dipole)"), (bflip + 1e-6, 80, "B (flipped dipole)")):
    b0 = brentq(f31, lo, hi, xtol=1e-12)
    t = coefs(b0); P = amp[2] / abs(t["q2"]); eps = amp[3] / (P * abs(t["o3"])); N = amp[0] - P * t["q0"]
    free.append(package("free eps " + lab, N, P / N, b0, eps))
ratio32 = amp[3] / amp[2]
kappa3 = 0.22806 / 0.49052 # |o3|/q2 per unit a
eps_min = ratio32 / (abs(coefs(0)["o3"]) / coefs(0)["q2"])

# ---------- 2. fixed eps, exact m0,m1,m2 -> predict m3 ----------
ratio12 = amp[1] / amp[2]
def exact012(eps):
    sols = []
    g = lambda b: eps * abs(coefs(b)["o1"]) / coefs(b)["q2"] - ratio12
    for lo, hi, br in ((0.0, bflip - 1e-9, "on-axis"), (bflip + 1e-9, 89.0, "flipped")):
        glo, ghi = g(lo), g(hi)
        if br == "on-axis":
            if glo < 0: continue  # max of eps|o1|/q2 on this branch is at b0=0
            if g(hi) > 0: continue
            b0 = brentq(g, lo, hi, xtol=1e-12)
        else:
            if glo * ghi > 0: continue
            b0 = brentq(g, lo, hi, xtol=1e-12)
        t = coefs(b0); P = amp[2] / t["q2"]; N = amp[0] - P * t["q0"]
        sols.append(package(f"exact m0,m1,m2 ({br})", N, P / N, b0, eps, dict(branch=br)))
    return sols

# ---------- 3. fixed eps, weighted LSQ on m1..m3 (m0 exact via N) ----------
def lsq(eps, S, bgrid=np.arange(0, 89.5, 0.05)):
    w = {m: 1 / S[m] ** 2 for m in (1, 2, 3)}
    def solveP(b0):
        t, cm = cm_of(b0, eps)
        P = sum(w[m] * amp[m] * cm[m] for m in cm) / sum(w[m] * cm[m] ** 2 for m in cm)
        return P, sum(w[m] * (amp[m] - P * cm[m]) ** 2 for m in cm), t
    ch = np.array([solveP(b)[1] for b in bgrid]); mins = []
    for i in range(len(bgrid)):
        if (i == 0 or ch[i] <= ch[i - 1]) and (i == len(bgrid) - 1 or ch[i] <= ch[i + 1]):
            lo_, hi_ = bgrid[max(i - 1, 0)], bgrid[min(i + 1, len(bgrid) - 1)]
            b = minimize_scalar(lambda b: solveP(b)[1], bounds=(lo_, hi_), method="bounded", options=dict(xatol=1e-9)).x
            mins.append((solveP(b)[1], b))
    mins.sort(); out = []
    for chi2w, b0 in mins:
        P, chi2w, t = solveP(b0); N = amp[0] - P * t["q0"]
        out.append(package("LSQ", N, P / N, float(b0), eps, dict(chi2_fit=float(chi2w), branch="on-axis" if t["o1"] > 0 else "flipped")))
    return out

DRS = [2.0, 3.0, 4.0]
res_exact = {Dr: exact012(1 / Dr) for Dr in DRS}
res_lsq_naive = {Dr: lsq(1 / Dr, sig_naive) for Dr in DRS}
res_lsq_corr = {Dr: lsq(1 / Dr, sig_corr) for Dr in DRS}

# compare with stored termmatch results
prev = json.load(open("/workspace/cosmology/figures/termmatch_Dr234_results.json"))["analytic"]
cmp = {}
for Dr in DRS:
    rows = []
    for p in prev[str(Dr)]:
        best = min(res_lsq_naive[Dr], key=lambda s: abs(s["b0"] - p["b0"]))
        rows.append(dict(b0_prev=p["b0"], b0_new=best["b0"], N_prev=p["N"], N_new=best["N"], K_prev=p["K"], K_new=best["K"],
                         R2_prev=p["R2_transformed"], R2_new=best["R2"], m3_prev=p["model_amp"]["3"], m3_new=best["model_amp"][3]))
    cmp[Dr] = rows

# ---------- 4. continuous curves vs D/r ----------
Drg = np.linspace(1.0, 6.0, 251)
cont = dict(Dr=list(Drg), m3_exact_on=[], m3_exact_flip=[], m3_lsq_naive=[], m3_lsq_corr=[], chi2_corr_profile=[], b0_lsq_naive=[])
for Dr in Drg:
    ex = {s["branch"]: s for s in exact012(1 / Dr)}
    cont["m3_exact_on"].append(ex["on-axis"]["model_amp"][3] if "on-axis" in ex else None)
    cont["m3_exact_flip"].append(ex["flipped"]["model_amp"][3] if "flipped" in ex else None)
    ln = lsq(1 / Dr, sig_naive, np.arange(0, 89.5, 0.25))[0]; lc = lsq(1 / Dr, sig_corr, np.arange(0, 89.5, 0.25))[0]
    cont["m3_lsq_naive"].append(ln["model_amp"][3]); cont["b0_lsq_naive"].append(ln["b0"])
    cont["m3_lsq_corr"].append(lc["model_amp"][3]); cont["chi2_corr_profile"].append(lc["chi2"]["corr"])
# D/r where exact on-axis m3 prediction crosses the data and its +-1 sigma band (corr sigma)
def cross(target):
    f = lambda Dr: exact012(1 / Dr)[0]["model_amp"][3] - target
    try: return brentq(f, 1.0, 2.0)
    except ValueError: return None
Dr_on_max = 1 / (ratio12 / (abs(coefs(0)["o1"]) / coefs(0)["q2"]))  # beyond this no on-axis exact solution
chi2p = np.array(cont["chi2_corr_profile"])
Dr_dchi1 = float(Drg[np.argmax(chi2p > chi2p.min() + 1)]) if np.any(chi2p > chi2p.min() + 1) else None

# ---------- JSON ----------
def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items() if k != "_curve"}
    if isinstance(o, list): return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    return o
R = dict(doc=__doc__, data=dict(amp=amp, peak=ph, l0=l0, coef=list(T.coef), R2_fourier=T.R2_fourier, resid_rms=resid_rms,
         mean_point_sigma=float(s_i.mean())), sigmas=SIGS, cov_corr=cov_corr.tolist(),
         ratios=dict(m3_over_m1=ratio31, m3_over_m2=ratio32, m1_over_m2=ratio12, o3_over_q2_at_a1=abs(coefs(0)["o3"]) / coefs(0)["q2"],
                     eps_min=eps_min, Dr_max=1 / eps_min, b0_o1_zero=bflip, Dr_on_axis_max=Dr_on_max),
         free=free, exact=res_exact, lsq_naive=res_lsq_naive, lsq_corr=res_lsq_corr, compare_prev=cmp, cont=cont,
         cross=dict(data=cross(amp[3]), lo1=cross(amp[3] - sig_corr[3]), lo2=cross(amp[3] - 2 * sig_corr[3])),
         chi2_profile_min=float(chi2p.min()), Dr_at_chi2min=float(Drg[np.argmin(chi2p)]), Dr_dchi2_1=Dr_dchi1)
json.dump(clean(R), open(OUT + "_results.json", "w"), indent=1)

# ---------- figure ----------
S = sig_corr
fig = plt.figure(figsize=(11, 8.2))
gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.22)
ax = fig.add_subplot(gs[0, 0])
cases = [("data", None), ("free $\\varepsilon$\n$D/r$=%.4f" % free[0]["D_over_r"], free[0])] + \
        [(f"$D/r$={int(Dr)}\n(LSQ)", res_lsq_naive[Dr][0]) for Dr in DRS]
cols = {1: "#1f77b4", 2: "#2ca02c", 3: "#d62728"}; wbar = 0.26
for i, (lab, s) in enumerate(cases):
    for k, m in enumerate((1, 2, 3)):
        v = amp[m] if s is None else s["model_amp"][m]
        ax.bar(i + (k - 1) * wbar, v, wbar * 0.92, color=cols[m], alpha=1.0 if s is None else 0.6,
               edgecolor="k", lw=0.6, label=(f"$|A_{m}|$" if i == 0 else None))
        if s is None: ax.errorbar(i + (k - 1) * wbar, v, yerr=S[m], color="k", capsize=3, lw=1)
for k, m in enumerate((1, 2, 3)):
    ax.axhline(amp[m], color=cols[m], lw=0.8, ls="--")
ax.set_xticks(range(len(cases))); ax.set_xticklabels([c[0] for c in cases], fontsize=8.5)
ax.set_ylabel("harmonic amplitude"); ax.set_ylim(0, 0.135)
ax.set_title("(a) amplitudes: model (pale) vs data (solid, $\\pm1\\sigma$ sector-correlated)", fontsize=9.5)
ax.legend(fontsize=8, ncol=3, loc="upper right")
for i, (lab, s) in enumerate(cases[1:], start=1):
    p = s["pull"]["corr"][3]; ax.text(i + wbar, s["model_amp"][3] + 0.004, f"{p:+.1f}$\\sigma$", ha="center", fontsize=7.5, color=cols[3])

ax = fig.add_subplot(gs[0, 1])
ax.axhspan(amp[3] - S[3], amp[3] + S[3], color="0.85", label="data $|A_3|\\pm1\\sigma$ (sector-correlated)")
ax.axhspan(amp[3] - sig_naive[3], amp[3] + sig_naive[3], color="0.7", label="data $\\pm1\\sigma$ (naive, 71 indep.)")
ax.axhline(amp[3], color="k", lw=1)
g = np.array(cont["Dr"])
on = np.array([np.nan if v is None else v for v in cont["m3_exact_on"]]); fl = np.array([np.nan if v is None else v for v in cont["m3_exact_flip"]])
ax.plot(g, on, color=cols[3], lw=2, label="exact $m_0,m_1,m_2$, on-axis branch")
ax.plot(g, fl, color=cols[3], lw=1.5, ls="--", label="exact $m_0,m_1,m_2$, flipped branch")
ax.plot(g, cont["m3_lsq_naive"], color="k", lw=1.2, ls=":", label="LSQ over $m_1..m_3$ (prev. weights)")
ax.plot(g, 0.46494 * amp[2] / g, color="0.4", lw=0.8, ls="-.", label="upper limit $0.465\\,|A_2|\\,\\varepsilon$ ($b_0=0$)")
for s in [free[0]] + [res_lsq_naive[Dr][0] for Dr in DRS]:
    ax.plot(s["D_over_r"], s["model_amp"][3], "o", color="k", ms=5, zorder=5)
    if s["D_over_r"] > 1.5:
        ax.annotate(f"{s['model_amp'][3]:.4f}\n({s['pull']['corr'][3]:+.1f}$\\sigma$)", (s["D_over_r"], s["model_amp"][3]),
                    xytext=(4, 12), textcoords="offset points", fontsize=7.5)
ax.annotate(f"free $\\varepsilon$: $D/r$={free[0]['D_over_r']:.4f}", (free[0]["D_over_r"], amp[3]), xytext=(12, 22), textcoords="offset points",
            fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7))
ax.set_xlim(1, 6); ax.set_ylim(0, 0.06); ax.set_xlabel("$D/r$ (held fixed)"); ax.set_ylabel("predicted $|A_3|$")
ax.set_title("(b) octupole predicted once $m_0,m_1,m_2$ are matched", fontsize=9.5)
ax.legend(fontsize=7, loc="upper right")

ax = fig.add_subplot(gs[1, 0])
ll = np.linspace(0, 360, 721)
ax.errorbar(l_d, y, yerr=[T.elo, T.ehi], fmt="o", ms=2.5, color="0.55", elinewidth=0.5, capsize=0, label="Migkas Fig.6 (71 overlapping sectors)")
ax.plot(ll, T.fourier(ll), color="k", lw=1.8, label="data Fourier fit ($m\\leq 3$)")
lc = {1: "#9467bd", 2: "#1f77b4", 3: "#ff7f0e", 4: "#d62728"}
curves = [(free[0], f"free $\\varepsilon$ ($D/r$=1.0075)", "#9467bd")] + [(res_lsq_naive[Dr][0], f"$D/r$={int(Dr)}", c) for Dr, c in zip(DRS, ("#1f77b4", "#2ca02c", "#d62728"))]
for s, lab, c in curves:
    ax.plot(ll, s["_curve"](ll), color=c, lw=1.1, label=f"{lab}: $R^2$={s['R2']:.3f}")
ax.set_xlim(0, 360); ax.set_xticks(range(0, 361, 60)); ax.set_xlabel("Galactic longitude $\\ell$ [deg]"); ax.set_ylabel("$A/A_{\\rm all}$")
ax.set_title("(c) $A(\\ell)$ with rotation + $2\\ell$ transforms (octupole at natural phase)", fontsize=9.5)
ax.legend(fontsize=6.8, loc="lower center", ncol=2, framealpha=0.95); ax.set_ylim(0.66, 1.26)

ax = fig.add_subplot(gs[1, 1])
for s, lab, c in curves:
    ax.plot(ll, s["_curve"](ll) - T.fourier(ll), color=c, lw=1.2, label=lab)
w3 = amp[3] * np.cos(3 * (ll - ph[3]) * d2r)
ax.plot(ll, -w3, color="k", lw=0.8, ls="--", label="$-$(data $3\\ell$ term)")
ax.axhspan(-S[3], S[3], color="0.9", zorder=0)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlim(0, 360); ax.set_xticks(range(0, 361, 60)); ax.set_xlabel("Galactic longitude $\\ell$ [deg]"); ax.set_ylabel("model $-$ data Fourier fit")
ax.set_title("(d) the difference is (almost) only the missing $3\\ell$ wiggle", fontsize=9.5)
ax.legend(fontsize=7, loc="upper center", ncol=3); ax.set_ylim(-0.045, 0.075)
fig.savefig(OUT + ".png", dpi=200, bbox_inches="tight"); fig.savefig(OUT + ".pdf", bbox_inches="tight")

# ---------- print ----------
print("sigmas", {k: {m: round(v, 4) for m, v in S_.items()} for k, S_ in SIGS.items()}, "resid rms", resid_rms)
print("ratios", R["ratios"]); print("cross", R["cross"], "chi2prof min", R["chi2_profile_min"], R["Dr_at_chi2min"], "dchi2=1 at", Dr_dchi1)
show = lambda s: print(f"  {s['label']:32s} D/r={s['D_over_r']:.4f} b0={s['b0']:.3f} N={s['N']:.4f} K={s['K']:.4f} m={[round(s['model_amp'][m],4) for m in range(4)]} "
                       f"pull_corr={[round(s['pull']['corr'][m],2) for m in (1,2,3)]} pull_naive={[round(s['pull']['naive'][m],2) for m in (1,2,3)]} "
                       f"pull16={[round(s['pull']['neff16'][m],2) for m in (1,2,3)]} chi2c={s['chi2']['corr']:.3f} R2={s['R2']:.4f} phi1={s['phi1']:.2f}")
for s in free: show(s)
for Dr in DRS:
    print("D/r", Dr)
    for s in res_exact[Dr]: show(s)
    for s in res_lsq_naive[Dr]: show(s)
    print("  corr-weighted LSQ:")
    for s in res_lsq_corr[Dr]: show(s)
print("compare", json.dumps(clean(cmp), indent=0))

# ---------- extra: where continuous m3 curves drop below data-1sigma; profile chi2 thresholds ----------
def first_below(vals, thr):
    v = np.array([np.nan if x is None else x for x in vals]); idx = np.where(v < thr)[0]
    return float(g[idx[0]]) if len(idx) else None
extra = {}
for key in ("m3_exact_on", "m3_exact_flip", "m3_lsq_naive", "m3_lsq_corr"):
    extra[key] = dict(below_1sig_corr=first_below(cont[key], amp[3] - sig_corr[3]), below_1sig_naive=first_below(cont[key], amp[3] - sig_naive[3]),
                      below_2sig_naive=first_below(cont[key], amp[3] - 2 * sig_naive[3]))
prof = dict(dchi2_1=float(g[np.argmax(chi2p > chi2p.min() + 1)]) if np.any(chi2p > chi2p.min() + 1) else None,
            dchi2_4=float(g[np.argmax(chi2p > chi2p.min() + 4)]) if np.any(chi2p > chi2p.min() + 4) else None,
            chi2_at={str(Dr): float(np.interp(Dr, g, chi2p)) for Dr in (1.0075, 2, 3, 4, 5, 6)})
# naive-sigma profile for comparison
chi2n = []
for Dr in g:
    chi2n.append(lsq(1 / Dr, sig_naive, np.arange(0, 89.5, 0.25))[0]["chi2"]["naive"])
chi2n = np.array(chi2n)
prof["naive_dchi2_1"] = float(g[np.argmax(chi2n > chi2n.min() + 1)]) if np.any(chi2n > chi2n.min() + 1) else None
prof["naive_dchi2_4"] = float(g[np.argmax(chi2n > chi2n.min() + 4)]) if np.any(chi2n > chi2n.min() + 4) else None
prof["naive_chi2_at"] = {str(Dr): float(np.interp(Dr, g, chi2n)) for Dr in (1.0075, 2, 3, 4, 5, 6)}
R2j = json.load(open(OUT + "_results.json")); R2j["below"] = extra; R2j["profile"] = prof
json.dump(R2j, open(OUT + "_results.json", "w"), indent=1)
print("below", extra); print("profile", prof)
