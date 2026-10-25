#!/usr/bin/env python3
"""Figures for the exact term-matching sheet: exact_termmatch.{png,pdf} (exact vs series) and exact_profile.{png,pdf} (continuous D/r fit)."""
import sys, json
import numpy as np
sys.path.insert(0, "/workspace/cosmology/figures")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from exact_termmatch import tau_pm, amp, ph, l0, SIG, package, T, d2r, E
from exact_profile import gcurve, r2fit, chi2

R = json.load(open("/workspace/cosmology/figures/exact_termmatch_results.json"))
P = json.load(open("/workspace/cosmology/figures/exact_profile_results.json"))
O = json.load(open("/workspace/cosmology/figures/fixedDr_octmiss_results.json"))
Sc, Sn = SIG["corr"], SIG["naive"]
C = R["cont"]; Dr = np.array([c["Dr"] for c in C])
cols = ["k", "#1f77b4", "#2ca02c", "#d62728"]

def tau_ser_an(b0, eps):
    c = T.coefs(b0); return np.array([c["q0"], eps * c["o1"], c["q2"], eps * c["o3"]])

# ---------------- Figure 1: exact vs series ----------------
fig, axs = plt.subplots(2, 2, figsize=(11, 8.0)); fig.subplots_adjust(hspace=0.33, wspace=0.22)
for ax, b0, tag in ((axs[0, 0], 0.0, "(a)"), (axs[0, 1], 19.1, "(b)")):
    rat = np.array([c[f"ratio_b{b0}"] for c in C])
    for m in range(4):
        ax.plot(Dr, rat[:, m], color=cols[m], lw=1.6, label=[r"$\tau_0/q_0$ (monopole)", r"$\tau_1/(\varepsilon o_1)$ (dipole)", r"$\tau_2/q_2$ (quadrupole)", r"$\tau_3/(\varepsilon o_3)$ (octupole)"][m])
    ax.axhline(1, color="0.5", lw=0.7, ls="--"); ax.set_xlim(1, 6); ax.set_xlabel(r"$D/r$")
    ax.set_ylabel("exact / series coefficient"); ax.set_title(f"{tag} exact/series ratio of sector-averaged coefficients, $b_0={b0:g}^\\circ$", fontsize=9.5)
    ax.legend(fontsize=8, loc="lower right"); ax.set_ylim(-0.2 if b0 == 0 else 0.5, 1.9 if b0 == 0 else 1.25)
ax = axs[1, 0]
ax.axhspan(amp[3] - Sc[3], amp[3] + Sc[3], color="0.88", label=r"data $A_3\pm1\sigma$ (sector-correlated)")
ax.axhspan(amp[3] - Sn[3], amp[3] + Sn[3], color="0.74", label=r"data $A_3\pm1\sigma$ (naive)")
ax.axhline(amp[3], color="k", lw=1)
ex_on = [(c["Dr"], [e[1] for e in c["exact012"] if e[2] == "on-axis"]) for c in C]
ex_fl = [(c["Dr"], [e[1] for e in c["exact012"] if e[2] == "flipped"]) for c in C]
ax.plot([d for d, v in ex_on if v], [v[0] for d, v in ex_on if v], color="#d62728", lw=2, label="exact: match $m_0$-$m_2$, on-axis")
ax.plot([d for d, v in ex_fl if v], [v[0] for d, v in ex_fl if v], color="#d62728", lw=1.4, ls="--", label="exact: match $m_0$-$m_2$, flipped")
oc = O["cont"]; g = np.array(oc["Dr"])
son = np.array([np.nan if v is None else v for v in oc["m3_exact_on"]]); sfl = np.array([np.nan if v is None else v for v in oc["m3_exact_flip"]])
ax.plot(g, son, color="#1f77b4", lw=2, label="series: match $m_0$-$m_2$, on-axis")
ax.plot(g, sfl, color="#1f77b4", lw=1.4, ls="--", label="series: match $m_0$-$m_2$, flipped")
for s in R["free_pm"]:
    ax.plot(s["D_over_r"], amp[3], "o", color="#d62728", ms=6, zorder=6)
ax.plot(O["free"][0]["D_over_r"], amp[3], "s", color="#1f77b4", ms=6, zorder=6)
ax.set_xlim(1, 6); ax.set_ylim(0, 0.06); ax.set_xlabel(r"$D/r$ (held fixed)"); ax.set_ylabel(r"predicted $A_3$")
ax.set_title("(c) predicted octupole once $m_0,m_1,m_2$ are matched: exact vs series", fontsize=9.5)
ax.legend(fontsize=7, loc="upper right")
ax = axs[1, 1]
ll = np.linspace(0, 360, 721)
ax.errorbar(T.l_d, T.A_d, yerr=[T.elo, T.ehi], fmt="o", ms=2.5, color="0.6", elinewidth=0.5, label="Migkas Fig. 6")
ax.plot(ll, T.fourier(ll), color="k", lw=1.8, label=r"data Fourier fit ($m\leq 3$)")
def mk(sd, eps):
    t = tau_pm(sd["b0"], eps); return package("", t, sd["N"], sd["K"], eps, sd["b0"])
curves = [(R["free_pm"][0], "#9467bd", "exact free, on-axis")]
for key, c in (("2.0", "#1f77b4"), ("3.0", "#2ca02c"), ("4.0", "#d62728")):
    curves.append((R["fixed_pm"][key]["lsq_naive"][0], c, f"exact $D/r$={float(key):g}"))
for sd, c, lab in curves:
    s = mk(sd, sd["eps"]); ax.plot(ll, s["_curve"](ll), color=c, lw=1.1, label=f"{lab}: $R^2$={s['R2']:.3f}")
ax.set_xlim(0, 360); ax.set_xticks(range(0, 361, 60)); ax.set_ylim(0.66, 1.26)
ax.set_xlabel(r"Galactic longitude $\ell$ [deg]"); ax.set_ylabel(r"$A/A_{\rm all}$")
ax.set_title(r"(d) exact point-mass curves with rotation + $2\ell$ transforms", fontsize=9.5)
ax.legend(fontsize=6.8, loc="lower center", ncol=2, framealpha=0.95)
fig.savefig("/workspace/cosmology/figures/exact_termmatch.png", dpi=200, bbox_inches="tight")
fig.savefig("/workspace/cosmology/figures/exact_termmatch.pdf", bbox_inches="tight")

# ---------------- Figure 2: continuous D/r profile ----------------
DRg = np.array(P["DR"])
def prof(kind, name): return np.array([P["profile"][kind][str(float(d))][name][0] for d in DRg])
# series profile (analytic coefficients), same objectives
ser = {"chi2_corr": [], "chi2_naive": [], "negR2": []}
bg = np.arange(0, 89.76, 0.25)
for d in DRg:
    eps = 1 / d; ts = [tau_ser_an(b, eps) for b in bg]
    ser["chi2_corr"].append(min(chi2(t, Sc)[0] for t in ts)); ser["chi2_naive"].append(min(chi2(t, Sn)[0] for t in ts))
    ser["negR2"].append(min(-r2fit(t)[0] for t in ts))
ser = {k: np.array(v) for k, v in ser.items()}
json.dump({k: list(map(float, v)) for k, v in ser.items()}, open("/workspace/cosmology/figures/series_profile.json", "w"))
fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.3)); fig.subplots_adjust(wspace=0.27)
ax = axs[0]
pmc, pmn = prof("pm", "chi2_corr"), prof("pm", "chi2_naive")
exc, exn = np.minimum(prof("ext", "chi2_corr"), pmc), np.minimum(prof("ext", "chi2_naive"), pmn)
ax.plot(DRg, pmc, color="#d62728", lw=2, label="exact point mass, sector-corr.")
ax.plot(DRg, pmn, color="#d62728", lw=1.3, ls="--", label="exact point mass, naive")
ax.plot(DRg, exc, color="#9467bd", lw=1.6, ls="-.", label="exact extended, sector-corr.")
ax.plot(DRg, ser["chi2_corr"], color="#1f77b4", lw=1.3, label="series, sector-corr.")
ax.plot(DRg, ser["chi2_naive"], color="#1f77b4", lw=1.0, ls="--", label="series, naive")
for y in (1, 4): ax.axhline(y, color="0.5", lw=0.7, ls=":")
ax.set_xlim(1, 6); ax.set_ylim(0, 5); ax.set_xlabel(r"$D/r$"); ax.set_ylabel(r"profile $\chi^2$ (min over $b_0$, $K$, $N$ [, $D/L$])")
ax.set_title(r"(a) amplitude $\chi^2$ profile", fontsize=10); ax.legend(fontsize=7.2, loc="upper left")
ax = axs[1]
ax.plot(DRg, -prof("pm", "negR2"), color="#d62728", lw=2, label="exact point mass")
ax.plot(DRg, -np.minimum(prof("ext", "negR2"), prof("pm", "negR2")), color="#9467bd", lw=1.6, ls="-.", label="exact extended")
ax.plot(DRg, -ser["negR2"], color="#1f77b4", lw=1.3, label="series")
ax.axhline(T.R2_fourier, color="k", lw=0.8, ls="--", label=r"data Fourier fit ($m\leq3$)")
ax.set_xlim(1, 6); ax.set_ylim(0.8, 0.93); ax.set_xlabel(r"$D/r$"); ax.set_ylabel(r"profile $R^2$ (max over $b_0$, $K$, $N$ [, $D/L$])")
ax.set_title(r"(b) curve-fit $R^2$ profile (transforms included)", fontsize=10); ax.legend(fontsize=7.5, loc="lower center")
ax = axs[2]
fam = R["ext_family"]
for br, c, mk_ in (("on-axis", "#d62728", "o-"), ("flipped", "#1f77b4", "s--")):
    pts = sorted([(s["D_over_r"], s["D_over_L"], s["b0"]) for s in fam if s["branch"] == br])
    if pts:
        ax.plot([p[0] for p in pts], [p[1] for p in pts], mk_, color=c, ms=4, label=f"exact 4x4 solutions, {br}")
for s in R["free_pm"]:
    ax.axvline(s["D_over_r"], color="#d62728" if s["tau"]["1"] > 0 else "#1f77b4", lw=0.8, ls=":")
ax.set_yscale("log"); ax.set_xlim(1.0, 1.2); ax.set_xlabel(r"$D/r$"); ax.set_ylabel(r"$D/L$ of the exact solution")
ax.set_title("(c) extended density: where exact solutions exist", fontsize=10); ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(0.0, 0.86)); ax.set_ylim(20, 900)
ax.text(0.97, 0.97, "dotted: point-mass roots ($D/L\\to\\infty$)\nno exact solution beyond them", transform=ax.transAxes, ha="right", va="top", fontsize=7)
fig.savefig("/workspace/cosmology/figures/exact_profile.png", dpi=200, bbox_inches="tight")
fig.savefig("/workspace/cosmology/figures/exact_profile.pdf", bbox_inches="tight")
print("figs done")
