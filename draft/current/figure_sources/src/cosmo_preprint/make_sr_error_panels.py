"""Regenerate separate SR error-reduction panels from existing run-2 fit (no new DE)."""
import json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

NOTES = Path("/workspace/cosmology-course/notes")
OUT = Path("/workspace/cosmo_preprint/figs")

src = (NOTES / "cmb-cancel.py").read_text().split("if __name__")[0]
g = {}
exec(src, g)
out = json.load(open(NOTES / "cmb-cancel.json"))
r = out["run2"]
P = [r[k] for k in ("beta1", "q", "mu", "s", "phi_deg")]
ELL, G, full = g["ELL"], g["G_FULL"], g["full"]
S, Dm, Dg = g["S_model"](*P)
phi = P[4]
th = np.radians(ELL - phi)
Sreq = (1 / G) / np.mean(1 / G)
R = g["resid"](G, S)
M = 12
c = np.array([2 * np.mean(R * np.cos(m * th)) for m in range(1, M + 1)])
modes = np.array([c[m - 1] * np.cos(m * th) for m in range(1, M + 1)])
series = modes.sum(0)
left = R - R.mean() - series
bins = np.loadtxt(NOTES / "migkas-A-vs-angle-bins.csv", delimiter=",", skiprows=1)
plt.rcParams.update({"font.size": 14, "axes.titlesize": 15})

def save(fig, name):
    fig.tight_layout()
    path = OUT / name
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)

fig, a = plt.subplots(figsize=(11, 4.2))
a.errorbar(bins[:, 0], bins[:, 1], yerr=[bins[:, 2], bins[:, 3]], fmt="o", color="gray", ms=4, capsize=2, label="Migkas Fig. 6 bins")
a.plot(ELL, full, "k", lw=2.4, label=r"shell geometry: $D$=58.2, $R$=14.3 Gpc, $\ell_0$=15.2°")
a.set_xlim(0, 360); a.grid(alpha=0.3); a.set_xlabel(r"Galactic longitude $\ell$ (deg)"); a.set_ylabel(r"$A/A_{\rm all}$")
a.set_title("1. Our Migkas fit (no sine), nothing removed"); a.legend(loc="best", frameon=False, fontsize=12)
save(fig, "Fig_SR_panel1_migkas_fit.png")

fig, a = plt.subplots(figsize=(11, 4.2))
a.plot(ELL, Sreq, "k", lw=3, label=r"required $S_{\rm req}$")
a.plot(ELL, S / np.mean(S), "C1--", lw=2.2, label=r"model $S=D_{\rm motion}D_{\rm grav}$")
a.plot(ELL, Dm / np.mean(Dm), "C0", lw=1.5, label=r"$D_{\rm motion}$")
a.plot(ELL, Dg / np.mean(Dg), "C2", lw=1.5, label=r"$D_{\rm grav}$")
a.set_xlim(0, 360); a.grid(alpha=0.3); a.set_xlabel(r"Galactic longitude $\ell$ (deg)"); a.set_ylabel("factor (mean 1)")
a.set_title(rf"2. SR + mass: $\beta_1$={P[0]:.3f}, $q$={P[1]:.3f}, $\mu$={P[2]:.4f}, $s$={P[3]:.0f}, $\varphi$={phi:.2f}°")
a.legend(loc="best", frameon=False, fontsize=11, ncol=2)
save(fig, "Fig_SR_panel2_S_correction.png")

fig, a = plt.subplots(figsize=(11, 4.2))
a.plot(ELL, R * 1e3, "C3", lw=2.4, label=rf"residual $R$: RMS {np.sqrt(np.mean(R**2)):.2e}")
a.plot(ELL, series * 1e3, "k--", lw=1.6, label=r"$\sum_{m=1}^{12} c_m\cos m(\ell-\varphi)$")
a.set_xlim(0, 360); a.grid(alpha=0.3); a.set_xlabel(r"Galactic longitude $\ell$ (deg)"); a.set_ylabel(r"$R\times10^{3}$")
a.set_title("3. After SR + mass, no modes removed"); a.legend(loc="best", frameon=False, fontsize=12)
save(fig, "Fig_SR_panel3_residual.png")

fig, a = plt.subplots(figsize=(11, 4.2))
cols = list(plt.cm.tab10.colors) + ["k", "gold"]
for m in range(1, M + 1):
    if abs(c[m - 1]) > 2e-5:
        a.plot(ELL, modes[m - 1] * 1e3, color=cols[m - 1], lw=2.0, label=f"m={m}")
a.set_xlim(0, 360); a.grid(alpha=0.3); a.set_xlabel(r"Galactic longitude $\ell$ (deg)"); a.set_ylabel(r"$c_m\cos m(\ell-\varphi)\times10^{3}$")
a.set_title(r"4. Modes $m=1$–$12$"); a.legend(loc="best", frameon=False, fontsize=10, ncol=6)
save(fig, "Fig_SR_panel4_modes.png")

fig, a = plt.subplots(figsize=(11, 4.2))
a.plot(ELL, left * 1e5, "C0", lw=2, label=rf"leftover RMS {np.sqrt(np.mean(left**2)):.2e}")
a.axhline(1, ls="--", color="gray"); a.axhline(-1, ls="--", color="gray", label=r"$\pm10^{-5}$")
a.set_xlim(0, 360); a.grid(alpha=0.3); a.set_xlabel(r"Galactic longitude $\ell$ (deg)"); a.set_ylabel(r"leftover $\times10^{5}$")
a.set_title("5. Full-precision, $m=1$–$12$ removed"); a.legend(loc="best", frameon=False, fontsize=12)
save(fig, "Fig_SR_panel5_leftover_full.png")

rr = out["run2_rounded"]; Pr = [rr[k] for k in ("beta1", "q", "mu", "s", "phi_deg")]
Rr = g["resid"](g["G_RND"], g["S_model"](*Pr)[0]); thr = np.radians(ELL - Pr[4]); Lr = Rr - Rr.mean()
for mm in range(1, M + 1):
    Lr = Lr - 2 * np.mean(Lr * np.cos(mm * thr)) * np.cos(mm * thr) - 2 * np.mean(Lr * np.sin(mm * thr)) * np.sin(mm * thr)
fig, a = plt.subplots(figsize=(11, 4.2))
a.plot(ELL, Lr * 1e5, "C3", lw=1.3, label=rf"4-decimal: RMS {np.sqrt(np.mean(Lr**2)):.2e}")
a.axhline(1, ls="--", color="gray"); a.axhline(-1, ls="--", color="gray", label=r"$\pm10^{-5}$"); a.axhline(0, color="k", lw=0.6)
a.set_xlim(0, 360); a.grid(alpha=0.3); a.set_xlabel(r"Galactic longitude $\ell$ (deg)"); a.set_ylabel(r"leftover $\times10^{5}$")
a.set_title("6. Rounded CSV leftover (rounding noise)"); a.legend(loc="best", frameon=False, fontsize=12)
save(fig, "Fig_SR_panel6_leftover_rounded.png")

rms_chain = [5.751e-02, 2.373e-03, 6.113e-04, 2.168e-04, 1.833e-06]; cmb = 1e-5
lab = ["Shell fit\n(start)", "+ SR +\nmass", "m=1–5\nremoved", "m=1–6\nremoved", "m=1–12\nremoved"]
x = [v / cmb for v in rms_chain]
fig, ax = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw=dict(hspace=0.85, height_ratios=[3, 1]))
a = ax[0]; a.bar(range(5), x, color=["k", "C3", "C1", "C2", "C0"]); a.set_yscale("log"); a.set_ylim(0.05, 3e4)
a.axhline(1, ls="--", color="purple", lw=2, label="CMB level (1×)")
for i, v in enumerate(x):
    a.text(i, v * 1.3, f"{v:,.0f}×" if v >= 1 else f"{v:.2f}×", ha="center", fontsize=14)
for i in range(1, 5):
    a.text(i - 0.5, 1.2e4, f"÷{x[i-1]/x[i]:.1f}", ha="center", fontsize=13, color="dimgray", fontweight="bold")
a.set_xticks(range(5)); a.set_xticklabels(lab); a.set_ylabel("RMS ΔT/T in CMB units (log)")
a.legend(loc="upper right", frameon=False); a.set_title("Anisotropy left after each step")
b = ax[1]; parts = [x[i - 1] - x[i] for i in range(1, 5)]; leftx = 0
for p, n, cc in zip(parts, ["SR + mass", "m=1–5", "m=6", "m=7–12"], ["C3", "C1", "C2", "C0"]):
    b.barh(0, p, left=leftx, color=cc, label=f"{n}: {p:,.0f}× ({100*p/x[0]:.3f}%)"); leftx += p
b.barh(0, x[4], left=leftx, color="purple", label=f"left: {x[4]:.2f}×"); b.set_xlim(0, x[0] * 1.02); b.set_yticks([])
b.set_xlabel("share of the starting 5,751× (linear)")
b.legend(loc="upper center", bbox_to_anchor=(0.5, 1.55), ncol=2, frameon=False, fontsize=11)
b.set_title("Where the starting anisotropy went", pad=8)
save(fig, "Fig_SR_waterfall.png")
print("done", P)
