import sys, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
sys.path.insert(0, "../spreadsheet_model/integrated")
import integ_core as IC
import analyze as an
plt.rcParams.update({"font.size": 8.5, "font.family": "serif"})
B = np.load("res_B.npz"); C = np.load("res_C.npz"); SB = json.load(open("res_B.json")); scan = json.load(open("scan.json"))
OBS = {1: 1.23e-3}
for l, Dl in ((2, 225.9), (3, 936.9), (4, 692.2)):
    Cl = 2 * np.pi * Dl / (l * (l + 1) * (an.T0K * 1e6)**2); OBS[l] = (2 * l + 1) * np.sqrt(Cl / (4 * np.pi))

def save(fig, name):
    fig.savefig(name + ".png", dpi=200, bbox_inches="tight"); plt.close(fig)
    im = Image.open(name + ".png").convert("RGB"); s = 1700 / max(im.size)
    if s < 1: im = im.resize((int(im.width * s), int(im.height * s)), Image.LANCZOS)
    im.save(name + "_web.jpg", quality=88)

def boosted_curve(mu, T, beta, cf_smooth):
    from scipy.interpolate import PchipInterpolator
    from numpy.polynomial.legendre import legfit, legval
    x = np.linspace(-1, 1, 721)
    ip = (lambda y: legval(y, legfit(mu, T, len(mu) - 1))) if cf_smooth else PchipInterpolator(np.r_[-1, mu, 1], np.r_[T[0], T, T[-1]])
    g = 1 / np.sqrt(1 - beta**2); m0 = (x - beta) / (1 - beta * x); Tb = ip(m0) / (g * (1 - beta * x))
    w = np.polynomial.legendre.leggauss(400); Tm = 0.5 * np.sum(w[1] * ((ip((w[0] - beta) / (1 - beta * w[0]))) / (g * (1 - beta * w[0]))))
    return np.degrees(np.arccos(x)), Tb / Tm - 1

# ---------- Fig 1: Delta T / T vs angle from the CoE
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.9))
psi = np.degrees(B["psi"])
a.plot(psi, B["dT_a"], "C3", lw=1.5, label="fixed $t_{\\rm rec}$, dust frame")
xb, yb = boosted_curve(B["mu"], 1 + B["dT_a"], SB["beta_a"], False)
a.plot(xb, yb, "C0", lw=1.2, ls="--", label="after dipole-nulling boost ($\\beta=%.2f$)" % SB["beta_a"])
ok = np.isfinite(B["z_b"]); Tb = 1 / (1 + B["z_b"][ok])
a.plot(psi[ok], Tb / Tb.mean() - 1, "k.", ms=3, label="fixed density (exists only here)")
a.axvspan(SB["psi_max_b_deg"], 180, color="0.92", zorder=0)
a.text(112, 2.2, "these rays never reach\nthe recombination density", ha="center", fontsize=7, color="0.35")
a.set_xlim(0, 180); a.set_xlabel("angle from CoE direction $\\psi$ (deg)"); a.set_ylabel("$\\Delta T/T$")
a.set_title("(a) best-fit LTB (e$^{-r/L}$ universe)", fontsize=8.5); a.legend(fontsize=6.3, frameon=False, loc="upper right")
psiC = np.degrees(C["psi"])
b.plot(psiC, 1e4 * C["dT_a"], "C3", lw=1.5, label="dust frame, fixed $t_{\\rm rec}$")
b.plot(psiC, 1e4 * C["dT_b"], "C1", lw=1.0, ls=":", label="dust frame, fixed density")
xb, yb = boosted_curve(C["mu"], 1 + C["dT_a"], float(C["beta_a"]), True)
b.plot(xb, 1e4 * yb, "C0", lw=1.2, ls="--", label="after dipole-nulling boost")
b.axhspan(-1e4 * OBS[2], 1e4 * OBS[2], color="C2", alpha=0.25, label="observed quadrupole level")
b.set_xlim(0, 180); b.set_xlabel("angle from CoE direction $\\psi$ (deg)"); b.set_ylabel("$10^4\\,\\Delta T/T$")
b.set_title("(b) $\\Omega=0.3$ background + lump, $\\delta M=0.1M$", fontsize=8.5); b.legend(fontsize=6.3, frameon=False)
fig.tight_layout(); save(fig, "fig_dT")

# ---------- Fig 2: multipoles
lin = [r for r in scan if r["label"] == "Om_bg=0.3, dM=0.01 Mfit"][0]
lin9 = [r for r in scan if r["label"] == "Om_bg=0.9, dM=0.01 Mfit"][0]
AF4 = an.fourier(np.arange(360.), IC.AFfour(np.arange(360.)))
fac, fac9 = AF4[1] / lin["A_m2"], AF4[1] / lin9["A_m2"]
ells = np.array([1, 2, 3, 4])
sets = [("best-fit LTB, dust frame", np.abs(SB["c_a"][1:5]), "C3"),
        ("best-fit LTB, boosted", np.r_[np.nan, np.abs(SB["cb_a"][2:5])], "C1"),
        ("$\\Omega=0.3$ bg+lump, scaled, boosted", np.r_[fac * abs(lin["c1_a"]), fac * np.abs([lin["c2_a_boost"], lin["c3_a_boost"]]), np.nan], "C0"),
        ("$\\Omega=0.9$ bg+lump, scaled, boosted", np.r_[fac9 * abs(lin9["c1_a"]), fac9 * np.abs([lin9["c2_a_boost"], lin9["c3_a_boost"]]), np.nan], "C4"),
        ("observed (Planck)", np.array([OBS[l] for l in ells]), "k")]
fig, ax = plt.subplots(figsize=(4.8, 3.0)); wdt = 0.16
for i, (lab, v, col) in enumerate(sets):
    ax.bar(ells + (i - 2) * wdt, v, wdt, color=col, label=lab, alpha=0.85 if col != "k" else 1)
ax.set_yscale("log"); ax.set_xticks(ells); ax.set_xticklabels(["dipole", "quadrupole", "octupole", "$\\ell=4$"])
ax.set_ylabel("$|c_\\ell|$ in $\\Delta T/T=\\sum c_\\ell P_\\ell(\\cos\\psi)$"); ax.set_ylim(1e-6, 3)
ax.legend(fontsize=6.3, frameon=False, loc="upper right", ncol=1); ax.grid(axis="y", alpha=0.25)
fig.tight_layout(); save(fig, "fig_multipoles")

# ---------- Fig 3: nearby A(l)
gl = json.load(open("../spreadsheet_model/integrated/integ_results.json"))["mig_l0"]
th = np.arange(0, 360, 1.0)
af, _ = IC.AF(th, gl["DR"], gl["LD"], gl["mc"] / IC.MCRIT_FAC)
ell = np.arange(0, 360, 1.0)
afl = np.interp((ell - gl["l0"]) % 360, np.append(th, 360), np.append(af, af[0]))
A_newt = IC.A_from(afl, ell - gl["l0"], gl["amp"], gl["phi"])
fig, (a, b) = plt.subplots(2, 1, figsize=(6.5, 4.8), sharex=True, gridspec_kw={"height_ratios": [1.6, 1], "hspace": 0.07})
a.errorbar(IC.ML, IC.MA, yerr=[IC.MLO, IC.MHI], fmt="o", ms=2.8, color="0.5", ecolor="0.8", lw=0.7, label="Migkas et al. 2020, Fig. 6")
a.plot(ell, IC.AFfour(ell), "C0", lw=1.3, label="$A_F$ (Fourier fit to data)")
a.plot(ell, A_newt, "k", lw=1.2, ls="--", label="Newtonian best fit (AF $\\times$ sine)")
nbA = 1.0 / B["nb_A"]; nbA = nbA / nbA.mean(axis=0); zl = B["nb_zlist"]   # A proportional to D_L^-2 (H0^2) at fixed z
for k, (z, col) in enumerate(zip(zl, ["C3", "C3", "C1", "C1"])):
    if z in (0.05, 0.3): a.plot(B["nb_ell"], nbA[:, k], color=col, lw=1.4 if z == 0.05 else 1.0, ls="-" if z == 0.05 else "-.", label="LTB, $A\\propto D_L^{-2}$, $z=%.2f$" % z)
a.set_ylabel("$A/\\langle A\\rangle$"); a.set_ylim(0.72, 1.36); a.legend(fontsize=6.4, frameon=False, ncol=3, loc="upper center")
b.axhline(0, color="k", lw=0.5)
b.plot(ell, 100 * (afl / afl.mean() - 1), "k", lw=1.2, ls="--", label="Newtonian AF (gravity part)")
dHp = np.interp(np.cos(np.radians(ell - an.L_COE)), B["nb_mu"], B["nb_dH"][:, 1])
b.plot(ell, 100 * dHp, "C3", lw=1.3, label="LTB $\\delta H_0/H_0$ at $z=0.1$")
b.plot(ell, 100 * (IC.AFfour(ell) ** 0.5 / np.mean(IC.AFfour(ell) ** 0.5) - 1), "C0", lw=1.1, label="data: $A_F^{1/2}-1$ ($H_0\\propto A^{1/2}$)")
b.set_ylabel("fractional $H_0$ (%)"); b.set_xlabel("Galactic longitude $\\ell$ (deg), $b=0$"); b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45))
b.legend(fontsize=6.4, frameon=False, ncol=3, loc="lower center"); b.set_ylim(-18, 14)
fig.savefig("tmp.png"); save(fig, "fig_nearby")

# ---------- numbers for the tables
def m2phase(ell, y):
    t = np.radians(ell); c = np.mean(y * np.cos(2 * t)); s = np.mean(y * np.sin(2 * t)); return (np.degrees(np.arctan2(s, c)) / 2) % 180
tab = dict(obs=OBS, AF4=AF4.tolist(), AF4_phase2=m2phase(ell, IC.AFfour(ell)),
           newt_AF=an.fourier(ell, afl / afl.mean()).tolist(), newt_AF_phase2=m2phase(ell, afl),
           newt_A=an.fourier(ell, A_newt / A_newt.mean()).tolist(), newt_A_phase2=m2phase(ell, A_newt),
           ltb_A={str(z): an.fourier(B["nb_ell"], nbA[:, k]).tolist() for k, z in enumerate(zl)},
           ltb_A_phase2={str(z): m2phase(B["nb_ell"], nbA[:, k]) for k, z in enumerate(zl)},
           ltb_dH={str(z): an.fourier(ell, np.interp(np.cos(np.radians(ell - an.L_COE)), B["nb_mu"], B["nb_dH"][:, k])).tolist() for k, z in enumerate(zl)},
           data_H=an.fourier(ell, IC.AFfour(ell) ** 0.5 / np.mean(IC.AFfour(ell) ** 0.5)).tolist(), data_H_phase2=m2phase(ell, IC.AFfour(ell) ** 0.5),
           scale_fac=[fac, fac9],
           corr_rows=dict(newt_A=float(np.corrcoef(np.interp(IC.ML, ell, A_newt), IC.MA)[0, 1]), newt_AF=float(np.corrcoef(np.interp(IC.ML, ell, afl), IC.MA)[0, 1]),
                          **{"ltb_" + str(z): float(np.corrcoef(np.interp(IC.ML, B["nb_ell"], nbA[:, k]), IC.MA)[0, 1]) for k, z in enumerate(zl)},
                          ltb_dH=float(np.corrcoef(np.interp(IC.ML, ell, dHp), IC.MA)[0, 1])),
           corr_ltbA_data=float(np.corrcoef(np.interp(IC.ML, B["nb_ell"], nbA[:, 1]), IC.MA)[0, 1]),
           corr_newtA_data=float(np.corrcoef(np.interp(IC.ML, ell, A_newt), IC.MA)[0, 1]))
json.dump(tab, open("tables.json", "w"), indent=1, default=float); print(json.dumps(tab, indent=1, default=float))
