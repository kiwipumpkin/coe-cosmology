import numpy as np, csv, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
z = np.load("primary_series.npz")
th, BH, Af, lmap = z["th"], z["BH"], z["Af"], z["lmap"]
mig = list(csv.DictReader(open("../../migkas/fig6_digitized.csv")))
ml = np.array([float(r["l_deg"]) for r in mig]); mA = np.array([float(r["A_over_Aall"]) for r in mig])
lo = np.array([float(r["A_err_low"]) for r in mig]); hi = np.array([float(r["A_err_high"]) for r in mig])
def AFf(l):
    l = np.radians(l); return (0.9965+0.0296*np.cos(l)+0.0434*np.sin(l)+0.0636*np.cos(2*l)+0.0366*np.sin(2*l)-0.0299*np.cos(3*l)-0.0114*np.sin(3*l))
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
fig, (a, b) = plt.subplots(2, 1, figsize=(6.5, 5.2), sharex=True, gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.06})
a.errorbar(ml, mA, yerr=[lo, hi], fmt="o", ms=3, color="0.45", ecolor="0.7", lw=0.8, label="Migkas et al. Fig. 6 (71 pts)")
a.plot(th, Af, color="C0", lw=1.6, label=r"$A_F(\ell)$, Fourier fit")
a.plot(th, BH, color="C3", lw=1.6, label=r"Chris col. BH, $\theta=\ell$")
o = np.argsort(lmap[:-1]); a.plot(lmap[:-1][o], BH[:-1][o], color="C3", lw=0.9, ls="--", label=r"Chris col. BH at sheet lookup $\ell(\theta)$")
a.set_ylabel(r"$A/A_{\rm all}$"); a.legend(fontsize=7.5, loc="lower left", ncol=1, frameon=False); a.set_ylim(0.74, 1.22)
b.axhline(0, color="k", lw=0.6)
b.plot(th, BH - Af, color="C3", lw=1.4, label=r"BH$-A_F$ ($\theta=\ell$)")
b.plot(ml, mA - AFf(ml), "o", ms=2.5, color="0.45", label=r"data$-A_F$")
b.set_ylabel("residual"); b.set_xlabel(r"angle (deg): sheet $\theta$ = Galactic $\ell$")
b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45)); b.legend(fontsize=7.5, frameon=False, loc="upper left", ncol=2)
fig.savefig("A_compare.png", dpi=200, bbox_inches="tight"); fig.savefig("A_compare.pdf", bbox_inches="tight")
