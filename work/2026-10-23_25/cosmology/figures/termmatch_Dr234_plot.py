#!/usr/bin/env python3
"""Figure for termmatch_Dr234.py results (reads termmatch_Dr234_results.json)."""
import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "/workspace/cosmology/figures/termmatch_Dr234"
R = json.load(open(OUT + "_results.json")); d2r = np.pi / 180
dat = np.genfromtxt("/workspace/cosmology/migkas/fig6_digitized.csv", delimiter=",", names=True)
l_d, A_d, elo, ehi = dat["l_deg"], dat["A_over_Aall"], dat["A_err_low"], dat["A_err_high"]
coef = np.array(R["data"]["coef"]); LG = np.array(R["LG"]); lf = np.linspace(0, 360, 721)
fourier = lambda l: np.column_stack([np.ones_like(l)] + [f(m * l * d2r) for m in (1, 2, 3) for f in (np.cos, np.sin)]) @ coef
def an_curve(s, l, tf=True):
    Ms = {int(k): v for k, v in s["model_signed"].items()}; nat = {int(k): v for k, v in s["nat_peak"].items()}
    p1, p2 = (s["phi1"], s["phi2"]) if tf else (0, 0); L = l * d2r
    return (Ms[0] + abs(Ms[1]) * np.cos(L - (nat[1] + p1) * d2r) + abs(Ms[2]) * np.cos(2 * (L - (nat[2] + p2) * d2r))
            + abs(Ms[3]) * np.cos(3 * (L - nat[3] * d2r)))
cyc = lambda y: np.append(y, y[0]); Lp = np.append(LG, 360)
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12.5, "legend.fontsize": 9.5})
fig = plt.figure(figsize=(17, 10.5))
gs = fig.add_gridspec(2, 3, height_ratios=[2.3, 1], hspace=0.34, wspace=0.08)
amp_d = R["data"]["amp"]; sig = R["data"]["sigma"]
for j, Dr in enumerate(["2.0", "3.0", "4.0"]):
    ax = fig.add_subplot(gs[0, j]); a = R["analytic"][Dr][0]; n = R["numeric"][Dr]
    ax.errorbar(l_d, A_d, yerr=[elo, ehi], fmt="o", ms=3.5, color="0.3", ecolor="0.72", elinewidth=0.9, zorder=1, label="Migkas+20 Fig. 6 (71 pts)")
    ax.plot(lf, fourier(lf), "-", color="0.45", lw=1.0, zorder=2, label=f"data Fourier fit m≤3  $R^2$={R['data']['R2_fourier']:.3f}")
    ax.plot(Lp, cyc(np.array(n["curve_natural"])), ":", color="C0", lw=2.0, zorder=3, label=f"numeric, no transforms  $R^2$={n['R2_notransform']:.3f}")
    ax.plot(Lp, cyc(np.array(n["curve_transformed"])), "-", color="C0", lw=2.4, zorder=4, label=f"numeric + rot. + 2ℓ  $R^2$={n['R2_transformed']:.3f}")
    ax.plot(lf, an_curve(a, lf), "-", color="C3", lw=1.9, zorder=5, label=f"analytic term-matched + rot. + 2ℓ  $R^2$={a['R2_transformed']:.3f}")
    ax.set_title(f"D/r = {float(Dr):g}   (ε = r/D = {1/float(Dr):.3f},  D = {float(Dr)*14.3:.1f} Gpc)")
    txt = (f"numeric: κ={n['kappa']:.3f}, D/L={n['D_over_L']:.3g}{'*' if n['DL_at_bound'] else ''}, f$_u$={n['f_u']:.2g}\n"
           f"   b$_0$={n['b0']:.0f}°{'*' if n['b0_at_bound'] else ''}, ℓ$_0$={n['l0']:.2f}°, φ$_1$={n['phi1']:.1f}°, φ$_2$={n['phi2']:.2f}°\n"
           f"analytic: N={a['N']:.4f}, K={a['K']:.4f}, b$_0$={a['b0']:.1f}°\n"
           f"   φ$_1$={a['phi1']:.1f}°, φ$_2$={a['phi2']:.0f}°, R$^2_{{no\\,tf}}$={a['R2_notransform']:.3f}\n"
           f"octupole mismatch δ$_3$ = {n['delta3']:.1f}° (both)")
    ax.text(0.02, 0.025, txt, transform=ax.transAxes, fontsize=9, va="bottom", family="DejaVu Sans",
            bbox=dict(fc="white", ec="0.75", alpha=0.9))
    ax.legend(loc="upper right", framealpha=0.9, handlelength=2.0)
    ax.set_xlim(0, 360); ax.set_xticks(range(0, 361, 60)); ax.grid(alpha=0.25); ax.set_ylim(0.62, 1.42)
    ax.set_xlabel("Galactic longitude ℓ [deg]")
    if j == 0: ax.set_ylabel("A / A$_{all}$")
    else: ax.set_yticklabels([])
    # harmonic bars
    bx = fig.add_subplot(gs[1, j]); ms = np.arange(1, 4); w = 0.26
    dv = [amp_d[str(m)] for m in ms]; nv = [n["harm_natural"][str(m)][0] for m in ms]; av = [a["model_amp"][str(m)] for m in ms]
    bx.bar(ms - w, dv, w, yerr=[sig[str(m)] for m in ms], color="0.6", capsize=3, label="data (±1σ)")
    bx.bar(ms, nv, w, color="C0", label="numeric")
    bx.bar(ms + w, av, w, color="C3", label="analytic")
    for x, v in zip(ms, nv): bx.text(x, v + 0.002, f"{v:.3f}", ha="center", fontsize=9, color="C0")
    for x, v in zip(ms + w, av): bx.text(x, v + 0.002, f"{v:.3f}", ha="center", fontsize=9, color="C3")
    for x, v, m in zip(ms - w, dv, ms): bx.text(x, v + sig[str(m)] + 0.002, f"{v:.3f}", ha="center", fontsize=9, color="k")
    bx.set_xticks(ms); bx.set_xticklabels([f"|m{m}|" for m in ms]); bx.set_ylim(0, 0.12); bx.grid(axis="y", alpha=0.25)
    bx.set_title(f"D/r = {float(Dr):g}: harmonic amplitudes\n(m0 = {amp_d['0']:.4f} data & analytic; 1.000 numeric, mean-normalised)", fontsize=10.5)
    if j == 0: bx.set_ylabel("amplitude"); bx.legend(loc="upper right", fontsize=9.5, ncol=3)
    else: bx.set_yticklabels([])
fig.suptitle("Sedov–Taylor CoE term matching at fixed D/r (sources on the vision shell, d = r): numeric ideal-drag model vs analytic point-mass form\n"
             "Rotation φ$_1$ shifts the m=1 harmonic, 2ℓ phase φ$_2$ shifts m=2; octupole left at its natural phase (ℓ$_0$ = 14.96° = data m=2 peak).  "
             "* = parameter at fit bound", fontsize=12.5)
fig.subplots_adjust(left=0.05, right=0.99, top=0.9, bottom=0.05)
fig.savefig(OUT + ".png", dpi=150); fig.savefig(OUT + ".pdf")
print("saved")
