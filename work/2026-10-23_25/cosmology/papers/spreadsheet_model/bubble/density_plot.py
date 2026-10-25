import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from prof_core import shape, RG
from scipy.integrate import cumulative_trapezoid
from PIL import Image
A = json.load(open("analysis_rows.json"))
def get(kind, mat):
    for k, l, m, d in A:
        if k == kind and m == mat and d: return d
sel = [("exp", "LCDM-equiv.", "exponential $e^{-r/L}$ (main fit)", "C3", "-", 2.4),
       ("pl", "LCDM-equiv.", "power law $(1+r/L)^{-n}$", "C0", "--", 1.6),
       ("gauss", "LCDM-equiv.", "Gaussian core", "C2", "-.", 1.6),
       ("cored", "LCDM-equiv.", "cored $1/(1+(r/L)^2)$", "C9", (0, (1, 1)), 1.8),
       ("shell", "LCDM-equiv.", "Gaussian shell", "C1", "-", 1.6),
       ("st53", "LCDM-equiv.", r"Sedov-Taylor $\gamma=5/3$", "C4", "--", 1.6),
       ("st43", "LCDM-equiv.", r"Sedov-Taylor $\gamma=4/3$", "C6", ":", 2.0)]
R = 14.3
r = np.geomspace(0.02, 400, 4000)       # Gpc
fig, (a, b) = plt.subplots(2, 1, figsize=(8.5, 13), sharex=True, gridspec_kw=dict(height_ratios=[1.15, 1], hspace=0.06))
rows = []
for kind, mat, lab, c, ls, lw in sel:
    d = get(kind, mat); D = d["D_Gpc"]; x = RG
    s = shape(kind, x, d); m = cumulative_trapezoid(4 * np.pi * x**2 * s, x, initial=0); m1 = np.interp(1.0, x, m)
    rho = d["mc"] * (4 * np.pi / 3) * np.interp(r / D, x, s / m1, right=0)          # rho/rho_crit
    M = d["mc"] * np.interp(r / D, x, m / m1)                                         # M(<r)/Mcrit(D)
    rho = np.where(rho > 0, rho, np.nan)
    Msun = d["M_sun"]
    if kind in ("exp", "gauss", "cored"): par = "L = %.2f Gpc" % (d["LD"] * D)
    elif kind == "pl": par = "L = %.2f Gpc, n = %.1f" % (d["LD"] * D, d["n"])
    elif kind == "shell": par = "$r_s$ = %.1f Gpc, w = %.1f Gpc" % (d["rs"] * D, d["w"] * D)
    else: par = "front $R_s$ = %.1f Gpc" % (d["rst"] * D)
    leg = "%s, $R^2$=%.3f\n   %s, D = %.1f Gpc, M = %.2g M$_\\odot$" % (lab, d["R2_mig"], par, D, Msun)
    a.plot(r, rho, color=c, ls=ls, lw=lw, label=leg); b.plot(r, M, color=c, ls=ls, lw=lw)
    for ax in (a, b): ax.plot([D], [ax is a and 3e4 or 0.24], marker="v", color=c, ms=7, clip_on=False)
    rows.append((lab, kind, d, par))
# point + ST shell, baryons only
d = get("ptshell", "baryons only"); D = d["D_Gpc"]; c = "k"
rin, rout = d["rin"] * D, d["rout"] * D
shell = np.where((r > rin) & (r < rout), d["Os"], np.nan)
leg = ("point + constant-density shell (baryons only), $R^2$=%.3f\n   shell $\\Omega_s$=0.05 from %.1f to %.0f Gpc, D = %.1f Gpc\n"
       "   point core M$_c$ = %.2g M$_\\odot$ (arrow, scale %.2f Gpc)") % (d["R2_mig"], rin, rout, D, d["Mc_sun"], d["Lc"] * D)
a.plot(r, shell, color=c, lw=1.4, ls=(0, (5, 2, 1, 2)), label=leg)
a.vlines([rin], 1e-6, d["Os"], color=c, lw=1.4, ls=(0, (5, 2, 1, 2)))
a.annotate("", xy=(0.021, 2e4), xytext=(0.021, 3e2), arrowprops=dict(arrowstyle="-|>", color=c, lw=1.8))
a.text(0.024, 1.5e3, "point core\n%.1e M$_\\odot$" % d["Mc_sun"], fontsize=7.5, va="center")
Mpt = d["mcc"] + np.clip(d["Os"] * ((np.minimum(r, rout) / D)**3 - (rin / D)**3), 0, None)
b.plot(r, Mpt, color=c, lw=1.4, ls=(0, (5, 2, 1, 2)))
b.plot([0.021], [d["mcc"]], marker="o", color=c, ms=6)
b.axhline(d["mcc"], xmax=0.25, color=c, lw=0.8, alpha=0.5)
for ax in (a, b): ax.plot([D], [ax is a and 3e4 or 0.24], marker="v", color=c, ms=7, clip_on=False)
rows.append(("point + shell (baryons only)", "ptshell", d, ""))
# background and visible sphere (main fit)
dm = get("exp", "LCDM-equiv."); Dm = dm["D_Gpc"]; Rb = 1.5 * (Dm + R)
a.plot(r, np.where(r < Rb, 0.3, np.nan), color="0.45", lw=2.2, alpha=0.7, label="uniform background $\\Omega_b$=0.3 (LCDM-equivalent), to %.0f Gpc" % Rb)
a.plot(r, np.where(r < Rb, 0.05, np.nan), color="0.45", lw=1.4, ls="--", alpha=0.7, label="uniform background $\\Omega_b$=0.05 (baryons only)")
b.plot(r, 0.3 * (np.minimum(r, Rb) / Dm)**3, color="0.45", lw=2.2, alpha=0.5, label="background $\\Omega_b$=0.3 enclosed (main fit's D)")
for ax in (a, b):
    ax.axvspan(Dm - R, Dm + R, color="C2", alpha=0.10, lw=0)
    ax.axvline(Dm, color="C3", lw=1.0, ls="--")
a.text(Dm * 1.5, 1.5e3, "Earth, D = %.1f Gpc\n(main fit)" % Dm, color="C3", fontsize=8)
a.text((Dm - R) * 1.03, 3e2, "our visible\nsphere\n(D$-$R to D+R)", color="darkgreen", fontsize=8)
a.set_xscale("log"); a.set_yscale("log"); a.set_ylim(1e-6, 3e4); a.set_xlim(0.02, 400)
a.set_ylabel("density / critical density today")
a.set_title("Best-fit CoE density profiles (contrast ≤ 2%, with $\\ell_0$); triangles mark each fit's Earth distance D", fontsize=10)
h, l = a.get_legend_handles_labels(); hb, lb = b.get_legend_handles_labels()
fig.legend(h + hb, l + lb, fontsize=7.2, loc="upper center", bbox_to_anchor=(0.5, 0.105), ncol=2, frameon=False, labelspacing=0.5, columnspacing=1.5)
b.set_yscale("log"); b.set_ylim(1e-4, 0.24)
b.set_ylabel("mass inside r / critical mass inside D")
b.set_xlabel("distance from the centre of expansion (Gpc), log scale")
fig.text(0.5, 0.112, "All compact profiles reach their full mass well before D, so by Gauss's law the fit\nsees only the enclosed mass (their curves are flat across the visible sphere).\n"
       "The cored profile keeps growing ($M\\propto r$) and fits worse ($R^2$=0.830). The point + shell mass includes the shell.", fontsize=8.5, ha="center", va="bottom")
for ax in (a, b): ax.grid(alpha=0.25, which="major")
fig.subplots_adjust(bottom=0.19)
fig.savefig("../density_profiles.png", dpi=200, bbox_inches="tight")
im = Image.open("../density_profiles.png").convert("RGB"); f = 1700 / max(im.size)
if f < 1: im = im.resize((round(im.size[0] * f), round(im.size[1] * f)), Image.LANCZOS)
im.save("../density_profiles_web.jpg", quality=88)
for lab, kind, d, par in rows:
    print(lab, "|", "%.4f" % d["R2_mig"], "|", "%.2f" % d["D_Gpc"], "|", par, "|", "%.3g" % d["M_sun"], "|", "%.4g" % d.get("mc", 0))
