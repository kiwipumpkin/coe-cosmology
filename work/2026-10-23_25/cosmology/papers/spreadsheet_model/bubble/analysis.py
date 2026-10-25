import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
import prof_core as P
import ptshell as S
from bubble_core import AFfour, ML, MA, MLO, MHI, TH1, r2, A_from, shifted
from scipy.integrate import cumulative_trapezoid
plt.rcParams.update({"font.size": 8.5, "font.family": "serif"})
B = json.load(open("bubble_results.json")); PR = json.load(open("prof_results.json")); PS = json.load(open("ptshell_results.json"))
def psi(kind, d):
    """central potential depth of the CoE mass alone, in units G*M_crit/D: mc * int_0^inf m(r)/r^2 dr"""
    if kind == "ptshell":
        a = dict(d); a["Os"] = 0.0; a["rin"] = 0.5; a["rout"] = 0.6
        core = np.exp(-P.RG / d["Lc"]); m = cumulative_trapezoid(4*np.pi*P.RG**2*core, P.RG, initial=0); m /= np.interp(1, P.RG, m); mc = d["mcc"]
    else:
        _, m = P.tables(kind, d, 100.0); mc = d["mc"]
    rg = P.RG[1:]; return float(mc * np.trapezoid(m[1:] / rg**2, rg))
def curve(kind, d, norm="mean"):
    if kind == "ptshell":
        _, _, m1 = S.build(d); rr = P.run("ptshell", d, m1, 0.0)
    else: rr = P.run(kind, d, d["mc"], d["Ob"])
    af = rr["af"]; af = af / af.mean() if norm == "mean" else af**2 / np.mean(af**2)
    l0 = d.get("l0", 0.0); afs = shifted(af, l0) if l0 else af
    return A_from(afs, TH1 - l0, d["amp"], d["phi"])
rows = [("exp", "e$^{-r/L}$ (main)", "LCDM-equiv.", B["mean|c=0.02|a|Ob=0.3|mig_l0"], "mean"),
        ("exp", "e$^{-r/L}$", "baryons only", PR["exp|c=0.02|mig_l0|mean|Ob=0.05"], "mean"),
        ("pl", "power law $(1+r/L)^{-n}$", "LCDM-equiv.", PR["pl|c=0.02|mig_l0|mean"], "mean"),
        ("gauss", "Gaussian core", "LCDM-equiv.", PR["gauss|c=0.02|mig_l0|mean"], "mean"),
        ("cored", "cored $1/(1+r^2/L^2)$", "LCDM-equiv.", PR["cored|c=0.02|mig_l0|mean"], "mean"),
        ("shell", "Gaussian shell", "LCDM-equiv.", PR["shell|c=0.02|mig_l0|mean"], "mean"),
        ("st53", r"Sedov--Taylor $\gamma=5/3$", "LCDM-equiv.", PR["st53|c=0.02|mig_l0|mean"], "mean"),
        ("st43", r"Sedov--Taylor $\gamma=4/3$", "LCDM-equiv.", PR["st43|c=0.02|mig_l0|mean"], "mean"),
        ("ptshell", r"point + ST shell, $\Omega_s$ free", "free", PS["inside|Os=free|mig_l0"], "mean"),
        ("ptshell", r"point + ST shell, $\Omega_s=0.3$", "LCDM-equiv.", PS.get("inside|Os=0.3|mig_l0|c=0.02"), "mean"),
        ("ptshell", r"point + ST shell, $\Omega_s=0.05$", "baryons only", PS.get("inside|Os=0.05|mig_l0|c=0.02"), "mean"),
        ("ptshell", r"point + ST shell, front free", "free", PS["rout_free|Os=free|mig_l0"], "mean"),
        ("exp", r"e$^{-r/L}$, $A\propto(H/H_0)^2$", "LCDM-equiv.", PR["exp|c=0.02|mig_l0|mean2"], "mean2")]
json.dump([(k, lab, mat, d and dict(d, psi=psi(k, d))) for k, lab, mat, d, n in rows], open("analysis_rows.json", "w"), indent=1)
for k, lab, mat, d, n in rows:
    if d: print(lab, mat, "psi=%.3f" % psi(k, d), "R2 %.4f" % d["R2_mig"])
# comparison A(l) plot
fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, 4.8), sharex=True, gridspec_kw={"height_ratios": [2.0, 1], "hspace": 0.06})
Af = AFfour(TH1)
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.5", ecolor="0.8", lw=0.8, label="Migkas Fig. 6")
a.plot(TH1, Af, color="k", lw=1.6, label=r"$A_F$")
sel = [0, 1, 2, 4, 5, 6, 8, 11]; cols = ["C3", "C1", "C0", "C9", "C2", "C4", "C5", "C7"]; lss = ["-", "--", ":", "-.", "-", "--", "-", ":"]
for i, c, ls in zip(sel, cols, lss):
    k, lab, mat, d, n = rows[i]
    if not d: continue
    A = curve(k, d, n)
    lab2 = "%s, %s ($R^2$=%.3f)" % (lab, mat, d["R2_mig"])
    a.plot(TH1, A, color=c, ls=ls, lw=1.2, label=lab2); b.plot(TH1, 100 * (A - Af) / Af, color=c, ls=ls, lw=1.1)
b.plot(ML, 100 * (MA - AFfour(ML)) / AFfour(ML), "o", ms=2.3, color="0.5")
b.axhline(0, color="k", lw=0.6); b.set_ylabel(r"$100\,(\cdot-A_F)/A_F$ [%]", fontsize=8)
a.set_ylabel(r"$A/A_{\rm all}$"); a.set_ylim(0.62, 1.2); a.legend(fontsize=6.3, loc="lower center", ncol=2, frameon=False)
b.set_xlabel(r"Galactic longitude $\ell$ (deg)"); b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45))
fig.savefig("../uniform_bubble_profiles.png", dpi=200, bbox_inches="tight")
im = Image.open("../uniform_bubble_profiles.png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../uniform_bubble_profiles_web.jpg", quality=88)
