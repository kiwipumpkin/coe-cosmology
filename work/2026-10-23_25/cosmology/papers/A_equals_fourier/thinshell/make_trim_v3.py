"""Trimmed sheet: Sec 1 (model), Sec 2 (term-by-term solution, was 3), Sec 3 (inputs vs Table 1, was 7; Newtonian columns only),
Sec 4 (model vs generally accepted values). Wording unchanged except numbering / cross-references / removed relativistic parts."""
import json, re
L = open("../thinshell_termmatch.tex").read().split("\n")
def idx(pfx): return next(i for i, l in enumerate(L) if l.startswith(pfx))
i1 = idx(r"\textbf{1. Model"); i2 = idx(r"\textbf{2. $A_F$ rotated"); i3 = idx(r"\textbf{3. Term-by-term"); i4 = idx(r"\textbf{4. Cross-check"); i7 = idx(r"\textbf{7. Input")
head = [r"\documentclass[10pt]{article}", r"\usepackage[letterpaper,margin=0.45in]{geometry}", r"\usepackage{array}"] + L[2:9]
sec1 = L[i1:i2]
sec2 = L[i3:i4]
# section 2: keep the table only (as in termmatch_A_equations.jpg): stop after the first \end{tabular}\end{center}
e = next(i for i, l in enumerate(sec2) if l.startswith(r"\end{tabular}\end{center}")); sec2 = sec2[:e + 1]
sec2[0] = sec2[0].replace(r"\textbf{3. Term-by-term solution}", r"\textbf{2. Term-by-term solution}")
# section 3 (was 7)
s7 = L[i7:]; e7 = next(i for i, l in enumerate(s7) if l.startswith(r"{\footnotesize $^a$")); s7 = s7[:e7 + 1]
t = s7[0]
t = t.replace(r"\textbf{7. Input parameters", r"\textbf{3. Input parameters")
t = t.replace(r"from the Sec.~3 equations", r"from the Sec.~2 equations")
s7[0] = t
out7 = []; intab = False
for l in s7[1:-1]:
    if l.startswith(r"\begin{center}"): l = l.replace("{lllllll}", "{llllll}").replace(r"\setlength{\tabcolsep}{2.6pt}", r"\setlength{\tabcolsep}{4.5pt}")
    elif l.startswith(" & \\multicolumn"): l = r" & \multicolumn{2}{c}{Newtonian, $\Omega_s$ free} & \multicolumn{2}{c}{Newtonian, $\Omega_s$=0.049} & Relativistic (LTB)\\"
    elif l.startswith("Quantity &"): l = r"Quantity & term-by-term & Table 1 & term-by-term & Table 1 & Table 1\\\midrule"
    elif " & " in l and l.endswith("\\\\"):
        c = l[:-2].split(" & "); l = " & ".join(c[:5] + [c[6]]) + "\\\\"
    l = l.replace("$^b$", "$^a$")
    out7.append(l)
fn = r"{\footnotesize $^a$ Held fixed (that column's shell assumption). The harmonic solve weights $C_m,S_m$ equally, not the 71 points, so its $R^2$ can sit slightly below the point fit.}"
sec3 = [s7[0]] + out7 + [fn]
# section 4 (v3): all rate-dependent model rows in units of h70 = <H>/(70 km/s/Mpc), <H> = model sky-mean rate
A = json.load(open("accepted.json")); a = A["accepted"]; Hm = json.load(open("hmean.json"))
BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble/"
TT = json.load(open(BUB + "universe_thin.json")); RL = json.load(open(BUB + "universe_rel_mcfree.json"))
c = 299792.458; G = 6.674e-11; MPC = 3.0857e19
def rhoc(H): h = H / MPC; return 3 * h * h / (8 * 3.141592653589793 * G)
M = {}
for key, tk in (("free", "free"), ("omb", "newt_omb")):
    F = A[key]; m = Hm[key]; t = TT[tk]; Ht = m["H_E/k (Earth radial rate, code Ht)"]; hs = m["R=14.3"]["sine"]; hg = m["R=14.3"]["geom"]
    k = 70 / hs; f = k / (70 / Ht)                                   # k at <H> = 70; f = rescale from the v2 convention (Earth rate = 70)
    M[key] = dict(Hk_s=hs, Hk_g=hg, Hk_loc=m["local (0.01 Gpc)"]["geom"], k=k, t0=977.8 / k, rho_s=F["rho_s"] * f * f, Mvis=F["Mvis"] * f * f, Mc=F["Mc"] * f * f, Mtot=F["Mtot"] * f * f,
                  VO=m["vE"] * k * t["D"] * 1000, bulk=F["bulk"] * f, Os=t["Os"], Os_H=t["Os"] * (Ht / hs)**2)
g = 70 / ((RL["H_perp"] + RL["H_r"]) / 2)        # LTB: in-plane sky mean of the local rate = (H_perp + H_r)/2 = 72.04 in the fit's units (3D mean = 70)
T3 = TT["rel"]; kL = 977.8 / RL["age_Gyr"]
M["rel"] = dict(Hk_s=None, Hk_loc=((RL["H_perp"] + RL["H_r"]) / 2) / kL, k=kL * g, t0=RL["age_Gyr"] / g, rho_s=0.049 * rhoc(70) * g * g, Mvis=T3["Mvis"] * g * g, Mc=T3["Mc"] * g * g,
                Mtot=T3["Mtot"] * g * g, VO=T3["vE"] * g * c, bulk=None, Os=0.049, Os_H=0.049 * g * g)
MG = json.load(open("migkas_anchor.json"))["summary"]
for key in M:                                      # anchor: <H> = Migkas 2021 all-combined (joint) fit; rescale from <H> = 70
    x = M[key]; h = MG[key]["joint"] / 70; x["Hm"] = MG[key]["joint"]; x["Hm_err"] = MG[key]["joint_err"]; x["Hm_mean"] = MG[key]["mean"]; x["Hm_sd"] = MG[key]["sd"]; x["h"] = h
    x["k"] *= h; x["t0"] /= h; x["VO"] *= h
    for q in ("rho_s", "Mvis", "Mc", "Mtot"): x[q] *= h * h
    if x["bulk"] is not None: x["bulk"] *= h
json.dump(M, open("section4_v3_values.json", "w"), indent=1, default=float)
def e(x, n=2):
    m_, ex = ("%.*e" % (n, x)).split("e"); return r"$%s{\times}10^{%d}$" % (m_, int(ex))
F, O, R = M["free"], M["omb"], M["rel"]
def trio(fn): return [fn(F, "free"), fn(O, "omb"), fn(R, "rel")]
ratio = lambda v: "{:,.0f}".format(v).replace(",", "{,}")
H2 = ""; H1 = ""
rows = [
 (r"$H_0$ [km\,s$^{-1}$\,Mpc$^{-1}$]", r"$67.36\pm0.54$ (CMB) [1]; $73.04\pm1.04$ (local) [2]",
  r"$\langle H\rangle$ = %.1f$\pm$%.1f$^a$ (= %.3f\,$k$)" % (F["Hm"], F["Hm_err"], F["Hk_s"]),
  r"$\langle H\rangle$ = %.1f$\pm$%.1f$^a$ (= %.3f\,$k$)" % (O["Hm"], O["Hm_err"], O["Hk_s"]),
  r"$\langle H\rangle$ = %.1f$\pm$%.1f$^{a,d}$ (= %.3f/$t_0$)" % (R["Hm"], R["Hm_err"], R["Hk_loc"])),
 (r"$k$ (LTB: $1/t_0$) [km\,s$^{-1}$\,Mpc$^{-1}$]", "--", r"%.1f%s" % (F["k"], H1), r"%.1f%s" % (O["k"], H1), r"%.1f%s" % (R["k"], H1)),
 (r"Age / $k^{-1}$", r"$13.797\pm0.023$ Gyr [1]", r"$k^{-1}=t_0$ = %.2f\ Gyr" % F["t0"], r"$k^{-1}=t_0$ = %.2f\ Gyr" % O["t0"], r"$t_0$ = %.2f\ Gyr" % R["t0"]),
 (r"Ordinary-matter mass within $R$", e(a["Mb"]) + r" kg [1]$^b$", e(F["Mvis"]) + " kg" + H2, e(O["Mvis"]) + " kg" + H2, e(R["Mvis"]) + " kg" + H2),
 (r"Total matter mass within $R$", e(a["Mm"]) + r" kg [1]$^b$", r"same (no dark matter)", r"same (no dark matter)", r"same (no dark matter)"),
 (r"Ordinary-matter density $\rho_b$", e(a["rho_b"]) + r" kg\,m$^{-3}$ [1]", e(F["rho_s"]) + r" kg\,m$^{-3}$" + H2, e(O["rho_s"]) + r" kg\,m$^{-3}$" + H2, e(R["rho_s"]) + r" kg\,m$^{-3}$" + H2),
 (r"$\Omega_b$ (= $\Omega_m$ in the model)", r"$\Omega_b=%.4f$, $\Omega_m=0.3153\pm0.0073$ [1]" % a["Ob"], r"%.3f (%.3f vs $\langle H\rangle$)" % (F["Os"], F["Os_H"]), r"%.3f (%.3f vs $\langle H\rangle$)" % (O["Os"], O["Os_H"]), r"%.3f (%.3f vs $\langle H\rangle$)" % (R["Os"], R["Os_H"])),
 (r"Critical density $\rho_c$ at $\langle H\rangle$", e(a["rhoc_planck"]) + r" kg\,m$^{-3}$ ($H_0$=67.36) [1]", e(rhoc(F["Hm"])) + r" kg\,m$^{-3}$", e(rhoc(O["Hm"])) + r" kg\,m$^{-3}$", e(rhoc(R["Hm"])) + r" kg\,m$^{-3}$"),
 (r"Earth (Sun) speed vs.\ CMB", r"$369.82\pm0.11$ km\,s$^{-1}$ [3]", r"$V_O$ = %s km\,s$^{-1}$ = %.1f\,$c$" % (e(F["VO"]), F["VO"] / c), r"$V_O$ = %s km\,s$^{-1}$ = %.1f\,$c$" % (e(O["VO"]), O["VO"] / c), r"$V_O$ = %s km\,s$^{-1}$ = %.1f\,$c$" % (e(R["VO"]), R["VO"] / c)),
 (r"\quad ratio $V_O$/(CMB dipole speed)", "1", ratio(F["VO"] / 369.82) + H1, ratio(O["VO"] / 369.82) + H1, ratio(R["VO"] / 369.82) + H1),
 (r"Local Group speed vs.\ CMB", r"$620\pm15$ km\,s$^{-1}$ [3]", "-- (not modelled)", "-- (not modelled)", "-- (not modelled)"),
 (r"Bulk flow within $150\,h^{-1}$ Mpc", r"$395\pm29$ km\,s$^{-1}$ [4]", r"%.1f km\,s$^{-1}$%s$^c$" % (F["bulk"], H1), r"%.2f km\,s$^{-1}$%s$^c$" % (O["bulk"], H1), "not computed"),
 (r"Core mass $M_c$ (vs.\ accepted $M_b$ in $R$)", "--", e(F["Mc"]) + " kg" + H2 + r" (%.0f$\times$)" % (F["Mc"] / a["Mb"]), e(O["Mc"]) + " kg" + H2 + r" (%.1f$\times$)" % (O["Mc"] / a["Mb"]), e(R["Mc"]) + " kg" + H2 + r" (%.0f$\times$)" % (R["Mc"] / a["Mb"])),
 (r"Mass inside $D$ (same comparison)", "--", e(F["Mtot"]) + " kg" + H2 + r" (%s$\times$)" % ratio(F["Mtot"] / a["Mb"]), e(O["Mtot"]) + " kg" + H2 + r" (%.0f$\times$)" % (O["Mtot"] / a["Mb"]), e(R["Mtot"]) + " kg" + H2 + r" (%.0f$\times$)" % (R["Mtot"] / a["Mb"])),
]
body = "\n".join(" & ".join(r) + r"\\" for r in rows)
sec4 = r"""
\textbf{4. Model vs.\ generally accepted values.} Model columns use the Table~1 fits. $\langle H\rangle$ is the model's sky-mean expansion rate (mean of $a_f=|\Delta\vec V|/$separation over the code's directions at $R$, with the fitted sine). The fits fix only its shape, so its absolute value is solved from Migkas et al.'s own $H_0$ values by fitting $H(\ell)=\langle H\rangle\,A(\ell)/\langle A\rangle$ to them (anchor: their all-combined joint result$^a$). All rate-dependent rows follow from $\langle H\rangle$: rates, $k$, $V_O\propto\langle H\rangle$; masses, densities $\propto\langle H\rangle^2$; $t_0\propto1/\langle H\rangle$. Lengths use $R=14.3$ Gpc.
\begin{center}\footnotesize\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.12}\begin{tabular}{>{\raggedright\arraybackslash}p{4.45cm}>{\raggedright\arraybackslash}p{3.2cm}>{\raggedright\arraybackslash}p{3.48cm}>{\raggedright\arraybackslash}p{3.48cm}>{\raggedright\arraybackslash}p{3.48cm}}\toprule
Quantity & Accepted value [source] & Newt.\ $\Omega_s$ free (0.909) & Newt.\ $\Omega_s$=0.049 (0.893) & Relativistic LTB (0.885)\\\midrule
BODY
\bottomrule
\end{tabular}\end{center}
\vspace{-3pt}
{\footnotesize $^a$ Weighted least-squares fit of $\langle H\rangle A(\ell)/\langle A\rangle$ to Migkas et al.\ (2021, A\&A, 649, A151) all-combined $H_0$: $66.5\pm1.0$ toward $\ell=273^\circ$ and $72.8\pm0.6$ for the rest of the sky (model value at the cone centre, and the model mean outside a $75^\circ$ cone). Across all 11 Migkas 2020/2021 choices (each relation, sample and joint result) $\langle H\rangle$ = PRANGE. Migkas's absolute $H_0$ values rest on their fiducial calibration ($H_0=70$ assumed to calibrate the scaling relations), so this anchor inherits it. $^b$ Derived from [1]: $\Omega_b h^2=0.02237$, $\Omega_m=0.3153$, $H_0=67.36$, times the volume of the $R=14.3$ Gpc sphere; ratios in brackets are model / this accepted ordinary-matter mass. $^c$ Mean of $\Delta\vec V$ over a $150\,h^{-1}$ Mpc ball around Earth, from the code's exact expansion-rate profile (second order in radius/$D$); the sine factor is not a velocity. $^d$ LTB: mean of the local expansion rate over the in-plane directions, $(H_\perp+H_r)/2$; the fit's $z=0.3$ curve is normalised and carries no absolute rate. LTB values rescaled from the fit (which set the 3D-mean local rate to 70) with the same power laws (approximate in GR; the fitted valley shows the shape does not depend on $HD/c$).

\textbf{Sources.} [1] Planck Collaboration 2020, A\&A, 641, A6 (Planck 2018 VI, TT,TE,EE+lowE+lensing). [2] Riess, A.~G., et al. 2022, ApJL, 934, L7 (SH0ES). [3] Planck Collaboration 2020, A\&A, 641, A1 (Planck 2018 I; Local Group value after Kogut et al. 1993, ApJ, 419, 1). [4] Watkins, R., et al. 2023, MNRAS, 524, 1885 (CosmicFlows-4). Sedov--Taylor: Taylor, G.~I. 1950, Proc.\ R.\ Soc.\ Lond.\ A, 201, 159; Sedov, L.~I. 1959, Similarity and Dimensional Methods in Mechanics (Academic Press).}
""".replace("BODY", body).replace("PRANGE", r"%.1f$\pm$%.1f / %.1f$\pm$%.1f / %.1f$\pm$%.1f (mean$\pm$s.d., the three fits)" % (F["Hm_mean"], F["Hm_sd"], O["Hm_mean"], O["Hm_sd"], R["Hm_mean"], R["Hm_sd"]))
def doc(parts): return "\n".join(head + parts + ["}", r"\end{document}"])
# A: title + sec1 + sec2 (image); B: sec3 + sec4 (image, with title-less header); full PDF: everything
open("../termmatch_A_v3.tex", "w").write(doc(sec1 + [""] + sec2))
hb = [l for l in head if not l.startswith(r"{\large") and not l.startswith(r"\rule")]
open("../termmatch_B_v3.tex", "w").write("\n".join(hb + sec3 + [""] + [sec4, "}", r"\end{document}"]))
open("../thinshell_termmatch_equations_table_v3.tex", "w").write(doc(sec1 + [""] + sec2 + [r"\vspace{6pt}", ""] + sec3 + [""] + [sec4]))
