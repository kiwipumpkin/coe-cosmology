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
t = t.replace(r" Relativistic: $R_{\rm in}=D-1.040\,R$, $R_{\rm out}=D+1.171\,R$.", "")
t = t.replace(r"from the Sec.~3 equations", r"from the Sec.~2 equations")
s7[0] = t
out7 = []; intab = False
for l in s7[1:-1]:
    if l.startswith(r"\begin{center}"): l = l.replace("{lllllll}", "{lllll}").replace(r"\setlength{\tabcolsep}{2.6pt}", r"\setlength{\tabcolsep}{5pt}")
    elif l.startswith(" & \\multicolumn"): l = r" & \multicolumn{2}{c}{Newtonian, $\Omega_s$ free} & \multicolumn{2}{c}{Newtonian, $\Omega_s$=0.049}\\"
    elif l.startswith("Quantity &"): l = r"Quantity & term-by-term & Table 1 & term-by-term & Table 1\\\midrule"
    elif " & " in l and l.endswith("\\\\"):
        l = " & ".join(l[:-2].split(" & ")[:5]) + "\\\\"
    l = l.replace("$^b$", "$^a$")
    out7.append(l)
fn = r"{\footnotesize $^a$ Held fixed (that column's shell assumption). The harmonic solve weights $C_m,S_m$ equally, not the 71 points, so its $R^2$ can sit slightly below the point fit.}"
sec3 = [s7[0]] + out7 + [fn]
# section 4
A = json.load(open("accepted.json")); a = A["accepted"]; F, O = A["free"], A["omb"]
def e(x, n=2):
    m, ex = ("%.*e" % (n, x)).split("e"); return r"$%s{\times}10^{%d}$" % (m, int(ex))
rows = [
 (r"Ordinary-matter mass within $R$", e(a["Mb"]) + r" kg", "[1]$^a$", e(F["Mvis"]) + r" kg (%.0f$\times$)" % F["Mvis_ratio"], e(O["Mvis"]) + r" kg (%.2f$\times$)" % O["Mvis_ratio"]),
 (r"Total matter mass within $R$", e(a["Mm"]) + r" kg", "[1]$^a$", e(F["Mvis"]) + r" kg$^b$", e(O["Mvis"]) + r" kg$^b$"),
 (r"Ordinary-matter density $\rho_b$", e(a["rho_b"]) + r" kg\,m$^{-3}$", "[1]", e(F["rho_s"]) + r" kg\,m$^{-3}$", e(O["rho_s"]) + r" kg\,m$^{-3}$"),
 (r"$\Omega_b$", "%.4f" % a["Ob"], "[1]", r"%.3f ($=\Omega_s$)" % F["Os"], r"%.3f ($=\Omega_s$)" % O["Os"]),
 (r"Total matter $\Omega_m$", r"$%.4f\pm0.0073$" % a["Om"], "[1]", r"%.3f ($=\Omega_s$, no dark matter)" % F["Os"], r"%.3f ($=\Omega_s$, no dark matter)" % O["Os"]),
 (r"Critical density $\rho_c$", e(a["rhoc_planck"]) + r" kg\,m$^{-3}$ at $H_0$=67.36", "[1]", e(a["rhoc70"]) + r" kg\,m$^{-3}$ at $H_0$=70", e(a["rhoc70"]) + r" kg\,m$^{-3}$ at $H_0$=70"),
 (r"$H_0$ [km\,s$^{-1}$\,Mpc$^{-1}$]", r"$67.36\pm0.54$ (CMB);\newline $73.04\pm1.04$ (local)", "[1]; [2]", "70 (input)", "70 (input)"),
 (r"Earth (Sun) speed relative to the CMB", r"$369.82\pm0.11$ km\,s$^{-1}$", "[3]", r"$V_O=k D$ = %s km\,s$^{-1}$ = %.1f\,$c$" % (e(F["VO"]), F["VO_c"]), r"$V_O=kD$ = %s km\,s$^{-1}$ = %.1f\,$c$" % (e(O["VO"]), O["VO_c"])),
 (r"\quad ratio $V_O$/(CMB dipole speed)", "1", "", r"%s" % "{:,.0f}".format(F["ratio_cmb"]).replace(",", "{,}"), r"%s" % "{:,.0f}".format(O["ratio_cmb"]).replace(",", "{,}")),
 (r"Local Group speed relative to the CMB", r"$620\pm15$ km\,s$^{-1}$", "[3]", r"-- (not modelled)", r"-- (not modelled)"),
 (r"Bulk flow within $150\,h^{-1}$ Mpc", r"$395\pm29$ km\,s$^{-1}$", "[4]", r"%.1f km\,s$^{-1}$$^c$" % F["bulk"], r"%.1f km\,s$^{-1}$$^c$" % O["bulk"]),
 (r"Core mass $M_c$ (vs.\ ordinary mass within $R$)", "--", "", e(F["Mc"]) + r" kg (%.0f$\times$)" % F["Mc_ratio"], e(O["Mc"]) + r" kg (%.1f$\times$)" % O["Mc_ratio"]),
 (r"Mass inside $D$ (vs.\ ordinary mass within $R$)", "--", "", e(F["Mtot"]) + r" kg ({:,.0f}$\times$)".format(F["Mtot_ratio"]).replace(",", "{,}"), e(O["Mtot"]) + r" kg (%.0f$\times$)" % O["Mtot_ratio"]),
 (r"Age / $k^{-1}$", r"$13.797\pm0.023$ Gyr", "[1]", r"$k^{-1}$ = %.2f Gyr; integration age %.2f Gyr" % (a["invH70"], F["age"]), r"$k^{-1}$ = %.2f Gyr; integration age %.2f Gyr" % (a["invH70"], O["age"])),
]
body = "\n".join(" & ".join(r) + r"\\" for r in rows)
sec4 = r"""
\textbf{4. Model vs.\ generally accepted values.} Model columns use the Table~1 values of the two Newtonian fits (their saved results and code; $H_0=70$). $R=14.3$ Gpc is the comoving radius of the observable universe.
\begin{center}\footnotesize\setlength{\tabcolsep}{3.5pt}\renewcommand{\arraystretch}{1.12}\begin{tabular}{>{\raggedright\arraybackslash}p{4.4cm}>{\raggedright\arraybackslash}p{4.1cm}l>{\raggedright\arraybackslash}p{4.25cm}>{\raggedright\arraybackslash}p{4.25cm}}\toprule
Quantity & Accepted value & Source & Newt.\ $\Omega_s$ free (0.909) & Newt.\ $\Omega_s$=0.049 (0.893)\\\midrule
BODY
\bottomrule
\end{tabular}\end{center}
\vspace{-3pt}
{\footnotesize $^a$ Derived from [1]: $\Omega_b h^2=0.02237$, $\Omega_m=0.3153$, $H_0=67.36$, times the volume of the $R=14.3$ Gpc sphere. Ratios in brackets: model / accepted ordinary-matter mass within $R$. $^b$ The model has no dark matter, so total = ordinary. $^c$ Mean of $\Delta\vec V$ over a $150\,h^{-1}$ Mpc ball around Earth ($h=0.70$), from the code's exact expansion-rate profile $h(q)$ (second order in radius/$D$; the first-order term averages to zero). The sine factor is not a velocity and is not included.

\textbf{Sources.} [1] Planck Collaboration 2020, A\&A, 641, A6 (Planck 2018 VI, TT,TE,EE+lowE+lensing). [2] Riess, A.~G., et al. 2022, ApJL, 934, L7 (SH0ES). [3] Planck Collaboration 2020, A\&A, 641, A1 (Planck 2018 I; Local Group value after Kogut et al. 1993, ApJ, 419, 1). [4] Watkins, R., et al. 2023, MNRAS, 524, 1885 (CosmicFlows-4). Sedov--Taylor: Taylor, G.~I. 1950, Proc.\ R.\ Soc.\ Lond.\ A, 201, 159; Sedov, L.~I. 1959, Similarity and Dimensional Methods in Mechanics (Academic Press).}
""".replace("BODY", body)
def doc(parts): return "\n".join(head + parts + ["}", r"\end{document}"])
# A: title + sec1 + sec2 (image); B: sec3 + sec4 (image, with title-less header); full PDF: everything
open("../termmatch_A_v2.tex", "w").write(doc(sec1 + [""] + sec2))
hb = [l for l in head if not l.startswith(r"{\large") and not l.startswith(r"\rule")]
open("../termmatch_B_v2.tex", "w").write("\n".join(hb + sec3 + [""] + [sec4, "}", r"\end{document}"]))
open("../thinshell_termmatch_equations_table_v2.tex", "w").write(doc(sec1 + [""] + sec2 + [r"\vspace{6pt}", ""] + sec3 + [""] + [sec4]))
