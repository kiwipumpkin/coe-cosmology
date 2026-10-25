import json
J = json.load(open("universe_thin.json")); A, B = J["free"], J["baryons"]
Msun = 1.989e30
def e(x, n=3):
    m, ex = ("%.*e" % (n, x)).split("e"); return r"$%s\times10^{%d}$" % (m, int(ex))
def c(x): return "{:,.0f}".format(x).replace(",", "{,}")
def mn(x, f="%.4f"): return (f % x).replace("-", "$-$")
rows = [("C23", "$D/R$", "%.3f" % A["DR"], "%.3f" % B["DR"]),
        ("D23", "$L_c/D$ (core size)", "0.0050 (bound)", "0.0050 (bound)"),
        ("B29", "$R$", r"$4.4127\times10^{23}$ km", "same"),
        ("B28", "$D$", e(A["Dkm"], 4) + " km", e(B["Dkm"], 4) + " km"),
        ("", "", "(%.2f Gpc)" % A["D"], "(%.2f Gpc)" % B["D"]),
        ("B30", "$L_c$", e(A["Lckm"], 2) + " km", e(B["Lckm"], 2) + " km"),
        ("B27", "$M(<D)$, total", e(A["Mtot"]) + " kg", e(B["Mtot"]) + " kg"),
        ("", "", e(A["Mtot"] / Msun, 2) + r" $M_\odot$", e(B["Mtot"] / Msun, 2) + r" $M_\odot$"),
        ("", "$M_c$, core", e(A["Mc"], 2) + " kg", e(B["Mc"], 2) + " kg"),
        ("", "$M_c/M_{\\rm crit}(D)$", "%.5f" % A["mcc"], "%.5f" % B["mcc"]),
        ("", "$M_{\\rm shell}(<D)$", e(A["Msh"], 2) + " kg", e(B["Msh"], 2) + " kg"),
        ("", "$M_{\\rm hole}$", e(A["Mhole"], 2) + " kg", e(B["Mhole"], 2) + " kg"),
        ("BA57", "Amp", mn(A["amp"]), mn(B["amp"])),
        ("BA59", "$\\varphi$", "%.3f rad" % A["phi"], "%.3f rad" % B["phi"]),
        ("", "$\\ell_0$", "%.2f$^\\circ$" % A["l0"], "%.2f$^\\circ$" % B["l0"]),
        ("MID",),
        ("", "$\\Omega_s$", "%.3f" % A["Os"], "0.05 (baryons)"),
        ("", "$\\rho_s$", e(A["Os"] * 9.2042e-27, 2) + r" kg\,m$^{-3}$", e(0.05 * 9.2042e-27, 2) + r" kg\,m$^{-3}$"),
        ("", "$R_{\\rm in}=D-R$", "%.1f Gpc" % A["Rin"], "%.1f Gpc" % B["Rin"]),
        ("", "$R_{\\rm out}=D+R$", "%.1f Gpc" % A["Rout"], "%.1f Gpc" % B["Rout"]),
        ("", "$r_s$, core", "%.0f Gpc (%.2f\\,$D$)" % (A["rs_c"], A["rs_c"] / A["D"]), "%.1f Gpc (%.2f\\,$D$)" % (B["rs_c"], B["rs_c"] / B["D"])),
        ("", "$r_s$, all mass $<D$", "%s Gpc (%.0f\\,$D$)" % (c(A["rs_D"]), A["rs_D"] / A["D"]), "%.0f Gpc (%.1f\\,$D$)" % (B["rs_D"], B["rs_D"] / B["D"])),
        ("MID",),
        ("", "$(v/v_{\\rm esc})^2$", "%.3f (bound)" % A["vesc2"], "%.2f (unbound)" % B["vesc2"]),
        ("", "Turnaround radius", "%.0f Gpc (%.1f\\,$D$)" % (A["rt"], A["rt"] / A["D"]), "none"),
        ("", "Age $t_0$", "%.1f Gyr" % A["age"], "%.1f Gyr" % B["age"]),
        ("", "$\\delta$ across sphere", "0", "0"),
        ("", "$R^2$ vs Migkas", "%.3f" % A["R2"], "%.3f" % B["R2"]),
        ("", "$R^2$ vs $A_F$", "%.3f" % A["R2F"], "%.3f" % B["R2F"])]
body = "\n".join(r"\midrule" if r[0] == "MID" else " & ".join(r) + r"\\" for r in rows)
tex = r"""\documentclass[10pt,letterpaper]{article}
\usepackage[margin=0.55in]{geometry}
\usepackage{amsmath,graphicx,booktabs,caption}
\captionsetup{font=footnotesize,labelfont=bf,skip=3pt}
\pagestyle{empty}
\begin{document}
\begin{center}{\large\bfseries Your universe at the best fit: thin constant-density shell hugging our visible sphere plus a CoE core}\end{center}
\vspace{-8pt}
\begin{figure}[h]\centering
\includegraphics[width=0.93\textwidth]{universe_thin_drawing.png}
\caption{To scale, in Gpc with equal axes, at the best fit with $\Omega_s$ free. The grey shell has constant density $\Omega_s=%(Os).2f$, about 20 times the baryon density, and runs from $R_{\rm in}=D-R=%(Rin).1f$~Gpc to $R_{\rm out}=D+R=%(Rout).1f$~Gpc. The interior is empty apart from the core at the CoE. (a) The dashed blue circle of radius $D$ runs through Earth, the green disc is the vision shell of radius $R$, and sky directions are Galactic longitude with $\theta=\ell-\ell_0$. The red arc is the Schwarzschild radius $r_s=2GM/c^2$ of the core alone, %(rsc).0f~Gpc, which already lies beyond Earth. (b) Zoomed out, showing the shell, the core's $r_s$ and Earth's turnaround radius (%(rt).0f~Gpc). The $r_s$ of all the mass inside $D$ is %(rsD)s~Gpc ($%(rsDD).0f\,D$), off the scale. Blue arrows show the outward velocity: Earth's recession away from the CoE in (a), and the outward, decelerating motion of the shell toward the turnaround circle in (b). The model is Newtonian, and Earth moves outward at $%(vE).1f\,c$. It is bound, with $(v/v_{\rm esc})^2=%(v2).2f$, and turns around at $%(rtD).1f\,D$ after %(tt).0f ages. The mass inside that radius has $r_s\approx%(rst)s$~Gpc, so turnaround lies inside $r_s$, and in general relativity the region would recollapse.}
\end{figure}
\vspace{-10pt}
\noindent\begin{minipage}[t]{0.50\textwidth}\vspace{0pt}\centering
\includegraphics[width=\textwidth]{universe_thin_Achart.png}
\captionof{figure}{Migkas points with error bars, $A_F$, and your model A for the thin shell with a core, with $\Omega_s$ free (solid) and at baryon density (dashed), normalised to the sky mean. The lower panel shows data and models as percent differences from $A_F$. The sources start at 1\%% of their present scale and coast through the empty interior, where only the core pulls, entering the shell only in the last $\sim R/D$ of the age. The core therefore sets the pattern ($\ell_0\approx13^\circ$) even though it is lighter than the mass missing from the hole.}
\end{minipage}\hfill
\begin{minipage}[t]{0.48\textwidth}\vspace{0pt}\centering\scriptsize\setlength{\tabcolsep}{2.2pt}\renewcommand{\arraystretch}{0.93}
\captionof{table}{Best-fit values in your sheet's cells: integrated model ($G$ real, from $0.01\,t_0$, fixed in space, $H_0=70$), shell from $D-R$ to $D+R$, $\ell_0$ free. B27 is all gravitating mass inside $D$. The second column is the same geometry with the shell at baryon density.}
\begin{tabular}{llll}\toprule
Cell & Quantity & $\Omega_s$ free & $\Omega_s=0.05$\\\midrule
BODY
\bottomrule
\end{tabular}
\end{minipage}
\end{document}
"""
v = dict(Os=A["Os"], Rin=A["Rin"], Rout=A["Rout"], rsc=A["rs_c"], rt=A["rt"], rsD=c(A["rs_D"]), rsDD=A["rs_D"] / A["D"], vE=A["vE"], v2=A["vesc2"],
         rtD=A["rt"] / A["D"], tt=A["t_turn"], rst=c(round(A["rs_t"], -3)))
open("../universe_bestfit_thinshell.tex", "w").write((tex % v).replace("BODY", body))
