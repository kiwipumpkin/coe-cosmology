import os
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
        ("", "", "(%.2f Gpc)" % A["R"], ""),
        ("B28", "$D$", e(A["Dkm"], 4) + " km", e(B["Dkm"], 4) + " km"),
        ("", "", "(%.2f Gpc)" % A["D"], "(%.2f Gpc)" % B["D"]),
        ("B30", "$L_c$", e(A["Lckm"], 2) + " km", e(B["Lckm"], 2) + " km"),
        ("", "", "(%.3f Gpc)" % A["Lc"], "(%.3f Gpc)" % B["Lc"]),
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
        ("", "$\\Omega_s$", "%.3f" % A["Os"], "0.05 (ordinary matter)"),
        ("", "$\\rho_s$", e(A["Os"] * 9.2042e-27, 2) + r" kg\,m$^{-3}$", e(0.05 * 9.2042e-27, 2) + r" kg\,m$^{-3}$"),
        ("", "$R_{\\rm in}=D-R$", "%.1f Gpc" % A["Rin"], "%.1f Gpc" % B["Rin"]),
        ("", "$R_{\\rm out}=D+R$", "%.1f Gpc" % A["Rout"], "%.1f Gpc" % B["Rout"]),
        ("", "$r_s$, core", "%.0f Gpc (%.2f\\,$D$)" % (A["rs_c"], A["rs_c"] / A["D"]), "%.1f Gpc (%.2f\\,$D$)" % (B["rs_c"], B["rs_c"] / B["D"])),
        ("", "$r_s$, all mass $<D$", "%s Gpc (%.0f\\,$D$)" % (c(A["rs_D"]), A["rs_D"] / A["D"]), "%.0f Gpc (%.1f\\,$D$)" % (B["rs_D"], B["rs_D"] / B["D"])),
        ("MID",),
        ("", "$(V_O/v_{\\rm esc})^2$", "%.3f (bound)" % A["vesc2"], "%.2f (unbound)" % B["vesc2"]),
        ("", "Turnaround radius", "%.0f Gpc (%.1f\\,$D$)" % (A["rt"], A["rt"] / A["D"]), "none"),
        ("", "Age $t_0$", "%.1f Gyr" % A["age"], "%.1f Gyr" % B["age"]),
        ("", "$\\delta$ across sphere", "0", "0"),
        ("", "fit $R^2$ vs Migkas", "%.3f" % A["R2"], "%.3f" % B["R2"]),
        ("", "fit $R^2$ vs $A_F$", "%.3f" % A["R2F"], "%.3f" % B["R2F"])]
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
\includegraphics[width=\textwidth]{universe_thin_drawing.png}
\caption{Three scales, left to right, at the best fit with $\Omega_s$ free; (b) and (c) are to scale in Gpc with equal axes. (a) Hubble velocity vectors. In a Sedov--Taylor blast wave (Taylor 1950; Sedov 1959) the velocity grows in proportion to distance from the centre, $\vec V = k\,\times$(position from the CoE) for some constant $k$. So we, at position $\vec D$, move with $\vec V_O = k\vec D$, and a galaxy at $\vec D+\vec r$ moves with $\vec V_P = k(\vec D+\vec r)$. Their relative velocity is $\vec V_H = \vec V_P-\vec V_O = k\vec r$, which is Hubble's law, so $k = H$. It depends only on the separation $\vec r$, so the CoE drops out and every observer sees the same isotropic Hubble flow. (b) Hubble bubble: the grey shell has constant density $\Omega_s=%(Os).2f$ (about 20 times the ordinary-matter density) from $R_{\rm in}=D-R=%(Rin).1f$ to $R_{\rm out}=D+R=%(Rout).1f$~Gpc, the interior is empty apart from the CoE core, the dashed circle of radius $D=%(Dg).1f$~Gpc runs through Earth, and the green disc is the vision sphere of radius $R=%(Rg).1f$~Gpc. The red arc is the core's own Schwarzschild radius $r_s=2GM/c^2=%(rsc).0f$~Gpc ($%(rscD).2f\,D$), already beyond Earth. Earth moves outward along $+x$, the anti-CoE direction $\ell_0=%(l0).1f^\circ$ (sky directions are Galactic longitude, $\theta=\ell-\ell_0$). (c) Cosmic shell: the shell, the core's $r_s$ and Earth's turnaround radius (%(rt).0f~Gpc, $%(rtD).1f\,D$); the $r_s$ of all mass inside $D$ is %(rsD)s~Gpc ($%(rsDD).0f\,D$), off scale. Blue arrows are velocities in all panels. In this Newtonian model Earth moves outward at $V_O=%(vE1).1f\,c$, is bound ($(V_O/v_{\rm esc})^2=%(v2).2f$) and turns around after %(tt).0f ages; the mass inside has $r_s\approx%(rst)s$~Gpc, so turnaround lies inside $r_s$ and in general relativity the region would recollapse.}
\end{figure}
\vspace{-10pt}
\noindent\begin{minipage}[t]{0.50\textwidth}\vspace{0pt}\centering
\includegraphics[width=\textwidth]{universe_thin_Achart.png}
\captionof{figure}{Migkas points with error bars, $A_F$, and your model A for the thin shell with a core, normalised to the sky mean. Solid: shell density free ($\Omega_s=%(Os).2f$). Dashed: the same model with the shell density fixed at 0.05 of critical, about the density of ordinary matter. The lower panel shows data and models as percent differences from $A_F$. The sources spend almost all their history coasting through the empty interior, where only the CoE core pulls on them, so the core sets the pattern ($\ell_0\approx13^\circ$).}
\end{minipage}\hfill
\begin{minipage}[t]{0.48\textwidth}\vspace{0pt}\centering\scriptsize\setlength{\tabcolsep}{2.2pt}\renewcommand{\arraystretch}{0.93}
\captionof{table}{Best-fit values in your sheet's cells: integrated model ($G$ real, from $0.01\,t_0$, fixed in space, $H_0=70$), shell from $D-R$ to $D+R$, $\ell_0$ free. B27 is all gravitating mass inside $D$. The second column is the same geometry with the shell at ordinary-matter density ($\Omega_s=0.05$).}
\begin{tabular}{llll}\toprule
Cell & Quantity & $\Omega_s$ free & $\Omega_s=0.05$\\\midrule
BODY
\bottomrule
\end{tabular}
\end{minipage}
\par\vspace{6pt}\noindent{\scriptsize\textbf{References.} Taylor, G.~I. 1950, Proc.\ R.\ Soc.\ Lond.\ A, 201, 159. Sedov, L.~I. 1959, Similarity and Dimensional Methods in Mechanics (New York: Academic Press).\par}
\end{document}
"""
v = dict(Dg=A['D'], Rg=A['R'], rscD=A['rs_c']/A['D'], l0=A['l0'], vE1=A['vE'], arrowsent=("Thin blue arrows show outward velocity: Earth's recession from the CoE (b) and the shell's decelerating expansion toward turnaround (c). " if os.environ.get('ARROWS','1')=='1' else ''), Os=A["Os"], Rin=A["Rin"], Rout=A["Rout"], rsc=A["rs_c"], rt=A["rt"], rsD=c(A["rs_D"]), rsDD=A["rs_D"] / A["D"], vE=A["vE"], v2=A["vesc2"],
         rtD=A["rt"] / A["D"], tt=A["t_turn"], rst=c(round(A["rs_t"], -3)))
open("../universe_bestfit_thinshell.tex", "w").write((tex % v).replace("BODY", body))
