import os
import json
J = json.load(open("universe_thin.json")); A, B, X = J["free"], J["newt_omb"], J["rel"]   # Newtonian (free), Newtonian Omega_s=0.049 forced, relativistic LTB (all free, z=0.3)
Msun = 1.989e30
def e(x, n=3):
    m, ex = ("%.*e" % (n, x)).split("e"); return r"$%s\times10^{%d}$" % (m, int(ex))
def c(x): return "{:,.0f}".format(x).replace(",", "{,}")
def sig(x): return "%.2g" % x if x < 0.01 else "%.3f" % x
def gp(r, D): return "%s Gpc (%.1f\\,$D$)" % (c(r), r / D) if r >= 10 else "%.2f Gpc (%.2f\\,$D$)" % (r, r / D)
def mn(x, f="%.4f"): return (f % x).replace("-", "$-$")
rows = [("C23", "$D/R$", "%.3f" % A["DR"], "%.3f" % B["DR"], "%.3f" % X["DR"]),
        ("D23", "$L_c/D$ (core size)", "0.0050 (bound)", "0.0050", "0.0050"),
        ("B29", "$R$", r"$4.4127\times10^{23}$ km", "same", "same"),
        ("", "", "(%.2f Gpc)" % A["R"], "", ""),
        ("B28", "$D$", e(A["Dkm"], 4) + " km", e(B["Dkm"], 4) + " km", e(X["Dkm"], 4) + " km"),
        ("", "", "(%.2f Gpc)" % A["D"], "(%.2f Gpc)" % B["D"], "(%.2f Gpc)" % X["D"]),
        ("B30", "$L_c$", e(A["Lckm"], 2) + " km", e(B["Lckm"], 2) + " km", e(X["Lckm"], 2) + " km"),
        ("", "", "(%.3f Gpc)" % A["Lc"], "(%.3f Gpc)" % B["Lc"], "(%.3f Gpc)" % X["Lc"]),
        ("B27", "$M(<D)$, total", e(A["Mtot"], 2) + " kg", e(B["Mtot"], 2) + " kg", e(X["Mtot"], 2) + " kg"),
        ("", "", e(A["Mtot"] / Msun, 2) + r" $M_\odot$", e(B["Mtot"] / Msun, 2) + r" $M_\odot$", e(X["Mtot"] / Msun, 2) + r" $M_\odot$"),
        ("", "$M_c$, core", e(A["Mc"], 2) + " kg", e(B["Mc"], 2) + " kg", e(X["Mc"], 2) + " kg"),
        ("", "", e(A["Mc"] / Msun, 2) + r" $M_\odot$", e(B["Mc"] / Msun, 2) + r" $M_\odot$", e(X["Mc"] / Msun, 2) + r" $M_\odot$"),
        ("", "$M_c/M_{\\rm crit}(D)$", "%.5f" % A["mcc"], "%.5f" % B["mcc"], sig(X["mcc"])),
        ("", "$M_{\\rm shell}(<D)$", e(A["Msh"], 2) + " kg", e(B["Msh"], 2) + " kg", e(X["Msh"], 2) + " kg"),
        ("", "$M_{\\rm shell}$, total", e(A["Mshtot"], 2) + " kg", e(B["Mshtot"], 2) + " kg", e(X["Mshtot"], 2) + " kg"),
        ("", "$M$ within $R$", e(A["Mvis"], 2) + " kg", e(B["Mvis"], 2) + " kg", e(X["Mvis"], 2) + " kg"),
        ("", "$M_{\\rm hole}$", e(A["Mhole"], 2) + " kg", e(B["Mhole"], 2) + " kg", e(X["Mhole"], 2) + " kg"),
        ("BA57", "Amp", mn(A["amp"]), mn(B["amp"]), mn(X["amp"])),
        ("BA59", "$\\varphi$", "%.3f rad" % A["phi"], "%.3f rad" % B["phi"], "%.3f rad" % X["phi"]),
        ("", "$\\ell_0$", "%.2f$^\\circ$" % A["l0"], "%.1f$^\\circ$" % B["l0"], "%.1f$^\\circ$" % X["l0"]),
        ("MID",),
        ("", "$\\Omega_s$", "%.3f" % A["Os"], "%.3f" % B["Os"], "%.3f" % X["Os"]),
        ("", "$\\rho_s$", e(A["Os"] * 9.2042e-27, 2) + r" kg\,m$^{-3}$", e(B["Os"] * 9.2042e-27, 2) + r" kg\,m$^{-3}$", e(X["Os"] * 9.2042e-27, 2) + r" kg\,m$^{-3}$"),
        ("", "$R_{\\rm in}$", "%.1f Gpc ($D-R$)" % A["Rin"], "%.1f Gpc" % B["Rin"], "%.1f Gpc" % X["Rin"]),
        ("", "$R_{\\rm out}$", "%.1f Gpc ($D+R$)" % A["Rout"], "%.1f Gpc" % B["Rout"], "%.1f Gpc" % X["Rout"]),
        ("", "$r_s$, core", "%.0f Gpc (%.2f\\,$D$)" % (A["rs_c"], A["rs_c"] / A["D"]), gp(B["rs_c"], B["D"]), gp(X["rs_c"], X["D"])),
        ("", "$r_s$, all mass $<D$", "%s Gpc (%.0f\\,$D$)" % (c(A["rs_D"]), A["rs_D"] / A["D"]), gp(B["rs_D"], B["D"]), gp(X["rs_D"], X["D"])),
        ("MID",),
        ("", "$(V_O/v_{\\rm esc})^2$ / $2E/c^2$", "%.3f (bound)" % A["vesc2"], "%.1f (unbound)" % B["vesc2"], "%+.0f (unbound)" % X["twoE"] if X["twoE"] > 0 else "%.2f (bound)" % X["twoE"]),
        ("", "Turnaround radius", "%.0f Gpc (%.1f\\,$D$)" % (A["rt"], A["rt"] / A["D"]), "none", "none"),
        ("", "Age $t_0$", "%.1f Gyr" % A["age"], "%.1f Gyr" % B["age"], "%.1f Gyr" % X["age"]),
        ("", "$\\delta$ across sphere", "0", "0", "0"),
        ("", "fit $R^2$ vs Migkas", "%.3f" % A["R2"], "%.3f (%.2f)" % (B["R2"], B["R2geo"]), "%.3f (%.2f)" % (X["R2"], X["R2ns"])),
        ("", "fit $R^2$ vs $A_F$", "%.3f" % A["R2F"], "%.3f" % B["R2F"], "%.3f" % X["R2F"])]
body = "\n".join(r"\midrule" if r[0] == "MID" else " & ".join(r) + r"\\" for r in rows)
tex = r"""\documentclass[10pt,letterpaper]{article}
\usepackage[margin=0.55in]{geometry}
\usepackage{amsmath,graphicx,booktabs,caption}
\captionsetup{font=footnotesize,labelfont=bf,skip=3pt}
\pagestyle{empty}
\begin{document}
\begin{center}{\large\bfseries Your universe at the best fit: thick ordinary-matter shell plus a CoE core (relativistic fit)}\end{center}
\vspace{-8pt}
\begin{figure}[h]\centering
\includegraphics[width=\textwidth]{universe_thin_drawing.png}
\caption{Three scales at the relativistic best fit; (b) and (c) are to scale in Gpc. (a) Hubble velocity vectors. In a Sedov--Taylor blast wave (Taylor 1950; Sedov 1959) $\vec V = k\,\times$(position from the CoE), so we move with $\vec V_O = k\vec D$ and a galaxy with $\vec V_P = k(\vec D+\vec r)$; then $\vec V_P-\vec V_O = k\vec r$ is Hubble's law, so $k = H$ and the CoE drops out. (b) Hubble bubble: ordinary-matter shell (grey, $\Omega_s=%(Os).3f$) from $R_{\rm in}=%(Rin).1f$ to $R_{\rm out}=%(Rout).1f$~Gpc around an empty interior and the CoE core; dashed circle $D=%(Dg).1f$~Gpc through Earth; green vision sphere $R=%(Rg).1f$~Gpc; Earth moves along $+x$ (anti-CoE, $\ell_0=%(l0).1f^\circ$) at $V_O=H_0D/c\approx%(vE).1f\,c$. (c) Cosmic shell: the core's Schwarzschild radius $r_s=%(rsc).0f$~Gpc ($%(rscD).1f\,D$); $r_s$ of all mass inside $D$ is %(rsD).0f~Gpc ($%(rsDD).1f\,D$), off scale. Blue arrows are velocities. The fit is exact relativistic (LTB), with ordinary matter only.}
\end{figure}
\vspace{-10pt}
\noindent\begin{minipage}[t]{0.40\textwidth}\vspace{0pt}\centering
\includegraphics[width=\textwidth]{universe_thin_Achart.png}
\captionof{figure}{Migkas points with error bars, $A_F$ and model A, normalised to the sky mean. Solid red: Newtonian thin shell with core ($\Omega_s=%(NOs).2f$, $R^2=%(NR2).3f$). Dotted green: Newtonian with $\Omega_s=%(BOs).3f$ forced, $R^2=%(BR2).3f$ (%(BR2g).2f geometry only). Dashed orange: relativistic (LTB) ordinary-matter shell with core ($\Omega_s=%(Os).3f$, all parameters free), sampled at the Migkas cluster distance ($z=0.3$, 1.19~Gpc), $R^2=%(R2).3f$ (%(R2ns).2f without the sine). Lower panel: percent differences from $A_F$.}
\end{minipage}\hfill
\begin{minipage}[t]{0.59\textwidth}\vspace{0pt}\centering\scriptsize\setlength{\tabcolsep}{1.5pt}\renewcommand{\arraystretch}{0.93}
\captionof{table}{Best-fit values in your sheet's cells. Newtonian: integrated model ($G$ real, from $0.01\,t_0$, fixed in space, $H_0=70$), thin shell $D-R$ to $D+R$. Newtonian $\Omega_s$=0.049: same model with the shell density forced to ordinary matter. Relativistic: exact LTB, ordinary-matter shell; $D$, $M_c$, $R_{\rm in}$, $R_{\rm out}$, $\ell_0$, Amp, $\varphi$ all free (shell $\ge R$ each side of Earth); $H=%(Hr).0f$/%(Hp).0f km\,s$^{-1}$\,Mpc$^{-1}$ along/across the axis. $R^2$ in brackets: geometry only (no sine). B27 is all mass inside $D$.}
\begin{tabular}{lllll}\toprule
Cell & Quantity & \begin{tabular}[b]{@{}l@{}}Newtonian,\\ $\Omega_s$ free\end{tabular} & \begin{tabular}[b]{@{}l@{}}Newtonian,\\ $\Omega_s$=0.049\end{tabular} & \begin{tabular}[b]{@{}l@{}}Relativistic,\\ $\Omega_s$=0.049,\\ $M_c$ free\end{tabular}\\\midrule
BODY
\bottomrule
\end{tabular}
\end{minipage}
\par\vspace{6pt}\noindent{\scriptsize\textbf{References.} Taylor, G.~I. 1950, Proc.\ R.\ Soc.\ Lond.\ A, 201, 159. Sedov, L.~I. 1959, Similarity and Dimensional Methods in Mechanics (New York: Academic Press).\par}
\end{document}
"""
v = dict(Dg=X["D"], Rg=X["R"], rscD=X["rs_c"] / X["D"], l0=X["l0"], Os=X["Os"], Rin=X["Rin"], Rout=X["Rout"], rsc=X["rs_c"], rsD=X["rs_D"], rsDD=X["rs_D"] / X["D"], vE=X["vE"],
         NOs=A["Os"], NR2=A["R2"], BOs=B["Os"], BR2=B["R2"], BR2g=B["R2geo"], R2=X["R2"], R2ns=X["R2ns"], Hr=X["H_r"], Hp=X["H_perp"])
open("../universe_bestfit_thinshell.tex", "w").write((tex % v).replace("BODY", body))
