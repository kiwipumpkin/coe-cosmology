import json
c = json.load(open("universe_real.json"))
def e(x, n=3):
    m, ex = ("%.*e" % (n, x)).split("e"); return r"$%s\times10^{%d}$" % (m, int(ex))
Msun = 1.989e30; G = 6.674e-11
tex = r"""\documentclass[10pt,letterpaper]{article}
\usepackage[margin=0.55in]{geometry}
\usepackage{amsmath,graphicx,booktabs,caption}
\captionsetup{font=footnotesize,labelfont=bf,skip=3pt}
\pagestyle{empty}
\begin{document}
\begin{center}{\large\bfseries Your universe at the best real-matter fit: point core plus a baryon-only Sedov--Taylor shell}\end{center}
\vspace{-8pt}
\begin{figure}[h]\centering
\includegraphics[width=0.93\textwidth]{universe_real_drawing.png}
\caption{To scale, in Gpc with equal axes. Grey is the shell of constant density $\Omega_s=0.05$ (baryons only), running from $R_{\rm in}=%(Rin).1f$~Gpc to the front at $R_{\rm out}=%(Rout).0f$~Gpc; white inside $R_{\rm in}$ is the hole, and the black dot is the point core (scale %(Lc).2f~Gpc). (a) The dashed blue circle of radius $D$ runs through Earth, the green disc is the vision shell of radius $R$, and sky directions are Galactic longitude with $\theta=\ell-\ell_0$. The red circle is the Schwarzschild radius $r_s=2GM/c^2$ of the core alone, %(rsc).1f~Gpc ($%(rscD).2f\,D$). (b) Zoomed out to show the front, $r_s$ of all the mass inside $D$ (%(rsD).0f~Gpc, $%(rsDD).1f\,D$), and Earth's turnaround radius (%(rt).0f~Gpc). The mass inside $D$ has $r_s>D$. The model is Newtonian, and Earth moves outward at $%(vE).1f\,c$, so it has not collapsed. The fit is bound, with $(v/v_{\rm esc})^2=%(vesc2).2f$, and Earth turns around at $%(rtD).1f\,D$ after %(tt).0f ages. The mass inside that radius has $r_s\approx%(rst)s$~Gpc, so turnaround falls inside $r_s$, and in general relativity the region would recollapse. The same holds for any uniform sphere with $\Omega=0.05$ larger than about 19~Gpc, including the matter of a standard cosmology.}
\end{figure}
\vspace{-10pt}
\noindent\begin{minipage}[t]{0.50\textwidth}\vspace{0pt}\centering
\includegraphics[width=\textwidth]{universe_real_Achart.png}
\captionof{figure}{Migkas points with error bars, $A_F$, and your model A for the point core plus baryon-only shell, normalised to the sky mean. The lower panel shows data and model as percent differences from $A_F$.}
\end{minipage}\hfill
\begin{minipage}[t]{0.48\textwidth}\vspace{0pt}\centering\scriptsize\setlength{\tabcolsep}{2.5pt}\renewcommand{\arraystretch}{0.95}
\captionof{table}{Best-fit values in your sheet's cells, integrated model ($G$ real, from $0.01\,t_0$, fixed in space, $H_0=70$), with density contrast across the visible sphere $\le2\%%$ and $\ell_0$ free. B27 is all the gravitating mass inside $D$. The CoE mass is unlabelled; only the shell is set to the baryon density.}
\begin{tabular}{lll}\toprule
Cell & Quantity & Best fit\\\midrule
C23 & $D/R$ & %(DR).3f\\
D23 & $L_c/D$ (core size) & 0.0050 (lower bound)\\
B29 & $R$ & $4.4127\times10^{23}$ km (%(R).2f Gpc)\\
B28 & $D$ & %(Dkm)s km (%(D).2f Gpc)\\
B30 & $L_c$ & %(Lckm)s km (%(Lc).3f Gpc)\\
B27 & $M(<D)$, total & %(Mtot)s kg = %(Mtots)s $M_\odot$\\
 & $M_c$, point core & %(Mc)s kg = %(Mcs)s $M_\odot$\\
 & $M_c/M_{\rm crit}(D)$ & %(mcc).5f\\
 & $M_{\rm shell}(<D)$ & %(Msh)s kg\\
BA57 & Amp & %(amp)s\\
BA59 & $\varphi$ & %(phi).3f rad\\
 & $\ell_0$ & %(l0).2f$^\circ$\\\midrule
 & $\Omega_s$ (baryons only) & 0.05\\
 & $\rho_s$ & %(rhos)s kg m$^{-3}$\\
 & $R_{\rm in}$ & %(Rin).2f Gpc ($%(rinD).3f\,D$)\\
 & $R_{\rm out}$ & %(Rout).1f Gpc ($%(routD).2f\,D$)\\
 & $R_b$ (matter edge) & $=R_{\rm out}$\\
 & $r_s$, point core & %(rsc).2f Gpc ($%(rscD).3f\,D$)\\
 & $r_s$, all mass $<D$ & %(rsD).0f Gpc ($%(rsDD).2f\,D$)\\\midrule
 & $(v/v_{\rm esc})^2$ & %(vesc2).3f (bound)\\
 & Age $t_0$ & %(age).1f Gyr\\
 & $\delta$ across sphere & 0 (shell uniform there)\\
 & $R^2$ vs 71 Migkas points & %(R2).3f\\
 & $R^2$ vs $A_F$ & %(R2F).3f\\\bottomrule
\end{tabular}
\end{minipage}
\end{document}
"""
v = dict(Rin=c["Rin_Gpc"], Rout=c["Rout_Gpc"], Lc=c["Lc_Gpc"], rsc=c["rs_c"], rscD=c["rs_c_D"], rsD=c["rs_D"], rsDD=c["rs_D_D"], rt=c["rt_Gpc"],
         vE=c["vE_over_c"], vesc2=c["vesc2"], rtD=c["rt_D"], tt=c["t_turn"], rst="{:,}".format(int(round(c['rs_t_Gpc'], -2))).replace(',', '{,}'), DR=c["DR"], R=c["R_Gpc"], Dkm=e(c["Dkm"], 4), D=c["D_Gpc"],
         Lckm=e(c["Lckm"]), Mtot=e(c["Mtot_kg"]), Mtots=e(c["Mtot_kg"] / Msun, 2), Mc=e(c["Mc_kg"]), Mcs=e(c["Mc_sun"], 2), mcc=c["Mc_kg"] / (c["Mtot_kg"] / 0.051635474886684206),
         Msh=e(c["Msh_kg"]), amp=str('%.4f' % c['amp']).replace('-', '$-$'), phi=c["phi"], l0=c["l0"], rhos=e(c["rhos"], 2), rinD=c["Rin_Gpc"] / c["D_Gpc"], routD=c["Rout_Gpc"] / c["D_Gpc"],
         age=c["age"], R2=c["R2"], R2F=c["R2F"])
out = tex % v
p2 = open("page2_real.tex").read().replace("ALTROWS", open("../alt_rows.tex").read().strip())
out = out.replace("\\end{document}", p2 + "\\end{document}")
open("../universe_bestfit_real.tex", "w").write(out)
