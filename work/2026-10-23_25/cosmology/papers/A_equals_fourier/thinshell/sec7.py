import json, numpy as np
T = json.load(open("/workspace/cosmology/papers/spreadsheet_model/bubble/universe_thin.json")); SI = json.load(open("solve_inputs.json")); TT = json.load(open("taui_test.json"))
MS = 1.989e30
def e(x, n=2):
    m, ex = ("%.*e" % (n, x)).split("e"); return r"$%s{\times}10^{%d}$" % (m, int(ex))
def mm(x, f="%.3f"): return "$" + (f % x) + "$"
A, B = SI["free"], SI["omb"]; Rk = sorted(SI["rel"].keys(), key=float); R0 = SI["rel"][Rk[0]]
ta, tb, tc = T["free"], T["newt_omb"], T["rel"]
mrel = [SI["rel"][k]["mcc"] for k in Rk]
nf = r"\multicolumn{1}{c}{not fixed}"
rows = [
 ("$D/R$", "%.3f" % A["DR"], "%.3f" % ta["DR"], "%.3f" % B["DR"], "%.3f" % tb["DR"], nf, "%.3f" % tc["DR"]),
 ("$D$", e(A["Dkm"], 3) + " km", e(ta["Dkm"], 3) + " km", e(B["Dkm"], 3) + " km", e(tb["Dkm"], 3) + " km", r"any $D$ (tried 33--330)", e(tc["Dkm"], 3) + " km"),
 ("", "(%.1f Gpc)" % A["D"], "(%.1f Gpc)" % ta["D"], "(%.1f Gpc)" % B["D"], "(%.1f Gpc)" % tb["D"], "", "(%.1f Gpc)" % tc["D"]),
 ("$R_{\\rm in}$", "%.1f Gpc" % A["R_in"], "%.1f Gpc" % ta["Rin"], "%.1f Gpc" % B["R_in"], "%.1f Gpc" % tb["Rin"], "$D-14.9$ Gpc", "%.1f Gpc" % tc["Rin"]),
 ("$R_{\\rm out}$", "%.1f Gpc" % A["R_out"], "%.1f Gpc" % ta["Rout"], "%.1f Gpc" % B["R_out"], "%.1f Gpc" % tb["Rout"], "$D+16.7$ Gpc", "%.1f Gpc" % tc["Rout"]),
 ("$M_c$", e(A["Mc"]) + " kg", e(ta["Mc"]) + " kg", e(B["Mc"]) + " kg", e(tb["Mc"]) + " kg", r"$\propto D^3$ (%s kg at 32.9 Gpc)" % e(R0["Mc"], 1), e(tc["Mc"]) + " kg"),
 ("", e(A["Mc"] / MS) + r" $M_\odot$", e(ta["Mc"] / MS) + r" $M_\odot$", e(B["Mc"] / MS) + r" $M_\odot$", e(tb["Mc"] / MS) + r" $M_\odot$", "", e(tc["Mc"] / MS) + r" $M_\odot$"),
 ("$M_c/M_{\\rm crit}(D)$", "%.5f" % A["mcc"], "%.5f" % ta["mcc"], "%.5f" % B["mcc"], "%.5f" % tb["mcc"], r"\textbf{%.4f--%.4f}$^a$" % (min(mrel), max(mrel)), "%.4f" % tc["mcc"]),
 ("$M(<D)$, total", e(A["Mtot"]) + " kg", e(ta["Mtot"]) + " kg", e(B["Mtot"]) + " kg", e(tb["Mtot"]) + " kg", r"$\propto D^3$", e(tc["Mtot"]) + " kg"),
 ("$M_{\\rm shell}$, total", e(A["Mshtot"]) + " kg", e(ta["Mshtot"]) + " kg", e(B["Mshtot"]) + " kg", e(tb["Mshtot"]) + " kg", r"$\propto D^2$", e(tc["Mshtot"]) + " kg"),
 ("$M$ within $R$", e(A["Mvis"]) + " kg", e(ta["Mvis"]) + " kg", e(B["Mvis"]) + " kg", e(tb["Mvis"]) + " kg", e(R0["Mvis"]) + " kg", e(tc["Mvis"]) + " kg"),
 ("$\\Omega_s$", "%.3f$^b$" % A["Os"], "%.3f" % ta["Os"], "0.049$^b$", "0.049", "0.049$^b$", "0.049"),
 ("$\\ell_0$", "%.1f$^\\circ$" % A["l0"], "%.1f$^\\circ$" % ta["l0"], "%.1f$^\\circ$" % B["l0"], "%.1f$^\\circ$" % tb["l0"], "%.1f$^\\circ$" % R0["l0"], "%.1f$^\\circ$" % tc["l0"]),
 ("Amp", mm(A["amp"], "%.4f"), mm(ta["amp"], "%.4f"), mm(B["amp"], "%.4f"), mm(tb["amp"], "%.4f"), mm(R0["amp"], "%.4f"), mm(tc["amp"], "%.4f")),
 ("$\\varphi$", "%.3f rad" % A["phi"], "%.3f rad" % ta["phi"], "%.3f rad" % B["phi"], "%.3f rad" % tb["phi"], "%.3f rad" % R0["phi"], "%.3f rad" % tc["phi"]),
 ("$R^2$ vs Migkas", "%.3f" % A["R2mig"], "%.3f" % ta["R2"], "%.3f" % B["R2mig"], "%.3f" % tb["R2"], "%.3f$^c$" % R0["R2mig"], "%.3f" % tc["R2"]),
 ("$R^2$ vs $A_F$", "%.3f" % A["R2F"], "%.3f" % ta["R2F"], "%.3f" % B["R2F"], "%.3f" % tb["R2F"], "%.3f$^c$" % R0["R2F"], "%.3f" % tc["R2F"]),
]
body = "\n".join(" & ".join(r) + r"\\" for r in rows)
t = {k: TT[k] for k in TT}
sec = r"""
\textbf{7. Input parameters implied by the term-by-term solution vs.\ Table~1} of universe\_bestfit\_thinshell.pdf. Each column keeps that fit's shell assumption. Newtonian $\Omega_s$ free: $\Omega_s=%(Osf).3f$, thin shell $R_{\rm in}=D-R$, $R_{\rm out}=D+R$. Newtonian $\Omega_s$=0.049: $R_{\rm in}=D-%(kino).3f\,R$, $R_{\rm out}=D+%(kouto).3f\,R$. Relativistic: $R_{\rm in}=D-%(kinr).3f\,R$, $R_{\rm out}=D+%(koutr).3f\,R$. The unknowns are solved from the Sec.~3 equations ($2\ell$ phase $\to\ell_0$, $2\ell$ amplitude $\to\beta\to M_c/M_{\rm crit}$, $1\ell\to$ Amp, $\varphi$, $3\ell$ cos $\to\varepsilon Q\to D$; $C_m,S_m$, $m=1$--3, by least squares). The $a_m(\beta,\gamma)$ come from the code's exact $h(q)$, not the linear formula. Masses use $H_0=70$, $\rho_s=\Omega_s\rho_c$.
\begin{center}\footnotesize\setlength{\tabcolsep}{2.6pt}\begin{tabular}{lllllll}\toprule
 & \multicolumn{2}{c}{Newtonian, $\Omega_s$ free} & \multicolumn{2}{c}{Newtonian, $\Omega_s$=0.049} & \multicolumn{2}{c}{Relativistic, all free}\\
Quantity & term-by-term & Table 1 & term-by-term & Table 1 & term-by-term (Newt., $r$=1.19) & Table 1 (LTB)\\\midrule
BODY
\bottomrule
\end{tabular}\end{center}
\vspace{-3pt}
{\footnotesize $^a$ Newtonian code, $\tau_i=0.01$, sampling $r=1.19$ Gpc: the $r/D$-odd terms vanish, so $D$ is not fixed; any $D$ works with $M_c/M_{\rm crit}(D)\approx0.0017$ (%(m1).4f, %(m2).4f, %(m3).4f at $D$ = 33, 100, 330 Gpc). The exact LTB curve instead gives $M_c/M_{\rm crit}(D)=0.061$--0.083 for any $D$ (Sec.~6). $^b$ Held fixed (that column's shell assumption). $^c$ At $D=32.9$ Gpc. The harmonic solve weights $C_m,S_m$ equally, not the 71 points, so its $R^2$ can sit slightly below the point fit.}

The \textbf{0.909 fit} is recovered: $D$ 92.6 vs 93.7 Gpc, the same $M_c/M_{\rm crit}$, and $\ell_0$, Amp, $\varphi$ to 3\%%. For \textbf{$\Omega_s$=0.049} the harmonic solve gives the same $M_c/M_{\rm crit}$ (%(mo).5f vs %(mto).5f), $\ell_0$, Amp and $\varphi$, but a smaller $D$ (%(Do).0f vs %(Dto).0f Gpc): the harmonic weighting pulls toward a stronger $3\ell$ term.
\textbf{$\tau_i$ test} (Chris's hypothesis that the code pins $M_c/M_{\rm crit}\cdot(1/\tau_i-1)$). The $\Omega_s$=0.049 fit was redone at other start times, with $D$, $M_c$, both shell edges, $\ell_0$, Amp and $\varphi$ free:
\begin{center}\footnotesize\begin{tabular}{lcccc}\toprule
$\tau_i$ & $M_c/M_{\rm crit}(D)$ & $M_c/M_{\rm crit}\cdot(1/\tau_i-1)$ & $D$ [Gpc] & $R^2$\\\midrule
0.1 & %(m01).4f & %(x01).3f & %(D01).1f & %(r01).3f\\
0.01 (Table 1) & %(m001).4f & %(x001).3f & %(D001).1f & %(r001).3f\\
0.001 & %(m0001).4f & %(x0001).3f & %(D0001).1f & %(r0001).3f\\\bottomrule
\end{tabular}\end{center}
\vspace{-3pt}
The hypothesis holds only partly. From $\tau_i$ = 0.01 to 0.1, $M_c/M_{\rm crit}$ rises %(f1).1f$\times$ (11$\times$ predicted). From 0.01 to 0.001 it falls only %(f2).0f\%% (10$\times$ predicted), because the core's own size ($L_c=0.005\,D$) cuts off the early pull. So the Newtonian $M_c/M_{\rm crit}$ is set by the code's start time and core size, and $M_c/M_{\rm crit}\cdot(1/\tau_i-1)$ is not constant (0.10--1.5). The $\sim$50$\times$ gap to the relativistic 0.083 is mostly this regularisation: at $\tau_i=0.1$ it shrinks to %(g01).0f$\times$. The LTB value (single bang time, no $\tau_i$) is the physically calibrated one.
"""
v = dict(Osf=A["Os"], kino=1.028, kouto=1.475, kinr=1.040, koutr=1.171, m1=mrel[0], m2=mrel[1], m3=mrel[2], mo=B["mcc"], mto=tb["mcc"], Do=B["D"], Dto=tb["D"],
         m01=t["0.1"]["full"]["mcc"], x01=t["0.1"]["full"]["mcc_x"], D01=t["0.1"]["full"]["D"], r01=t["0.1"]["full"]["R2"],
         m001=t["0.01"]["full"]["mcc"], x001=t["0.01"]["full"]["mcc_x"], D001=t["0.01"]["full"]["D"], r001=t["0.01"]["full"]["R2"],
         m0001=t["0.001"]["full"]["mcc"], x0001=t["0.001"]["full"]["mcc_x"], D0001=t["0.001"]["full"]["D"], r0001=t["0.001"]["full"]["R2"],
         f1=t["0.1"]["full"]["mcc"] / t["0.01"]["full"]["mcc"], f2=100 * (1 - t["0.001"]["full"]["mcc"] / t["0.01"]["full"]["mcc"]), g01=0.0835 / t["0.1"]["full"]["mcc"])
SEC7 = (sec % v).replace("BODY", body)
