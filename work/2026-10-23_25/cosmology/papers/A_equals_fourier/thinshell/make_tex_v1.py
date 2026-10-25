import json, numpy as np
J = json.load(open("termmatch.json")); V = json.load(open("relvalley.json")); F = J["fits"]; H = J["hand"]; S = J["solve"]
def ap(h, m):
    a, b = h[str(m)]; return np.hypot(a, b), (np.degrees(np.arctan2(b, a)) / m) % (360 / m)
def f3(x): return "%.3f" % x
def f4(x): return "%.4f" % x
hp = H["pointmass"]; C = hp["C"]
AF = J["AF_amp"]
# fits' full-curve harmonics
def fullrow(k):
    f = F[k]; out = []
    for m in (1, 2, 3):
        A, l = ap(f["full_num"], m); out.append("%.4f @ %.0f$^\\circ$" % (A, l))
    return out
Fr = {k: fullrow(k) for k in F}
tex = r"""\documentclass[10pt]{article}
\usepackage[paperwidth=8.5in,paperheight=12.9in,margin=0.45in]{geometry}
\usepackage{amsmath,amssymb,booktabs}
\pagestyle{empty}\setlength{\parindent}{0pt}\setlength{\parskip}{2pt}
\renewcommand{\arraystretch}{1.05}
\begin{document}
{\large\bfseries Current model (ordinary-matter shell + CoE core) vs.\ $A_F$, term by term}\hfill{\footnotesize worked sheet, 27 Sep 2026}\\[-2pt]
\rule{\textwidth}{0.4pt}
{\small
\textbf{1. Model and where each term comes from.} Sedov--Taylor flow $\vec V=k\,\vec x$ (position from the CoE; Taylor 1950; Sedov 1959), decelerated by a point core $M_c$ at the CoE and a uniform shell $\Omega_s$ between $R_{\rm in}$ and $R_{\rm out}$ (fixed in space; Earth and all sources inside it). Spherical symmetry makes the field radial, $\vec V=h(q)\,\vec q$, with $q=s/D$ ($D$ = CoE--Earth distance). Linear perturbation of the code's integration (from $\tau_i=0.01$, $H t\simeq1$):
\[
\frac{h(q)}{k}\simeq1-\tfrac12\Big[\underbrace{\tfrac{M_c}{M_{\rm crit}(D)}\big(\tfrac1{\tau_i}-1\big)q^{-3}}_{\text{core tidal term}}+\Omega_s\Big(\underbrace{\tfrac12}_{\text{uniform: }m=0\text{ only}}\ \underbrace{-\tfrac32 r_{\rm in}^2q^{-2}+r_{\rm in}^3q^{-3}}_{\text{hole }(r_{\rm in}=R_{\rm in}/D)}\Big)\Big],\qquad
\beta\equiv\frac{h'}{h}\Big|_{q=1},\ \ \gamma\equiv\frac{h''}{h}\Big|_{q=1}.
\]
Chris's $A$: $a_f=|\Delta\vec V|/{\rm separation}$ with sources at $R=14.3$~Gpc, $\varepsilon\equiv R/D$, $\mu=\cos b_0\cos\psi$, $\psi=\ell-\ell_0$ ($\ell_0$ = anti-CoE direction). Exactly $a_f\propto|h(q_s)(\hat x+\varepsilon\hat n)-h(1)\hat x|/\varepsilon$; to first order in $\varepsilon$ (the higher orders are the $r/D$ terms):
\[
a_f^2\propto1+B\mu^2+\varepsilon\,(P\mu+Q\mu^3),\qquad B=2\beta+\beta^2,\ \ P=\beta(3+\beta),\ \ Q=(1+\beta)(\gamma+\beta)-2\beta .
\]
Take the square root, then use $\cos^2\psi=\tfrac12(1+\cos2\psi)$, $\cos^3\psi=\tfrac14(3\cos\psi+\cos3\psi)$, $\cos^4\psi=\tfrac38+\tfrac12\cos2\psi+\tfrac18\cos4\psi$, $\cos^5\psi=\tfrac1{16}(10\cos\psi+5\cos3\psi+\cos5\psi)$ ($b_0=0$ as in the code):
\[
\frac{a_f}{\langle a_f\rangle}=1+a_1\cos\psi+a_2\cos2\psi+a_3\cos3\psi+\dots,\quad
a_2=\frac{\frac B4-\frac{B^2}{16}}{m_0},\ \
a_1=\frac{\varepsilon}{2m_0}\Big[P+\tfrac34\big(Q-\tfrac{PB}2\big)-\tfrac{5QB}{16}\Big],\ \
a_3=\frac{\varepsilon}{2m_0}\Big[\tfrac14\big(Q-\tfrac{PB}2\big)-\tfrac{5QB}{32}\Big],
\]
with $m_0=1+\frac B4-\frac{3B^2}{64}$. For a pure point mass ($\gamma=-4\beta$, $\beta\to3K$) this is the earlier line-of-sight form $K[(3\mu^2-1)+\varepsilon(\frac92\mu-\frac{15}2\mu^3)]$. The geometry gives only $\cos m\psi$ terms: even harmonics set by $\beta$, odd ones by $\varepsilon Q$. \emph{Rotation} (the sine factor, $A=a_f[1-{\rm Amp}\sin(\psi+\varphi)]$) adds
$C_1\!=\!a_1-{\rm Amp}\sin\varphi(1+\tfrac{a_2}2)$, $S_1\!=\!-{\rm Amp}\cos\varphi(1-\tfrac{a_2}2)$, $C_3\!=\!a_3-\tfrac{{\rm Amp}\,a_2}2\sin\varphi$, $S_3\!=\!-\tfrac{{\rm Amp}\,a_2}2\cos\varphi$ ($C_m,S_m$: $\cos m\psi,\sin m\psi$ parts).
\textbf{Latitude:} the code evaluates the model in the Galactic plane only (one direction per degree, no $|b|$ cut, no sector average), so the match is done the same way. A tilt $b_0$ would give $B\to B\cos^2b_0$ and $\varepsilon$-terms $\propto\cos b_0,\cos^3b_0$, degenerate with $\beta$ and $\varepsilon Q$; therefore $b_0=0$.

\textbf{2. $A_F$ rotated onto the model axis} ($\ell_0=14.96^\circ$, the $2\ell$ peak), $c_m\cos m\ell+s_m\sin m\ell=C_m\cos m\psi+S_m\sin m\psi$:
\begin{center}\begin{tabular}{lcccccccc}\toprule
 & $A_0$ & $A_1$, $\ell_1$ & $A_2$, $\ell_2$ & $A_3$, $\ell_3$ & $C_1$ & $S_1$ & $C_3$ & $S_3$\\\midrule
$A_F$ & 0.9965 & %(A1).4f, %(l1).1f$^\circ$ & %(A2).4f, %(l2).2f$^\circ$ & %(A3).4f, %(l3).1f$^\circ$ & $%(C1)s$ & $%(S1)s$ & $%(C3)s$ & $%(S3)s$\\\bottomrule
\end{tabular}\end{center}
\vspace{-4pt}
\textbf{3. Term-by-term solution} (Newtonian, series above, point-core shape $\gamma=-4\beta$; exact least-squares in brackets):
\begin{center}\footnotesize\setlength{\tabcolsep}{3pt}\begin{tabular}{llll}\toprule
Term & $A_F$ coefficient & Analytic expression & Solved value\\\midrule
$2\ell$ phase & $\ell_2=14.96^\circ$ & $\ell_0$ (axis; $\beta>0$: fast along it) & $\ell_0=14.96^\circ$ \ (%(xl0).1f$^\circ$)\\
$2\ell$ amplitude & $C_2=0.0734$ & $a_2(\beta)$ & $\beta=%(beta).3f$ \ (%(xbeta).3f)\\
$1\ell$ (rotation) & $C_1=%(C1)s$, $S_1=%(S1)s$ & $C_1,S_1$ of Sec.~1 & Amp $=%(amp)s$, $\varphi=%(phi).2f$ \ (%(xamp)s, %(xphi).2f)\\
$3\ell$ cos & $C_3=%(C3)s$ & $a_3(\varepsilon Q)-\frac{{\rm Amp}\,a_2}2\sin\varphi$ & $\varepsilon Q=%(epsQ)s$ $\Rightarrow$ $\varepsilon=%(eps).3f$, $D=%(D).0f$~Gpc \ (%(xD).0f)\\
$3\ell$ sin & $S_3=%(S3)s$ & $-\frac{{\rm Amp}\,a_2}2\cos\varphi$ (check) & predicted $%(S3p).4f$: \textbf{not produced}\\
$m=0$ & $A_0=0.9965$ & $N[1-\frac12{\rm Amp}\,a_1\sin\varphi]$ & $N=%(N).4f$\\\midrule
mass & -- & $\beta\simeq\frac{3K}{1-K}$, $K=\frac12\frac{M_c}{M_{\rm crit}(D)}(\frac1{\tau_i}-1)$ & $M_c/M_{\rm crit}(D)\approx%(mcc_lin).4f$ (linear; $\approx$0.0017 nonlinear)\\\bottomrule
\end{tabular}\end{center}
\vspace{-4pt}
The $3\ell$ step fixes only $\varepsilon Q$: with the $0.909$-fit shape ($\gamma/\beta=%(gbf).1f$) it gives $D=%(Df).0f$~Gpc (exact %(xDf).0f), with the $\Omega_s=0.049$ shape ($\gamma/\beta=%(gbo).1f$) $D=%(Do).0f$~Gpc (exact %(xDo).0f). Every Newtonian solution reaches $R^2=%(xR2F).3f$ vs $A_F$ ($%(xR2m).3f$ vs Migkas).

\textbf{4. Cross-check: analytic series vs exact numeric projection of Chris's code} (geometric part, $a_m$ of $a_f/\langle a_f\rangle$; the code's own $\beta,\gamma$ from $h(q)$):
\begin{center}\footnotesize\setlength{\tabcolsep}{3pt}\begin{tabular}{lcccccccc}\toprule
Fit & $D$ [Gpc] & $\varepsilon$ & $\beta$ ($\beta_{\rm lin}$) & $\gamma$ & $a_1$ series / exact & $a_2$ series / exact & $a_3$ series / exact & $a_4$ exact\\\midrule
%(rowsX)s\bottomrule
\end{tabular}\end{center}
\vspace{-4pt}
The linear formula for $h(q)$ overestimates $\beta$ by about 1.5--1.7 (the core pull at $\tau\sim\tau_i$ is not small), so $\beta,\gamma$ are read from the code; with them the series agrees with the exact projection to $\lesssim$30\%% ($a_2$ to 6\%%). Using the exact radial field $h(q)$ reproduces the code to $10^{-4}$. For the relativistic fit $A\propto D_L^{-1}$ at $z=0.3$ (sampling $r=1.19$~Gpc); the local shear gives $a_2\simeq\frac{\beta/2}{1+\beta/2}$ with $\beta=H_r/H_\perp-1$, and $r/D=0.036$ would make $a_1,a_3\approx0$, but the exact LTB curve has $a_3=-0.018$ (light-cone/shell-edge structure, not $r/D$).

\textbf{5. The three numerical fits, harmonic by harmonic} (full model with sine; amplitude @ peak longitude):
\begin{center}\footnotesize\setlength{\tabcolsep}{4pt}\begin{tabular}{lllllll}\toprule
 & $R^2$ Migkas & $R^2$ vs $A_F$ & $A_1$ & $A_2$ & $A_3$ & $M_c/M_{\rm crit}(D)$, $D$\\\midrule
$A_F$ & -- & 1 & %(A1).4f @ %(l1).0f$^\circ$ & %(A2).4f @ %(l2).0f$^\circ$ & %(A3).4f @ %(l3).0f$^\circ$ & --\\
Newtonian, $\Omega_s$ free & %(R2free).3f & %(R2Ffree).3f & %(Ffree1)s & %(Ffree2)s & %(Ffree3)s & %(mccfree).4f, 93.7 Gpc\\
Newtonian, $\Omega_s$=0.049 & %(R2omb).3f & %(R2Fomb).3f & %(Fomb1)s & %(Fomb2)s & %(Fomb3)s & %(mccomb).4f, 49.2 Gpc\\
Relativistic, all free & %(R2rel).3f & %(R2Frel).3f & %(Frel1)s & %(Frel2)s & %(Frel3)s & %(mccrel).3f, 32.9 Gpc\\\bottomrule
\end{tabular}\end{center}
\vspace{-4pt}
All three get $A_1$ (via the rotation) and $A_2$ (via $\beta$) right; they differ only in $A_3$, which ranks them (0.909 $>$ 0.893 $>$ 0.885). The hand solution ($\beta\approx0.16$, Amp $\approx-0.07$, $\varphi\approx1.0$, $\ell_0\approx13$--15$^\circ$, $\varepsilon Q\approx-0.31$) matches the $\Omega_s$=0.049 fit ($\beta=%(bomb).3f$, Amp $=%(aomb)s$, $\varphi=%(pomb).2f$, $\ell_0=13^\circ$, $D=49$ vs 41--45~Gpc).

\textbf{6. Degeneracies and what the geometry cannot do.}
(i) \emph{Only $\beta$ is fixed by $2\ell$.} $\beta$ depends on $M_c/M_{\rm crit}(D)$, $\Omega_s$ and $R_{\rm in}/D$, all ratios: the core's tidal field in units of $H^2$ is $GM_c/(H^2D^3)=\frac12M_c/M_{\rm crit}(D)$. $D$ enters only through $\varepsilon=r/D$ in the odd terms.
(ii) \emph{Newtonian:} $r=R=14.3$~Gpc, so $3\ell$ fixes $\varepsilon Q$. $D$ trades against the curvature $\gamma$ of $h(q)$, i.e. the shell structure ($D\approx41$--47, 39--45 and 114--135~Gpc for the point-core, $\Omega_s$=0.049 and 0.909 shapes, same $R^2$).
(iii) \emph{Relativistic:} $r=1.19$~Gpc, so $r/D\le0.036$ and the odd terms do not depend on $D$. Along the fitted valley $a_2$ = %(va2)s and $a_3$ = %(va3)s stay the same for $D$ = %(vD)s~Gpc while $M_c/M_{\rm crit}(D)$ = %(vm)s. Hence the LTB fit pins $M_c/M_{\rm crit}(D)$, not $D$.
(iv) \emph{Not produced by the geometry:} the dipole is along the axis ($a_1\cos\psi$, peak at $\ell_0$ or $\ell_0+180^\circ$), but $A_F$ peaks at $\ell_1=55.7^\circ$, so the off-axis dipole needs the rotation (sine). The $3\ell$ phase offset ($S_3=+0.013$ vs $\le0.002$ from rotation) cannot be made: the model's $3\ell$ peak is locked to $\ell_0+60^\circ$ (data $\ell_3=67^\circ$ vs $\ell_0+60\approx73$--$75^\circ$). This missing ``$3\ell$ transform'' is the main residual of every fit.

{\footnotesize\textbf{References.} Taylor, G.~I. 1950, Proc.\ R.\ Soc.\ Lond.\ A, 201, 159. Sedov, L.~I. 1959, Similarity and Dimensional Methods in Mechanics (Academic Press). Migkas, K., et al. 2020, A\&A, 636, A15. Scripts: A\_equals\_fourier/thinshell/tm\_core.py, termmatch.py, relvalley.py, make\_tex.py; numbers in termmatch.json.}
}
\end{document}
"""
def row(k, name):
    f = F[k]; s = f["ser"]; n = f["num"]
    return "%s & %.1f & %.3f & %.3f (%.2f) & $%s$ & $%s$ / $%s$ & %.4f / %.4f & $%s$ / $%s$ & %.4f\\\\" % (name, f["D"], f["eps"], f["beta"], f["beta_lin"], f3(f["gamma"]), f4(s["a1"]), f4(n["1"][0]), s["a2"], n["2"][0], f4(s["a3"]), f4(n["3"][0]), n["4"][0])
rel = F["rel"]
rowsX = row("free", "Newt.\\ $\\Omega_s$ free") + "\n" + row("omb", "Newt.\\ $\\Omega_s$=0.049") + "\n" + \
    "Relativistic & 32.9 & 0.036$^a$ & %.3f$^b$ & -- & 0 / $%s$ & %.4f / %.4f & 0 / $%s$ & %.4f\\\\" % (rel["beta"], f4(rel["num"]["1"][0]), rel["ser"]["a2"], rel["num"]["2"][0], f4(rel["num"]["3"][0]), rel["num"]["4"][0])
Kc = hp["beta"] / (3 + hp["beta"]); mcc_lin = 2 * Kc / 99
v = dict(A1=AF["1"][0], l1=AF["1"][1], A2=AF["2"][0], l2=AF["2"][1], A3=AF["3"][0], l3=AF["3"][1], C1=f4(C["1"][0]), S1=f4(C["1"][1]), C3=f4(C["3"][0]), S3=f4(C["3"][1]),
         beta=hp["beta"], amp=f3(hp["amp"]), phi=hp["phi"], epsQ=f3(hp["epsQ"]), eps=hp["eps"], D=hp["D"], S3p=hp["S3_pred"], N=hp["N"], mcc_lin=mcc_lin,
         xl0=S["pointmass"]["p"][0], xbeta=S["pointmass"]["p"][1], xamp=f3(S["pointmass"]["p"][3]), xphi=S["pointmass"]["p"][4] % (2 * np.pi), xD=S["pointmass"]["D"],
         gbf=F["free"]["gamma"] / F["free"]["beta"], Df=H["shape_free"]["D"], xDf=S["shape_free"]["D"], gbo=F["omb"]["gamma"] / F["omb"]["beta"], Do=H["shape_omb"]["D"], xDo=S["shape_omb"]["D"],
         xR2F=S["pointmass"]["R2F"], xR2m=S["pointmass"]["R2mig"], rowsX=rowsX,
         R2free=F["free"]["R2"], R2omb=F["omb"]["R2"], R2rel=F["rel"]["R2"], R2Ffree=F["free"]["R2F"], R2Fomb=F["omb"]["R2F"], R2Frel=F["rel"]["R2F"],
         mccfree=F["free"]["mcc"], mccomb=F["omb"]["mcc"], mccrel=F["rel"]["mcc"], bomb=F["omb"]["beta"], aomb=f3(F["omb"]["amp"]), pomb=F["omb"]["phi"],
         va2="%.3f--%.3f" % (min(x["a"][2] for x in V), max(x["a"][2] for x in V)), va3="$-$%.3f--$-$%.3f" % (min(-x["a"][3] for x in V), max(-x["a"][3] for x in V)),
         vD="33--$7\\times10^6$", vm="%.3f--%.3f" % (min(x["mcc"] for x in V), max(x["mcc"] for x in V)))
for k in ("free", "omb", "rel"):
    for i in range(3): v["F%s%d" % (k, i + 1)] = Fr[k][i]
tex = tex % v
tex = tex.replace(r"\end{tabular}\end{center}" + "\n" + r"\vspace{-4pt}" + "\n" + "The linear formula", r"\end{tabular}\end{center}" + "\n" + r"{\footnotesize $^a$ $r/D$, $r=1.19$ Gpc. $^b$ $H_r/H_\perp-1$ at Earth.}\par\vspace{-2pt}" + "\n" + "The linear formula", 1)
open("../thinshell_termmatch.tex", "w").write(tex)
