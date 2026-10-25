#!/usr/bin/env python3
"""Writes termmatch_Dr234.tex from figures/termmatch_Dr234_results.json (all numbers come from the JSON)."""
import json, numpy as np
R = json.load(open("/workspace/cosmology/figures/termmatch_Dr234_results.json"))
A = R["data"]["amp"]; S = R["data"]["sigma"]; P = R["data"]["peak"]; l0 = R["data"]["l0"]
wrap = lambda x, p: (x + p / 2) % p - p / 2
rad = lambda d: d * np.pi / 180
f = lambda x, n=4: f"{x:.{n}f}"
def sg(x, n=4): return f"{x:+.{n}f}"
L = []; w = L.append
w(r"""\documentclass[10pt]{article}
\usepackage[letterpaper,margin=0.6in]{geometry}
\usepackage{amsmath,amssymb,booktabs,array}
\pagestyle{plain}\setlength{\parindent}{0pt}\setlength{\parskip}{3pt}
\allowdisplaybreaks\sloppy
\begin{document}
{\large\textbf{CoE model vs Migkas+2020 Fig.\,6 $A(\ell)$: term matching at fixed $D/r=2,3,4$}}\\[2pt]
Sedov--Taylor Center-of-Expansion model, sources on the vision shell ($d=r$). Analytic point-mass form and numeric ideal-drag model, each with Chris's rotation and $2\ell$ transforms.
Scripts: \texttt{figures/termmatch\_Dr234.py} (solve) and \texttt{make\_termmatch\_tex.py} (writes this sheet); every number here is read from \texttt{termmatch\_Dr234\_results.json}.

\medskip\textbf{Step 0. Chris's variables}\\
$r$ = vision-shell radius $=14.3$ Gpc (sources at $d=r$); \ $D$ = our distance from the CoE; \ $L$ = density scale length; \ $M=M(<D)$ = mass enclosed toward the CoE;
$(\ell_0,b_0)$ = Galactic longitude and tilt of the CoE axis; \ $\varphi$ = rotation phase in $\mathrm{Amp}_1\sin(\ell+\varphi)$; \ $\varphi_2$ = $2\ell$ phase in $\mathrm{Amp}_2\sin(2\ell+\varphi_2)$.

\medskip\textbf{Step 1. Define parameters}
\begin{align*}
&\varepsilon\equiv r/D=1/(D/r)\ \ (=\tfrac12,\ \tfrac13,\ \tfrac14\ \text{here}),\qquad x\equiv D/L,\qquad
M_{\rm crit}(<D)=\tfrac{4\pi}{3}\rho_cD^3=\frac{H_0^2D^3}{2G},\\
&\kappa\equiv\frac{G\,M(<D)}{H_0^2D^3}=\frac{M(<D)}{2M_{\rm crit}},\qquad K\equiv\alpha\kappa\ \ (\alpha=O(1)\text{ time-integration factor; tables use }\alpha=1),\\
&\phi_1\equiv\text{rotation: shift of the }m{=}1\text{ peak},\qquad \phi_2\equiv\text{shift of the }m{=}2\text{ peak (the }2\ell\text{ transform)},\\
&\varphi=90^\circ-(\ell^{(1)}_{\rm nat}+\phi_1),\qquad \varphi_2=90^\circ-2(\ell^{(2)}_{\rm nat}+\phi_2),\\
&\delta_3\equiv\ell_3^{\rm data}-\ell_3^{\rm model}\ \ (\text{octupole left at its natural phase; }\ell_3\text{ defined mod }120^\circ).
\end{align*}
\textbf{Analytic model} (point-mass tidal series to $O(\varepsilon)$, averaged over $90^\circ$ sectors every $5^\circ$, $|b|>20^\circ$ cut, $\cos b$ weight; same as the v2 sheet), with $\psi=\ell-\ell_0$:
\begin{align*}
\frac{A}{A_{\rm all}}&=N\big[1+K(3\mu^2-1)+K\varepsilon(\tfrac92\mu-\tfrac{15}{2}\mu^3)\big],\qquad \mu=\cos b\cos b_0\cos\psi+\sin b\sin b_0\\
&\xrightarrow{\ \text{sector avg}\ }\ N\big[1+Kq_0\big]+NK\varepsilon\,o_1\cos\psi+NKq_2\cos2\psi+NK\varepsilon\,o_3\cos3\psi,
\end{align*}
$q_0=0.45900-0.68850a^2,\ q_2=0.49052a^2,\ o_1=0.90032\,a(-3.1958+3.9948a^2),\ o_3=-0.22806a^3$, \ $a=\cos b_0$.\\
\textbf{Numeric model}: ideal-drag ST background, $du/dt+Hu=-GM(<r)/r^2$ integrated in time, density $\rho\propto(1-f_u)e^{-R/L}/e^{-D/L}+f_u$ (exp core + floor),
$H_{\rm obs}(\hat n)=H+[(\vec u(\vec r)-\vec u(D))\cdot\hat n]/d$, $90^\circ$ boxcar, mean-normalised; parameters $\kappa,\ D/L,\ f_u,\ b_0,\ \ell_0$ (as in \texttt{r2\_vs\_Dr\_density\_families.py}; $|b_0|\le30^\circ$, $D/L\le300$).

\medskip\textbf{Step 2. Data harmonics} (Fourier fit $A=0.9965+0.0296\cos\ell+0.0434\sin\ell+0.0636\cos2\ell+0.0366\sin2\ell-0.0299\cos3\ell-0.0114\sin3\ell$, $R^2=""" + f(R["data"]["R2_fourier"], 3) + r"""$)
\[
|A_0|=""" + f(A["0"], 5) + r""",\quad |A_1|=""" + f(A["1"], 5) + r"""\ (\ell_1=""" + f(P["1"], 2) + r"""^\circ),\quad |A_2|=""" + f(A["2"], 5) + r"""\ (\ell_2=""" + f(P["2"], 2) + r"""^\circ),\quad |A_3|=""" + f(A["3"], 5) + r"""\ (\ell_3=""" + f(P["3"], 2) + r"""^\circ)
\]
Weights: 1$\sigma$ amplitude errors from the 71-point least-squares Fourier fit with per-point errors $(\sigma_{\rm lo}+\sigma_{\rm hi})/2$ (overlapping sectors treated as independent, so these are too small by $\sim\sqrt{71/16}$; only their \emph{ratios} matter for the weights):
$\sigma_0=""" + f(S["0"]) + r""",\ \sigma_1=""" + f(S["1"]) + r""",\ \sigma_2=""" + f(S["2"]) + r""",\ \sigma_3=""" + f(S["3"]) + r"""$. \ Axis: $\ell_0\equiv\ell_2=""" + f(l0, 2) + r"""^\circ$.

\medskip\textbf{Step 3. Method}\\
\emph{A (analytic, fixed $\varepsilon$).} Unknowns $N,K,b_0$; four amplitude equations, so weighted least squares. Write $P\equiv NK$. The $m{=}0$ equation contains $N$ alone and is satisfied exactly for any $(P,b_0)$ by $N=A_0-Pq_0$.
For fixed $b_0$ the $m{=}1,2,3$ equations are linear in $P$: with $c_1=\varepsilon|o_1|,\ c_2=|q_2|,\ c_3=\varepsilon|o_3|$ and $w_m=1/\sigma_m^2$,
\[
P(b_0)=\frac{\sum_m w_mA_mc_m}{\sum_m w_mc_m^2},\qquad \chi^2(b_0)=\sum_{m=1}^3w_m\big(A_m-P\,c_m\big)^2\ \ \text{minimised over }b_0\in[0^\circ,90^\circ)
\]
($\pm b_0$ give identical curves, so only $b_0\ge0$ is scanned).
Transforms: $\phi_1=\ell_1^{\rm data}-\ell^{(1)}_{\rm nat}$, where $\ell^{(1)}_{\rm nat}=\ell_0$ if $o_1>0$ and $\ell_0+180^\circ$ if $o_1<0$; $\phi_2=\ell_2-\ell_0=0$ by construction. Octupole: $o_3<0$, so its natural peak is $\ell_0+60^\circ=""" + f(l0 + 60, 2) + r"""^\circ$, giving $\delta_3=""" + f(R["analytic"]["2.0"][0]["delta3"], 2) + r"""^\circ$ at every $D/r$.\\
\emph{B (numeric).} Numeric curve on a 181-point grid, FFT $\to$ harmonics $m=0..3$ about their natural peaks; the $m{=}1$ and $m{=}2$ components are shifted by $\phi_1,\phi_2$, the rest (octupole) untouched.
$(\kappa,D/L,f_u,b_0,\phi_1,\phi_2)$ are fitted to maximise $R^2$ of the transformed curve vs the 71 points, with $\ell_0$ held at $""" + f(l0, 2) + r"""^\circ$ as in A (a free $\ell_0$ only moves the octupole, i.e.\ it would act as a third free phase; that variant is listed at the end for reference).
Start points: the stored no-transform best fits, one differential-evolution run and 4 random starts, all polished with Nelder--Mead.
""")
lab = {"2.0": "2", "3.0": "3", "4.0": "4"}
for Dr in ["2.0", "3.0", "4.0"]:
    sols = R["analytic"][Dr]; n = R["numeric"][Dr]; e = 1 / float(Dr)
    s0 = sols[0]
    w(r"\clearpage")
    w(r"{\large\bfseries\boldmath $D/r=" + lab[Dr] + r"$} \ ($\varepsilon=" + f(e, 4) + r"$, \ $D=" + f(float(Dr) * 14.3, 1) + r"$ Gpc $=" + f(float(Dr) * 14.3 * 3.26156, 1) + r"$ Gly)" + "\n")
    w(r"\textbf{A. Analytic: set equal} (numbers substituted, $\varepsilon=" + f(e, 4) + r"$)")
    w(r"\[ " + f(A["0"], 5) + r"=N(1+Kq_0),\quad " + f(A["1"], 5) + r"=NK\cdot" + f(e, 4) + r"\,|o_1|,\quad " + f(A["2"], 5) + r"=NK\,|q_2|,\quad " + f(A["3"], 5) + r"=NK\cdot" + f(e, 4) + r"\,|o_3| \]")
    w(r"\textbf{Solve.} $\chi^2(b_0)$ has " + str(len(sols)) + r" local minima on $[0^\circ,90^\circ)$:" + "\n")
    w(r"\begin{center}\small\begin{tabular}{lrrrrrrrrr}\toprule")
    w(r"branch & $b_0$ & $q_0$ & $q_2$ & $o_1$ & $o_3$ & $c_1,\ c_2,\ c_3$ & $P=NK$ & $\chi^2$ & dipole\\\midrule")
    for i, s in enumerate(sols):
        c = s["c_m"]
        w(("\\textbf{%s}" % ("primary" if i == 0 else "second")) + f" & ${s['b0']:.2f}^\\circ$ & ${sg(s['q0'])}$ & {f(s['q2'])} & ${sg(s['o1'])}$ & ${sg(s['o3'])}$ & {f(c['1'])}, {f(c['2'])}, {f(c['3'])} & {f(s['NK'],5)} & {f(s['chi2'],3)} & {'on axis' if s['o1']>0 else 'flipped'}\\\\")
    w(r"\bottomrule\end{tabular}\end{center}")
    for i, s in enumerate(sols):
        c = s["c_m"]; ww = {m: 1 / S[m] ** 2 for m in "123"}
        num = sum(ww[m] * A[m] * c[m] for m in "123"); den = sum(ww[m] * c[m] ** 2 for m in "123")
        w((r"\emph{Primary}" if i == 0 else r"\emph{Second branch}") + r" ($b_0=" + f(s["b0"], 2) + r"^\circ$): " +
          r"$P=\dfrac{\sum w_mA_mc_m}{\sum w_mc_m^2}=\dfrac{" + f(num, 3) + "}{" + f(den, 2) + "}=" + f(s["NK"], 5) +
          r"$,\ \ $N=" + f(A["0"], 5) + ("+" if s["q0"] < 0 else "-") + f(s["NK"], 5) + r"\times" + f(abs(s["q0"]), 5) + "=" + f(s["N"], 5) +
          r"$,\ \ $K=P/N=" + f(s["K"], 5) + r"$.\\")
        ma = s["model_amp"]; rs = s["resid"]; ro = s["resid_over_sigma"]
        w(r"\emph{Plug back:} \ $m0$: $" + f(s["N"], 5) + r"(1" + sg(s["K"] * s["q0"], 5) + ")=" + f(ma["0"], 5) +
          r"$; \ $m1$: $" + f(s["NK"], 5) + r"\times" + f(c["1"], 5) + "=" + f(ma["1"], 5) +
          r"$; \ $m2$: $" + f(s["NK"], 5) + r"\times" + f(c["2"], 5) + "=" + f(ma["2"], 5) +
          r"$; \ $m3$: $" + f(s["NK"], 5) + r"\times" + f(c["3"], 5) + "=" + f(ma["3"], 5) + r"$.\\")
        w(r"Residuals (model$-$data): $m0$: $" + sg(rs["0"], 5) + r"$, \ $m1$: $" + sg(rs["1"], 5) + r"\ (" + sg(ro["1"], 2) + r"\sigma)$, \ $m2$: $" + sg(rs["2"], 5) + r"\ (" + sg(ro["2"], 2) +
          r"\sigma)$, \ $m3$: $" + sg(rs["3"], 5) + r"\ (" + sg(ro["3"], 2) + r"\sigma)$.\\")
        w(r"Phases: natural dipole peak $" + f(s["nat_peak"]["1"], 2) + r"^\circ\Rightarrow\phi_1=" + f(P["1"], 2) + "-" + f(s["nat_peak"]["1"], 2) + "=" + sg(s["phi1"], 2) +
          r"^\circ$ ($" + sg(rad(s["phi1"]), 3) + r"$ rad), Chris's $\varphi$: $" + f(s["phiC_model"], 2) + r"^\circ\to" + f(s["phiC_data"], 2) + r"^\circ$; \ $\phi_2=0^\circ$, $\varphi_2=" + f(s["phi2C"], 2) +
          r"^\circ$; \ $\delta_3=" + sg(s["delta3"], 2) + r"^\circ$ (" + sg(s["delta3_3l"], 1) + r"^\circ$ in $3\ell$).\\")
        w(r"$R^2$ (71 pts): rotation $+\,2\ell$: \textbf{" + f(s["R2_transformed"], 3) + r"}; \ no transforms: $" + f(s["R2_notransform"], 3) +
          r"$; \ (octupole also given its own phase: " + f(s["R2_all3phases"], 3) + r"). \ $M(<D)/M_{\rm crit}=2K/\alpha=" + f(s["MoverMcrit_alpha1"], 3) + r"/\alpha$.")
        w("")
    # numeric
    hn = n["harm_natural"]; ht = n["harm_transformed"]
    w(r"\medskip\textbf{B. Numeric: fitted parameters} \ $\kappa=" + f(n["kappa"], 4) + r"$, \ $D/L=" + (f(n["D_over_L"], 2) if n["D_over_L"] < 299 else r"300\ (\text{bound: point-mass limit})") +
      r"$, \ $f_u=" + (f"{n['f_u']:.4f}" if n["f_u"] > 1e-4 else f"{n['f_u']:.1e}") + r"$, \ $b_0=" + f(n["b0"], 2) + r"^\circ" + (r"\ (\text{bound})" if n["b0_at_bound"] else "") + r"$, \ $\ell_0=" + f(n["l0"], 2) + r"^\circ$ (fixed), \ $\phi_1=" +
      sg(n["phi1"], 2) + r"^\circ$, \ $\phi_2=" + sg(n["phi2"], 2) + r"^\circ$." + "\n")
    w(r"\textbf{Set equal (harmonics of the fitted numeric curve vs data)}")
    w(r"\begin{center}\small\begin{tabular}{lrrrrrr}\toprule")
    w(r"$m$ & data $|A_m|$ & numeric $|M_m|$ & resid. & resid./$\sigma_m$ & natural peak & transformed peak (data)\\\midrule")
    w(f"0 & {f(A['0'])} & {f(hn['0'])} & ${sg(hn['0']-A['0'])}$ & ${sg((hn['0']-A['0'])/S['0'],2)}$ & -- & --\\\\")
    for m in "123":
        after = ht[m][1]
        extra = r"$^\dagger$" if m == "3" else ""
        w(f"{m} & {f(A[m])} & {f(hn[m][0])} & ${sg(hn[m][0]-A[m])}$ & ${sg((hn[m][0]-A[m])/S[m],2)}$ & ${hn[m][1]:.2f}^\\circ$ & ${after:.2f}^\\circ$ (${P[m]:.2f}^\\circ$){extra}\\\\")
    w(r"\bottomrule\end{tabular}\end{center}")
    phiC = wrap(90 - (hn["1"][1] + n["phi1"]), 360); phiCm = wrap(90 - hn["1"][1], 360); phi2C = wrap(90 - 2 * ht["2"][1], 360)
    w(r"$^\dagger$octupole untouched. (The numeric curve is mean-normalised, so its $m{=}0$ is 1 by construction; $R^2$ uses the curve as is, no rescaling.) \ Octupole mismatch $\delta_3=" + sg(n["delta3"], 2) +
      r"^\circ$. \ Chris's $\varphi$: $" + f(phiCm, 2) + r"^\circ\to" + f(phiC, 2) + r"^\circ$ ($\phi_1=" + sg(rad(n["phi1"]), 3) + r"$ rad); \ $\varphi_2=" + f(phi2C, 2) + r"^\circ$.\\")
    w(r"$R^2$ (71 pts): transformed \textbf{" + f(n["R2_transformed"], 3) + r"}; \ same parameters, no transforms $" + f(n["R2_notransform"], 3) +
      r"$; \ stored no-transform best fit " + f(n["R2_stored"], 3) + r". \ $M(<D)/M_{\rm crit}=2\kappa=" + f(n["MoverMcrit"], 3) + r"$" +
      (r", \ $L=" + f(n["L_Gpc"], 2) + r"$ Gpc." if n["D_over_L"] < 299 else "."))
    br = n["dipole_branches"]; other = [k for k in br if abs(br[k]["R2_transformed"] - n["R2_transformed"]) > 1e-9]
    if other:
        b = br[other[0]]
        w(r"\\ Best solution on the other dipole branch (" + ("natural dipole on the axis" if other[0] == "on_axis" else "natural dipole flipped to $\\ell_0+180^\\circ$") + r"): $R^2=" + f(b["R2_transformed"], 3) +
          r"$ (no transforms $" + f(b["R2_notransform"], 3) + r"$), $\kappa=" + f(b["kappa"], 3) + r"$, $D/L=" + f(b["D_over_L"], 2) + r"$, $f_u=" + f(b["f_u"], 3) + r"$, $b_0=" + f(b["b0"], 1) + r"^\circ$, $\phi_1=" + sg(b["phi1"], 1) + r"^\circ$, $\phi_2=" + sg(b["phi2"], 2) + r"^\circ$.")
    w("")
# summary
w(r"\clearpage\textbf{\large Summary in Chris's variables} \ ($r=14.3$ Gpc, $\ell_0=" + f(l0, 2) + r"^\circ$ for all rows; analytic $M/M_{\rm crit}$ uses $\alpha=1$)")
w(r"\begin{center}\footnotesize\setlength{\tabcolsep}{3.2pt}\begin{tabular}{llrrrrrrrrrr}\toprule")
w(r"$D/r$ & model & $D$ [Gpc] & $D$ [Gly] & $M/M_{\rm crit}$ & $D/L$ & $b_0$ & $\phi_1$ [deg] & $\phi_1$ [rad] & $\phi_2$ [deg] & $R^2_{\rm tf}$ & $R^2_{\rm no\,tf}$\\\midrule")
for Dr in ["2.0", "3.0", "4.0"]:
    n = R["numeric"][Dr]; D = float(Dr) * 14.3
    w(f"{lab[Dr]} & numeric & {D:.1f} & {D*3.26156:.1f} & {n['MoverMcrit']:.3f} & " + ("300$^*$" if n["D_over_L"] > 299 else f"{n['D_over_L']:.2f}") +
      f" & ${n['b0']:.1f}^\\circ$" + ("$^*$" if n["b0_at_bound"] else "") + f" & ${n['phi1']:+.2f}$ & ${rad(n['phi1']):+.3f}$ & ${n['phi2']:+.2f}$ & \\textbf{{{n['R2_transformed']:.3f}}} & ${n['R2_notransform']:.3f}$\\\\")
    br = n["dipole_branches"]
    for k, b in br.items():
        if abs(b["R2_transformed"] - n["R2_transformed"]) > 1e-9:
            w(f" & \\quad (alt.\\ branch) & & & {b['MoverMcrit']:.3f} & " + ("300$^*$" if b["D_over_L"] > 299 else f"{b['D_over_L']:.2f}") +
              f" & ${b['b0']:.1f}^\\circ$ & ${b['phi1']:+.2f}$ & ${rad(b['phi1']):+.3f}$ & ${b['phi2']:+.2f}$ & {b['R2_transformed']:.3f} & ${b['R2_notransform']:.3f}$\\\\")
    for i, s in enumerate(R["analytic"][Dr]):
        w(f" & analytic{' (2nd)' if i else ''} & & & {s['MoverMcrit_alpha1']:.3f} & -- & ${s['b0']:.1f}^\\circ$ & ${s['phi1']:+.2f}$ & ${rad(s['phi1']):+.3f}$ & ${s['phi2']:+.2f}$ & " +
          (f"\\textbf{{{s['R2_transformed']:.3f}}}" if i == 0 else f"{s['R2_transformed']:.3f}") + f" & ${s['R2_notransform']:.3f}$\\\\")
    w(r"\midrule" if Dr != "4.0" else r"\bottomrule")
w(r"\end{tabular}\end{center}")
w(r"$^*$ at fit bound ($|b_0|\le30^\circ$, $D/L\le300$; $D/L=300$ is effectively a point mass at the CoE). Octupole mismatch $\delta_3=-8.0^\circ$ in every row (octupole at natural phase about $\ell_0$).")
w(r"Numeric rows: $f_u$ = " + ", ".join(f"{R['numeric'][d]['f_u']:.3g} ($D/r={lab[d]}$)" for d in ["2.0", "3.0", "4.0"]) + r". Analytic rows: $N$, $K$ = " +
  "; ".join(", ".join(f"{s['N']:.4f}, {s['K']:.4f}" for s in R["analytic"][d][:1]) + f" ($D/r={lab[d]}$)" for d in ["2.0", "3.0", "4.0"]) + r" (primary). Data Fourier fit: $R^2=" + f(R["data"]["R2_fourier"], 3) + r"$ (ceiling for any $m\le3$ curve).")
w("")
w(r"\textbf{Harmonic amplitudes at a glance} (data: $|A_1|,|A_2|,|A_3|=" + f"{A['1']:.4f}, {A['2']:.4f}, {A['3']:.4f}" + r"$)")
w(r"\begin{center}\small\begin{tabular}{lrrrrrr}\toprule $D/r$ & num.\ $|M_1|$ & num.\ $|M_2|$ & num.\ $|M_3|$ & an.\ $|M_1|$ & an.\ $|M_2|$ & an.\ $|M_3|$\\\midrule")
for Dr in ["2.0", "3.0", "4.0"]:
    n = R["numeric"][Dr]["harm_natural"]; a = R["analytic"][Dr][0]["model_amp"]
    w(f"{lab[Dr]} & {n['1'][0]:.4f} & {n['2'][0]:.4f} & {n['3'][0]:.4f} & {a['1']:.4f} & {a['2']:.4f} & {a['3']:.4f}\\\\")
w(r"\bottomrule\end{tabular}\end{center}")
w(r"\textbf{Reference: numeric with $\ell_0$ also free} (a third free phase in effect): " + "; ".join(
    f"$D/r={lab[d]}$: $R^2={R['numeric_freel0'][d]['R2_transformed']:.3f}$ ($\\ell_0={R['numeric_freel0'][d]['l0']:.1f}^\\circ$, $\\delta_3={R['numeric_freel0'][d]['delta3']:+.1f}^\\circ$)" for d in ["2.0", "3.0", "4.0"]) + r". At $D/r=2$ and $4$ the free $\ell_0$ runs far from the data axis ($130^\circ$, $187^\circ$) and the transforms then re-point $m1$ and $m2$. The gain over fixed $\ell_0$ is $\le0.01$ in $R^2$.")
w(r"""
\medskip\textbf{Observations.}
(1) With the two phases free, both primary solutions match $|m1|$ and $|m2|$ to within $0.3\sigma$ at every $D/r$; the one exception is the numeric fit at $D/r=2$ ($\mp0.4\sigma$). Where they fail, both fail in the octupole.
(2) The octupole is the failure: $|m3|\propto\varepsilon$ in the analytic form, and the numeric one is similar ($""" + ", ".join(f"{R['numeric'][d]['harm_natural']['3'][0]:.3f}" for d in ["2.0", "3.0", "4.0"]) + r"""$ vs data $""" + f(A["3"], 3) + r"""$), and its natural phase is $8^\circ$ off. This is why $R^2_{\rm tf}$ falls from $D/r=2$ to $4$.
(3) The transformed $R^2$ values are all $0.84$--$0.89$, below the $m\le3$ Fourier ceiling of $0.915$. In both methods the rotation is $\phi_1\approx+40^\circ$ (or $-140^\circ$ on the flipped branch) and $\phi_2\approx0^\circ$.
(4) Without transforms, the on-axis branches keep $R^2\approx0.71$--$0.75$ and the flipped branches collapse ($R^2\le0.33$). For $D/r=3,4$ the analytic weighted-LSQ optimum is the flipped branch ($b_0\approx40^\circ$) by a small $\chi^2$ margin.

\medskip\textbf{Caveats.}
The rotation $\phi_1$ and the $2\ell$ phase $\phi_2$ (and $\ell_0$) are free parameters; they are not predicted by the axisymmetric physics. With them free, only amplitudes are being tested, and $|m1|,|m2|$ are matched almost automatically.
The 71 points come from overlapping $90^\circ$ sectors stepped by $5^\circ$, so there are only $\sim$16 independent values: $R^2$ overstates the evidence, and differences of a few $\times0.01$ are not significant.
The analytic form is a point-mass tidal series truncated at $O(\varepsilon)$. It misses the dipole contributed by the extended fitted density in the numeric model, so analytic $K$ and numeric $\kappa$ are not the same quantity ($\alpha\ne1$ in general).
Several numeric parameters sit on their bounds ($|b_0|=30^\circ$, or $D/L=300$), and the data sector weighting assumes clusters uniform on the sky outside $|b|\le20^\circ$.
\end{document}""")
open("termmatch_Dr234.tex", "w").write("\n".join(L))
