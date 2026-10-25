#!/usr/bin/env python3
"""Writes exact_termmatch.tex (journal style) from the JSON outputs of exact_termmatch.py, exact_termmatch_ext2.py,
exact_profile.py and exact_termmatch_plot.py. Every number in the sheet is read from those files."""
import json, numpy as np
F = "/workspace/cosmology/figures/"
R = json.load(open(F + "exact_termmatch_results.json")); P = json.load(open(F + "exact_profile_results.json"))
O = json.load(open(F + "fixedDr_octmiss_results.json")); SP = json.load(open(F + "series_profile.json"))
A = R["data"]["amp"]; Sc = R["sigmas"]["corr"]; Sn = R["sigmas"]["naive"]; pk = R["data"]["peak"]; l0 = R["data"]["l0"]
f = lambda x, n=4: f"{x:.{n}f}"; sg = lambda x, n=2: f"{x:+.{n}f}"
DRg = np.array(P["DR"])
L = []; w = L.append

def interval(y, ymin, i, dth):
    inside = y <= ymin + dth; lo_i = i
    while lo_i > 0 and inside[lo_i - 1]: lo_i -= 1
    hi_i = i
    while hi_i < len(y) - 1 and inside[hi_i + 1]: hi_i += 1
    lo = "1 (bound)" if lo_i == 0 else f(float(np.interp(ymin + dth, [y[lo_i], y[lo_i - 1]], [DRg[lo_i], DRg[lo_i - 1]])), 2)
    hi = "$>6$" if hi_i == len(y) - 1 else f(float(np.interp(ymin + dth, [y[hi_i], y[hi_i + 1]], [DRg[hi_i], DRg[hi_i + 1]])), 2)
    return lo, hi
def prof(kind, name): return np.array([P["profile"][kind][str(float(d))][name][0] for d in DRg])

w(r"""\documentclass[10pt]{article}
\usepackage[letterpaper,margin=0.75in]{geometry}
\usepackage{amsmath,amssymb,booktabs,array,graphicx,caption,float}
\captionsetup{font=small,labelfont=bf}
\setlength{\parskip}{3pt}\sloppy
\begin{document}
\begin{center}
{\Large\textbf{Exact line-of-sight tidal term for the Sedov--Taylor Center-of-Expansion model and its fit to the Migkas et al.\ (2020) cluster anisotropy}}\\[6pt]
{\small Worked sheet. Scripts: \texttt{figures/exact\_tidal\_core.py}, \texttt{exact\_termmatch.py}, \texttt{exact\_termmatch\_ext2.py}, \texttt{exact\_profile.py}, \texttt{exact\_termmatch\_plot.py}; sheet written by \texttt{papers/A\_equals\_fourier/make\_exact\_termmatch\_tex.py} from the result files.}
\end{center}

\section{Introduction}
Earlier sheets compared the Center-of-Expansion (CoE) model with the longitude dependence of the Migkas et al.\ cluster scaling-relation anisotropy $A(\ell)$. They used an analytic form of the tidal term expanded to first order in the ratio of the vision-shell radius to the distance from the CoE. That expansion was chosen only because it gives the sector-averaged harmonics in closed form. It is not a physical approximation, and the preferred solutions sit where the expansion parameter is close to one. In this sheet the line-of-sight tidal term is evaluated exactly, with no expansion. It is computed first for a point mass and then for an extended, decreasing density profile, averaged over the Migkas sectors, and projected numerically onto longitude harmonics. The harmonic-matching procedure is then repeated, and the distance ratio is treated as a continuous parameter that is fitted rather than sampled.

\section{Chris's variables and the model}
The model is written in Chris's variables. These are the vision-shell radius $r$ on which the sources sit, the observer's distance $D$ from the CoE, the density scale length $L$, the enclosed mass $M$ (the spreadsheet's $M/M_0$), the Galactic longitude $\ell_0$ and tilt $b_0$ of the CoE axis, the rotation phase $\varphi$ applied to the dipole, the phase $\varphi_2$ of the $2\ell$ term, and the harmonic amplitudes $\mathrm{Amp}_m$. The first-order analytic form used previously is
\begin{equation}
\frac{A}{A_{\rm all}}=N\Big[1+K\,(3\mu^2-1)+K\varepsilon\big(\tfrac92\mu-\tfrac{15}{2}\mu^3\big)\Big],
\label{eq:series}
\end{equation}
where $N$ is an overall normalisation, $K$ is the tidal strength defined in Sec.~\ref{sec:define}, $\varepsilon$ is the distance ratio defined there, and $\mu$ is the cosine of the angle between the line of sight and the CoE axis.

\section{Define parameters}\label{sec:define}
The compact parameters used below are
\begin{align}
\varepsilon&\equiv\frac{r}{D},\qquad x\equiv\frac{D}{L},\qquad \mu=\cos b\cos b_0\cos(\ell-\ell_0)+\sin b\sin b_0,\label{eq:def1}\\
\kappa&\equiv\frac{G\,M(<D)}{H^2D^3}=\frac{M(<D)}{2M_{\rm crit}(<D)},\qquad M_{\rm crit}(<D)=\frac{4\pi}{3}\rho_cD^3,\qquad K\equiv-\frac{u(D)}{H\,D}=\alpha\,\kappa,\label{eq:def2}
\end{align}
where $(\ell,b)$ is the direction of a source, $H$ is the Sedov--Taylor background expansion rate, $\rho_c$ is the critical density, $u(R)$ is the radial peculiar velocity at distance $R$ from the CoE, and $\alpha$ is the ratio between the tidal strength and the enclosed-mass parameter. In the earlier sheets $\alpha$ was left as an order-unity factor, and tables quoted $\alpha=1$. Here it is computed from the ideal-drag integration of the numerical code,
\begin{equation}
u(R,t_0)=-\int_{t_i}^{t_0}\exp\!\Big[-\!\int_{t'}^{t_0}H\,dt''\Big]\,\frac{G\,M\big(<R(t')\big)}{R(t')^2}\,dt',\qquad R(t')=R\,(t'/t_0)^{2/5},
\label{eq:drag}
\end{equation}
which for a point mass gives $\alpha=2/3$ analytically. Both conventions are tabulated.

\section{Exact line-of-sight tidal term}
With the CoE at the origin, the observer at $\vec x_o=D\hat x$ and a source at $\vec x_s=\vec x_o+r\hat n$, the fractional change of the observed expansion rate is
\begin{equation}
\frac{\delta H}{H}(\hat n)=\frac{\big[\vec u(\vec x_s)-\vec u(\vec x_o)\big]\cdot\hat n}{r\,H}=K\,\mathcal T(\mu,\varepsilon),\qquad
\mathcal T=\frac{1}{\varepsilon}\Big[\mu-f(\rho_s)\,\frac{\mu+\varepsilon}{\rho_s}\Big],\qquad \rho_s=\sqrt{1+2\varepsilon\mu+\varepsilon^2},
\label{eq:exact}
\end{equation}
where $f(\rho)=u(\rho D)/u(D)$ is the velocity shape and $\rho_s$ is the source's distance from the CoE in units of $D$. For a point mass the field is exactly inverse-square, so $f=\rho^{-2}$ and
\begin{equation}
\mathcal T_{\rm pm}=\frac{\mu-(\mu+\varepsilon)\,(1+2\varepsilon\mu+\varepsilon^2)^{-3/2}}{\varepsilon}
=(3\mu^2-1)+\varepsilon\big(\tfrac92\mu-\tfrac{15}{2}\mu^3\big)+O(\varepsilon^2).
\label{eq:pm}
\end{equation}
The expansion on the right is exactly Eq.~\eqref{eq:series} with $A/A_{\rm all}=N(1+K\mathcal T)$, so the normalisation and constant are unchanged, and the series is the small-$\varepsilon$ limit of the exact form. For the extended model the density is the pure exponential of the numerical code, $\rho\propto e^{-R/L}$ with $x$ free. The shape $f$ is obtained from Eq.~\eqref{eq:drag} with the enclosed mass computed by quadrature, and it tends to the point-mass shape as $x\to\infty$.

\section{Sector averaging and harmonic projection}
Each Migkas value is an average over a longitude sector of width $\Delta\ell$ centred at $\ell_c$, over all latitudes outside the Galactic-plane cut, with area weight $\cos b$ (Table~\ref{tab:setup}). Because the average is linear,
\begin{equation}
\Big\langle\frac{A}{A_{\rm all}}\Big\rangle(\ell_c)=N\big[1+K\tau_0\big]+NK\sum_{m\ge1}\tau_m\cos m(\ell_c-\ell_0),\qquad
\tau_m=W_m\,\frac{1}{\pi}\!\int_0^{2\pi}\!\bar{\mathcal T}(\psi)\cos m\psi\,d\psi,\qquad W_m=\frac{\sin(m\Delta\ell/2)}{m\Delta\ell/2},
\label{eq:tau}
\end{equation}
where $\bar{\mathcal T}(\psi)$ is the latitude average of $\mathcal T$ at longitude offset $\psi=\ell-\ell_0$, $W_m$ is the exact response of the sector boxcar, and $\tau_0$ is the mean with the $1/\pi$ replaced by $1/2\pi$. The sine terms vanish because $\mathcal T$ depends on $\psi$ only through $\cos\psi$. The first-order form gives $\tau_0=q_0$, $\tau_1=\varepsilon o_1$, $\tau_2=q_2$ and $\tau_3=\varepsilon o_3$, where $q_0,q_2,o_1,o_3$ are the closed-form latitude coefficients of the earlier sheets. The exact form also produces harmonics above $m=3$. The latitude integral uses Gauss--Legendre nodes and the longitude integral a uniform periodic grid; the numerical settings and checks are listed in Table~\ref{tab:setup}.
""")
lim = R["limit"]
w(r"""
\begin{table}[!htb]\centering\small
\caption{Analysis settings and numerical checks.}\label{tab:setup}
\begin{tabular}{ll}\toprule
Quantity & Value\\\midrule
Vision-shell radius $r$ & $14.3$ Gpc\\
Sector width $\Delta\ell$; sector spacing; number of data points & $90^\circ$; $5^\circ$; $71$\\
Latitude cut; area weight & $|b|>20^\circ$; $\cos b$\\
CoE axis longitude $\ell_0$ (quadrupole peak of the data) & $""" + f(l0, 2) + r"""^\circ$\\
Quadrature (scans; convergence check) & $300\times2$ latitude nodes, $720$ longitudes; $800\times2$, $2880$\\
Harmonics retained in model curves & $m\le30$\\
Max.\ difference, harmonic route vs brute-force $72$-sector average & $""" + (lambda v: f"{v/10**np.floor(np.log10(v)):.1f}\\times10^{{{int(np.floor(np.log10(v)))}}}")(R['sector_check_maxdiff']) + r"""$\\
$\max|\mathcal T_{\rm pm}-\mathcal T_{\rm series}|/\varepsilon^2$ at $\varepsilon=0.1,\,0.01,\,0.001$ & """ + ", ".join(f(lim["pointwise_maxdiff_over_eps2"][k], 3) for k in ("0.1", "0.01", "0.001")) + r"""\\
Point-mass $\alpha$: analytic; numerical (Euler) integration & $2/3$; $""" + f(R["alpha_point_mass_numeric"], 4) + r"""$\\\bottomrule
\end{tabular}
\end{table}
""")
w(r"""
\section{Data harmonics and uncertainties}
The data are summarised by the least-squares Fourier fit to the digitised points,
\begin{equation}
A_{\rm fit}(\ell)=A_0+\sum_{m=1}^{3}A_m\cos m(\ell-\ell_m),
\label{eq:datafit}
\end{equation}
where $A_m$ and $\ell_m$ are the amplitude and peak longitude of harmonic $m$. Their uncertainties are propagated from the coefficient covariance in two ways. The adopted ``sector-correlated'' covariance accounts for the overlap of neighbouring sectors, which share a fraction $\max(0,1-|\ell_i-\ell_j|/\Delta\ell)$ of their clusters. The ``naive'' covariance treats all points as independent. Both are given in Table~\ref{tab:data}.

\begin{table}[!htb]\centering\small
\caption{Data harmonics of the Migkas $A(\ell)$ curve and their $1\sigma$ uncertainties.}\label{tab:data}
\begin{tabular}{lcccc}\toprule
 & $A_0$ & $A_1$ & $A_2$ & $A_3$\\\midrule
Amplitude & $""" + "$ & $".join(f(A[str(m)]) for m in range(4)) + r"""$\\
Peak longitude $\ell_m$ & -- & $""" + f(pk["1"], 2) + r"""^\circ$ & $""" + f(pk["2"], 2) + r"""^\circ$ & $""" + f(pk["3"], 2) + r"""^\circ$\\
$\sigma$, sector-correlated (adopted) & $""" + "$ & $".join(f(Sc[str(m)]) for m in range(4)) + r"""$\\
$\sigma$, naive & $""" + "$ & $".join(f(Sn[str(m)]) for m in range(4)) + r"""$\\\bottomrule
\end{tabular}
\end{table}

\section{Matching procedure}
Setting model and data amplitudes equal gives four equations,
\begin{equation}
N(1+K\tau_0)=A_0,\qquad NK|\tau_1|=A_1,\qquad NK|\tau_2|=A_2,\qquad NK|\tau_3|=A_3,
\label{eq:match}
\end{equation}
where every $\tau_m$ is a function of $(b_0,\varepsilon)$ for the point mass and of $(b_0,\varepsilon,x)$ for the extended density. With $\varepsilon$ free, the point-mass problem has four unknowns $(N,K,b_0,\varepsilon)$ and can be solved exactly. With $\varepsilon$ fixed it is overdetermined and is solved by least squares. Because $N$ appears only in the first equation, it is eliminated through $N=A_0-P\tau_0$ with $P\equiv NK$, and the objective is
\begin{equation}
\chi^2(b_0,\varepsilon\,[,x])=\min_P\sum_{m=1}^{3}\frac{\big(A_m-P|\tau_m|\big)^2}{\sigma_m^2},\qquad
P=\frac{\sum_m A_m|\tau_m|/\sigma_m^2}{\sum_m\tau_m^2/\sigma_m^2}.
\label{eq:chi2}
\end{equation}
The transforms are applied as before. The rotation sets the dipole peak to the data, $\phi_1=\ell_1-\ell^{(1)}_{\rm nat}$, where $\ell^{(1)}_{\rm nat}$ is $\ell_0$ for $\tau_1>0$ (``on-axis'' dipole) and $\ell_0+180^\circ$ for $\tau_1<0$ (``flipped''). The $2\ell$ phase $\phi_2$ sets the quadrupole peak to the data. The octupole and all higher harmonics keep their natural phases. The model curve compared with the points is then
\begin{equation}
A_{\rm mod}(\ell)=N(1+K\tau_0)+\sum_{m=1}^{2}NK|\tau_m|\cos m(\ell-\ell_m)+NK|\tau_3|\cos 3(\ell-\ell^{(3)}_{\rm nat})+\sum_{m\ge4}NK\tau_m\cos m(\ell-\ell_0),
\label{eq:curve}
\end{equation}
and in Chris's notation $\varphi=90^\circ-\ell_1$ and $\varphi_2=90^\circ-2\ell_2$. Two objectives are used when $D/r$ is fitted as a continuous parameter. The first is the amplitude $\chi^2$ of Eq.~\eqref{eq:chi2}. The second is the coefficient of determination $R^2$ of Eq.~\eqref{eq:curve} against the points, with $N$ and $P$ fitted linearly. Each is profiled over the remaining parameters,
\begin{equation}
\chi^2_{\rm prof}(D/r)=\min_{b_0\,[,x]}\chi^2,\qquad R^2_{\rm prof}(D/r)=\max_{N,P,b_0\,[,x]}R^2,\qquad 1<D/r\le6,
\label{eq:profile}
\end{equation}
and the intervals quoted are the ranges where $\chi^2_{\rm prof}$ lies within one and four of its minimum.
""")
# Table 3: limit + exact/series ratios
C = R["cont"]
def crow(Dr):
    c = min(C, key=lambda c: abs(c["Dr"] - Dr)); return c
w(r"""
\section{Results}
\subsection{Exact versus first-order coefficients}
Table~\ref{tab:limit} checks the small-$\varepsilon$ limit. The exact coefficients converge to the first-order ones, and the residual falls as $\varepsilon^2$, as Eq.~\eqref{eq:pm} requires. Table~\ref{tab:ratio} and Fig.~\ref{fig:exact}a,b show that the two forms separate strongly when $D/r$ approaches one. The largest change is a suppression of the quadrupole and octupole coefficients when the CoE axis lies near the plane. When the tilt is closer to the free-solution value the changes are milder but still noticeable.

\begin{table}[!htb]\centering\small
\caption{Ratio of exact to first-order sector-averaged coefficients, $\tau_m^{\rm exact}/\tau_m^{\rm series}$, in the small-$\varepsilon$ limit.}\label{tab:limit}
\begin{tabular}{cccccc}\toprule
$b_0$ & $\varepsilon$ & monopole & dipole & quadrupole & octupole\\\midrule
""")
for l in lim["tau_ratios"]:
    w(f"${f(l['b0'],1)}^\\circ$ & ${l['eps']:g}$ & " + " & ".join(f"${f(v,5)}$" for v in l["ratio"]) + r"\\")
w(r"""\bottomrule
\end{tabular}
\end{table}

\begin{table}[!htb]\centering\small
\caption{Ratio of exact to first-order coefficients at finite distance ratio, for an in-plane axis and for a tilt near the free solution.}\label{tab:ratio}
\begin{tabular}{c cccc cccc}\toprule
 & \multicolumn{4}{c}{$b_0=0^\circ$} & \multicolumn{4}{c}{$b_0=19.1^\circ$}\\\cmidrule(lr){2-5}\cmidrule(lr){6-9}
$D/r$ & $\tau_0$ & $\tau_1$ & $\tau_2$ & $\tau_3$ & $\tau_0$ & $\tau_1$ & $\tau_2$ & $\tau_3$\\\midrule
""")
for Dr in (1.05, 1.25, 1.45, 2.1, 2.9, 4.1, 6.0):
    c = crow(Dr)
    w(f"${f(c['Dr'],2)}$ & " + " & ".join(f"${f(v,3)}$" for v in c["ratio_b0.0"]) + " & " + " & ".join(f"${f(v,3)}$" for v in c["ratio_b19.1"]) + r"\\")
w(r"""\bottomrule
\end{tabular}
\end{table}
""")
# free solutions
fp = R["free_pm"]; fs = O["free"]
w(r"""
\subsection{Free distance ratio: exact solutions}
Solving Eq.~\eqref{eq:match} with all four unknowns free gives the roots in Table~\ref{tab:free}. The first-order form has one root with $D$ just above $r$ and one with $D<r$. With the exact point-mass form, both roots lie at $D>r$, so an exact four-harmonic match with the observer outside the vision shell exists. The roots were located on a $(b_0,\varepsilon)$ grid and refined, and they are unchanged on the finer quadrature. No other root exists in the scanned range. Close to $D/r=1$, however, the exact model also carries sizeable harmonics above the octupole, which the four equations do not constrain. These harmonics and the octupole phase offset lower $R^2$ at the exact roots.

\begin{table}[!htb]\centering\small
\caption{Exact four-harmonic solutions with free $D/r$ (amplitudes matched exactly, so all residuals vanish). The last column is the largest amplitude among the harmonics with $m\ge4$.}\label{tab:free}
\begin{tabular}{lcccccccc}\toprule
Model, branch & $D/r$ & $b_0$ & $N$ & $K$ & $\phi_1$ & $\delta_3$ & $R^2$ & max amp.\ $m\ge4$\\\midrule
""")
for s in fs:
    w(f"First order, {'on-axis' if s['o1']>0 else 'flipped'} & ${f(s['D_over_r'],4)}$ & ${f(s['b0'],2)}^\\circ$ & ${f(s['N'])}$ & ${f(s['K'])}$ & ${sg(s['phi1'],2)}^\\circ$ & ${f(s['delta3'],2)}^\\circ$ & ${f(s['R2'],3)}$ & $0$\\\\")
for s in fp:
    w(f"Exact point mass, {'on-axis' if s['tau']['1']>0 else 'flipped'} & ${f(s['D_over_r'],4)}$ & ${f(s['b0'],2)}^\\circ$ & ${f(s['N'])}$ & ${f(s['K'])}$ & ${sg(s['phi1'],2)}^\\circ$ & ${f(s['delta3'],2)}^\\circ$ & ${f(s['R2'],3)}$ & ${f(s['tau_hi_max']*s['NK'],4)}$\\\\")
w(r"""\bottomrule
\end{tabular}
\end{table}
""")
# continuous fit
pmc, pmn, pmr = prof("pm", "chi2_corr"), prof("pm", "chi2_naive"), -prof("pm", "negR2")
exc, exn, exr = np.minimum(prof("ext", "chi2_corr"), pmc), np.minimum(prof("ext", "chi2_naive"), pmn), np.maximum(-prof("ext", "negR2"), pmr)
sc, sn, sr = np.array(SP["chi2_corr"]), np.array(SP["chi2_naive"]), -np.array(SP["negR2"])
def rowfit(name, yc, yn, yr, exact_note):
    out = []
    for lab, y in (("$\\chi^2$, sector-corr.", yc), ("$\\chi^2$, naive", yn)):
        i = int(np.argmin(y)); ymin = y[i]
        zero = DRg[y < 1e-8]
        best = exact_note if exact_note else (f"{f(zero.min(),2)}--{f(zero.max(),2)} (family)" if len(zero) > 2 else f(DRg[i], 2))
        lo1, hi1 = interval(y, ymin, i, 1.0); lo4, hi4 = interval(y, ymin, i, 4.0)
        out.append(f"{name} & {lab} & {best} & ${f(max(ymin,0),3)}$ & {lo1} -- {hi1} & {lo4} -- {hi4}\\\\")
    i = int(np.argmax(yr))
    out.append(f"{name} & $R^2$ (max) & {f(DRg[i],2)} & ${f(yr[i],4)}$ & -- & --\\\\")
    return out
w(r"""
\subsection{Continuous fit of the distance ratio}
The profiles of Eq.~\eqref{eq:profile} are shown in Fig.~\ref{fig:profile}, and the best-fitting distance ratios with their intervals are collected in Table~\ref{tab:bestDr}. For the point mass the amplitude $\chi^2$ reaches zero at the exact roots of Table~\ref{tab:free}. It then rises slowly with $D/r$, and the rise is slower with the sector-correlated errors than with the naive ones. With either error model the four-unit threshold is not reached anywhere in the allowed range, because the octupole is itself only a marginal detection. For the extended density the $\chi^2$ is zero along a whole family of exact solutions, discussed below, so the extended fit does not select a unique $D/r$. The $R^2$ profile peaks slightly above $D/r=1$, falls steeply towards $D/r=1$ for the point mass because of the higher harmonics, and declines slowly at large $D/r$. The first-order profiles are included for comparison. The best-fit parameters are listed in Table~\ref{tab:bestpar} and translated into Chris's variables in Table~\ref{tab:chris}.

\begin{table}[!htb]\centering\small
\caption{Best-fit $D/r$ from the continuous profile fit, $1<D/r\le6$. For $\chi^2$ the best value is where the profile vanishes (an exact solution) or is smallest; the intervals are where $\Delta\chi^2\le1$ and $\le4$. ``1 (bound)'' means the interval extends to the constraint $D/r>1$.}\label{tab:bestDr}
\begin{tabular}{llcccc}\toprule
Model & Objective & best $D/r$ & min.\ $\chi^2$ or max.\ $R^2$ & $\Delta\chi^2\le1$ & $\Delta\chi^2\le4$\\\midrule
""")
pmroots = ", ".join(f(s["D_over_r"], 3) for s in sorted(fp, key=lambda s: s["D_over_r"])) + " (roots)"
for line in rowfit("Exact point mass", pmc, pmn, pmr, pmroots): w(line)
w(r"\midrule")
for line in rowfit("Exact extended", exc, exn, exr, None): w(line)
w(r"\midrule")
for line in rowfit("First order", sc, sn, sr, f(fs[0]["D_over_r"], 4) + " (root)"): w(line)
w(r"""\bottomrule
\end{tabular}
\end{table}
""")
# profile values table
w(r"""
\begin{table}[!htb]\centering\small
\caption{Profile values at selected distance ratios (minimum $\chi^2$ and maximum $R^2$ over the remaining parameters, transforms included).}\label{tab:profile}
\begin{tabular}{c ccc ccc ccc}\toprule
 & \multicolumn{3}{c}{exact point mass} & \multicolumn{3}{c}{exact extended} & \multicolumn{3}{c}{first order}\\\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-10}
$D/r$ & $\chi^2_{\rm corr}$ & $\chi^2_{\rm naive}$ & $R^2$ & $\chi^2_{\rm corr}$ & $\chi^2_{\rm naive}$ & $R^2$ & $\chi^2_{\rm corr}$ & $\chi^2_{\rm naive}$ & $R^2$\\\midrule
""")
for Dr in (1.05, 1.1, 1.15, 1.2, 1.3, 1.44, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0):
    i = int(np.argmin(np.abs(DRg - Dr)))
    w(f"${f(DRg[i],2)}$ & " + " & ".join(f"${f(max(v,0),3)}$" for v in (pmc[i], pmn[i])) + f" & ${f(pmr[i],3)}$ & " + " & ".join(f"${f(max(v,0),3)}$" for v in (exc[i], exn[i])) + f" & ${f(exr[i],3)}$ & " + " & ".join(f"${f(max(v,0),3)}$" for v in (sc[i], sn[i])) + f" & ${f(sr[i],3)}$\\\\")
w(r"""\bottomrule
\end{tabular}
\end{table}
""")
w(r"""
\begin{figure}[!t]\centering
\includegraphics[width=\textwidth]{../../figures/exact_profile.pdf}
\caption{Continuous fit of the distance ratio. (a) Profile amplitude $\chi^2$ for the exact point mass, the exact extended density and the first-order form, with sector-correlated and naive errors; dotted lines mark $\Delta\chi^2=1$ and $4$. (b) Profile $R^2$ of the transformed curve against the digitised points; the dashed line is the data's own Fourier fit. (c) Exact four-harmonic solutions of the extended model: $D/L$ against $D/r$ for the two dipole branches; the vertical lines are the point-mass roots, which the extended family approaches as $D/L\to\infty$.}\label{fig:profile}
\end{figure}
""")
# best-fit parameter table
B = P["best"]
bests = [("Exact pm, $\\chi^2$ root (on-axis)", fp[0] if fp[0]["tau"]["1"] > 0 else fp[1]), ("Exact pm, $\\chi^2$ root (flipped)", fp[1] if fp[0]["tau"]["1"] > 0 else fp[0]),
         ("Exact pm, max $R^2$", B["negR2"]), ("Exact extended, max $R^2$", B["ext_negR2"])]
w(r"""
\begin{table}[!htb]\centering\small
\caption{Best-fit parameters and predicted harmonics. Pulls are $(A_m^{\rm mod}-A_m)/\sigma_m$ with sector-correlated (c) and naive (n) errors; the $\chi^2$ roots match all amplitudes exactly.}\label{tab:bestpar}
\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{lcccccccccc}\toprule
Solution & $D/r$ & $b_0$ & $D/L$ & $A_1$ & $A_2$ & $A_3$ & pull $A_3$ (c) & pull $A_3$ (n) & $\chi^2_{\rm corr}$ & $R^2$\\\midrule
""")
for lab, s in bests:
    dl = f"${f(s['D_over_L'],1)}$" if s.get("D_over_L") else "$\\infty$"
    w(f"{lab} & ${f(s['D_over_r'],3)}$ & ${f(s['b0'],2)}^\\circ$ & {dl} & ${f(s['model_amp']['1'])}$ & ${f(s['model_amp']['2'])}$ & ${f(s['model_amp']['3'])}$ & ${sg(s['pull']['corr']['3'])}$ & ${sg(s['pull']['naive']['3'])}$ & ${f(s['chi2']['corr'],3)}$ & ${f(s['R2'],3)}$\\\\")
w(r"""\bottomrule
\end{tabular}
\end{table}
""")
# fixed D/r reference rows
w(r"""
\subsection{Reference rows at fixed distance ratio}
Table~\ref{tab:fixed} repeats the fixed-$D/r$ least-squares fits of the earlier sheet with the exact point-mass form, using the earlier weights. The dipole and quadrupole are still reproduced almost exactly, and the whole misfit remains in the octupole. From $D/r=2$ upward the exact and first-order octupole predictions agree closely, so the octupole shortfall found before is a property of the physics and not of the expansion. Only close to $D/r=1$ does the exact form raise the predicted octupole appreciably. The comparison is shown in Fig.~\ref{fig:exact}c.

\begin{table}[!htb]\centering\small
\caption{Fixed-$D/r$ least-squares fits with the exact point-mass form (weights from the naive errors, as in the earlier sheets; best branch). Pulls use the sector-correlated (c) and naive (n) errors. The last column is the first-order octupole at the same $D/r$.}\label{tab:fixed}
\setlength{\tabcolsep}{3.2pt}
\begin{tabular}{ccccccccccccc}\toprule
$D/r$ & branch & $b_0$ & $N$ & $K$ & $A_1$ & $A_2$ & $A_3$ & pulls $A_{1,2,3}$ (c) & pull $A_3$ (n) & $R^2$ & $A_3$ first order\\\midrule
data & & & & & $""" + f(A["1"]) + "$ & $" + f(A["2"]) + "$ & $" + f(A["3"]) + r"""$ & & & $""" + f(0.915048, 3) + r"""$ & \\
""")
for Dr in ("1.3", "2.0", "3.0", "4.0"):
    s = R["fixed_pm"][Dr]["lsq_naive"][0]; ss = R["fixed_pm"][Dr]["series_lsq_naive"][0]
    w(f"${float(Dr):g}$ & {s['branch']} & ${f(s['b0'],2)}^\\circ$ & ${f(s['N'])}$ & ${f(s['K'])}$ & ${f(s['model_amp']['1'])}$ & ${f(s['model_amp']['2'])}$ & ${f(s['model_amp']['3'])}$ & "
      f"${sg(s['pull']['corr']['1'])},{sg(s['pull']['corr']['2'])},{sg(s['pull']['corr']['3'])}$ & ${sg(s['pull']['naive']['3'])}$ & ${f(s['R2'],3)}$ & ${f(ss['model_amp']['3'])}$\\\\")
w(r"""\bottomrule
\end{tabular}
\end{table}
""")
w(r"""
\begin{figure}[!t]\centering
\includegraphics[width=\textwidth]{../../figures/exact_termmatch.pdf}
\caption{Exact versus first-order point-mass model. (a, b) Ratio of exact to first-order sector-averaged coefficients against $D/r$ for an in-plane axis and for a tilt near the free solution. (c) Octupole predicted after the monopole, dipole and quadrupole are matched exactly, for the exact form (red) and the first-order form (blue), on-axis (solid) and flipped (dashed) branches; grey bands show the data octupole with its sector-correlated and naive $1\sigma$ ranges; markers show the free-$D/r$ roots. (d) Digitised points, the data Fourier fit and exact point-mass curves with the transforms, for the exact on-axis root and the fixed-$D/r$ reference fits.}\label{fig:exact}
\end{figure}
""")
# extended
fam = R["ext_family"]
w(r"""
\subsection{Extended density}
With the exponential density the extended mass restores a dipole that is larger relative to the quadrupole than for a point mass, and the scale length adds one unknown. At fixed $D/r$ the extended problem therefore has as many unknowns as equations. Exact solutions exist, but only on the family shown in Fig.~\ref{fig:profile}c and listed in Table~\ref{tab:ext}. That family runs from the point-mass roots, reached as $D/L\to\infty$, towards smaller $D/r$. No exact extended solution was found at larger distance ratios. There the least-squares optimum sits at the upper bound of the scale-length range, which is the point-mass limit, except for a narrow range where a finite scale length gives a marginally lower $\chi^2$. A finite scale length therefore does not relieve the octupole shortfall at large $D/r$. Along the exact family the curve-fit $R^2$ is lower than at the $R^2$-optimised solutions, again because of the higher harmonics and the octupole phase.

\begin{table}[!htb]\centering\small
\caption{Extended (exponential) density. Upper block: exact four-harmonic solutions ($N,K,b_0,D/L$ at fixed $D/r$). Lower block: least-squares optima where no exact solution exists ($D/L$ restricted to the scanned range, whose upper end is effectively a point mass).}\label{tab:ext}
\begin{tabular}{lcccccccc}\toprule
Type & $D/r$ & branch & $b_0$ & $D/L$ & $K$ & $\alpha$ & $A_3$ & $R^2$\\\midrule
""")
for s in sorted(fam, key=lambda s: (s["branch"], s["D_over_r"])):
    if s["D_over_r"] in (1.01, 1.05, 1.09, 1.11, 1.13, 1.15) or (s["branch"] == "flipped" and s["D_over_r"] in (1.01, 1.05, 1.09, 1.11)):
        w(f"exact & ${f(s['D_over_r'],2)}$ & {s['branch']} & ${f(s['b0'],2)}^\\circ$ & ${f(s['D_over_L'],1)}$ & ${f(s['K'])}$ & ${f(s['alpha'],3)}$ & ${f(s['model_amp']['3'])}$ & ${f(s['R2'],3)}$\\\\")
w(r"\midrule")
for Dr in ("1.3", "1.5", "2.0", "3.0", "4.0", "6.0"):
    s = R["fixed_ext"][Dr]["lsq"]
    w(f"LSQ & ${float(Dr):g}$ & {'on-axis' if s['tau']['1']>0 else 'flipped'} & ${f(s['b0'],2)}^\\circ$ & ${f(s['D_over_L'],1)}$ & ${f(s['K'])}$ & ${f(s['alpha'],3)}$ & ${f(s['model_amp']['3'])}$ & ${f(s['R2'],3)}$\\\\")
w(r"""\bottomrule
\end{tabular}
\end{table}
""")
# Chris's variables
w(r"""
\section{Translation to Chris's variables}
The fitted parameters are translated back with
\begin{equation}
D=r\,(D/r),\qquad \frac{M(<D)}{M_{\rm crit}}=2\kappa=\frac{2K}{\alpha},\qquad L=\frac{D}{x},\qquad \varphi=90^\circ-\ell_1,\qquad \varphi_2=90^\circ-2\ell_2,
\label{eq:chris}
\end{equation}
where the mass ratio is given both for $\alpha=1$ (the convention of the earlier sheets) and for the ideal-drag value of Eq.~\eqref{eq:drag}. Because the transforms place the dipole and quadrupole peaks on the data, $\varphi$ and $\varphi_2$ are the same for every solution. The results are in Table~\ref{tab:chris}.

\begin{table}[!htb]\centering\small
\caption{Best-fit solutions and fixed-$D/r$ reference fits in Chris's variables ($r=14.3$ Gpc). For all rows $\ell_0=""" + f(l0, 2) + r"""^\circ$, $\varphi=""" + f(fp[0]["phiC"], 2) + r"""^\circ$ and $\varphi_2=""" + f(fp[0]["phi2C"], 2) + r"""^\circ$.}\label{tab:chris}
\setlength{\tabcolsep}{3.3pt}
\begin{tabular}{lccccccccc}\toprule
Solution & $D/r$ & $D$ [Gpc] & $D$ [Gly] & $K$ & $\frac{M}{M_{\rm crit}}$ ($\alpha{=}1$) & $\frac{M}{M_{\rm crit}}$ (drag) & $b_0$ & $\phi_1$ & $L$ [Gpc]\\\midrule
""")
ALPM = R["alpha_point_mass_numeric"]
rows = bests + [(f"Exact pm, fixed $D/r={float(Dr):g}$", R["fixed_pm"][Dr]["lsq_naive"][0]) for Dr in ("2.0", "3.0", "4.0")]
for lab, s in rows:
    al = s.get("alpha", ALPM); Lg = f"${f(s['D_Gpc']/s['D_over_L'],3)}$" if s.get("D_over_L") else "--"
    w(f"{lab} & ${f(s['D_over_r'],3)}$ & ${f(s['D_Gpc'],2)}$ & ${f(s['D_Gly'],1)}$ & ${f(s['K'])}$ & ${f(2*s['K'],3)}$ & ${f(2*s['K']/al,3)}$ & ${f(s['b0'],2)}^\\circ$ & ${sg(s['phi1'],2)}^\\circ$ & {Lg}\\\\")
w(r"""\bottomrule
\end{tabular}
\end{table}

\section{Conclusions}
Evaluating the tidal term exactly changes the harmonic coefficients substantially only when the observer is within a few vision-shell radii of the CoE (Table~\ref{tab:ratio}). Farther out, the first-order form was an adequate description. With four unknowns, an exact match of all four measured harmonics now exists with the observer outside the vision shell, on both dipole branches (Table~\ref{tab:free}). The preferred distance ratio therefore moves from essentially unity to slightly above it, but it stays close to one.

Treating $D/r$ as a continuous parameter does not tighten this much (Table~\ref{tab:bestDr}, Fig.~\ref{fig:profile}). The amplitude $\chi^2$ grows only slowly with distance, and even the four-unit threshold is not crossed within the allowed range, so the octupole alone cannot exclude large distance ratios even at the level of the wider interval. The curve-fit $R^2$ is highest slightly above unity and degrades gradually, and the differences are small compared with what overlapping sectors can support.

The octupole tension at fixed $D/r\ge2$ is essentially the same with the exact form as with the series (Table~\ref{tab:fixed}, Fig.~\ref{fig:exact}c). It is not an artefact of the expansion. An extended exponential density does not rescue large distance ratios either: its extra freedom produces exact solutions only between the point-mass roots and $D/r=1$ (Table~\ref{tab:ext}).

Close to $D/r=1$ the exact model develops harmonics above the octupole, which the four-harmonic matching ignores. As a result, the exactly matching solutions have a lower $R^2$ than the $R^2$-optimised ones. Caveats carried over from the earlier sheets remain. The rotation and $2\ell$ phases are free rather than predicted. The error model assumes uniform clusters and statistical, overlap-correlated error bars read from a digitised figure. The sector averages assume a uniform sky outside the plane cut. Only a single-exponential density family was tried for the extended mass.
\end{document}
""")
open("/workspace/cosmology/papers/A_equals_fourier/exact_termmatch.tex", "w").write("\n".join(L))
print("written")
