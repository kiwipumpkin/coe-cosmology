# Two-column equation card for the SR + central-mass treatment, sized to fit one preprint page.
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif"})
left = [
 r"\textbf{Geometry and residual}",
 r"$c=1,\quad d_1=1,\quad \theta=\ell-\phi$",
 r"$G(\ell)=\dfrac{A/A_{\rm all}(\ell)}{\langle A/A_{\rm all}\rangle}$",
 r"$S_{\rm req}=\dfrac{1/G}{\langle 1/G\rangle}$",
 r"$R(\ell)=\dfrac{G\,S_{\rm model}}{\langle G\,S_{\rm model}\rangle}-1$",
 r"\textbf{Positions and speeds}",
 r"$\mathbf{x}_{e0}=(1+q\cos\theta,\ q\sin\theta)$",
 r"$r_{e0}=\sqrt{1+q^2+2q\cos\theta}$",
 r"$r_o=s,\qquad r_e=s\,r_{e0}$",
 r"$v^2=\beta_1^2 r_0^2+2\mu(\dfrac{1}{r}-\dfrac{1}{r_0})$",
 r"$\boldsymbol{\beta}_o=v_o\,\hat{\mathbf{x}},\qquad \boldsymbol{\beta}_e=v_e\,\mathbf{x}_e/r_e$",
]
right = [
 r"\textbf{Light path}",
 r"$\mathbf{d}=\mathbf{x}_o-\mathbf{x}_e$",
 r"$(v_o^2-1)\,t_f^2+2d_x v_o\,t_f+|\mathbf{d}|^2=0,\ \ t_f>0$",
 r"$\hat{\mathbf{n}}=(\mathbf{d}+v_o t_f\,\hat{\mathbf{x}})/t_f$",
 r"\textbf{Redshift factors}",
 r"$D_{\rm motion}=\dfrac{\gamma_o(1-\boldsymbol{\beta}_o\cdot\hat{\mathbf{n}})}{\gamma_e(1-\boldsymbol{\beta}_e\cdot\hat{\mathbf{n}})}$",
 r"$D_{\rm grav}=\sqrt{\dfrac{1-2\mu/r_e}{1-2\mu/r_o}}$",
 r"$S_{\rm model}=D_{\rm motion}\,D_{\rm grav}$",
 r"\textbf{Run 2} ($s\leq 3$)",
 r"$\beta_1=0.5866,\ q=0.2303,\ \mu=0.09332$",
 r"$s=3.000\ \mathrm{(bound)},\ \phi=15.16^\circ$",
]
bottom = (r"$\mathrm{RMS}(G-1)=5.75\times10^{-2}\ \rightarrow\ \mathrm{RMS}(R)=2.37\times10^{-3}$"
          r"  (95.9% of starting RMS, 99.83% of variance)")
bottom2 = (r"$R(\ell)\approx\sum_m c_m\cos m(\ell-\phi)\ \ (|s_m|\sim10^{-8})$;"
           r"  after solving $m=1$–$12$, leftover RMS $=1.83\times10^{-6}$")
fig = plt.figure(figsize=(7.5, 6.6))
fig.text(0.5, 0.975, "Special relativity + central mass (CoE shell)", ha="center", va="top", fontsize=15, fontweight="bold")
def col(lines, x):
    y = 0.90
    for L in lines:
        tall = r"\dfrac" in L
        if L.startswith(r"\textbf") and y < 0.89: y -= 0.012
        bold = L.startswith(r"\textbf")
        txt = re.sub(r"\\textbf\{(.*?)\}", r"\1", L) if bold else L
        fig.text(x, y, txt, ha="left", va="top", fontsize=13, fontweight="bold" if bold else "normal")
        y -= 0.095 if tall else 0.058
    return y
y1 = col(left, 0.03); y2 = col(right, 0.52)
yb = min(y1, y2) - 0.01
fig.text(0.03, yb, bottom, ha="left", va="top", fontsize=12)
fig.text(0.03, yb-0.05, bottom2, ha="left", va="top", fontsize=12)
fig.savefig("sr-treatment-eqs-2col.png", dpi=200, bbox_inches="tight", pad_inches=0.12)
