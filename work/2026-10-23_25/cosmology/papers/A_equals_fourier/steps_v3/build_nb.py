import nbformat as nbf
C = []
def md(s): C.append(nbf.v4.new_markdown_cell(s.strip()))
def py(s): C.append(nbf.v4.new_code_cell(s.strip()))

md(r"""
# CoE shell model vs. Migkas cluster anisotropy — every step (v3)

This notebook reproduces the v3 term-by-term sheet (`thinshell_termmatch_equations_table_v3.pdf`, Sections 1–4) and the Migkas $\langle H\rangle$ anchor, one step at a time. Each step has a short explanation and then the code that does it.

**Rules used throughout**
* Gravity only, $\Lambda=0$, no dark matter. All parameters are free, with no age, boundedness or ΛCDM limits.
* Chris's $A$: $a_f=|\Delta\vec V|/\text{separation}$, used linearly. It is never squared and never turned into $D_L^{-2}$.
* Symbols: $D$ = CoE–Earth distance; $R$ = the 14.3 Gpc vision sphere; shell edges $R_{\rm in}$, $R_{\rm out}$. Lengths are in Gpc.
* **$H_0=70$ is not an input.** It appears in only two labelled places: (i) Section 3's Table 1 masses, which were originally computed with Earth's rate set to 70, and are shown unrescaled with a label; (ii) the unit convention of the LTB integration, which is then rescaled to the Migkas-anchored $\langle H\rangle$.

**Sources.** Taylor, G. I. 1950, Proc. R. Soc. Lond. A 201, 159; Sedov, L. I. 1959, *Similarity and Dimensional Methods in Mechanics*; Migkas, K. et al. 2020, A&A 636, A15 (arXiv:2004.03305); Migkas, K. et al. 2021, A&A 649, A151 (arXiv:2103.13904); Planck Collaboration 2020, A&A 641, A6 and A1 (Planck 2018 VI and I); Riess, A. G. et al. 2022, ApJL 934, L7 (SH0ES); Watkins, R. et al. 2023, MNRAS 524, 1885 (CosmicFlows-4 bulk flow).

The heavy integrators are Chris's code, used unchanged: `ptshell.py`/`prof_core.py` (Newtonian) and `relfree.py`/`rm_core.py` (LTB). Everything else is written out in the cells.
""")
py(r"""
# Step 0 - paths, Chris's integrators (unchanged), plotting
import sys, os, json, csv, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
%matplotlib inline
BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"
THS = "/workspace/cosmology/papers/A_equals_fourier/thinshell"
sys.path.insert(0, BUB); sys.path.insert(0, BUB + "/gr_ltb/refit_mass"); os.chdir(BUB)
import ptshell as S, prof_core as P                 # Newtonian shell + core integrator (Chris's code)
np.set_printoptions(precision=4, suppress=True)
TH1 = np.arange(0, 360, 1.0)                        # Galactic longitude grid, 1 deg (the code's 360 directions)
RV = 14.30048287260589                              # R, vision sphere [Gpc]
G = 6.674e-11; c = 299792.458; MPC = 3.0857e19; GPCm = 3.0857e25; MSUN = 1.989e30
def rhoc(H_kms_Mpc):                                # critical density 3H^2/(8 pi G) [kg/m^3]
    h = H_kms_Mpc / MPC; return 3 * h * h / (8 * np.pi * G)
print("R =", round(RV, 2), "Gpc;", len(TH1), "directions")
""")

md(r"""
## Step 1 — Model setup (Section 1)

**Flow.** Sedov–Taylor self-similar expansion (Taylor 1950; Sedov 1959). Every point moves away from the Centre of Expansion (CoE) as
$$\vec V=k\,\vec r\qquad(\vec r=\text{position from the CoE}).$$
**Vector subtraction.** An observer O (Earth) sees a source P move at
$$\vec V_P-\vec V_O=k(\vec r_P-\vec r_O)=k\,\vec s,$$
so with no mass present, $|\Delta\vec V|/s=k$ in every direction and the sky is isotropic. The anisotropy comes only from gravity.

**Mass.** A point-like core $M_c$ at the CoE, plus a **uniform spherical shell** of ordinary matter with density $\rho_s=\Omega_s\rho_c$ between $R_{\rm in}$ and $R_{\rm out}$ (both centred on the CoE and fixed). Earth sits at distance $D$ from the CoE, and all sources are on the $R=14.3$ Gpc sphere around Earth, inside the shell. Each fit keeps its own shell rule: Newtonian $\Omega_s$ free uses $R_{\rm in}=D-R$, $R_{\rm out}=D+R$; Newtonian $\Omega_s=0.049$ uses $D-1.028R$, $D+1.475R$; LTB uses $D-1.040R$, $D+1.171R$.
""")
py(r"""
# Step 1a - with no mass, V = k r gives |dV|/s = k in every direction (vector subtraction)
k = 1.0; D = 3.0; rO = np.array([D, 0.0])                     # Earth at distance D from the CoE (units arbitrary)
th = np.radians(TH1); rP = rO + RV / 14.3 * np.c_[np.cos(th), np.sin(th)]
dV = k * rP - k * rO                                          # V_P - V_O
af = np.linalg.norm(dV, axis=1) / np.linalg.norm(rP - rO, axis=1)
print("pure Sedov-Taylor: af/k min, max =", af.min(), af.max())
""")
py(r"""
# Step 1b - the three fitted models (parameters from the fit files; all free, nothing imposed)
T1 = json.load(open(BUB + "/universe_thin.json"))            # Table 1 of universe_bestfit_thinshell.pdf
hug = json.load(open(BUB + "/thin_b_free_hug.json"))          # Newtonian, Omega_s free
omb = json.load(open(BUB + "/universe_newt_omb.json"))["Rs14.3|sine"]   # Newtonian, Omega_s = 0.049
REL = json.load(open(BUB + "/universe_rel_mcfree.json"))      # relativistic LTB, all free
NEWT = {
 "free": dict(DR=hug["DR"], Os=hug["Os"], mcc=hug["mcc"], rin=hug["rin"], rout=hug["rout"], Lc=0.005, mc=0.0, Ob=0.0, l0=hug["l0"], amp=hug["amp"], phi=hug["phi"]),
 "omb":  dict(DR=omb["D"] / RV, Os=0.049, mcc=omb["mcc"], rin=omb["R_in"] / omb["D"], rout=omb["R_out"] / omb["D"], Lc=0.005, mc=0.0, Ob=0.0, l0=omb["l0"], amp=omb["amp"], phi=omb["phi"]),
}
for kk, a in NEWT.items():
    D = a["DR"] * RV
    print(f"{kk:5s} D = {D:6.2f} Gpc  D/R = {a['DR']:.3f}  R_in = {a['rin']*D:6.1f}  R_out = {a['rout']*D:6.1f} Gpc  "
          f"Omega_s = {a['Os']:.3f}  Mc/Mcrit(D) = {a['mcc']:.5f}")
print(f"LTB   D = {REL['D']:6.2f} Gpc  D/R = {REL['D']/RV:.3f}  R_in = {REL['R_in']:6.1f}  R_out = {REL['R_out']:6.1f} Gpc  "
      f"Omega_s = 0.049  Mc/Mcrit(D) = {REL['mcc']:.4f}  (k_in {REL['kin']:.3f}, k_out {REL['kout']:.3f})")
""")

md(r"""
## Step 2 — Chris's $A$, exactly as his code computes it

The code integrates the flow from $\tau_i=0.01$ to today (time unit $t_0$, length unit $D$), with Earth and 360 sources on the $R$ sphere (one per degree of Galactic longitude, in the Galactic plane). For each source,
$$a_f(\ell)=\frac{|\vec V_P-\vec V_O|}{s_{\rm nominal}},\qquad A(\ell)=\frac{a_f(\ell)}{\langle a_f\rangle_{\rm sky}} .$$
This is compared **linearly** with Migkas's $A/A_{\rm all}$: it is not squared, and not converted to $D_L^{-2}$ or $H^2$. In code units $1/t_0=k$, so $\langle a_f\rangle/k$ is dimensionless. The code also returns Earth's radial rate $H_E/k$ (`Ht`) and speed $v_E$ (units $D/t_0$).
""")
py(r"""
# Step 2 - run Chris's integrator (inlined from prof_core.run) and form af = |dV| / separation
def run_newt(a, taui=0.01, nsteps=200, iters=3):
    _, _, m1 = S.build(a)                                   # core + uniform shell enclosed-mass table (ptshell.build)
    fac = max(1.5, 1.02 * a["rout"] / (1 + 1 / a["DR"]))
    q, rr = P.geometry(TH1, a["DR"]); Rb = P.rb_of(a["DR"], fac); Ht = 0.95
    _, MT = P.tables("ptshell", a, Rb)
    for _ in range(iters):                                   # self-consistent Earth rate H_E t0
        mu = m1 * Ht**2 / 2
        X, V, VR = P.integrate(q, mu, P.RG, MT, 0.0, Rb, taui, nsteps)
        Ht = np.hypot(*V[0]) / np.hypot(*X[0])
    dV = V[1:] - V[0]                                        # V_P - V_O
    af = np.hypot(dV[:, 0], dV[:, 1]) / rr                   # / nominal separation (R in units of D)
    return dict(af=af, Ht=Ht, vE=np.hypot(*V[0]), rE=np.hypot(*X[0]), mu=mu, MT=MT, Rb=Rb)
RUN = {kk: run_newt(a) for kk, a in NEWT.items()}
for kk, r in RUN.items():
    A = r["af"] / r["af"].mean()
    print(f"{kk:5s} <af>/k = {r['af'].mean():.4f}  A min/max = {A.min():.3f}/{A.max():.3f}  H_E/k = {r['Ht']:.4f}  v_E = {r['vE']:.4f} D/t0")
""")

md(r"""
## Step 3 — Migkas data: $A/A_{\rm all}$ vs. $\ell$ and the Fourier form $A_F$

The points are the digitised Migkas et al. 2020 Fig. 6 profile, $A/A_{\rm all}$ against Galactic longitude (71 points). $A_F$ is its 3-harmonic Fourier fit:
$$A_F(\ell)=0.9965+0.0296\cos\ell+0.0434\sin\ell+0.0636\cos2\ell+0.0366\sin2\ell-0.0299\cos3\ell-0.0114\sin3\ell .$$
""")
py(r"""
# Step 3 - load the Migkas profile and define A_F
mig = list(csv.DictReader(open("/workspace/cosmology/migkas/fig6_digitized.csv")))
ML = np.array([float(r["l_deg"]) for r in mig]); MA = np.array([float(r["A_over_Aall"]) for r in mig]); IMIG = ML.astype(int)
AFc = {0: (0.9965, 0.0), 1: (0.0296, 0.0434), 2: (0.0636, 0.0366), 3: (-0.0299, -0.0114)}      # (cos, sin) coefficients
def AFfour(l):
    l = np.radians(l); return sum(cs * np.cos(m * l) + sn * np.sin(m * l) for m, (cs, sn) in AFc.items())
def r2(y, m): return 1 - np.sum((y - m) ** 2) / np.sum((y - y.mean()) ** 2)
print(len(ML), "points; l range", ML.min(), "-", ML.max(), "; R2(A_F vs points) =", round(r2(MA, AFfour(ML)), 3))
""")

md(r"""
## Step 4 — The diagnostic sine term and the headline $R^2$

The geometry alone gives only $\cos m\psi$ terms about the CoE axis ($\psi=\ell-\ell_0$). A single rotation/tilt factor absorbs the sine parts:
$$A(\ell)=a_f(\ell-\ell_0)\,\big[1-\mathrm{Amp}\,\sin(\psi+\varphi)\big].$$
It is a diagnostic factor, not a velocity. The headline number is $R^2$ of this full curve against the 71 Migkas points. Each fit's $\ell_0$, Amp and $\varphi$ are its Table 1 values. The LTB curve comes from the relativistic ray-tracing code (heavy cell below, about 1 s).
""")
py(r"""
# Step 4a - Newtonian: shift to l0, apply the sine factor, R^2 vs Migkas
def shifted(af, l0): return np.interp((TH1 - l0) % 360, np.append(TH1, 360), np.append(af, af[0]))
def A_from(af, th, amp, phi): return (1 - amp * np.sin(np.radians(th) + phi)) * af
CURVE, GEO = {}, {}
for kk, a in NEWT.items():
    g = RUN[kk]["af"] / RUN[kk]["af"].mean(); GEO[kk] = shifted(g, a["l0"])
    CURVE[kk] = A_from(GEO[kk], TH1 - a["l0"], a["amp"], a["phi"])
    print(f"{kk:5s} l0 = {a['l0']:.1f}  Amp = {a['amp']:.4f}  phi = {a['phi']:.3f}   R2 vs Migkas = {r2(MA, CURVE[kk][IMIG]):.3f} "
          f"(geometry only {r2(MA, GEO[kk][IMIG]):.3f})   R2 vs A_F = {r2(AFfour(TH1), CURVE[kk]):.3f}")
""")
py(r"""
# Step 4b - HEAVY CELL (relativistic LTB): Chris's LTB solver + ray tracing at the Table-1 parameters (~1 s)
# Unit note: the LTB code works in units where the 3D-mean local expansion rate is 70 km/s/Mpc (its internal convention);
# only dimensionless ratios are used from it here, and Section 4 rescales to the Migkas-anchored <H>.
import relfree as RF, rm_core as RC
ltb = RF.ev(REL["Mc_kg"], REL["D"], REL["kin"], REL["kout"], full=True)
Ad = np.array(ltb["Ad"])
def ltb_curve(Ad, l0, amp=0.0, phi=0.0):                     # inlined rm_core.curve: 0.5-deg ray table -> 1-deg grid, normalise, sine
    idx = (((TH1 - l0 - 180.0) % 360) / 0.5).astype(int); A = Ad[idx]; A = A / A.mean()
    return A_from(A, TH1 - l0, amp, phi) if amp else A
GEO["rel"] = ltb_curve(Ad, ltb["l0_s"]); CURVE["rel"] = ltb_curve(Ad, ltb["l0_s"], ltb["amp"], ltb["phi"])
json.dump(dict(ltb={k: v for k, v in ltb.items() if k != "Ad"}, A=CURVE["rel"].tolist()), open("/workspace/cosmology/papers/A_equals_fourier/steps_v3/ltb_cache.json", "w"), default=float)
print(f"LTB   l0 = {ltb['l0_s']:.1f}  Amp = {ltb['amp']:.4f}  phi = {ltb['phi']:.3f}   R2 vs Migkas = {r2(MA, CURVE['rel'][IMIG]):.3f} "
      f"(geometry only {r2(MA, GEO['rel'][IMIG]):.3f})   R2 vs A_F = {r2(AFfour(TH1), CURVE['rel']):.3f}")
print(f"      H_perp = {ltb['H_perp']:.2f}, H_r = {ltb['H_r']:.2f} (fit units), age = {ltb['age_Gyr']:.2f} Gyr (fit units)")
""")
py(r"""
# Step 4c - figure: the three full curves vs the Migkas points
fig, ax = plt.subplots(figsize=(9, 4.2))
ax.plot(ML, MA, "ko", ms=3.5, label="Migkas 2020 Fig. 6 (A/A_all)")
ax.plot(TH1, AFfour(TH1), "k--", lw=1, label="A_F (Fourier)")
for kk, col, lab in (("free", "C0", "Newt. Ω_s free"), ("omb", "C1", "Newt. Ω_s=0.049"), ("rel", "C2", "LTB")):
    ax.plot(TH1, CURVE[kk], col, label=f"{lab}  R²={r2(MA, CURVE[kk][IMIG]):.3f}")
ax.set_xlabel("Galactic longitude ℓ (deg)"); ax.set_ylabel("A / ⟨A⟩"); ax.set_xlim(0, 360); ax.legend(fontsize=8, ncol=2); plt.tight_layout(); plt.show()
""")

md(r"""
## Step 5 — Section 2: term-by-term solution against $A_F$

To first order in $\varepsilon=R/D$, with $h(q)=V/(kq)$ the radial rate profile, $\beta=h'/h$ and $\gamma=h''/h$ at Earth:
$$a_f^2\propto1+B\mu^2+\varepsilon(P\mu+Q\mu^3),\quad B=2\beta+\beta^2,\;P=\beta(3+\beta),\;Q=(1+\beta)(\gamma+\beta)-2\beta .$$
Taking the square root and projecting onto harmonics gives
$$\frac{a_f}{\langle a_f\rangle}=1+a_1\cos\psi+a_2\cos2\psi+a_3\cos3\psi,\quad a_2=\frac{B/4-B^2/16}{m_0},\; m_0=1+\tfrac B4-\tfrac{3B^2}{64},$$
$$a_1=\frac{\varepsilon}{2m_0}\Big[P+\tfrac34\big(Q-\tfrac{PB}{2}\big)-\tfrac{5QB}{16}\Big],\qquad a_3=\frac{\varepsilon}{2m_0}\Big[\tfrac14\big(Q-\tfrac{PB}{2}\big)-\tfrac{5QB}{32}\Big].$$
The sine factor adds $C_1=a_1-\mathrm{Amp}\sin\varphi(1+a_2/2)$, $S_1=-\mathrm{Amp}\cos\varphi(1-a_2/2)$, $C_3=a_3-\tfrac{\mathrm{Amp}\,a_2}{2}\sin\varphi$, $S_3=-\tfrac{\mathrm{Amp}\,a_2}{2}\cos\varphi$.
Solve order: 2ℓ phase → $\ell_0$; 2ℓ amplitude → $\beta$; $C_1,S_1,C_3$ → Amp, $\varphi$, $\varepsilon$ (so $D=R/\varepsilon$); $S_3$ is then a check. Point-core shape: $\gamma=-4\beta$. Mass: $\beta\simeq3K/(1-K)$, $K=\tfrac12\frac{M_c}{M_{\rm crit}(D)}(\tfrac1{\tau_i}-1)$.
""")
py(r"""
# Step 5 - sequential (hand) solve with the series, then exact least squares on C_m, S_m
from scipy.optimize import brentq, fsolve, least_squares
TH = np.radians(TH1)
def harm(y, m): return (y.mean(), 0.0) if m == 0 else (2 * np.mean(y * np.cos(m * TH)), 2 * np.mean(y * np.sin(m * TH)))
def rot(cs, sn, m, l0):                 # coefficients about l -> about psi = l - l0
    t = np.radians(m * l0); return cs * np.cos(t) + sn * np.sin(t), -cs * np.sin(t) + sn * np.cos(t)
def series(beta, gam, eps):
    B = 2 * beta + beta**2; Pp = beta * (3 + beta); Q = (1 + beta) * (gam + beta) - 2 * beta
    m0 = 1 + B / 4 - 3 * B**2 / 64; m2 = B / 4 - B**2 / 16
    m1 = eps / 2 * (Pp + 0.75 * (Q - Pp * B / 2) - Q * B / 2 * 10 / 16); m3 = eps / 2 * (0.25 * (Q - Pp * B / 2) - Q * B / 2 * 5 / 16)
    return dict(a1=m1 / m0, a2=m2 / m0, a3=m3 / m0, Q=Q)
C2, S2 = AFc[2]; l0 = (np.degrees(np.arctan2(S2, C2)) / 2) % 180                     # 2-l phase -> axis
Dt = {m: rot(*AFc[m], m, l0) for m in range(4)}; C0 = Dt[0][0]; (C1, S1), (C2r, _), (C3, S3) = Dt[1], Dt[2], Dt[3]
gb = -4.0                                                                             # point core: gamma = -4 beta
beta = brentq(lambda b: series(b, gb * b, 0)["a2"] - C2r / C0, 1e-4, 1.0)
s1 = series(beta, gb * beta, 1.0); a2, g1, g3 = s1["a2"], s1["a1"], s1["a3"]
def F(x):
    amp, phi, eps = x
    return [C0 * (eps * g1 - amp * np.sin(phi) * (1 + a2 / 2)) - C1, -C0 * amp * np.cos(phi) * (1 - a2 / 2) - S1, C0 * (eps * g3 - amp * a2 / 2 * np.sin(phi)) - C3]
amp, phi, eps = fsolve(F, [-0.06, 0.9, 0.3])
S3_pred = -C0 * amp * a2 / 2 * np.cos(phi); N = C0 / (1 - amp * eps * g1 * np.sin(phi) / 2)
K = beta / (3 + beta); mcc_lin = 2 * K / (1 / 0.01 - 1)
print(f"2l phase  : l0 = {l0:.2f} deg      2l amplitude: C2 = {np.hypot(*AFc[2]):.4f} -> beta = {beta:.3f}")
print(f"1l, 3l cos: Amp = {amp:.3f}, phi = {phi % (2*np.pi):.2f}, eps = {eps:.3f} -> D = R/eps = {RV/eps:.0f} Gpc, eps*Q = {eps*s1['Q']:.3f}")
print(f"3l sin    : predicted S3 = {S3_pred:.4f} vs data {S3:.4f} (not produced)     m=0: N = {N:.4f}")
print(f"mass      : K = {K:.4f} -> Mc/Mcrit(D) = {mcc_lin:.4f} (linear)")
# exact least squares: exact radial field with h(q) = 1 + beta (q-1) + (gamma/2)(q-1)^2, gamma = -4 beta
def A_radial(hfun, eps):
    x = 1 + eps * np.cos(TH); y = eps * np.sin(TH); qs = np.hypot(x, y); hs = hfun(qs); h1 = hfun(np.array([1.0]))[0]
    af = np.hypot(hs * x - h1, hs * y) / eps; return af / af.mean()
def model(p):
    l0_, b, e, am, ph = p; g = A_radial(lambda q: 1 + b * (q - 1) + 0.5 * gb * b * (q - 1)**2, e); return A_from(shifted(g, l0_), TH1 - l0_, am, ph)
def res(p):
    A = model(p); out = []
    for m in range(4):
        a_, b_ = harm(A, m); cs, sn = AFc[m]; out += [a_ - cs, b_ - sn] if m else [a_ - cs]
    return np.array(out)
ls = least_squares(res, [15, 0.16, 0.3, -0.05, 0.9], bounds=([-30, 0, 0.02, -0.3, -10], [60, 1, 0.95, 0.3, 10])).x
print(f"exact LSQ : l0 = {ls[0]:.1f}, beta = {ls[1]:.3f}, Amp = {ls[3]:.3f}, phi = {ls[4]:.2f}, D = {RV/ls[2]:.0f} Gpc")
""")

md(r"""
## Step 6 — Section 3: input parameters, term-by-term solve vs. Table 1

For each Newtonian fit, the unknowns $M_c/M_{\rm crit}(D)$, $D$, $\ell_0$, Amp and $\varphi$ are solved from $C_m,S_m$ ($m=1$–3) of $A_F$ by least squares. Each solve uses the code's **exact** $h(q)$ from the integrator and keeps that fit's shell rule. The Table 1 columns are the full point fits; for LTB only the Table 1 column is shown.

**Label: the masses in this section were computed with Earth's rate set to $H_E=70$ km s⁻¹ Mpc⁻¹, as in the original sheet** ($\rho_s=\Omega_s\rho_c(70)$, $M_{\rm crit}=\tfrac{4\pi}{3}\rho_c(70)D^3$). They are shown unrescaled; Section 4 (Step 8) derives the masses from the anchored $\langle H\rangle$.
""")
py(r"""
# Step 6a - term-by-term solve (exact h(q)) for the two Newtonian cases  [masses at H_E = 70: LABELLED, original-sheet convention]
RHOC70 = rhoc(70.0)                                            # <-- only for the Section 3 mass columns (original convention)
def af_code(D, mcc, Os, kin, kout):
    a = dict(DR=D / RV, Os=Os, mcc=mcc, rin=1 - kin * RV / D, rout=1 + kout * RV / D, Lc=0.005, mc=0.0, Ob=0.0)
    if not (0.01 < a["rin"] < 1): return None
    af = run_newt(a)["af"]; return af / af.mean() if np.all(np.isfinite(af)) else None
def hres(A):
    out = []
    for m in (1, 2, 3):
        a_, b_ = harm(A, m); cs, sn = AFc[m]; out += [a_ - cs, b_ - sn]
    return np.array(out)
def masses(D, mcc, Os, kin, kout, rho_c):
    Dm = D * GPCm; Mc = mcc * 4 * np.pi / 3 * rho_c * Dm**3; rhos = Os * rho_c; Rin, Rout = D - kin * RV, D + kout * RV
    Msh = lambda Rr: 4 * np.pi / 3 * rhos * ((min(max(Rr, Rin), Rout) * GPCm)**3 - (Rin * GPCm)**3)
    return dict(Mc=Mc, Mtot=Mc + Msh(D), Mshtot=Msh(Rout), Mvis=4 * np.pi / 3 * rhos * (RV * GPCm)**3, R_in=Rin, R_out=Rout)
def tbt(Os, kin, kout, x0):
    def res6(x):
        lm, lD, l0_, am, ph = x; A = af_code(10**lD, 10**lm, Os, kin, kout)
        return np.ones(6) if A is None else hres(A_from(shifted(A, l0_), TH1 - l0_, am, ph))
    x = least_squares(res6, x0, diff_step=1e-3, x_scale=[0.3, 0.1, 3, 0.03, 0.2]).x
    lm, lD, l0_, am, ph = x; D = 10**lD; A = af_code(D, 10**lm, Os, kin, kout); A = A_from(shifted(A, l0_), TH1 - l0_, am, ph)
    if am > 0: am, ph = -am, ph + np.pi
    return dict(D=D, mcc=10**lm, Os=Os, l0=l0_ % 360, amp=am, phi=ph % (2 * np.pi), R2=r2(MA, A[IMIG]), R2F=r2(AFfour(TH1), A), **masses(D, 10**lm, Os, kin, kout, RHOC70))
TBT = {"free": tbt(T1["free"]["Os"], 1.0, 1.0, [np.log10(T1["free"]["mcc"]), np.log10(T1["free"]["D"]), 13, -0.12, 1.26]),
       "omb":  tbt(0.049, omb["kin"], omb["kout"], [np.log10(omb["mcc"]), np.log10(omb["D"]), 13, -0.06, 0.93])}
for kk, v in TBT.items(): print(kk, {q: float("%.4g" % v[q]) for q in ("D", "mcc", "l0", "amp", "phi", "R2", "R2F", "Mc", "Mtot")})
""")
py(r"""
# Step 6b - Section 3 table (Table-1 columns from the fit files; masses at H_E = 70 as originally computed - NOT rescaled)
def col_t1(t):
    return dict(D=t["D"], R_in=t["Rin"], R_out=t["Rout"], Mc=t["Mc"], mcc=t["mcc"], Mtot=t["Mtot"], Mshtot=t["Mshtot"], Mvis=t["Mvis"],
                Os=t["Os"], l0=t["l0"], amp=t["amp"], phi=t["phi"], R2=t["R2"], R2F=t["R2F"])
COLS = [("Newt Ω_s free, term-by-term", TBT["free"]), ("Newt Ω_s free, Table 1", col_t1(T1["free"])),
        ("Newt Ω_s=0.049, term-by-term", TBT["omb"]), ("Newt Ω_s=0.049, Table 1", col_t1(T1["newt_omb"])), ("LTB, Table 1", col_t1(T1["rel"]))]
rows = [("D/R", lambda v: f"{v['D']/RV:.3f}"), ("D [Gpc]", lambda v: f"{v['D']:.1f}"), ("R_in [Gpc]", lambda v: f"{v['R_in']:.1f}"),
        ("R_out [Gpc]", lambda v: f"{v['R_out']:.1f}"), ("Mc [kg] (H_E=70)", lambda v: f"{v['Mc']:.2e}"), ("Mc [M_sun] (H_E=70)", lambda v: f"{v['Mc']/MSUN:.2e}"),
        ("Mc/Mcrit(D)", lambda v: f"{v['mcc']:.3g}"), ("M(<D) total [kg] (H_E=70)", lambda v: f"{v['Mtot']:.2e}"),
        ("M_shell total [kg] (H_E=70)", lambda v: f"{v['Mshtot']:.2e}"), ("M within R [kg] (H_E=70)", lambda v: f"{v['Mvis']:.2e}"),
        ("Ω_s", lambda v: f"{v['Os']:.3f}"), ("ℓ0 [deg]", lambda v: f"{v['l0']:.1f}"), ("Amp", lambda v: f"{v['amp']:.4f}"),
        ("φ [rad]", lambda v: f"{v['phi']:.3f}"), ("R² vs Migkas", lambda v: f"{v['R2']:.3f}"), ("R² vs A_F", lambda v: f"{v['R2F']:.3f}")]
from IPython.display import Markdown, display
tab = "| Quantity | " + " | ".join(n for n, _ in COLS) + " |\n|" + "---|" * (len(COLS) + 1) + "\n"
tab += "\n".join("| " + lab + " | " + " | ".join(f(v) for _, v in COLS) + " |" for lab, f in rows)
display(Markdown(tab))
print("\nThe relativistic term-by-term column (D not fixed; R2 0.824/0.899) is left out, as in the v3 sheet; it can be restored.")
""")

md(r"""
## Step 7 — Anchoring $\langle H\rangle$ on Migkas's own $H_0$ values

The fits fix only the **shape** $A(\ell)/\langle A\rangle$. The absolute sky-mean rate comes from Migkas et al.'s own $H_0$ numbers by fitting
$$H(\ell)=\langle H\rangle\,\frac{A(\ell)}{\langle A\rangle}$$
by weighted least squares, $\langle H\rangle=\sum_i w_iH_im_i/\sum_iw_im_i^2$ with $w_i=1/\sigma_i^2$. Here $m_i$ is the model at the region's cone centre, or the model mean outside the cone for "rest of sky" (75° cones for LX–T and the joint result, 60° for YSZ–T, 90° for "opposite hemisphere"). Only $\ell$ is used because the model is planar. The values below were checked against the arXiv full texts (2004.03305 §7, §8.1–8.2; 2103.13904 §6.1–6.4, App. B). The **representative value** is the 2021 all-combined joint result.

*Caveat:* Migkas's absolute $H_0$ values rest on their fiducial calibration ($H_0=70$ assumed when calibrating the scaling relations; M21 §3), so the anchor inherits that calibration. The H0 maps are published only as figures, not tables, so they are not used.
""")
py(r"""
# Step 7 - all 11 Migkas pairs -> <H> for each fit
CH = [  # (label, [(H0, sigma, kind, l, cone)])   kind: 'pt' = cone centre, 'out' = outside cone / opposite hemisphere
 ("M20 our sample LX-T, max/min",           [(75.17, 1.81, "pt", 17, 0), (66.20, 1.72, "pt", 274, 0)]),
 ("M20 ACC LX-T, max/min",                  [(78.76, 4.15, "pt", 77, 0), (58.12, 2.68, "pt", 317, 0)]),
 ("M20 XCS-DR1 LX-T, max/min",              [(77.91, 2.20, "pt", 31, 0), (63.56, 2.32, "pt", 281, 0)]),
 ("M20 joint (3 samples), max/min",         [(76.64, 1.41, "pt", 34, 0), (65.20, 1.48, "pt", 303, 0)]),
 ("M21 our LX-T, region/rest",              [(66.0, 1.7, "pt", 274, 0), (72.1, 1.4, "out", 274, 75)]),
 ("M21 ours+ACC LX-T, region/rest",         [(66.2, 1.6, "pt", 284, 0), (72.7, 1.5, "out", 284, 75)]),
 ("M21 our YSZ-T, region/rest",             [(62.3, 2.2, "pt", 268, 0), (72.1, 0.9, "out", 268, 60)]),
 ("M21 LBCG-T, region/opp. hemisphere",     [(66.2, 2.5, "pt", 257, 0), (72.6, 2.3, "out", 257, 90)]),
 ("M21 ACC LX-T, region/rest",              [(61.3, 4.2, "pt", 318, 0), (72.2, 2.4, "out", 318, 75)]),
 ("M21 ACC YSZ-T, region/opp. hemisphere",  [(60.6, 3.6, "pt", 311, 0), (73.4, 1.9, "out", 311, 90)]),
 ("M21 ALL COMBINED (joint), region/rest",  [(66.5, 1.0, "pt", 273, 0), (72.8, 0.6, "out", 273, 75)]),
]
SHAPE = {kk: CURVE[kk] / CURVE[kk].mean() for kk in ("free", "omb", "rel")}
def m_of(s, kind, l, th):
    if kind == "pt": return s[int(round(l)) % 360]
    dl = np.abs((TH1 - l + 180) % 360 - 180); return s[dl > th].mean()
HM = {kk: [] for kk in SHAPE}
print("choice".ljust(42) + "Newt free    Newt 0.049   LTB")
for lab, pts in CH:
    line = lab.ljust(42)
    for kk, s in SHAPE.items():
        m = np.array([m_of(s, kd, l, th) for _, _, kd, l, th in pts]); H = np.array([p[0] for p in pts]); wt = 1 / np.array([p[1] for p in pts])**2
        Hm = (wt * H * m).sum() / (wt * m * m).sum(); err = 1 / np.sqrt((wt * m * m).sum()); HM[kk].append((Hm, err)); line += f"{Hm:6.2f}±{err:.2f}  "
    print(line)
ANCH = {kk: dict(H=v[-1][0], err=v[-1][1], mean=np.mean([x[0] for x in v]), sd=np.std([x[0] for x in v], ddof=1)) for kk, v in HM.items()}
print("\nrepresentative (M21 joint):", {kk: f"{a['H']:.1f}±{a['err']:.1f}" for kk, a in ANCH.items()})
print("mean ± s.d. of all 11     :", {kk: f"{a['mean']:.1f}±{a['sd']:.1f}" for kk, a in ANCH.items()})
""")

md(r"""
## Step 8 — Section 4: everything that follows from $\langle H\rangle$

**Newtonian.** $\langle H\rangle$ is the sky mean of the full curve $a_f[1-\mathrm{Amp}\sin(\psi+\varphi)]$ in code units ($1/t_0=k$), so $\langle H\rangle/k$ is a pure number and
$$k=\frac{\langle H\rangle}{\langle H\rangle/k},\quad t_0=k^{-1},\quad H_E=(H_E/k)\,k,\quad \rho_s=\Omega_s\rho_c(H_E),\quad M_{\rm crit}=\tfrac{4\pi}{3}\rho_c(H_E)D^3,$$
$$V_O=v_E\,k\,D,\qquad \text{bulk flow}=\Big|\tfrac{\epsilon_b^2}{10}(h''+4h')\Big|V_O,\;\epsilon_b=r_b/D,\; r_b=150\,h^{-1}\,\text{Mpc with } h=\langle H\rangle/100 .$$
**LTB.** $\langle H\rangle$ is the in-plane mean of the local rate, $(H_\perp+H_r)/2$. From the LTB run, $\langle H\rangle t_0=0.956$ (dimensionless). Masses are rescaled from the fit's unit convention (3D-mean local rate 70) by $(\langle H\rangle/72.04)^2$. This rescaling is approximate in GR, and the bulk flow is not computed.

Accepted values: Planck 2018 VI ($H_0=67.36\pm0.54$, $\Omega_bh^2=0.02237$, $\Omega_m=0.3153$, age 13.797 Gyr); SH0ES ($73.04\pm1.04$); Planck 2018 I (CMB dipole $369.82\pm0.11$ km/s; Local Group $620\pm15$ km/s); Watkins 2023 (bulk flow $395\pm29$ km/s at $150\,h^{-1}$ Mpc).
""")
py(r"""
# Step 8a - accepted values (derived masses use Planck's own H0 = 67.36, not 70)
Pl = dict(H0=67.36, Obh2=0.02237, Om=0.3153, age=13.797); VOL = 4 * np.pi / 3 * (RV * GPCm)**3
Ob = Pl["Obh2"] / (Pl["H0"] / 100)**2; rcP = rhoc(Pl["H0"])
ACC = dict(Ob=Ob, rho_b=Ob * rcP, rhoc=rcP, Mb=Ob * rcP * VOL, Mm=Pl["Om"] * rcP * VOL, vcmb=369.82, vlg=620.0, bulk=395.0)
print({k: float("%.4g" % v) for k, v in ACC.items()})
""")
py(r"""
# Step 8b - model columns from the anchored <H> (no H0 = 70 input for the Newtonian columns)
def radial_h(a, r, qs):                                     # h(q) = V(q)/q along the CoE-Earth axis, code units
    q = np.vstack([[1.0, 0.0], np.c_[qs, np.zeros_like(qs)]])
    X, V, VR = P.integrate(q, r["mu"], P.RG, r["MT"], 0.0, r["Rb"], 0.01, 200); return V[1:, 0] / qs
SEC4 = {}
for kk, a in NEWT.items():
    r = RUN[kk]
    Hk = A_from(shifted(r["af"], a["l0"]), TH1 - a["l0"], a["amp"], a["phi"]).mean()   # <H>/k: sky mean of af x sine factor, code units (1/t0 = k)
    H = ANCH[kk]["H"]; k = H / Hk; HE = r["Ht"] * k; D = a["DR"] * RV
    kin, kout = (1 - a["rin"]) * a["DR"], (a["rout"] - 1) * a["DR"]           # shell edges in units of R: R_in = D - kin R, R_out = D + kout R
    ms = masses(D, a["mcc"], a["Os"], kin, kout, rhoc(HE))                    # rho_c at Earth's rate H_E (the code's Omega reference)
    qs = np.linspace(0.95, 1.05, 41); h = radial_h(a, r, qs); cf = np.polyfit(qs - 1, h, 6); h1 = np.polyval(cf, 0)
    hp = np.polyval(np.polyder(cf), 0) / h1; hpp = np.polyval(np.polyder(cf, 2), 0) / h1
    VO = r["vE"] * k * D * 1000; eb = 0.150 / (H / 100) / D
    SEC4[kk] = dict(H=H, Hk=Hk, k=k, t0=977.8 / k, rho_s=a["Os"] * rhoc(HE), Mvis=ms["Mvis"], Mc=ms["Mc"], Mtot=ms["Mtot"], Os=a["Os"], Os_H=a["Os"] * (HE / H)**2,
                    rhoc=rhoc(H), VO=VO, bulk=abs(eb**2 / 10 * (hpp + 4 * hp)) * VO)
# LTB: dimensionless <H> t0 from the run; masses from the fit's unit convention (3D local mean = 70) rescaled by (<H>/<H>_fit)^2
Hloc_fit = (ltb["H_perp"] + ltb["H_r"]) / 2; Ht0 = Hloc_fit * ltb["age_Gyr"] / 977.8; H = ANCH["rel"]["H"]; g = H / Hloc_fit; t3 = T1["rel"]
SEC4["rel"] = dict(H=H, Hk=Ht0, k=H / Ht0, t0=977.8 * Ht0 / H, rho_s=0.049 * rhoc(70.0) * g**2, Mvis=t3["Mvis"] * g**2, Mc=t3["Mc"] * g**2, Mtot=t3["Mtot"] * g**2,
                   Os=0.049, Os_H=0.049 * (70.0 / Hloc_fit)**2,   # 70 = LTB fit's unit convention (3D local mean), labelled
                   rhoc=rhoc(H), VO=t3["vE"] * c * g, bulk=None)
print(f"LTB: <H>_fit = {Hloc_fit:.2f} (fit units), <H> t0 = {Ht0:.3f}, rescale g = {g:.4f}")
""")
py(r"""
# Step 8c - Section 4 table
F4 = [("⟨H⟩ [km/s/Mpc]", "67.36±0.54 (CMB); 73.04±1.04 (local)", lambda v: f"{v['H']:.1f}  (= {v['Hk']:.3f} k)" if v is not SEC4['rel'] else f"{v['H']:.1f}  (= {v['Hk']:.3f}/t0)"),
      ("k (LTB: 1/t0) [km/s/Mpc]", "-", lambda v: f"{v['k']:.1f}"),
      ("Age / k^-1 = t0 [Gyr]", "13.797±0.023", lambda v: f"{v['t0']:.2f}"),
      ("Ordinary mass within R [kg]", f"{ACC['Mb']:.2e}", lambda v: f"{v['Mvis']:.2e}"),
      ("Total matter within R [kg]", f"{ACC['Mm']:.2e}", lambda v: "same (no DM)"),
      ("rho_b [kg/m^3]", f"{ACC['rho_b']:.2e}", lambda v: f"{v['rho_s']:.2e}"),
      ("Omega_b (= Omega_m in model)", f"{ACC['Ob']:.4f}; Om 0.3153", lambda v: f"{v['Os']:.3f} ({v['Os_H']:.3f} vs ⟨H⟩)"),
      ("rho_c at ⟨H⟩ [kg/m^3]", f"{ACC['rhoc']:.2e} (67.36)", lambda v: f"{v['rhoc']:.2e}"),
      ("Earth speed vs CMB V_O", "369.82±0.11 km/s", lambda v: f"{v['VO']:.2e} km/s = {v['VO']/c:.1f}c"),
      ("V_O / CMB dipole", "1", lambda v: f"{v['VO']/369.82:,.0f}"),
      ("Local Group vs CMB", "620±15 km/s", lambda v: "not modelled"),
      ("Bulk flow, 150/h Mpc [km/s]", "395±29", lambda v: "not computed" if v["bulk"] is None else f"{v['bulk']:.2g}"),
      ("Core mass Mc [kg] (x M_b)", "-", lambda v: f"{v['Mc']:.2e} ({v['Mc']/ACC['Mb']:.1f}x)"),
      ("Mass inside D [kg] (x M_b)", "-", lambda v: f"{v['Mtot']:.2e} ({v['Mtot']/ACC['Mb']:,.0f}x)")]
from IPython.display import Markdown, display
tab = "| Quantity | Accepted | Newt Ω_s free (0.909) | Newt Ω_s=0.049 (0.893) | LTB (0.885) |\n|---|---|---|---|---|\n"
tab += "\n".join("| " + lab + " | " + accv + " | " + " | ".join(f(SEC4[kk]) for kk in ("free", "omb", "rel")) + " |" for lab, accv, f in F4)
display(Markdown(tab))
json.dump(dict(anchor=ANCH, sec4=SEC4, tbt=TBT, R2={kk: r2(MA, CURVE[kk][IMIG]) for kk in CURVE}), open("/workspace/cosmology/papers/A_equals_fourier/steps_v3/notebook_results.json", "w"), indent=1, default=float)
""")
md(r"""
## Summary

* Headline fits (full curve with sine, vs 71 Migkas points): $R^2$ = 0.909 (Newtonian, $\Omega_s$ free), 0.893 (Newtonian, $\Omega_s=0.049$), 0.885 (LTB).
* Migkas-anchored sky-mean rate (2021 joint result): $\langle H\rangle\approx71.8/71.6/71.6$ km s⁻¹ Mpc⁻¹; the mean over all 11 Migkas pairs is $\approx70.4\pm2$.
* With that anchor, $t_0=k^{-1}\approx13.2/13.8/13.05$ Gyr. Earth's speed about the CoE is $V_O\approx21c/11c/7.6c$, thousands of times the CMB dipole speed. This is the model's main tension.
* The absolute anchor inherits Migkas's fiducial calibration.
""")
import build_figs; build_figs.add(md, py)
nb = nbf.v4.new_notebook(); nb["cells"] = C
nb["metadata"]["kernelspec"] = dict(name="coe", display_name="coe", language="python")
nbf.write(nb, "coe_steps_v3.ipynb"); print(len(C), "cells")
