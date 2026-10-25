"""Arithmetic / derivation checks for dq_equality.tex."""
import json, numpy as np
J = json.load(open('/workspace/cosmology/figures/dq_model_equality.json'))
c_kms, z, H0 = 2.99792458e5, 0.075, 70.0
MLY = 3.2616  # Mly per Mpc
d = c_kms*z/H0
print(f"d = cz/H0 = {d:.1f} Mpc = {d*MLY/1e3:.3f} Gly")
for q in (0.8117, 0.790, 0.759, 0.755):
    print(f"  d/D={q}: D={d/q:.0f} Mpc = {d/q*MLY/1e3:.2f} Gly, D/r={d/q/14300:.4f}")
print(f"r=14.3 Gpc = {14.3*MLY:.1f} Gly; 46.5 Gly = {46.5/MLY:.2f} Gpc")
# rho_crit
G = 4.30091e-9  # Mpc (km/s)^2 / Msun
rc = 3*H0**2/(8*np.pi*G); print(f"rho_crit(h=0.7) = {rc:.4e} Msun/Mpc^3 (2.775e11 h^2 = {2.775e11*0.49:.4e})")
for D in (d/0.8117, d/0.759):
    Mc = 4*np.pi/3*rc*D**3
    print(f"  D={D:.0f}: Mcrit(<D)={Mc:.3e}; kappa*Mcrit(k=.45)={.45*Mc:.3e}; 2*kappa*Mcrit(k=.45..0.47)={0.9*Mc:.3e}..{0.94*Mc:.3e}; 0.3*Mcrit={0.3*Mc:.3e}; L={D/7.62:.1f} Mpc={D/7.62*MLY:.0f} Mly; L(D/L=9..6)={D/9:.0f}..{D/6:.0f}")
# kappa definition in code: G M/(D^3 H^2) = M/(2 Mcrit)
print("code kappa = GM/(H^2 D^3); Mcrit = H^2 D^3/(2G) -> M/Mcrit = 2 kappa")
# windows
W = {m: np.sin(m*np.pi/4)/(m*np.pi/4) for m in (1,2,3)}; print("W:", {k: round(v,4) for k,v in W.items()})
# check window numerically
x = np.linspace(-np.pi/4, np.pi/4, 200001)
for m in (1,2,3): print(f"  mean cos(m x') over 90deg: m={m}: {np.trapezoid(np.cos(m*x),x)/(np.pi/2):.4f}")
# ratio solve
p = J['projections_l0_15']['parallel']; m1,m2,m3,m0 = p['m1'],p['m2'],p['m3'],p['m0']
R = m1/m3
c2 = 4.5*W[1]/(45/8*W[1] - R*15/8*W[3]); c=np.sqrt(c2)
dD = (-m3/m2)/(1.25*c*W[3]/W[2]); Nak = m2/(1.5*c2*W[2]); N = m0 - Nak*(1.5*c2-1)
print(f"c2={c2:.4f} b0={np.degrees(np.arccos(c)):.2f} d/D={dD:.4f} Nak={Nak:.4f} N={N:.4f} ak={Nak/N:.4f}")
co = J['data_fit_m0_3']['coef']; l0=np.radians(15)
for m in (1,2,3):
    par = co[f'c{m}']*np.cos(m*l0)+co[f's{m}']*np.sin(m*l0); ort = -co[f'c{m}']*np.sin(m*l0)+co[f's{m}']*np.cos(m*l0)
    print(f"  m={m}: parallel {par:+.4f} orth {ort:+.4f}")
print("dipole phase (m0-3 fit):", J['data_fit_m0_3']['amp_phase']['m1'])
# tidal expansion: exact point-mass (u = -k Rvec/R^3) vs series
def exact(mu, e):  # delta H * D^3/k
    Rv = np.array([1+e*mu, e*np.sqrt(1-mu**2)]); R = np.hypot(*Rv)
    uS = -Rv/R**3; uO = -np.array([1.0,0.0]); n = np.array([mu, np.sqrt(1-mu**2)])
    return (uS-uO)@n/e
for e in (1e-3, 1e-2):
    mus = np.linspace(-1,1,41)
    ex = np.array([exact(m_, e) for m_ in mus]); lead = 3*mus**2-1
    nxt = (ex-lead)/e
    A = np.stack([mus, mus**3],1); coef,*_ = np.linalg.lstsq(A, nxt, rcond=None)
    print(f"  eps={e}: leading resid {np.max(np.abs(ex-lead))/e:.3f}/eps; next-order coeffs mu,mu^3 = {coef}")
e=0.76
mus=np.linspace(-1,1,5)
print("  eps=0.76 exact vs 2-term series:", [(round(exact(m_,e),3), round(3*m_**2-1+e*(4.5*m_-7.5*m_**3),3)) for m_ in mus])
# beyond-horizon implication
print(f"D>r=14300 with d={d:.0f} -> d/D < {d/14300:.4f}; m3 would be ~{m3*(d/14300)/dD:.5f} vs measured {m3:.4f}")
print(f"d/D=0.76-0.81 with D>r needs d > {0.76*14300:.0f} Mpc")

# ---- d < r < D: odd terms predicted by Eqs (5),(7) at d/D = d/r
e = d/14300.0
print(f"\nd<r<D: d/D < d/r = {e:.4f}")
for lab, c2 in (("b0=33.5 (sheet)", 0.6949), ("b0=0", 1.0), ("b0=15", np.cos(np.radians(15))**2)):
    cc = np.sqrt(c2); Nak2 = 0.0733/(1.5*c2*W[2]); N2 = 0.9965 - Nak2*(1.5*c2-1)
    m3p = -Nak2*e*15/8*cc**3*W[3]; m1p = Nak2*e*cc*(4.5-45/8*c2)*W[1]
    print(f"  {lab}: Nak={Nak2:.4f} N={N2:.4f} ak={Nak2/N2:.4f} m3={m3p:+.5f} (meas -0.0292, x{-0.0292/m3p:.0f}) m1={m1p:+.5f} (meas +0.0398)")
