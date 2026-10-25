"""Best-fit LTB models built from Chris's variables."""
import numpy as np, ltb
from scipy.optimize import brentq, fsolve
H0 = ltb.H0_70
DR, LD, MC = 2.904374576719836, 0.07027609488588729, 0.023946715374238915
Rvis = 14.3                       # Gpc, vision-shell radius
D = DR * Rvis; L = LD * D
MCRIT = H0**2 * D**3 / 2          # geometric mass (Gpc) of critical density inside D
M = MC * MCRIT                    # 46.8 Gpc = 1.944e54 kg
KG_PER_GPC = 3.0857e25 * 2.99792458e8**2 / 6.674e-11

def Hloc(p, r, t):
    Hp, Hr = ltb.local_rates(r, t, p); return (Hr + 2 * Hp) / 3

def model_A(M=M, D=D, L=L, tb=(0.0, 1.0)):
    """Profile on the areal radius today, uniform (or given) bang time, E(r) from the age; observer at r = D."""
    p = np.zeros(10); p[0] = 1; p[1] = M / ltb.Fexp(D / L)[0]; p[2] = L; p[3] = 1; p[6], p[7] = tb
    def f(t0):
        p[5] = t0; return Hloc(p, D, t0) / H0 - 1
    t0 = brentq(f, 0.5 / H0, 1.5 / H0, xtol=1e-14); p[5] = t0
    return p, t0, D

def model_B(M=M, D=D, L=L, tb=(0.0, 1.0), efac=1.0):
    """Profile on the Lagrangian label r, uniform bang, E = efac r^2/(2 t0^2) (Milne energy of the Newtonian
    initial condition v = x/t); r_obs and t0 fixed by R(r_obs,t0) = D and local H = H0."""
    p = np.zeros(10); p[0] = 1; p[2] = L; p[3] = 0; p[6], p[7] = tb
    def set_(ro, t0):
        p[1] = M / ltb.Fexp(ro / L)[0]; p[4] = efac / (2 * t0**2); p[5] = t0
    def eqs(v):
        ro, t0 = v; set_(ro, t0)
        return [ltb.geom(ro, t0, p)[0] / D - 1, Hloc(p, ro, t0) / H0 - 1]
    ro, t0 = fsolve(eqs, [D, 1 / H0], xtol=1e-13); set_(ro, t0)
    return p, t0, ro

def describe(p, t0, ro, name):
    R, Rr, Rrr, Rt, Rtr, E, E1, Mo, M1 = ltb.geom(ro, t0, p)
    rhoc = 3 * H0**2 / (8 * np.pi)
    Hp, Hr = ltb.local_rates(ro, t0, p)
    d = dict(name=name, t0_Gyr=t0 * ltb.GYR, r_obs=ro, R_obs=R, M_obs_kg=Mo * KG_PER_GPC, M_over_Mcrit=Mo / (H0**2 * R**3 / 2),
             rho_loc_over_crit=ltb.density(ro, t0, p) / rhoc, Hperp=Hp / H0, Hr=Hr / H0,
             shear=(Hr - Hp) / (Hr + 2 * Hp) * 3, twoM_over_R=2 * Mo / R, E_obs=E, Rdot=Rt)
    # centre and shell crossing scan
    rr = np.geomspace(1e-3 * p[2], 60 * p[2] if p[3] == 1 else 3 * ro, 400)
    tt = t0 * np.geomspace(1e-6, 1, 300)
    minRr = np.inf; coll = []; crunch = False
    for r in rr:
        for t in tt:
            try:
                g = ltb.geom(r, t, p)
            except Exception:
                crunch = True; continue
            if not np.isfinite(g[1]): crunch = True; continue
            minRr = min(minRr, g[1])
        g = ltb.geom(r, t0, p)
        coll.append(g[3] < 0)
    coll = np.array(coll)
    d["collapsing_today_r_max"] = float(rr[coll].max()) if coll.any() else 0.0
    d["min_Rprime"] = float(minRr)
    d["rho_c_centre_today"] = ltb.density(1e-3 * p[2], t0, p) / rhoc
    return d

def model_C(dM, Om=0.3, D=D, L=L, tb=(0.0, 1.0)):
    """Open FRW background (Omega_bg, areal labels) plus an e^{-r/L} lump of mass dM centred at the CoE,
    uniform bang, E = (1-Om) H0^2 r^2/2 (pure growing mode); observer at label r = D, t0 from local H = H0."""
    p = np.zeros(10); p[0] = 2; p[1] = dM; p[2] = L; p[3] = 0; p[6], p[7] = tb
    p[8] = Om * H0**2 / 2; p[4] = (1 - Om) * H0**2 / 2
    def f(t0):
        p[5] = t0; return Hloc(p, D, t0) / H0 - 1
    t0 = brentq(f, 0.3 / H0, 1.5 / H0, xtol=1e-14); p[5] = t0
    return p, t0, D
