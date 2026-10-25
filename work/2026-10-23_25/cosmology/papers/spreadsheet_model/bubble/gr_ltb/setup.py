"""Configurations and helpers. Imports rays/analyze from papers/cmb_redshift and points them at ltb_thin."""
import json, os, sys, numpy as np
from scipy.optimize import brentq, fsolve
HERE = os.path.dirname(os.path.abspath(__file__))
import ltb_thin as LT
sys.path.insert(0, os.path.join(HERE, "../../../cmb_redshift"))
import rays, analyze
rays.ltb = LT; analyze.ltb = LT
H0 = LT.H0_70; GYR = LT.GYR
KG_PER_GPC = 3.0857e25 * 2.99792458e8**2 / 6.674e-11
N = json.load(open(os.path.join(HERE, "../universe_thin.json")))["free"]
D, RV, Rin, Rout = N["D"], N["R"], N["Rin"], N["Rout"]
Mc = N["Mc"] / KG_PER_GPC
RHOC = 3 * H0**2 / (8 * np.pi)
rho_s = N["Os"] * RHOC                         # geometric (Gpc^-2)
M_D = N["Mtot"] / KG_PER_GPC                   # M(<D)
Lc = 0.005 * D; W = 0.5                        # core scale; shell-edge smoothing half-width (Gpc)

def Hloc(p, r, t):
    Hp, Hr = LT.local_rates(r, t, p); return (Hr + 2 * Hp) / 3

def base():
    p = np.zeros(12); p[1] = Mc; p[2] = Lc; p[8] = rho_s; p[9] = Rin; p[10] = Rout; p[11] = W; return p

def model_areal(t0=None, H=H0):
    """GR-A: label = areal radius today, so today's mass distribution is exactly the Newtonian one; uniform bang t_B = 0;
    E(r) fixed by R(r,t0) = r. One free number t0: solved for local H = H (default) or given."""
    p = base(); p[3] = 1
    if t0 is None:
        def f(t):
            p[5] = t; return Hloc(p, D, t) / H - 1
        t0 = brentq(f, 2.0, 6.0, xtol=1e-13)
    p[5] = t0; return p, t0, D

def model_lagr(H=H0):
    """GR-B (earlier convention, case B of ltb_cmb): Newtonian profile shape on the Lagrangian label, scaled by ro/D
    (core Lc = 0.005 ro, edges ro(1 -+ R/D)), M(<ro) = M(<D) and total mass as Newtonian; E = r^2/(2 t0^2)
    (energy of the Newtonian initial condition v = x/t); uniform bang. ro, t0 from R(ro,t0) = D and local H = H."""
    p = base(); p[3] = 0
    def set_(ro, t0):
        s = ro / D; p[2] = Lc * s; p[9] = Rin * s; p[10] = Rout * s; p[11] = W * s; p[8] = rho_s / s**3
        p[4] = 1 / (2 * t0**2); p[5] = t0
    def eqs(v):
        ro, t0 = v; set_(ro, t0); return [LT.geom(ro, t0, p)[0] / D - 1, Hloc(p, ro, t0) / H - 1]
    ro, t0 = fsolve(eqs, [D, 4.0], xtol=1e-13); set_(ro, t0)
    return p, t0, ro
