"""Radial blast-wave shells about the CoE, present ONLY inside our visible sphere (|x| < 14.3 Gpc).
delta = A*sqrt(2)*cos(2 pi r_CoE / lam + phi)  (rms A).  lam = 120 h^-1 Mpc = 0.1714 Gpc (Einasto 1997, h=0.7) or 0.150 Gpc (BAO-like).
Same ring kernel / sampling as ../checkerboard/cb_core.py (imported read-only)."""
import sys; sys.path.insert(0,"/workspace/clump_scan/checkerboard")
from cb_core import *
def delta_r(P,A,lam,phi=0.0):
    r=np.linalg.norm(P-C_G,axis=1); return A*np.sqrt(2)*np.cos(2*np.pi*r/lam+phi)
def pattern_r(A,lam,phi=0.0,DT=1.0,h=LAM/6,nth=360,pts=None):
    if pts is None: pts=points(DT,h,sphere_only=True)
    P,h=pts; M=delta_r(P,A,lam,phi)*seq_core.rho_gpc(P)*GPC**3*h**3
    return kern(P,M,h,nth)
