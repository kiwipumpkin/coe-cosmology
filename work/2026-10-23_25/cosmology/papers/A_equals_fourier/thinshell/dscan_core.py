"""D scan at fixed R = 14.3 Gpc: refit Mc (+ sine, + shell/Os as each fit frees them) at each D/R. Conventions as Table 1 (H = 70 for masses/V_O/t0)."""
from tm_core import *
import prof_core as PC
from bubble_core import best_sine
from scipy.optimize import minimize
G = 6.674e-11; c = 2.99792458e5; GPCm = 3.0857e25; MPC_PER_GPC = 1000.0
RHOC70 = 3 * (70 / 3.0857e19)**2 / (8 * np.pi * G)
def newt_eval(D, mcc, Os, kin, kout, l0, full=False):
    a = dict(DR=D / RV, Os=Os, mcc=mcc, rin=1 - kin * RV / D, rout=1 + kout * RV / D, Lc=0.005, mc=0.0, Ob=0.0)
    if not (0.005 < a["rin"] < 1): return (-5.0, None) if full else -5.0
    try:
        _, _, m1 = S.build(a); rr = PC.run("ptshell", a, m1, 0.0)
    except Exception: return (-5.0, None) if full else -5.0
    af = rr["af"]
    if not np.all(np.isfinite(af)) or rr["vrmin"] <= 0: return (-5.0, None) if full else -5.0
    af = shifted(af / af.mean(), l0); th = TH1 - l0
    amp, ph = best_sine(af[IMIG], th[IMIG], MA); v = r2(MA, A_from(af[IMIG], th[IMIG], amp, ph))
    if not full: return v
    geo = r2(MA, af[IMIG])
    return v, dict(amp=amp, phi=ph, R2geo=geo, Ht=rr["Ht"], rE=rr["rE"], vE=rr["vE"])
def masses(D, mcc, Os, kin, kout):
    Dm = D * GPCm; Mc = mcc * 4 * np.pi / 3 * RHOC70 * Dm**3; rs = Os * RHOC70
    Rin, Rout = D - kin * RV, D + kout * RV
    Msh = 4 * np.pi / 3 * rs * ((Rout * GPCm)**3 - (Rin * GPCm)**3)
    return dict(Mc=Mc, Mshell=Msh, share=Mc / (Mc + Msh), R_in=Rin, R_out=Rout)
