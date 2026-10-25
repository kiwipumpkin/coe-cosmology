"""Term matching for the current (thin/thick ordinary-matter shell + CoE core) model. Uses Chris's code unchanged."""
import sys, json, numpy as np
BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"
sys.path.insert(0, BUB); sys.path.insert(0, BUB + "/gr_ltb/refit_mass"); import os; os.chdir(BUB)
import ptshell as S, prof_core as P
from bubble_core import AFfour, TH1, MA, IMIG, r2, A_from, shifted
TH = np.radians(TH1)
def harm(y, m):            # y on TH1 grid; returns (a_m, b_m): y = a0 + sum a_m cos m th + b_m sin m th
    if m == 0: return y.mean(), 0.0
    return 2 * np.mean(y * np.cos(m * TH)), 2 * np.mean(y * np.sin(m * TH))
def amp_ph(a, b, m): return np.hypot(a, b), (np.degrees(np.arctan2(b, a)) / m) % (360 / m)
RV = 14.30048287260589
# --- data: A_F ---
AFc = {0: (0.9965, 0), 1: (0.0296, 0.0434), 2: (0.0636, 0.0366), 3: (-0.0299, -0.0114)}
def rot(c, s, m, l0):      # c cos ml + s sin ml = C cos m psi + S sin m psi, psi = l - l0
    t = np.radians(m * l0); return c * np.cos(t) + s * np.sin(t), -c * np.sin(t) + s * np.cos(t)
# --- Newtonian fits ---
def newt_params(key):
    if key == "free":
        d = json.load(open("thin_b_free_hug.json")); return dict(DR=d["DR"], Os=d["Os"], mcc=d["mcc"], rin=d["rin"], rout=d["rout"], Lc=0.005, mc=0.0, Ob=0.0), d
    o = json.load(open("universe_newt_omb.json"))["Rs14.3|sine"]; D = o["D"]
    a = dict(DR=D / RV, Os=0.049, mcc=o["mcc"], rin=o["R_in"] / D, rout=o["R_out"] / D, Lc=0.005, mc=0.0, Ob=0.0)
    return a, dict(l0=o["l0"], amp=o["amp"], phi=o["phi"], R2=o["R2"], D_Gpc=D)
def newt_run(a):
    _, _, m1 = S.build(a); rr = P.run("ptshell", a, m1, 0.0); return rr, m1
def radial_h(a, rr, qs):
    """numeric expansion-rate profile h(q) = V(q)/q along the axis (labels q, as the code uses nominal positions)"""
    q = np.c_[qs, np.zeros_like(qs)]; q = np.vstack([[1.0, 0.0], q])
    X, V, VR = P.integrate(q, rr["mu"], P.RG, rr["MT"], 0.0, rr["Rb"], 0.01, 200)
    return V[1:, 0] / qs, V[0, 0]
def A_radial(hfun, eps):
    """exact |dV|/separation for a radial field V = q h(q) (units D, Earth at q=1), sources at distance eps, angle th from anti-CoE"""
    x = 1 + eps * np.cos(TH); y = eps * np.sin(TH); qs = np.hypot(x, y); hs = hfun(qs); h1 = hfun(np.array([1.0]))[0]
    af = np.hypot(hs * x - h1, hs * y) / eps; return af / af.mean()
