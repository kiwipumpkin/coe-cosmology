import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
import prof_core as P
import ptshell as S
from bubble_core import AFfour, ML, MA, MLO, MHI, TH1, r2, A_from, shifted
from scipy.integrate import cumulative_trapezoid
B = json.load(open("bubble_results.json")); PR = json.load(open("prof_results.json")); PS = json.load(open("ptshell_results.json"))
def psi(kind, d):
    """central potential depth of the CoE mass alone, in units G*M_crit/D: mc * int_0^inf m(r)/r^2 dr"""
    if kind == "ptshell":
        a = dict(d); a["Os"] = 0.0; a["rin"] = 0.5; a["rout"] = 0.6
        core = np.exp(-P.RG / d["Lc"]); m = cumulative_trapezoid(4*np.pi*P.RG**2*core, P.RG, initial=0); m /= np.interp(1, P.RG, m); mc = d["mcc"]
    else:
        _, m = P.tables(kind, d, 100.0); mc = d["mc"]
    rg = P.RG[1:]; return float(mc * np.trapezoid(m[1:] / rg**2, rg))
def curve(kind, d, norm="mean"):
    if kind == "ptshell":
        _, _, m1 = S.build(d); rr = P.run("ptshell", d, m1, 0.0)
    else: rr = P.run(kind, d, d["mc"], d["Ob"])
    af = rr["af"]; af = af / af.mean() if norm == "mean" else af**2 / np.mean(af**2)
    l0 = d.get("l0", 0.0); afs = shifted(af, l0) if l0 else af
    return A_from(afs, TH1 - l0, d["amp"], d["phi"])
