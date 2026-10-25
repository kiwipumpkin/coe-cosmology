"""Shared exact-LTB helpers for cycle 2 (reuse Chris's LTB code unchanged: gr_ltb/refit_mass/ltb_fill.py, rm_core.py,
papers/cmb_redshift/rays.py + analyze.py). G=c=1, Gpc. LTB convention: single bang t_B=0, E(r) from R(r,t0)=r, t0 from local mean H=70 at Earth
(the code's unit convention; only dimensionless ratios are used)."""
import sys, os, json, numpy as np
RM = "/workspace/cosmology/papers/spreadsheet_model/bubble/gr_ltb/refit_mass"
sys.path.insert(0, RM)
import rm_core as C                    # sets rays.ltb = analyze.ltb = ltb_fill
import ltb_fill as LT, rays, analyze
from numpy.polynomial.legendre import leggauss
S3 = "/workspace/cosmology/papers/A_equals_fourier/steps_v3"
RV = 14.30048287260589
def fitparams():
    d = json.load(open(S3 + "/dscan/dscan.json"))["points"]; rel = json.load(open("/workspace/cosmology/papers/spreadsheet_model/bubble/universe_rel_mcfree.json"))
    pick = lambda f, DR: [p for p in d[f] if abs(p["DR"] - DR) < 1e-6][0]
    fr = pick("free", 6.553635010401956); om = pick("omb", 3.4386644382849716)
    return {"free": dict(Mc=fr["Mc"], D=fr["D"], Os=fr["Os"], Rin=fr["R_in"], Rout=fr["R_out"], mcc=fr["mcc"]),
            "omb": dict(Mc=om["Mc"], D=om["D"], Os=0.049, Rin=om["R_in"], Rout=om["R_out"], mcc=om["mcc"]),
            "rel": dict(Mc=rel["Mc_kg"], D=rel["D"], Os=0.049, Rin=rel["R_in"], Rout=rel["R_out"], mcc=rel["mcc"])}
def make_p(Mc_kg, D, Os, Rin, Rout, cavity=True, Lc=None):
    p = C.params(Mc_kg, D, RV, Os); p[9] = Rin; p[10] = Rout
    if Lc is not None: p[2] = Lc
    p[11] = min(0.5, 0.05 * min(RV, Rin if Rin > 0 else RV))
    if not cavity: p[12] = p[8]                           # fill the interior: uniform dust + core only (no empty cavity)
    return p
def setup(p, D):
    t0 = C.solve_t0(p, D); return t0
