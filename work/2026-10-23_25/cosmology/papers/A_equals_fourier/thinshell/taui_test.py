import sys, json, time, numpy as np
sys.path.insert(0, "/workspace/cosmology/papers/spreadsheet_model/bubble"); import os; os.chdir("/workspace/cosmology/papers/spreadsheet_model/bubble")
import prof_core as P, ptshell as S, newt_omb as NO
from scipy.optimize import minimize, minimize_scalar
base = S.run   # ptshell-patched run
TAU = {"v": 0.01}
def run_t(kind, a, mc, Ob, fac=1.5, **kw):
    kw.setdefault("taui", TAU["v"]); kw.setdefault("nsteps", 200 if TAU["v"] >= 0.01 else 320); return base(kind, a, mc, Ob, fac, **kw)
P.run = run_t
o = json.load(open("universe_newt_omb.json"))["Rs14.3|sine"]
out = {}
t = time.time(); print("obj", NO.objective([np.log10(o["D"]), np.log10(o["mcc"]), np.log(o["kin"] - 1), np.log(o["kout"] - 1)], NO.RV, True), time.time() - t)
for tau in (0.01, 0.001, 0.1):
    TAU["v"] = tau
    # (a) D, kin, kout fixed at the Table-1 fit; refit Mc/Mcrit (+ l0, Amp, phi)
    f = lambda lm: NO.objective([np.log10(o["D"]), lm, np.log(o["kin"] - 1), np.log(o["kout"] - 1)], NO.RV, True)
    r = minimize_scalar(f, bounds=(np.log10(o["mcc"]) - 2.5, np.log10(o["mcc"]) + 2.5), method="bounded", options=dict(xatol=1e-3))
    # (b) full refit, started from the Table-1 optimum with Mc rescaled by 99/(1/tau-1)
    u0 = [np.log10(o["D"]), np.log10(o["mcc"] * 99 / (1 / tau - 1)), np.log(o["kin"] - 1), np.log(o["kout"] - 1)]
    rf = minimize(NO.objective, u0, args=(NO.RV, True), method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-5, maxiter=250))
    D, mcc, kin, kout = NO.unpack(rf.x)
    out[tau] = dict(fixedgeom_mcc=10**r.x, fixedgeom_R2=1 - r.fun, full=dict(D=D, mcc=mcc, kin=kin, kout=kout, R2=1 - rf.fun))
    out[tau]["fixedgeom_mcc_x"] = 10**r.x * (1 / tau - 1); out[tau]["full"]["mcc_x"] = mcc * (1 / tau - 1)
    print(tau, out[tau], "%.0fs" % (time.time() - t), flush=True)
json.dump(out, open("/workspace/cosmology/papers/A_equals_fourier/thinshell/taui_test.json", "w"), indent=1, default=float)
