"""Term-by-term (harmonic) solve for the input parameters, per Table-1 column, using the code's exact geometry (exact h(q) -> a_m).
Unknowns: Mc/Mcrit(D), D (via eps = R/D), l0, Amp, phi. Equations: C_m, S_m of A_F for m = 1..3 (m=0 is the normalisation)."""
from tm_core import *
from scipy.optimize import least_squares
import prof_core as P
H0 = 70 / 3.0857e19; G = 6.674e-11; RHOC = 3 * H0**2 / (8 * np.pi * G); GM = 3.0857e25; MSUN = 1.989e30
def af_code(D, mcc, Os, kin, kout, Rs):
    a = dict(DR=D / Rs, Os=Os, mcc=mcc, rin=1 - kin * RV / D, rout=1 + kout * RV / D, Lc=0.005, mc=0.0, Ob=0.0)
    if not (0.01 < a["rin"] < 1): return None
    _, _, m1 = S.build(a); rr = P.run("ptshell", a, m1, 0.0); af = rr["af"]
    if not np.all(np.isfinite(af)) or rr["vrmin"] <= 0: return None
    return af / af.mean()
def curveA(af, l0, amp, phi): return A_from(shifted(af, l0), TH1 - l0, amp, phi)
def hres(A):
    out = []
    for m in (1, 2, 3):
        a, b = harm(A, m); c, s = AFc[m]; out += [a - c, b - s]
    return np.array(out)
def solve(Os, kin, kout, Rs, x0, D_fixed=None):
    def res(x):
        lm, lD, l0, amp, phi = x; D = D_fixed if D_fixed else 10**lD
        af = af_code(D, 10**lm, Os, kin, kout, Rs)
        return np.ones(6) if af is None else hres(curveA(af, l0, amp, phi))
    r = least_squares(res, x0, diff_step=1e-3, x_scale=[0.3, 0.1, 3, 0.03, 0.2])
    lm, lD, l0, amp, phi = r.x; D = D_fixed if D_fixed else 10**lD; mcc = 10**lm
    af = af_code(D, mcc, Os, kin, kout, Rs); A = curveA(af, l0, amp, phi)
    if amp > 0: amp, phi = -amp, phi + np.pi
    Dm = D * GM; Mcrit = 4 * np.pi / 3 * RHOC * Dm**3; Mc = mcc * Mcrit; rhos = Os * RHOC; Rin, Rout = D - kin * RV, D + kout * RV
    Msh = lambda R: 4 * np.pi / 3 * rhos * ((min(max(R, Rin), Rout) * GM)**3 - (Rin * GM)**3)
    return dict(D=D, DR=D / RV, Dkm=D * 3.0857e22, R_in=Rin, R_out=Rout, Mc=Mc, Mc_sun=Mc / MSUN, mcc=mcc, Mtot=Mc + Msh(D), Mshtot=Msh(Rout), Mvis=4 * np.pi / 3 * rhos * (RV * GM)**3,
                Os=Os, l0=l0 % 360, amp=amp, phi=phi % (2 * np.pi), R2mig=r2(MA, A[IMIG]), R2F=r2(AFfour(TH1), A), res=list(r.fun), cost=float(np.sum(r.fun**2)), eps=Rs / D)
T = json.load(open("universe_thin.json"))
out = {}
fr = T["free"]; om = json.load(open("universe_newt_omb.json"))["Rs14.3|sine"]; rl = json.load(open("universe_rel_mcfree.json"))
out["free"] = solve(fr["Os"], 1.0, 1.0, RV, [np.log10(fr["mcc"]), np.log10(fr["D"]), 13, -0.12, 1.26])
out["omb"] = solve(0.049, om["kin"], om["kout"], RV, [np.log10(om["mcc"]), np.log10(om["D"]), 13, -0.06, 0.93])
out["rel"] = {}
for D in (rl["D"], 100.0, 330.0):
    out["rel"][D] = solve(0.049, rl["kin"], rl["kout"], 1.19, [np.log10(0.002), 0, 13.5, -0.055, 0.85], D_fixed=D)
for k, v in out.items():
    if k == "rel":
        for D, w in v.items(): print("rel D=%.1f" % D, {kk: (float("%.4g" % vv) if isinstance(vv, float) else vv) for kk, vv in w.items() if kk != "res"})
    else: print(k, {kk: (float("%.4g" % vv) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "res"})
json.dump(out, open("/workspace/cosmology/papers/A_equals_fourier/thinshell/solve_inputs.json", "w"), indent=1, default=float)
