"""Way-out B: CMB-specific terms, sign and size, on the headline point core + Omega_s=0.049 shell (mcc re-solved per D so the
cluster-scale quadrupole at Rs = 0.31 Gpc equals the 'all' slice's 0.084). Legendre a_l about the CoE axis on the R = 14.3 Gpc sphere."""
import json, numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.optimize import brentq
import lib4 as L
T = L.T; RV = L.RV; C = L.CKM
DS = [14.7, 15.0, 16.0, 18.0, 20.0, 21.0, 22.0, 22.3, 23.0, 25.0, 30.0, 40.0, 49.1746, 70.0, 94.0, 140.0, 200.0]
def mcc_for(D, target=0.084, Rs=0.31):
    f = lambda lm: np.log(T.dq(L.run(D, 10**lm, "point", {}, [Rs])[0][0]["af"])["a2"] / target)
    return 10**brentq(f, -3.5, -2.3, xtol=1e-6)
def point_only_stat(D, mcc, Hconv):
    """static term of the core alone (no shell/cavity): -G Mc / r, Mc/Mcrit = mcc  ->  in code units mu_c m(r)/r with mu_c = mcc*Ht^2/2"""
    th = np.radians(L.TH1); rs = RV / D; rr = np.hypot(1 + rs * np.cos(th), rs * np.sin(th))
    return rr
out = []
for D in DS:
    mcc = mcc_for(D); c = L.cmb(D, mcc, "point", {})
    res, info = L.run(D, mcc, "point", {}, [0.31, RV]); Hconv = L.HM / res[0]["af"].mean(); vunit = Hconv * D * 1e3
    th = np.radians(L.TH1); rs = RV / D; rr = np.hypot(1 + rs * np.cos(th), rs * np.sin(th))
    rin = 1 - T.KIN * RV / D
    sc = L.legendre(L.static_pattern(D, mcc)); st = c["stat"]; kin = c["kin"]
    out.append(dict(D=D, mcc=mcc, Mhole=0.049 * rin**3, stat=st, stat_core=sc, stat_cavity=[a - b for a, b in zip(st, sc)], kin=kin, cl_a2=c["cl_a2"], far_a2=c["far_a2"]))
    print("D %6.1f mcc %.5f Mhole %.5f | stat a2 %+.3e (core %+.3e cavity %+.3e) a3 %+.3e | kin a2 %+.3e a3 %+.3e a4 %+.3e"
          % (D, mcc, 0.049 * rin**3, st[2], sc[2], st[2] - sc[2], st[3], kin[2], kin[3], kin[4]), flush=True)
# potential-history multipliers for the static tidal pattern (Phi_2 ~ chi^2 today)
Om = 0.049
def Dplus(a):
    x = (1 / Om - 1) * a
    return 1 + 3 / x + 3 * np.sqrt(1 + x) / x**1.5 * np.log(np.sqrt(1 + x) - np.sqrt(x))
E = lambda a: np.sqrt(Om / a**3 + (1 - Om) / a**2)
ae = 1 / 1090.0; ag = np.geomspace(ae, 1, 4000)
g = np.array([Dplus(a) / a for a in ag]); g /= g[-1]
chi = np.array([quad(lambda b: 1 / (b * b * E(b)), a, 1)[0] for a in ag]); chi /= chi[0]      # chi/chi_e, chi_e == R
dg = np.gradient(g, ag); isw = 2 * np.trapezoid(dg * chi**2, ag)
mult = dict(static_newtonian=1.0, SW_adiabatic_third=1 / 3, growing_SW_full=float(g[0]), growing_ISW=float(isw), growing_net_linear=float(g[0] + isw),
            growing_net_exactGR_cycle2=34 / 5.6)
print("multipliers", mult)
# ionised-shell optical depth (today's density, path R) and kSZ-type dilution
rhoc = T.rhoc(71.6); ne = 0.049 * rhoc / 1.6726e-27 * 0.88; tau = 6.652e-29 * ne * RV * T.GPCm
print("shell tau (fully ionised, today's density, path R): %.4f" % tau)
# SW + Doppler cancellation search: m * stat_l + kin_l = 0 at l = 2
canc = []
for o in out:
    m_star = -o["kin"][2] / o["stat"][2]
    res = [m_star * s + k for s, k in zip(o["stat"], o["kin"])]
    canc.append(dict(D=o["D"], m_star=m_star, a3=res[3], a4=res[4]))
    print("D %6.1f  multiplier needed to cancel l=2: %+.2f  -> residual a3 %+.2e a4 %+.2e" % (o["D"], m_star, res[3], res[4]))
# tuning of D for the static null (M_in = M_hole): stat a2 per unit Meff and dMeff/dD
Dz = brentq(lambda D: mcc_for(D) - 0.049 * (1 - T.KIN * RV / D)**3, 20.0, 25.0, xtol=1e-5)
a2u = L.legendre(L.static_pattern(Dz, 1.0))[2]; h = 0.01
dM = ((mcc_for(Dz + h) - 0.049 * (1 - T.KIN * RV / (Dz + h))**3) - (mcc_for(Dz - h) - 0.049 * (1 - T.KIN * RV / (Dz - h))**3)) / (2 * h)
tolD = 1e-5 / abs(a2u * dM)
cz = L.cmb(Dz, mcc_for(Dz), "point", {})
print("static null at D = %.4f Gpc (mcc %.5f); a2 per unit Meff %.3f; dMeff/dD %.2e /Gpc -> |a2|<1e-5 needs |dD| < %.4f Gpc; kin there a2 %+.3e a3 %+.3e a4 %+.3e"
      % (Dz, mcc_for(Dz), a2u, dM, tolD, cz["kin"][2], cz["kin"][3], cz["kin"][4]))
# emission-radius needed for each term to pass (stat ~ Re^2, kin ~ Re in the tidal limit), at D = 49.17
o49 = [o for o in out if abs(o["D"] - 49.1746) < 1e-3][0]
Re_kin = RV * 1e-5 / abs(o49["kin"][2]); Re_stat = RV * np.sqrt(1e-5 / abs(o49["stat_core"][2]))
print("emission radius needed at D=49: Doppler %.2e Gpc, core static %.3f Gpc" % (Re_kin, Re_stat))
json.dump(dict(static_null=dict(D=Dz, mcc=mcc_for(Dz), a2_per_Meff=a2u, dMeff_dD=dM, tol_D_Gpc=tolD, kin=cz["kin"]), Re_kin=Re_kin, Re_stat=Re_stat, rows=out, multipliers=mult, tau_shell=tau, cancel=canc), open("partB.json", "w"), default=float)
