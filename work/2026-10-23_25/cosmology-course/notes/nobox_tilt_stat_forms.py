"""Displayed (hand-grouped) solutions of d(ell)/d(input)=0, checked against sympy's own
solutions at random points. w = 1 - f_u. S, SD, SL = bracket and its D, L derivatives."""
import numpy as np, sympy as sp
from nobox_tilt_stat_sym import *
w = 1 - f; cph, sph = sp.cos(ph), sp.sin(ph); g = al*G*M
def forms(sv):
    F = {}
    # l_max, d/dd: same for both maxima
    F["max_d_M"] = H**2/(al*G)*S*(w*D**7 + 4*D**6*L - 4*D**3*L*S)/(w*D**7 + 8*D**6*L - 12*D**3*L*S + 4*L*S**2)
    F["max_L_M"] = D**3*H**2/(al*G)*(4*D*L**2*SL - sv*3*d*((w*D + 4*L)*L*SL + w*D*S)) / \
                   (4*D*L**2*SL - sv*3*d*(w*D**4 + 4*L**2*SL) - 9*w*d**2*(D**3 - L*SL - S))
    F["max_D_M"] = D**3*H**2/(al*G)*(4*D*L*(D**4*SD - 3*S**2) - sv*3*d*(w*D**5*SD + 4*D**4*L*SD + 4*D**3*L*S - 16*L*S**2)) / \
                   (4*D**4*L*(D*SD - 3*S) - sv*12*d*L*(2*D**6 + D**4*SD - 6*D**3*S + S**2) + 9*w*d**2*(D**6 + D**4*SD - 4*D**3*S))
    F["min_d_M"] = D**3*H**2/(al*G)*sv*a*cph*S/(sv*a*cph*S + 6*(S - D**3))
    F["min_M_d"] = -sv*4*D/3*sp.tan(ph)
    F["min_D_d"] = -4*a*D**4*sph*(D*(D**3*H**2 - g)*SD - 3*S*(H**2*S - g)) / \
        (3*(sv*a*cph*(D**4*(D**3*H**2 - g)*SD - D**3*(D**3*H**2 - 4*g)*S - (2*D**3*H**2 + g)*S**2) - 6*g*(S - D**3)**2))
    F["dMtheta"] = a*(4*D*sph + sv*3*d*cph)*D**2*S*H**2/(24*g*M*(S - D**3))
    F["dLtheta"] = a*D**2*SL*(4*D*sph + sv*3*d*cph)*(D**3*H**2 - g)/(24*g*(S - D**3)**2)
    F["dMdelta"] = a*cph*H**2/(g*M)/Q.subs(s, sv)
    return F
def check():
    rng = np.random.default_rng(5); worst = {}
    NUMV = lambda: {D: rng.uniform(.8, 1.3), L: rng.uniform(.1, .4), d: rng.uniform(.03, .2), M: rng.uniform(.1, .5),
                    al: rng.uniform(.5, 2), G: 1.3, H: 0.9, f: 0.94, a: -0.06, ph: 0.35}
    for sv in (1, -1):
        Fm = forms(sv)
        ref = {}
        for X, key in ((d, "max_d_M"), (L, "max_L_M"), (D, "max_D_M")):
            n = sp.fraction(sp.together(to_ig(Fmax(X).subs(s, sv))))[0]
            ref[key] = H**2/(al*G*sp.solve(n, ig)[0])
        n = sp.fraction(sp.together(to_ig(Fmin(d).subs(s, sv))))[0]; ref["min_d_M"] = H**2/(al*G*sp.solve(n, ig)[0])
        n = sp.fraction(sp.together(Fmin(M).subs(s, sv)))[0]; ref["min_M_d"] = sp.solve(n, d)[0]
        n = sp.fraction(sp.together(Fmin(D).subs(s, sv)))[0]; ref["min_D_d"] = sp.solve(n, d)[0]
        ref["dMtheta"] = Fmin(M).subs(s, sv); ref["dLtheta"] = Fmin(L).subs(s, sv)
        ref["dMdelta"] = pd(delta, M).subs(s, sv)
        for _ in range(4):
            v = NUMV()
            for k in ref:
                r1 = float(full(ref[k]).subs(v).evalf(30)); r2 = float(full(Fm[k]).subs(v).evalf(30))
                worst[k] = max(worst.get(k, 0), abs(r1 - r2)/abs(r1))
    return worst
if __name__ == "__main__":
    print(check())
    # physicality of the l_max d-condition: alpha G M/(D^3 H^2) - 1, factored
    kap = forms(1)["max_d_M"]*al*G/(D**3*H**2)
    print("kappa*-1 =", sp.factor(sp.together(kap - 1)))
