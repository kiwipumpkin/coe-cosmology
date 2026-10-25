"""Hand-grouped explicit partials of the four angles (no delta names), checked against sympy."""
import numpy as np, sympy as sp
from nobox_dl_dd_sym import *
P_or = lambda sv: (4*L*D*(S - D**3) - sv*3*d*(4*L*(S - D**3) - w*D**4))      # Qt/3
def forms(nm):
    sv = 1 if nm.endswith("1") else -1; Qb = P_or(sv)
    if nm.startswith("lmax"):
        Pd = D**3*H**2*S*(w*D**4 - 4*L*(S - D**3)) - g*(w*D**7 + 4*L*(S - D**3)*(S - 2*D**3))
        PL = (-4*D*L**2*SL*(D**3*H**2 - g) - 9*w*g*d**2*(D**3 - L*SL - S)
              + sv*3*d*(w*D**4*(H**2*(L*SL + S) - g) + 4*L**2*SL*(D**3*H**2 - g)))
        PD = (4*D**4*L*(D*SD*(D**3*H**2 - g) - 3*S*(H**2*S - g)) - 9*w*g*d**2*(D**6 + D**4*SD - 4*D**3*S)
              - sv*3*d*(w*D**8*H**2*SD + 4*D**4*L*SD*(D**3*H**2 - g) + 4*D**3*H**2*L*S*(D**3 - 4*S)
                        - 4*g*L*(2*D**6 - 6*D**3*S + S**2)))
        return {d: deg*2*a*cph*D*L*Pd/(g*Qb**2),
                L: -sv*deg*sp.Rational(2, 3)*a*cph*D**4*PL/(g*Qb**2),
                D: sv*deg*sp.Rational(2, 3)*a*cph*L*PD/(g*Qb**2),
                M: sv*deg*2*a*cph*L*D**4*H**2*S/(3*g*M*Qb)}
    Pm = 4*D*sph + sv*3*d*cph
    return {d: sv*135/(sp.pi*D)*(1 - sv*a*cph*(D**3*H**2 - g)*S/(6*g*(S - D**3))),
            M: sv*deg*a*Pm*D**2*H**2*S/(24*g*M*(S - D**3)),
            L: sv*deg*a*D**2*SL*Pm*(D**3*H**2 - g)/(24*g*(S - D**3)**2),
            D: sv*deg*(-sp.Rational(3, 4)*d/D**2 - a/(24*g)*((16*D**3*H**2*sph - 4*g*sph + sv*9*d*D**2*H**2*cph)*S/(D*(S - D**3))
                        + Pm*(D**3*H**2 - g)*(4*D**3*S - S**2 - D**4*SD)/(D**2*(S - D**3)**2)))}
def ell_explicit(nm):
    """explicit angle (for the page): l_max via Nt/Qt, l_min via theta."""
    return ell(nm)
if __name__ == "__main__":
    rng = np.random.default_rng(7); worst = {}
    for nm in NAMES4:
        F = forms(nm); e = ell(nm)
        for _ in range(4):
            v = {D: rng.uniform(.8, 1.3), L: rng.uniform(.1, .4), d: rng.uniform(.03, .2), M: rng.uniform(.1, .5),
                 al: rng.uniform(.5, 2), G: 1.2, H: .9, f: .94, a: -.06, ph: .35}
            for X in (D, L, d, M):
                r1 = float(full(pd(e, X)).subs(v).evalf(30)); r2 = float(full(F[X]).subs(v).evalf(30))
                worst[(nm, str(X))] = max(worst.get((nm, str(X)), 0), abs(r1 - r2)/abs(r1))
    print(max(worst.values())); print(worst)
