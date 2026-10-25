"""Explicit (delta-free) forms of the four tilt-card angles and their partials (sympy).
S, SD, SL: code placeholders for the bracket and its D, L derivatives (written out in full on pages)."""
import sympy as sp
from nobox_tilt_deriv_sym import D, L, d, M, al, G, H, f, a, ph, s, S, SD, SL, pd, full, quantities, collect_terms, Sfun
w = 1 - f; g = al*G*M; cph, sph = sp.cos(ph), sp.sin(ph); deg = 180/sp.pi
# l_max = -deg*Amp cos(phi) Nt/(g Qt)  (index 1), 180 + deg*Amp cos(phi) Nt/(g Qt) (index 2, s=-1)
Nt = 2*L*D**4*H**2*S + g*(4*L*D*S - 6*L*D**4) - s*3*d*g*(2*L*(S - D**3) - w*D**4)
Qt = 12*L*D*(S - D**3) - s*9*d*(4*L*(S - D**3) - w*D**4)
delta_x = -a*cph*Nt/(g*Qt)
theta_x = sp.Rational(3, 4)*d/D - a*(4*D*sph + s*3*d*cph)*(D**3*H**2 - g)*S/(24*D*g*(S - D**3))
def ell(name):
    sv = 1 if name.endswith("1") else -1
    if name == "lmax1": return deg*delta_x.subs(s, 1)
    if name == "lmax2": return 180 - deg*delta_x.subs(s, -1)
    if name == "lmin1": return 90 + deg*theta_x.subs(s, 1)
    if name == "lmin2": return 270 - deg*theta_x.subs(s, -1)
NAMES4 = ["lmax1", "lmax2", "lmin1", "lmin2"]
if __name__ == "__main__":
    import numpy as np
    q = quantities(); rng = np.random.default_rng(3)
    for nm in NAMES4:
        for _ in range(3):
            v = {D: rng.uniform(.8, 1.3), L: rng.uniform(.1, .4), d: rng.uniform(.03, .2), M: rng.uniform(.1, .5),
                 al: rng.uniform(.5, 2), G: 1.2, H: .9, f: .94, a: -.06, ph: .35}
            print(nm, float(full(ell(nm) - q[nm]).subs(v)))
    # numerators of d delta/dX
    for X in [d, L, D]:
        num = sp.expand(pd(Nt, X)*Qt - Nt*pd(Qt, X))
        num = sp.expand(num.subs(s**2, 1))
        print("X =", X, "factor:", sp.factor(num))
