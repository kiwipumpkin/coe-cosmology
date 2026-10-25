"""Symbolic derivatives of the no-boxcar tilt card formulas (sympy).
Placeholders S, SD, SL stand for the bracket and its D, L derivatives ONLY inside
the code; the rendered pages write them out in full."""
import sympy as sp
D, L, d, M, al, G, H = sp.symbols("D L d M alpha G H", positive=True)
f, a, ph = sp.symbols("f_u Amp_A phi_A", real=True)
s = sp.Symbol("s")          # +1 -> index 1 (upper sign), -1 -> index 2
S, SD, SL = sp.symbols("S S_D S_L", positive=True)
Sfun = (1-f)*(6*L**3*sp.exp(D/L) - 6*L**3 - 6*D*L**2 - 3*D**2*L) + f*D**3
SDfun = sp.diff(Sfun, D); SLfun = sp.diff(Sfun, L)
INPUTS = [D, L, d, M]
ig = H**2/(al*G*M); g = al*G*M/H**2
N = ig + 2/D**3 - 3/S - s*3*d/D*(1/D**3 - (2*L + (1-f)*D)/(2*L*S))
Q = 6/D**3 - 6/S - s*sp.Rational(9, 2)*d/D*(4/D**3 - (4*L + (1-f)*D)/(L*S))
delta = -a*sp.cos(ph)*N/Q
C = sp.sin(ph) + s*sp.Rational(3, 4)*d/D*sp.cos(ph)
K = ig - 1/D**3
R = 6/D**3 - 6/S
deltap = -a*C*K/R
Aall = 1 + g/2*(1/D**3 - 3/S)
hpk = 1 + g*(2/D**3 - 3/S - s*3*d/D*(1/D**3 - (2*L + (1-f)*D)/(2*L*S)))
hmin = 1 - g/D**3
Tmax = 1 - s*a*sp.sin(ph) - a*sp.cos(ph)*delta/2
Tmin = 1 - s*a*sp.cos(ph) + sp.Rational(3, 4)*d/D*a*sp.sin(ph) + a*sp.sin(ph)*deltap/2
deg = 180/sp.pi
def quantities():
    """8 card quantities (placeholder form)."""
    return {
        "lmax1": deg*delta.subs(s, 1), "lmax2": 180 - deg*delta.subs(s, -1),
        "lmin1": 90 + deg*(sp.Rational(3, 4)*d/D + deltap.subs(s, 1)),
        "lmin2": 270 - deg*(sp.Rational(3, 4)*d/D + deltap.subs(s, -1)),
        "Amax1": (Tmax*hpk/Aall).subs(s, 1), "Amax2": (Tmax*hpk/Aall).subs(s, -1),
        "Amin1": (Tmin*hmin/Aall).subs(s, 1), "Amin2": (Tmin*hmin/Aall).subs(s, -1)}
def pd(expr, X):
    """Partial derivative wrt input X with S = S(D, L)."""
    out = sp.diff(expr, X)
    if X == D: out += sp.diff(expr, S)*SD
    if X == L: out += sp.diff(expr, S)*SL
    return out
def full(expr):
    return expr.subs({SD: SDfun, SL: SLfun}).subs(S, Sfun)
def collect_terms(expr):
    """Expand and group by monomials in 1/S, SD, SL; combine each coefficient."""
    e = sp.expand(expr)
    groups = {}
    for t in sp.Add.make_args(e):
        mono = sp.Integer(1); coef = sp.Integer(1)
        for fac in sp.Mul.make_args(t):
            b, ex = fac.as_base_exp()
            if b in (S, SD, SL): mono *= fac
            else: coef *= fac
        groups[mono] = groups.get(mono, 0) + coef
    return [(sp.factor(sp.together(c)), m) for m, c in groups.items() if sp.simplify(c) != 0]
if __name__ == "__main__":
    for nm, ex in (("N", N), ("Q", Q), ("R", R), ("Aall", Aall), ("hpk", hpk), ("C", C), ("K", K)):
        for X in INPUTS:
            print(f"d{nm}/d{X}:", collect_terms(pd(ex, X)))
    print("SD:", sp.factor(SDfun)); print("SL:", sp.factor(SLfun))
    print("Euler:", sp.simplify(D*SDfun + L*SLfun - 3*Sfun))
