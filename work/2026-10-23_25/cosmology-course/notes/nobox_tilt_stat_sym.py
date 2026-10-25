"""Stationarity conditions d(ell)/d(input) = 0 for the four tilt-card angles (sympy).
S, SD, SL are code placeholders for the bracket and its derivatives only."""
import sympy as sp
from nobox_tilt_deriv_sym import *
ig = sp.Symbol("ig", positive=True)        # stands for H^2/(alpha G M) in code only
def to_ig(e): return e.subs(M, H**2/(al*G*ig))
Nig, Qig = to_ig(N), to_ig(Q)
b = sp.Rational(3, 4)*d/D
theta = b + deltap                          # ell_min = 90 + deg*theta (s=+1), 270 - deg*theta (s=-1)
def Fmax(X):
    """numerator of d(delta)/dX (zero <=> d ell_max/dX = 0): N_X Q - N Q_X."""
    return sp.expand(pd(N, X)*Q - N*pd(Q, X))
def Fmin(X):
    return pd(theta, X)
if __name__ == "__main__":
    for X in INPUTS:
        print("max", X, sp.factor(sp.expand(to_ig(Fmax(X)))))
