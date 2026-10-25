"""Numerical checks of the analytic derivatives of the no-boxcar tilt card formulas.
(1) displayed structures == direct sympy derivative (random points);
(2) analytic == central FD of the card formulas (tilt_formulas) at baseline;
(3) analytic vs central FD of the EXACT no-boxcar-with-tilt extrema.
Units: code units D = 1, G = H = 1 (so M is in D^3 H^2 / G). No fitting."""
import numpy as np, sympy as sp
from scipy.optimize import minimize_scalar
from nobox_tilt_deriv_sym import *
from nobox_tilt_check import tilt_formulas
from nobox_check import h_exact, K0
FU, AMP, PH = 0.94, -0.06, 0.35
BASE = {D: 1.0, L: 1/5.9, d: 0.1, M: K0, al: 1.0, G: 1.0, H: 1.0, f: FU, a: AMP, ph: PH}
ORDER = ["lmax1", "lmax2", "lmin1", "lmin2", "Amax1", "Amax2", "Amin1", "Amin2"]
q = quantities()
X5 = [D, L, d, M, al]

# ---- displayed structures (placeholder form), per input X and index sign sv ----
def disp(name, X):
    sv = 1 if name[-1] == "1" else -1
    sub = lambda e: e.subs(s, sv)
    dd = sub(-(a*sp.cos(ph)*pd(N, X) + delta*pd(Q, X))/Q)                 # d delta
    ddp = sub(deltap*(pd(C, X)/C + pd(K, X)/K - pd(R, X)/R))              # d delta'
    b = sp.Rational(3, 4)*d/D
    if name == "lmax1": return deg*dd
    if name == "lmax2": return -deg*dd
    if name == "lmin1": return deg*(pd(b, X) + ddp)
    if name == "lmin2": return -deg*(pd(b, X) + ddp)
    hX = g*pd(N, X) if X != M else (hpk - 1)/M                             # stated identity
    if name.startswith("Amax"):
        T = sub(Tmax); dT = -a*sp.cos(ph)*dd/2; h = sub(hpk); dh = sub(hX)
    else:
        T = sub(Tmin); dT = a*sp.sin(ph)*(pd(b, X) + ddp/2); h = hmin; dh = pd(hmin, X)
    return sub(T*h/Aall)*(dT/T + dh/h - pd(Aall, X)/Aall)

def num(e, vals):
    return float(full(e).subs(vals).evalf(20))

def card_vals(v):
    l, A, *_ = tilt_formulas(v[D], v[L], v[d], v[M], al=v[al], G=v[G], H=v[H])
    return dict(zip(ORDER, list(l) + list(A)))

def exact_vals(v):
    x = v[D]/v[L]; e = v[d]/v[D]; ka = v[G]*v[M]/(v[D]**3*v[H]**2); alp = v[al]
    th = np.linspace(0, 2*np.pi, 36000, endpoint=False)
    h, _ = h_exact(th, x, e, alp, ka); hbar = h.mean()
    Af = lambda t: (1 - AMP*np.sin(t + PH))*h_exact(np.atleast_1d(t), x, e, alp, ka)[0][0]/hbar
    A = (1 - AMP*np.sin(th + PH))*h/hbar
    out = []
    n = len(A)
    for i in range(n):
        a0, a1, a2 = A[i-1], A[i], A[(i+1) % n]
        kind = "max" if (a1 > a0 and a1 >= a2) else ("min" if (a1 < a0 and a1 <= a2) else None)
        if kind is None: continue
        sg = -1 if kind == "max" else 1
        r = minimize_scalar(lambda t: sg*Af(t), bounds=(th[i] - 2e-3, th[i] + 2e-3), method="bounded",
                            options={"xatol": 1e-11})
        out.append((kind, np.degrees(r.x) % 360, sg*r.fun))
    mx = [o for o in out if o[0] == "max"]; mn = [o for o in out if o[0] == "min"]
    assert len(mx) == 2 and len(mn) == 2, out
    m1 = min(mx, key=lambda o: min(o[1], 360 - o[1])); m2 = [o for o in mx if o is not m1][0]
    n1 = min(mn, key=lambda o: abs(o[1] - 90)); n2 = [o for o in mn if o is not n1][0]
    l = [m1[1] if m1[1] < 180 else m1[1] - 360, m2[1], n1[1], n2[1]]
    return dict(zip(ORDER, l + [m1[2], m2[2], n1[2], n2[2]]))

def fd(fun, X, rel):
    vp = dict(BASE); vm = dict(BASE); hstep = rel*BASE[X]
    vp[X] += hstep; vm[X] -= hstep
    fp, fm = fun(vp), fun(vm)
    return {k: (fp[k] - fm[k])/(2*hstep) for k in ORDER}

if __name__ == "__main__":
    import json
    rng = np.random.default_rng(1)
    worst = 0
    for trial in range(3):   # (1) displayed structure vs direct derivative at random points
        v = dict(BASE); v.update({D: rng.uniform(.8, 1.3), L: rng.uniform(.1, .4), d: rng.uniform(.03, .2),
                                  M: rng.uniform(.1, .5), al: rng.uniform(.5, 2)})
        for nm in ORDER:
            for X in INPUTS:
                a1 = num(pd(q[nm], X), v); a2 = num(disp(nm, X), v)
                worst = max(worst, abs(a1 - a2)/max(1e-12, abs(a1)))
    print("displayed-structure vs direct sympy, worst rel diff:", worst)
    ana = {nm: {str(X): num(pd(q[nm], X), BASE) for X in X5} for nm in ORDER}
    val = {nm: num(q[nm], BASE) for nm in ORDER}
    cfd = {str(X): fd(card_vals, X, 1e-5) for X in X5}
    efd = {str(X): fd(exact_vals, X, 1e-3) for X in X5}
    efd2 = {str(X): fd(exact_vals, X, 3e-4) for X in X5}
    ev = exact_vals(BASE)
    worst_c = 0
    for nm in ORDER:
        for X in X5:
            k = str(X); worst_c = max(worst_c, abs(ana[nm][k] - cfd[k][nm])/max(1e-9, abs(ana[nm][k])))
    print("analytic vs card FD, worst rel diff:", worst_c)
    for nm in ORDER:
        print(nm, "val %.5f exact %.5f" % (val[nm], ev[nm]))
        for X in X5:
            k = str(X)
            print("   d/d%-5s ana %+.5e  cardFD %+.5e  exactFD %+.5e (step2 %+.5e)" % (k, ana[nm][k], cfd[k][nm], efd[k][nm], efd2[k][nm]))
    # scaling identity check
    for nm in ORDER:
        sres = sum(BASE[X]*ana[nm][str(X)]*w for X, w in ((D, 1), (L, 1), (d, 1), (M, 3)))
        print("scaling identity", nm, sres)
    json.dump(dict(ana=ana, val=val, cfd=cfd, efd=efd, efd2=efd2, exact=ev,
                   base={str(k): v for k, v in BASE.items()}), open("nobox_tilt_deriv.json", "w"), indent=1)

# ---- hand-simplified delta' forms shown on the pages ----
def deltap_forms(X, sv):
    dp = deltap.subs(s, sv); c, sn = sp.cos(ph), sp.sin(ph)
    if X == M: return -dp*D**3*H**2/(M*(D**3*H**2 - al*G*M))
    if X == d: return sv*3*dp*c/(4*D*sn + sv*3*d*c)
    if X == D: return dp*(-sv*3*d*c/(D*(4*D*sn + sv*3*d*c)) + 3*al*G*M/(D*(D**3*H**2 - al*G*M))
                          + L*SL/(D*(S - D**3)) + SD/S)
    if X == L: return dp*(-SL/(S - D**3) + SL/S)
def check_deltap_forms():
    rng = np.random.default_rng(2); worst = 0
    for _ in range(3):
        v = dict(BASE); v.update({D: rng.uniform(.8, 1.3), L: rng.uniform(.1, .4), d: rng.uniform(.03, .2),
                                  M: rng.uniform(.1, .5), al: rng.uniform(.5, 2)})
        for sv in (1, -1):
            for X in INPUTS:
                a1 = num(pd(deltap.subs(s, sv), X), v); a2 = num(deltap_forms(X, sv), v)
                worst = max(worst, abs(a1 - a2)/abs(a1))
    return worst
