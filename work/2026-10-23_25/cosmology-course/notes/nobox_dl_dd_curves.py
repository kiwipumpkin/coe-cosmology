"""Angles vs d/D (others at baseline): exact no-boxcar curve with tilt (and without tilt,
for the geometric check), first-order card formula, and the arccos form for the minima. No fits."""
import json, numpy as np
from scipy.optimize import minimize_scalar
from nobox_check import h_exact, K0
FU, AMP, PH = 0.94, -0.06, 0.35
X0, KA, AL = 5.9, K0, 1.0
def exact4(e, amp=AMP, n=36000):
    th = np.linspace(0, 2*np.pi, n, endpoint=False)
    h, _ = h_exact(th, X0, e, AL, KA); hbar = h.mean()
    Af = lambda t: (1 - amp*np.sin(t + PH))*h_exact(np.atleast_1d(t), X0, e, AL, KA)[0][0]/hbar
    A = (1 - amp*np.sin(th + PH))*h/hbar; out = []
    for i in range(n):
        a0, a1, a2 = A[i-1], A[i], A[(i+1) % n]
        kind = "max" if (a1 > a0 and a1 >= a2) else ("min" if (a1 < a0 and a1 <= a2) else None)
        if kind is None: continue
        sg = -1 if kind == "max" else 1
        r = minimize_scalar(lambda t: sg*Af(t), bounds=(th[i] - 2e-3, th[i] + 2e-3), method="bounded", options={"xatol": 1e-11})
        out.append((kind, np.degrees(r.x) % 360))
    mx = [o[1] for o in out if o[0] == "max"]; mn = [o[1] for o in out if o[0] == "min"]
    if len(mx) != 2 or len(mn) != 2: return [np.nan]*4
    m1 = min(mx, key=lambda o: min(o, 360 - o)); m2 = [o for o in mx if o != m1][0]
    n1 = min(mn, key=lambda o: abs(o - 90)); n2 = [o for o in mn if o != n1][0]
    return [m1 if m1 < 180 else m1 - 360, m2, n1, n2]
def formulas(e):
    D, L, d, M = 1.0, 1/X0, e, KA
    S = (1-FU)*(6*L**3*np.exp(D/L) - 6*L**3 - 6*D*L**2 - 3*D**2*L) + FU*D**3; w = 1 - FU; g = AL*M
    c, sn = np.cos(PH), np.sin(PH)
    def dl(sv):   # l_max explicit
        Nt = 2*L*D**4*S + g*(4*L*D*S - 6*L*D**4) - sv*3*d*g*(2*L*(S - D**3) - w*D**4)
        Qt = 12*L*D*(S - D**3) - sv*9*d*(4*L*(S - D**3) - w*D**4)
        return -AMP*c*Nt/(g*Qt)
    def th_t(sv): return -AMP*(4*D*sn + sv*3*d*c)*(D**3 - g)*S/(24*D*g*(S - D**3))   # delta'
    lin = [np.degrees(dl(1)), 180 - np.degrees(dl(-1)), 90 + np.degrees(3*d/(4*D) + th_t(1)), 270 - np.degrees(3*d/(4*D) + th_t(-1))]
    ac = np.degrees(np.arccos(-3*d/(4*D)))
    arc = [ac + np.degrees(th_t(1)), 360 - ac - np.degrees(th_t(-1))]
    corr = AMP*c*(D**3 - g)*S/(6*g*(S - D**3))
    dlin = [None, None, 135/(np.pi*D)*(1 - corr), -135/(np.pi*D)*(1 + corr)]
    root = np.sqrt(1 - 9*d**2/(16*D**2))
    darc = [135/(np.pi*D)*(1/root - corr), -135/(np.pi*D)*(1/root + corr)]
    return lin, arc, dlin, darc
if __name__ == "__main__":
    es = np.round(np.arange(0.005, 0.9001, 0.005), 4)
    ex = np.array([exact4(e) for e in es]); ex0 = np.array([exact4(e, amp=0.0) for e in es])
    lin = np.array([formulas(e)[0] for e in es]); arc = np.array([formulas(e)[1] for e in es])
    np.savez("nobox_dl_dd_curves.npz", es=es, ex=ex, ex0=ex0, lin=lin, arc=arc)
    print("nan exact rows:", int(np.isnan(ex).any(1).sum()), "first nan at", es[np.isnan(ex).any(1)][:3])
    for e in (0.1, 0.2, 0.3, 0.5, 0.7, 0.9):
        i = int(np.argmin(abs(es - e)))
        print(e, "exact", np.round(ex[i], 2), "lin", np.round(lin[i], 2), "arc", np.round(arc[i], 2), "notilt exact", np.round(ex0[i], 2))
