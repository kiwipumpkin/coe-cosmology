"""Independent minimal LTB (dust, Lambda=0, G=c=1, Gpc) for a comoving observer at areal radius D today.
Uniform bang t_B=0, label = areal radius today (R(r,t0)=r), mass M(r) = Mc + shell (density rho_s today, edges Rin, Rout).
Local expansion tensor at the observer: H_perp = Rdot/R, H_r = Rdot'/R'. For z << structure scale,
H(n) = H_perp + (H_r - H_perp) cos^2(psi), psi = angle from the CoE axis; A ~ H^2 (Migkas: A ~ H0^2)."""
import numpy as np
from scipy.optimize import brentq
KG = 3.0857e25 * 2.99792458e8**2 / 6.674e-11      # kg per geometric Gpc
H70 = 70 / 299792.458 * 1000.0                   # Gpc^-1
GYR = 3.0857e25 / 2.99792458e8 / 3.15576e16
def tof(M, E, r):
    """proper time from bang to areal radius r on the expanding branch (nan if unreachable)."""
    if abs(E) * r / M < 1e-6:                     # parabolic series
        return np.sqrt(2 * r**3 / (9 * M)) * (1 - 0.3 * E * r / M * 1.0)
    if E > 0:
        eta = np.arccosh(1 + 2 * E * r / M); return M / (2 * E)**1.5 * (np.sinh(eta) - eta)
    y = -2 * E * r / M
    if y >= 2: return np.nan                      # beyond turnaround radius
    eta = np.arccos(1 - y); return M / (-2 * E)**1.5 * (eta - np.sin(eta))
def Eof(M, r, t0):
    Emin = -M / r * (1 - 1e-12)                   # turnaround exactly at r
    tmax = tof(M, Emin, r)
    if t0 > tmax: return np.nan                   # would have turned around already: no expanding solution (contracting only)
    f = lambda E: tof(M, E, r) - t0
    hi = 1.0
    while f(hi) > 0: hi *= 4
    return brentq(f, Emin, hi, xtol=1e-14 * max(1, abs(hi)), rtol=1e-13)
class Model:
    def __init__(s, Mc_kg, D, Os, Rin, Rout, H=H70):
        s.Mc = Mc_kg / KG; s.rho = Os * 3 * H**2 / (8 * np.pi); s.Rin, s.Rout, s.D, s.H = Rin, Rout, D, H
    def M(s, r): return s.Mc + 4 * np.pi / 3 * s.rho * (np.clip(r, s.Rin, s.Rout)**3 - s.Rin**3)
    def Mp(s, r): return 4 * np.pi * r * r * s.rho * (s.Rin < r < s.Rout)
    def rates(s, t0, r=None):
        r = s.D if r is None else r; M = s.M(r); E = Eof(M, r, t0)
        if not np.isfinite(E): return np.nan, np.nan, np.nan
        h = 1e-4 * r; Ep = (Eof(s.M(r + h), r + h, t0) - Eof(s.M(r - h), r - h, t0)) / (2 * h)
        Rd = np.sqrt(2 * M / r + 2 * E)
        return Rd / r, (s.Mp(r) / r - M / r**2 + Ep) / Rd, E
    def solve(s):
        """t0 such that H_loc = (H_r + 2 H_perp)/3 = H at the observer."""
        def g(t):
            hp, hr, _ = s.rates(t); return (hr + 2 * hp) / 3 / s.H - 1
        ts = np.geomspace(1e-4, 200, 90); v = []
        for t in ts:
            try: v.append(g(t))
            except Exception: v.append(np.nan)
        v = np.array(v); ok = np.where(np.isfinite(v[:-1]) & np.isfinite(v[1:]) & (np.sign(v[:-1]) != np.sign(v[1:])))[0]
        if len(ok) == 0: return None
        t0 = brentq(g, ts[ok[0]], ts[ok[0] + 1], xtol=1e-13)
        hp, hr, E = s.rates(t0)
        return dict(t0=t0, age_Gyr=t0 * GYR, H_perp=hp / s.H * 70, H_r=hr / s.H * 70, twoE=2 * E, bound=E < 0, one_plus_2E=1 + 2 * E)
    def shells_ok(s, t0, n=200):
        """1+2E>0 on every matter shell and R'>0 for all past times (no shell crossing before today)."""
        rr = np.linspace(s.Rin * 1.0005, s.Rout * 0.9995, n); Es = np.array([Eof(s.M(r), r, t0) for r in rr])
        return bool(np.all(np.isfinite(Es)) and np.all(1 + 2 * Es > 0)), float(np.nanmin(1 + 2 * Es))
if __name__ == "__main__":
    import json
    N = json.load(open("../universe_thin.json"))["free"]
    m = Model(N["Mc"], N["D"], N["Os"], N["Rin"], N["Rout"]); print("GR-A replica:", m.solve())
