"""Isolation tests on the Newtonian side (model A exactly as prof_core: af = |dV|/separation, mean-normalised, linear).
 N0 baseline: static (fixed-in-space) profile, sources on the R = 14.3 Gpc vision sphere.
 N1 co-moving profile only: enclosed-mass table evaluated at r/tau (profile expands with the v = x/t flow), else identical.
 N2 static profile, small sampling sphere (R_sample = 0.43 Gpc ~ z = 0.1), else identical."""
import os, sys, json, numpy as np, numba as nb
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, "../.."); os.chdir(B)
sys.path.insert(0, B); sys.path.insert(0, os.path.join(B, "../integrated"))
import ptshell as S, prof_core as P
from bubble_core import MA, IMIG, TH1, r2, A_from, best_sine, geometry, RKM, GPC
old = json.load(open("thin_b_free_hug.json"))
geo = {k: old[k] for k in ("DR", "Os", "mcc", "l0", "rin", "rout", "Lc")}

@nb.njit(cache=True)
def acc(x, y, tau, mu, RGt, MT, como):
    r = max(np.sqrt(x * x + y * y), 1e-12)
    m = np.interp(r / tau if como else r, RGt, MT)
    f = -mu * m / (r * r * r); return f * x, f * y

@nb.njit(cache=True)
def integ(q, mu, RGt, MT, como, taui, nsteps):
    N = q.shape[0]; X = np.empty((N, 2)); V = np.empty((N, 2)); s0 = np.log(taui); h = -s0 / nsteps
    for i in range(N):
        x = q[i, 0] * taui; y = q[i, 1] * taui; vx = q[i, 0]; vy = q[i, 1]; s = s0
        for n in range(nsteps):
            t1 = np.exp(s); t2 = np.exp(s + 0.5 * h); t4 = np.exp(s + h)
            ax, ay = acc(x, y, t1, mu, RGt, MT, como); k1x = t1*vx; k1y = t1*vy; l1x = t1*ax; l1y = t1*ay
            ax, ay = acc(x+0.5*h*k1x, y+0.5*h*k1y, t2, mu, RGt, MT, como); k2x = t2*(vx+0.5*h*l1x); k2y = t2*(vy+0.5*h*l1y); l2x = t2*ax; l2y = t2*ay
            ax, ay = acc(x+0.5*h*k2x, y+0.5*h*k2y, t2, mu, RGt, MT, como); k3x = t2*(vx+0.5*h*l2x); k3y = t2*(vy+0.5*h*l2y); l3x = t2*ax; l3y = t2*ay
            ax, ay = acc(x+h*k3x, y+h*k3y, t4, mu, RGt, MT, como); k4x = t4*(vx+h*l3x); k4y = t4*(vy+h*l3y); l4x = t4*ax; l4y = t4*ay
            x += h/6*(k1x+2*k2x+2*k3x+k4x); y += h/6*(k1y+2*k2y+2*k3y+k4y); vx += h/6*(l1x+2*l2x+2*l3x+l4x); vy += h/6*(l1y+2*l2y+2*l3y+l4y); s += h
        X[i, 0] = x; X[i, 1] = y; V[i, 0] = vx; V[i, 1] = vy
    return X, V

def af_of(como=False, DRs=None):
    Sd, MT, m1 = S.build(geo)
    q, rr = geometry(TH1, DRs or geo["DR"]); Ht = 0.95
    for _ in range(3):
        mu = m1 * Ht**2 / 2; X, V = integ(q, mu, P.RG, MT, como, 0.01, 200); Ht = np.hypot(*V[0]) / np.hypot(*X[0])
    af = np.hypot(*(V[1:] - V[0]).T) / rr; return af / af.mean(), Ht

def fitall(af):
    """orientation (af along the CoE axis / across), R2 at the Newtonian l0 (+sine), and with l0 refit (+sine)."""
    along = np.mean(af[[0, 180]]); across = np.mean(af[[90, 270]])
    def r2s(l0):
        A = np.interp((TH1 - l0) % 360, np.append(TH1, 360), np.append(af, af[0])); x = A[IMIG]; th = (TH1 - l0)[IMIG]
        am, ph = best_sine(x, th, MA); return r2(MA, x), r2(MA, A_from(x, th, am, ph))
    ns0, s0 = r2s(geo["l0"]); G = [(l0,) + r2s(l0) for l0 in np.arange(0, 360, 1.0)]
    bns = max(G, key=lambda g: g[1]); bs = max(G, key=lambda g: g[2])
    return dict(along_over_across=float(along / across), orientation="along" if along > across else "across",
                R2_ns_l0fixed=float(ns0), R2_s_l0fixed=float(s0), R2_ns_l0free=float(bns[1]), l0_ns=float(bns[0]), R2_s_l0free=float(bs[2]), l0_s=float(bs[0]))
out = {}
for name, kw in (("N0 Newtonian baseline (static profile, 14.3 Gpc sphere)", {}), ("N1 co-moving profile only", dict(como=True)),
                 ("N2 static profile, 0.43 Gpc sphere (~z=0.1)", dict(DRs=geo["DR"] * RKM / GPC / 0.43))):
    af, Ht = af_of(**kw); d = fitall(af); d["age_Gyr"] = Ht * 13.968; out[name] = d; print(name, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in d.items()})
json.dump(out, open(os.path.join(HERE, "iso_newton.json"), "w"), indent=1)
