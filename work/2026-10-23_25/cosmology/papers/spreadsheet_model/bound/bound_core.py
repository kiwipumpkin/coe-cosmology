"""Bound-expansion fits of the integrated fixed-profile model (units D, t0=1/k).
CoE mass: e^{-r/L} profile, mass mu inside r=1 (mu = G M t0^2/D^3).
Background: uniform sphere of radius D (=1) fixed in space, mass mub = Ob/MCRIT_FAC (Ob = M_b/M_crit)."""
import sys, numpy as np, numba as nb
sys.path.insert(0, "../integrated")
from integ_core import Fx, geometry, AFfour, ML, MA, MLO, MHI, r2, A_from, best_sine, MCRIT_FAC, K, H0, GSI, RKM, GPC
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
TH1 = np.arange(0, 360, 1.0); IMIG = ML.astype(int)

@nb.njit(cache=True)
def acc(x, y, mu, LD, FD, mub):
    r = max(np.sqrt(x * x + y * y), 1e-12)
    f = -mu * Fx(r / LD) / FD / (r * r * r)
    f -= mub if r < 1.0 else mub / (r * r * r)
    return f * x, f * y

@nb.njit(cache=True, parallel=True)
def integrate(q, mu, LD, mub, taui, nsteps):
    N = q.shape[0]; FD = Fx(1.0 / LD)
    X = np.empty((N, 2)); V = np.empty((N, 2)); VR = np.empty(N)
    s0 = np.log(taui); h = -s0 / nsteps
    for i in nb.prange(N):
        x = q[i, 0] * taui; y = q[i, 1] * taui; vx = q[i, 0]; vy = q[i, 1]; s = s0; vrmin = 1e30
        for n in range(nsteps):
            t1 = np.exp(s); t2 = np.exp(s + 0.5 * h); t4 = np.exp(s + h)
            ax, ay = acc(x, y, mu, LD, FD, mub)
            k1x = t1 * vx; k1y = t1 * vy; l1x = t1 * ax; l1y = t1 * ay
            ax, ay = acc(x + 0.5*h*k1x, y + 0.5*h*k1y, mu, LD, FD, mub)
            k2x = t2 * (vx + 0.5*h*l1x); k2y = t2 * (vy + 0.5*h*l1y); l2x = t2 * ax; l2y = t2 * ay
            ax, ay = acc(x + 0.5*h*k2x, y + 0.5*h*k2y, mu, LD, FD, mub)
            k3x = t2 * (vx + 0.5*h*l2x); k3y = t2 * (vy + 0.5*h*l2y); l3x = t2 * ax; l3y = t2 * ay
            ax, ay = acc(x + h*k3x, y + h*k3y, mu, LD, FD, mub)
            k4x = t4 * (vx + h*l3x); k4y = t4 * (vy + h*l3y); l4x = t4 * ax; l4y = t4 * ay
            x += h/6*(k1x + 2*k2x + 2*k3x + k4x); y += h/6*(k1y + 2*k2y + 2*k3y + k4y)
            vx += h/6*(l1x + 2*l2x + 2*l3x + l4x); vy += h/6*(l1y + 2*l2y + 2*l3y + l4y)
            s += h
            vr = (x*vx + y*vy) / np.sqrt(x*x + y*y)
            if vr < vrmin: vrmin = vr
        X[i, 0] = x; X[i, 1] = y; V[i, 0] = vx; V[i, 1] = vy; VR[i] = vrmin
    return X, V, VR

def phi_pot(r, mu, LD, mub):
    """exact potential of the e^{-r/L} profile (total mass mu/F(1/LD)) plus the uniform sphere of radius 1."""
    u = r / LD; Mt = mu / Fx(1.0 / LD)
    p = -Mt * (Fx(u) / r + (1 + u) * np.exp(-u) / (2 * LD))
    p += -mub * (3 - r * r) / 2 if r < 1 else -mub / r
    return p

def run(DR, LD, mc, ob=0.0, taui=0.01, nsteps=800):
    q, rho = geometry(TH1, DR); mu = mc / MCRIT_FAC; mub = ob / MCRIT_FAC
    X, V, VR = integrate(q, mu, LD, mub, taui, nsteps)
    vrel = V[1:] - V[0]; af = np.hypot(vrel[:, 0], vrel[:, 1]) / rho
    rE = np.hypot(*X[0]); vE = np.hypot(*V[0])
    E = 0.5 * vE**2 + phi_pot(rE, mu, LD, mub)
    esc2 = vE**2 / (-2 * phi_pot(rE, mu, LD, mub))      # (v/v_esc)^2 ; bound if <= 1
    return dict(af=af, vrmin=float(VR.min()), E=float(E), vesc2=float(esc2), rE=float(rE), vE=float(vE), mu=mu, mub=mub)

def t_turn(DR, LD, mc, ob=0.0, rr=None):
    """time (units t0) at which Earth's outward motion stops; inf if unbound."""
    rr = rr or run(DR, LD, mc, ob)
    if rr["E"] >= 0: return np.inf
    mu, mub = rr["mu"], rr["mub"]; FD = Fx(1 / LD)
    def f(t, y):
        r = y[0]; g = mu * Fx(r / LD) / FD / r**2 + (mub * r if r < 1 else mub / r**2); return [y[1], -g]
    ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(f, (1.0, 1e4), [rr["rE"], rr["vE"]], events=ev, rtol=1e-10, atol=1e-12)
    return float(s.t_events[0][0]) if len(s.t_events[0]) else np.inf

def shifted(af, l0): return np.interp((TH1 - l0) % 360, np.append(TH1, 360), np.append(af, af[0]))

def score(p, target, norm="kR", need_bound=False):
    rr = run(p["DR"], p["LD"], p["mc"], p.get("ob", 0.0))
    af = rr["af"]
    if not np.all(np.isfinite(af)) or rr["vrmin"] <= 0: return -np.inf, None, None, rr      # (b) violated
    if need_bound and rr["E"] > 0: return -np.inf, None, None, rr                               # (a) violated
    if norm == "mean": af = af / af.mean()
    l0 = p.get("l0", 0.0); af = shifted(af, l0) if l0 else af; th = TH1 - l0
    if target == "mig": afx, thx, y = af[IMIG], th[IMIG], MA
    else: afx, thx, y = af, th, AFfour(TH1)
    a, ph = best_sine(afx, thx, y)
    return r2(y, A_from(afx, thx, a, ph)), a, ph, rr

LB = dict(DR=(1.05, 10.0), LD=(0.02, 3.0), mc=(1e-3, 30.0), ob=(0.0, 3.0), l0=(0, 360))
LOGV = dict(DR=False, LD=True, mc=True, ob=False, l0=False)
rng = np.random.default_rng(5)
def unpack(names, u):
    out = {}
    for n, ui in zip(names, u):
        lo, hi = LB[n]
        if n == "l0": out[n] = (ui % 1) * 360; continue
        ui = min(max(ui, 0), 1)
        out[n] = float(np.exp(np.log(lo) + ui*(np.log(hi)-np.log(lo)))) if LOGV[n] else lo + ui*(hi-lo)
    return out
def pack(names, p):
    u = []
    for n in names:
        lo, hi = LB[n]; v = p[n]
        if n == "l0": u.append(v/360); continue
        u.append((np.log(v)-np.log(lo))/(np.log(hi)-np.log(lo)) if LOGV[n] else (v-lo)/(hi-lo))
    return np.array(u)
def fit(target, free, fixed, norm="kR", need_bound=False, starts=6, nscreen=150, seeds=()):
    def full(u): p = dict(fixed); p.update(unpack(free, u)); return p
    def obj(u):
        pen = sum(10*(max(0, -x) + max(0, x-1))**2 for n, x in zip(free, u) if n != "l0")
        v = score(full(u), target, norm, need_bound)[0]
        return (1 - v if np.isfinite(v) else 50.0) + pen
    U = [rng.random(len(free)) for _ in range(nscreen)] + [pack(free, s) for s in seeds if all(n in s for n in free)]
    vals = [obj(u) for u in U]; order = np.argsort(vals)
    cand = [U[i] for i in order[:starts]]
    best = (np.inf, None)
    for u0 in cand:
        if vals[[id(x) for x in U].index(id(u0))] >= 50: continue
        r = minimize(obj, u0, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9, maxiter=1500))
        if r.fun < best[0]: best = (r.fun, r.x)
    if best[1] is None: return None
    p = full(best[1]); v, a, ph, rr = score(p, target, norm, need_bound)
    p.update(amp=a, phi=ph, R2=v, E=rr["E"], vesc2=rr["vesc2"], vrmin=rr["vrmin"])
    return p
