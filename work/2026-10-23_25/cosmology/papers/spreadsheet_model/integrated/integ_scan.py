import sys, json, time, numpy as np
sys.path.insert(0, "."); from integ_core import *
from scipy.optimize import minimize, minimize_scalar
rng = np.random.default_rng(3)
TH1 = np.arange(0, 360, 1.0)                      # 1-degree grid (contains the Migkas l values)
IMIG = ML.astype(int)
CFG = dict(como=False, taui=0.01, nsteps=800)
def af_grid(DR, LD, mu, cfg=None):
    c = CFG if cfg is None else cfg
    return AF(TH1, DR, LD, mu, **c)[0]
def shifted(af, l0):                               # value of AF at theta = l - l0, l on TH1
    return np.interp((TH1 - l0) % 360, np.append(TH1, 360), np.append(af, af[0]))
def score(DR, LD, mu, target, l0=0.0, cfg=None, amp=None, phi=None):
    af = af_grid(DR, LD, mu, cfg)
    if not np.all(np.isfinite(af)): return -np.inf, np.nan, np.nan
    if l0: af = shifted(af, l0)
    th = TH1 - l0
    if target == "mig": afx, thx, y = af[IMIG], th[IMIG], MA
    else: afx, thx, y = af, th, AFfour(TH1)
    if amp is None and phi is None: a, p = best_sine(afx, thx, y)
    elif phi is None:                               # amp fixed: phase on a grid then refine
        ph = np.linspace(0, 2*np.pi, 361); v = [np.sum((y - A_from(afx, thx, amp, q))**2) for q in ph]
        p = minimize_scalar(lambda q: np.sum((y - A_from(afx, thx, amp, q))**2), bracket=(ph[np.argmin(v)]-0.02, ph[np.argmin(v)]+0.02)).x; a = amp
    elif amp is None:                               # phase fixed: amp is linear
        t = np.radians(thx); b = afx * np.sin(t + phi); a = float(np.clip(np.dot(b, afx - y) / np.dot(b, b), -0.15, 0.15)); p = phi
    else: a, p = amp, phi
    return r2(y, A_from(afx, thx, a, p)), a, p % (2*np.pi)
LB = dict(DR=(1.05, 10.0), LD=(0.02, 3.0), mc=(1e-3, 10.0), l0=(0, 360))
LOGV = dict(DR=False, LD=True, mc=True, l0=False)
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
def fit(target, free, fixed, starts=10, nscreen=300, seeds=(), cfg=None, amp=None, phi=None):
    def full(u):
        p = dict(fixed); p.update(unpack(free, u)); return p
    def obj(u):
        pen = sum(10*(max(0, -x) + max(0, x-1))**2 for n, x in zip(free, u) if n != "l0")
        p = full(u)
        v = score(p["DR"], p["LD"], p["mc"]/MCRIT_FAC, target, p.get("l0", 0.0), cfg, amp, phi)[0]
        return (1 - v if np.isfinite(v) else 50.0) + pen
    U = [rng.random(len(free)) for _ in range(nscreen)] + [pack(free, s) for s in seeds]
    vals = [obj(u) for u in U]; order = np.argsort(vals)
    cand = [U[i] for i in order[:starts]] + [pack(free, s) for s in seeds]
    best = (np.inf, None)
    for u0 in cand:
        r = minimize(obj, u0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-11, maxiter=3000))
        if r.fun < best[0]: best = (r.fun, r.x)
    p = full(best[1]); v, a, ph = score(p["DR"], p["LD"], p["mc"]/MCRIT_FAC, target, p.get("l0", 0.0), cfg, amp, phi)
    p.update(amp=a, phi=ph, R2=v); return p
if __name__ == "__main__":
    t0 = time.time(); out = {}
    for target in ("mig", "four"):
        for l0 in (False, True):
            free = ["DR", "LD", "mc"] + (["l0"] if l0 else [])
            p = fit(target, free, {} if l0 else {"l0": 0.0}, starts=20, nscreen=600)
            p["R2_mig"] = score(p["DR"], p["LD"], p["mc"]/MCRIT_FAC, "mig", p["l0"], amp=p["amp"], phi=p["phi"])[0]
            p["R2_four"] = score(p["DR"], p["LD"], p["mc"]/MCRIT_FAC, "four", p["l0"], amp=p["amp"], phi=p["phi"])[0]
            key = target + ("_l0" if l0 else ""); out[key] = p; print(key, p, round(time.time()-t0))
    # his geometry + sine, M refitted
    p = fit("mig", ["mc"], {"DR": 5.0, "LD": 0.18, "l0": 0.0}, starts=5, amp=-0.06, phi=0.35)
    p["R2_four"] = score(5.0, 0.18, p["mc"]/MCRIT_FAC, "four", 0.0, amp=-0.06, phi=0.35)[0]
    out["his_geom_Mfit"] = p; print("his geom", p)
    json.dump(out, open("integ_results.json", "w"), indent=1)

def vrmin_of(p, cfg=None):
    c = dict(CFG if cfg is None else cfg); c.pop("shoot", None)
    q0, _ = geometry(TH1, p["DR"]); X, V, VR = integrate(q0, p["mc"]/MCRIT_FAC, p["LD"], c["como"], c["taui"], c["nsteps"])
    return float(VR.min() / 1.0)
