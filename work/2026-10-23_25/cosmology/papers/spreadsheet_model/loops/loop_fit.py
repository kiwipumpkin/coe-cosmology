import json, sys, numpy as np
from scipy.optimize import minimize, minimize_scalar
sys.path.insert(0, "/workspace/cosmology/figures"); sys.path.insert(0, "../bubble")
from coe_half_loop_geometry import lb_to_n, n_to_lb, ring_basis, MIGKAS_MAXES, MIGKAS_MINS
import os
_here = os.getcwd(); os.chdir('../bubble')
import analysis_curve as AC
os.chdir(_here)
from bubble_core import TH1
d = json.load(open("../bubble/thin_b_free_hug.json")); A = AC.curve("ptshell", d); l0 = d["l0"]
ext = []
n = len(A)
for i in range(n):
    a, b, c = A[i - 1], A[i], A[(i + 1) % n]
    if (b > a and b >= c) or (b < a and b <= c):
        off = 0.5 * (a - c) / (a - 2 * b + c)          # parabolic refinement
        ext.append(("MAX" if b > a else "MIN", float((TH1[i] + off) % 360), float(b)))
print(ext)
maxs = [e for e in ext if e[0] == "MAX"]; mins = [e for e in ext if e[0] == "MIN"]
# pair with Migkas markers by nearest longitude
def pair(lst, marks):
    out = []
    for m in marks:
        e = min(lst, key=lambda e: abs(((e[1] - m[0] + 180) % 360) - 180)); out.append((e, m))
    return out
pairs = [("MAX",) + p for p in pair(maxs, MIGKAS_MAXES)] + [("MIN",) + p for p in pair(mins, MIGKAS_MINS)]
TH = [(p[1][1] - l0) % 360 for p in pairs]          # loop angle from anti-CoE
P = [lb_to_n(*p[2]) for p in pairs]
ang = lambda u, v: np.degrees(np.arccos(np.clip(np.dot(u, v), -1, 1)))
def ringres(D):
    return np.array([ang(D, Pi) - np.degrees(np.arccos(np.cos(np.radians(t)))) for Pi, t in zip(P, TH)])
def theta_pt(D, Y, t): return np.cos(np.radians(t)) * D + np.sin(np.radians(t)) * Y
def rayres(D, roll): Y, _ = ring_basis(D, roll); return np.array([ang(theta_pt(D, Y, t), Pi) for Pi, t in zip(P, TH)])
res = {}
D0 = lb_to_n(l0, 0.0); res["untilted"] = dict(l=l0, b=0.0, ring=ringres(D0).tolist())
f = lambda b: np.sum(ringres(lb_to_n(l0, b))**2); r = minimize_scalar(f, bounds=(-60, 60), method="bounded")
res["b0_only"] = dict(l=l0, b=float(r.x), ring=ringres(lb_to_n(l0, r.x)).tolist())
best = None
for lg in range(-40, 41, 10):
    for bg in range(-40, 41, 10):
        rr = minimize(lambda x: np.sum(ringres(lb_to_n(x[0], x[1]))**2), [l0 + lg, bg], method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9))
        if best is None or rr.fun < best.fun: best = rr
lf, bf = best.x; lf %= 360
res["free"] = dict(l=float(lf), b=float(bf), ring=ringres(lb_to_n(lf, bf)).tolist())
for key in ("untilted", "b0_only", "free"):
    D = lb_to_n(res[key]["l"], res[key]["b"])
    rolls = np.arange(0, 360, 0.25); cost = [np.sum(rayres(D, q)**2) for q in rolls]
    q = float(rolls[int(np.argmin(cost))]); res[key]["roll"] = q; res[key]["ray"] = rayres(D, q).tolist()
    res[key]["rms_ring"] = float(np.sqrt(np.mean(np.square(res[key]["ring"])))); res[key]["rms_ray"] = float(np.sqrt(np.mean(np.square(res[key]["ray"]))))
out = dict(l0=l0, amp=d["amp"], phi=d["phi"], R2=d["R2_mig"], extrema=ext,
           pairs=[dict(kind=p[0], l_model=p[1][1], A=p[1][2], theta=t, migkas=list(p[2])) for p, t in zip(pairs, TH)], fits=res)
json.dump(out, open("loop_fit.json", "w"), indent=1)
for k, v in res.items(): print(k, {kk: (np.round(vv, 2).tolist() if isinstance(vv, list) else round(vv, 3)) for kk, vv in v.items()})
for p, t in zip(pairs, TH): print(p[0], "l=%.1f theta=%.1f" % (p[1][1], t), "Migkas", p[2])
