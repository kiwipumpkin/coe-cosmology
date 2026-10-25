"""Thin shell + core ("b|free|hug") refit with NO added sine term.
Model A = af (mean-normalised over 360 deg, shifted by l0), compared directly with the Migkas 2020 A(l) points.
A_from(af, th, amp, phi) = (1 - amp*sin(th+phi)) * af  (pure multiplicative sine, no scale factor), so amp=0 <=> A = af exactly.
Same bounds (thin_shell overrides), same screening (rng seed 7, 160 random + 24 seeds, 4 Nelder-Mead starts), plus the old best as an extra seed.
Originals (thin_shell.py, thin_b_free_hug.json, universe_thin.json) are only read."""
import json, time, numpy as np
import thin_shell as T                        # sets LB/LOGV exactly as the original fit
import ptshell as S, prof_core as P
from bubble_core import MA, IMIG, TH1, r2, shifted, A_from
from scipy.optimize import minimize

def model_nosine(p, norm="mean"):
    q = dict(p, mc=0.0, Ob=0.0)
    if p["rin"] >= p["rout"]: return None                      # same guard as S.score(inside=False)
    rr = P.run("ptshell", q, 0.0, 0.0)
    af = rr["af"]
    if (not np.all(np.isfinite(af))) or rr["vrmin"] <= 0: return None   # same guard as P.score
    af = af / af.mean() if norm == "mean" else af**2 / np.mean(af**2)
    l0 = p.get("l0", 0.0); return shifted(af, l0) if l0 else af

def score_nosine(p, norm="mean"):
    A = model_nosine(p, norm)
    return -np.inf if A is None else r2(MA, A[IMIG])

def fit(free, fixed, kin, kout, norm="mean", seeds=(), starts=4, nscreen=160):
    def full(u): p = dict(fixed); p.update(P.unpack(free, u)); return T.geom(p, kin, kout)
    def obj(u):
        pen = sum(10*(max(0, -x) + max(0, x-1))**2 for n, x in zip(free, u) if n != "l0")
        v = score_nosine(full(u), norm)
        return (1 - v if np.isfinite(v) else 50.0) + pen
    rng = np.random.default_rng(7)
    U = [rng.random(len(free)) for _ in range(nscreen)] + [P.pack(free, dict(fixed, **s)) for s in seeds]
    vals = np.array([obj(u) for u in U]); order = list(np.argsort(vals)[:starts])
    if len(U) - 1 not in order: order.append(len(U) - 1)     # always also start from the old best (last seed)
    runs = []
    for i in order:
        if vals[i] >= 50: continue
        r = minimize(obj, U[i], method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-7, maxiter=700))
        runs.append(dict(start=int(i), R2_start=float(1 - vals[i]), R2_end=float(1 - r.fun), x=r.x.tolist(), p={k: float(v) for k, v in P.unpack(free, r.x).items()}))
        print("  start", i, "R2 %.4f -> %.4f" % (1 - vals[i], 1 - r.fun), runs[-1]["p"], flush=True)
    best = min(runs, key=lambda z: -z["R2_end"])
    p = full(np.array(best["x"])); p.update(amp=0.0, phi=0.0, R2=score_nosine(p, norm))
    d = S.describe(p, norm); d.update(kin=kin, kout=kout, norm=norm); return d, runs

def scale_r2(A):   # diagnostic only: best pure overall scale c*A (not used in fit)
    x = A[IMIG]; c = np.dot(x, MA) / np.dot(x, x); return float(c), float(r2(MA, c * x))

if __name__ == "__main__":
    t0 = time.time()
    old = json.load(open("thin_b_free_hug.json"))
    geo = {k: old[k] for k in ("DR", "Os", "mcc", "l0", "rin", "rout", "Lc")}
    # check 1: reproduce old R2 with sine
    q = dict(geo, mc=0.0, Ob=0.0); v_old = S.score(q, "mig", "mean", inside=False)[0]
    A0 = model_nosine(geo); v_rm = float(r2(MA, A0[IMIG]))
    Aw = A_from(A0, TH1 - geo["l0"], old["amp"], old["phi"]); v_chk = float(r2(MA, Aw[IMIG]))
    print("old R2 with sine (recomputed) %.6f  json %.6f  via A_from %.6f" % (v_old, old["R2"], v_chk))
    print("old geometry, sine removed: R2 = %.6f" % v_rm, " (diag: best overall scale c=%.4f -> R2 %.4f)" % scale_r2(A0), flush=True)
    free, fixed, kin, kout, norm = T.JOBS["b|free|hug"]
    seeds = []
    for DR in (1.5, 2.9, 5.0):
        for l0 in (13.0, 103.0, 193.0, 283.0):
            seeds.append(dict(DR=DR, Os=0.3, mcc=0.002, l0=l0)); seeds.append(dict(DR=DR, Os=0.05, mcc=0.02, l0=l0))
    seeds.append(dict(DR=old["DR"], Os=old["Os"], mcc=old["mcc"], l0=old["l0"]))
    d, runs = fit(free, fixed, kin, kout, norm, seeds=seeds)
    A1 = model_nosine(d)
    out = dict(refit=d, starts=runs,
               old_with_sine=dict(geo, amp=old["amp"], phi=old["phi"], R2=float(v_old), D_Gpc=old["D_Gpc"]),
               old_sine_removed=dict(geo, R2=v_rm, D_Gpc=old["D_Gpc"], diag_scale=scale_r2(A0)),
               refit_diag_scale=scale_r2(A1),
               note="A_from = (1 - amp sin(th+phi)) * af; no scale factor. No-sine model: A = shifted(af/mean(af), l0); no extra normalisation.")
    json.dump(out, open("thin_nosine.json", "w"), indent=1, default=float)
    print("REFIT", round(time.time() - t0), {k: round(d[k], 5) for k in ("DR", "Os", "mcc", "l0", "R2", "D_Gpc")}, flush=True)
