"""Way-out A: interior-profile ladder. For each profile and D, solve (mcc, l0, sine amp, phi) exactly on the slice's dipole+quadrupole terms
(same equations as termslice/solve.py), then score R2 against the rebuilt slice profile and compute CMB-scale diagnostics at R = 14.3 Gpc."""
import json, sys, numpy as np, multiprocessing as mp
from scipy.optimize import least_squares
import lib4 as L
T = L.T; sys.path.insert(0, "../cycle3"); import migkas_A as M
DQ = json.load(open("../migkas_zslices_dipquad.json")); S2 = json.load(open("../cycle3/slices2.json"))
z = M.C["z"]; rcl = 299792.458 * z / 71.6 / 1e3
SL = {"z>0.07": (0.07, 1.0), "all": (0, 1.0)}
PROF = [("point", {}), ("pl", {"n": 2.0}), ("pl", {"n": 1.0}), ("pl", {"n": 0.0}), ("st", {}), ("hollow", {"rh": 0.5}),
        ("remnant", {"f": 0.5}), ("remnant", {"f": 0.2}), ("remnant", {"f": 0.05}), ("comp_st", {}), ("comp_pl2", {})]
def r2(y, m): return 1 - np.sum((y - m)**2) / np.sum((y - y.mean())**2)
def solve(s, D, kind, a):
    d = DQ[s]; tg = T.target(d["a1"], d["p1"], d["a2"], d["p2"]); lo, hi = SL[s]; sel = (z >= lo) & (z < hi); Rs = float(np.median(rcl[sel]))
    cache = {}
    def afc(mcc):
        k = round(np.log10(mcc), 9)
        if k not in cache: cache[k] = L.run(D, mcc, kind, a, [Rs])[0][0]["af"]
        return cache[k]
    def res(x):
        lm, l0, amp, phi = x; A = T.A_from(T.shifted(afc(10**lm), l0), T.TH1 - l0, amp, phi); A = A / A.mean()
        c1, s1 = T.harm(A, 1); c2, s2 = T.harm(A, 2); return np.array([c1 - tg[1][0], s1 - tg[1][1], c2 - tg[2][0], s2 - tg[2][1]])
    # scan mcc for the quad amplitude first
    grid = np.geomspace(1e-5, 1.0, 46); a2s = []
    for m in grid:
        try: af = afc(m); a2s.append(T.dq(af)["a2"] if np.all(np.isfinite(af)) else np.nan)
        except Exception: a2s.append(np.nan)
    a2s = np.array(a2s); best = None
    ok = np.isfinite(a2s)
    if not ok.any(): return dict(slice=s, D=D, kind=kind, a=a, err="no finite")
    up = np.where(ok & (a2s >= d["a2"]))[0]                       # first crossing of the target quad amplitude (lowest-mass branch)
    i0 = int(up[0]) if len(up) else int(np.nanargmin(np.abs(np.log(np.where(ok, a2s, 1e-9) / d["a2"]))))
    if len(up) and up[0] > 0: i0 = int(up[0])
    p2m = T.dq(afc(grid[i0]))["p2"]
    for l0g in (d["p2"] - p2m, d["p2"] - p2m + 180):
        for pg in np.linspace(-np.pi, np.pi, 3, endpoint=False):
            try: r = least_squares(res, [np.log10(grid[i0]), l0g, -d["a1"], pg], x_scale=[0.1, 3, 0.02, 0.3], diff_step=1e-4)
            except Exception: continue
            if best is None or r.cost < best.cost: best = r
    lm, l0, amp, phi = best.x; mcc = 10**lm
    A = T.A_from(T.shifted(afc(mcc), l0), T.TH1 - l0, amp, phi); A = A / A.mean()
    prof = np.array(S2[s]["A"], float); okp = np.isfinite(prof); y = prof[okp] / np.nanmean(prof)
    c = L.cmb(D, mcc, kind, a, Rs_cl=Rs)
    rin = 1 - T.KIN * L.RV / D; Mh = 0.049 * rin**3
    if kind.startswith("comp_") and mcc > Mh: return dict(slice=s, D=D, kind=kind, a=a, err="compensation impossible: cavity mass %.5f < core %.5f" % (Mh, mcc))
    return dict(slice=s, D=D, kind=kind, a=a, Rs=Rs, resid=float(np.sqrt(2 * best.cost)), mcc=mcc, Mhole=Mh, Meff=(Mh if kind.startswith("comp") and mcc < Mh else mcc) - Mh,
                l0=float(l0 % 360), CoE_axis_offset=float(((l0 - d["p2"]) + 90) % 180 - 90), amp=float(amp), phi=float(phi % (2 * np.pi)),
                R2=float(r2(y, np.interp(M.LC[okp], T.TH1, A, period=360))), a2grid=[grid.tolist(), a2s.tolist()], curve=A.tolist(), **{k: c[k] for k in ("cl_a2", "cl_p2", "far_a2", "ratio", "kin", "stat", "vrmin")})
JOBS = [(s, D, k, a) for s in SL for D in (14.7, 16.0, 20.0, 22.3, 30.0, 49.1746, 94.0) for (k, a) in PROF]
def job(x):
    try: return solve(*x)
    except Exception as e: return dict(slice=x[0], D=x[1], kind=x[2], a=x[3], err=repr(e))
if __name__ == "__main__":
    with mp.Pool(7) as Pp: R = Pp.map(job, JOBS, chunksize=1)
    json.dump(R, open("solveA.json", "w"), default=float)
    for r in R:
        if "err" in r: print(r["slice"], r["D"], r["kind"], r["a"], "ERR", r["err"][:150]); continue
        print("%-6s D%6.1f %-9s %-14s resid %.1e mcc %.5f Meff %+.4f l0 %6.1f off %+5.1f R2 %.3f | cl a2 %.4f far/cl %.2f vrmin %.2f | kin a1 %+.2e a2 %+.2e a3 %+.2e | stat a1 %+.2e a2 %+.2e a3 %+.2e"
              % (r["slice"], r["D"], r["kind"], json.dumps(r["a"]), r["resid"], r["mcc"], r["Meff"], r["l0"], r["CoE_axis_offset"], r["R2"], r["cl_a2"], r["ratio"], r["vrmin"],
                 r["kin"][1], r["kin"][2], r["kin"][3], r["stat"][1], r["stat"][2], r["stat"][3]), flush=True)
