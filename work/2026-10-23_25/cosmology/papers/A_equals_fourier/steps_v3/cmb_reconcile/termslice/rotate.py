"""Orange (0.03-0.07) as the black (all) / red (z>0.07) solved model rotated by one constant longitude offset dl (no other refit).
Uncertainty of dl from 500 cluster bootstraps of the orange profile (L-T slope fixed at the slice fit), refitting dl each time."""
import json, sys, numpy as np, multiprocessing as mp
sys.path.insert(0, "../cycle3"); import migkas_A as M
R = json.load(open("solve.json")); S2 = json.load(open("../cycle3/slices2.json")); DQ = json.load(open("../migkas_zslices_dipquad.json"))
TH1 = np.arange(0, 360, 1.0); LC = M.LC; z = M.C["z"]; sel = (z >= 0.03) & (z < 0.07)
pick = lambda s: [r for r in R if r["slice"] == s and abs(r["D"] - 49.1746) < 1e-3 and r["mode"] == "median" and not r["cone"] and r["l0_branch"] == 0][0]
MOD = {s: np.array(pick(s)["curve"]) for s in ("all", "z>0.07")}
cone = lambda A: np.convolve(np.r_[A[-45:], A, A[:45]], np.ones(91) / 91, mode="valid")
r2 = lambda y, m: 1 - np.sum((y - m)**2) / np.sum((y - y.mean())**2)
def fit_dl(prof, A, binned=False):
    ok = np.isfinite(prof); y = prof[ok] / np.nanmean(prof); A2 = cone(A) if binned else A; best = (-np.inf, 0)
    for dl in np.arange(-180, 180, 0.5):
        m = np.interp((LC[ok] - dl) % 360, TH1, A2, period=360); v = r2(y, m)
        if v > best[0]: best = (v, dl)
    return best
B = M.fitAB_w(sel, sel.astype(int))[1]
def boot(seed):
    rng = np.random.default_rng(seed); idx0 = np.where(sel)[0]; out = []
    for _ in range(50):
        cnt = np.bincount(rng.choice(idx0, len(idx0)), minlength=len(sel)); p = M.profile_w(cnt, B=B)
        out.append([fit_dl(p, MOD[s], b)[1] for s in ("all", "z>0.07") for b in (False, True)])
    return out
if __name__ == "__main__":
    prof = np.array(S2["0.03-0.07"]["A"], float); res = {}
    for s in ("all", "z>0.07"):
        for b in (False, True):
            v, dl = fit_dl(prof, MOD[s], b); v0 = fit_dl.__wrapped__ if False else None
            ok = np.isfinite(prof); y = prof[ok] / np.nanmean(prof); A2 = cone(MOD[s]) if b else MOD[s]
            r0 = r2(y, np.interp(LC[ok] % 360, TH1, A2, period=360))
            res["%s|%s" % (s, "cone" if b else "point")] = dict(dl=dl, R2=v, R2_unrotated=r0)
    with mp.Pool(7) as P: Bs = np.array(sum(P.map(boot, range(10)), []))
    keys = ["all|point", "all|cone", "z>0.07|point", "z>0.07|cone"]
    for j, k in enumerate(keys):
        d = (Bs[:, j] - res[k]["dl"] + 180) % 360 - 180 + res[k]["dl"]
        res[k].update(dl_boot_median=float(np.median(d)), dl_err68=float(0.5 * (np.percentile(d, 84) - np.percentile(d, 16))), dl_err_std=float(np.std(d)),
                      frac_boot_dl_ge0=float(np.mean(d >= 0)))
        print(k, {q: round(v, 3) for q, v in res[k].items()})
    # local quadrupole needed to rotate the black/red quad into orange's (vector sum in 2l space)
    for s, (q0, p0) in {"all": (DQ["all"]["a2"], DQ["all"]["p2"]), "z>0.07": (DQ["z>0.07"]["a2"], DQ["z>0.07"]["p2"])}.items():
        zo = DQ["0.03-0.07"]["a2"] * np.exp(2j * np.radians(DQ["0.03-0.07"]["p2"])); z0 = q0 * np.exp(2j * np.radians(p0)); zl = zo - z0
        res["local_quad_needed_vs_" + s] = dict(amp=float(abs(zl)), axis_l=float((np.degrees(np.angle(zl)) / 2) % 180))
        print("local quad needed vs", s, res["local_quad_needed_vs_" + s])
    json.dump(res, open("rotate.json", "w"), indent=1)
