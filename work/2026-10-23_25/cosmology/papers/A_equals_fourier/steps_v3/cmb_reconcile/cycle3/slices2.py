"""Part 1b: robust version of the redshift split. Changes vs slices.py:
- L-T slope B fixed per slice at that slice's all-sky best fit inside the 90-deg bins (free B in 5-cluster bins made the z<0.03
  bootstrap blow up); also a run with B fixed at the full-sample value.
- Robust bootstrap errors (16-84 percentile half-width, 300 resamples).
- Isotropy test: 500 shuffles of sky positions among the slice's clusters -> p-value of amp1, amp2 (probability that an isotropic
  sky gives a 2l amplitude at least as large).
- Coherent 2l component: amp2 * cos(2 (phase2 - phase2_all)), i.e. the part of each slice's 2l that lines up with the full-sample pattern.
- Common-orientation model tests: l0 (and sine) taken from the full-sample fit, only applied to each slice (a single physical CoE
  direction must work at all depths). Near core D = 0.2/0.3/0.5 Gpc (mcc 0.08), far core headline (Newtonian omb, sampled at slice depth)."""
import os, json, numpy as np, multiprocessing as mp
import migkas_A as M
from slices import fit_curve, far_af, near_af, SL, TH1
HERE = os.path.dirname(os.path.abspath(__file__))
z = M.C["z"]
def prof_fixedB(sel, cnt=None, B=None):
    if cnt is None: cnt = sel.astype(int)
    if B is None: B = M.fitAB_w(cnt > 0, cnt)[1]
    return M.profile_w(cnt, B=B), B
def boot(sel, B, n=300, seed=1):
    rng = np.random.default_rng(seed); idx0 = np.where(sel)[0]; H = []; P = []
    for _ in range(n):
        cnt = np.bincount(rng.choice(idx0, len(idx0)), minlength=len(sel)); p = M.profile_w(cnt, B=B); P.append(p); H.append(M.harm(p))
    q = lambda a: float(0.5 * (np.nanpercentile(a, 84) - np.nanpercentile(a, 16)))
    return {k: q([h[k] for h in H]) for k in ("amp1", "amp2", "amp3")}, np.array([0.5 * (np.nanpercentile(c, 84) - np.nanpercentile(c, 16)) for c in np.array(P).T])
def perm(sel, B, n=500, seed=2):
    rng = np.random.default_rng(seed); idx = np.where(sel)[0]; g0 = M.C["GLON"].copy(); A1 = []; A2 = []
    for _ in range(n):
        g = g0.copy(); g[idx] = g0[rng.permutation(idx)]; M.C["GLON"] = g
        h = M.harm(M.profile_w(sel.astype(int), B=B)); A1.append(h["amp1"]); A2.append(h["amp2"])
    M.C["GLON"] = g0; return np.array(A1), np.array(A2)
if __name__ == "__main__":
    out = {}; Ball = M.fitAB(np.ones(len(z), bool))[1]
    zms = {s: float(np.median(z[(z >= lo) & (z < hi)])) for s, (lo, hi) in SL.items()}
    Ds = (0.2, 0.3, 0.5); jobs = [(0.08, D, zms[s]) for s in SL for D in Ds]
    with mp.Pool(6) as P: naf = dict(zip(jobs, P.map(near_af, jobs)))
    far = {s: far_af(299792.458 * zms[s] / 71.6 / 1e3) for s in SL}
    # common orientation from the full sample
    selA = np.ones(len(z), bool); pA, _ = prof_fixedB(selA); okA = np.isfinite(pA)
    commonF, _ = fit_curve(far["all"], M.LC[okA], pA[okA]); commonN = {D: fit_curve(naf[(0.08, D, zms["all"])], M.LC[okA], pA[okA])[0] for D in Ds}
    hA = M.harm(pA)
    for s, (lo, hi) in SL.items():
        sel = (z >= lo) & (z < hi); p, B = prof_fixedB(sel); pB, _ = prof_fixedB(sel, B=Ball)
        h = M.harm(p); hB = M.harm(pB); he, perr = boot(sel, B); a1, a2 = perm(sel, B)
        ok = np.isfinite(p); lc, y = M.LC[ok], p[ok]
        e = dict(n=int(sel.sum()), z_median=zms[s], B_slice=float(B), A=p.tolist(), A_err=perr.tolist(), harm=h, harm_err=he, harm_Bfull=hB,
                 p_iso_amp1=float(np.mean(a1 >= h["amp1"])), p_iso_amp2=float(np.mean(a2 >= h["amp2"])), iso_amp2_median=float(np.median(a2)),
                 iso_amp2_p95=float(np.percentile(a2, 95)), coherent2=float(h["amp2"] * np.cos(np.radians(2 * (h["phase2_deg"] - hA["phase2_deg"])))))
        e["far_refit"], _ = fit_curve(far[s], lc, y)
        e["far_common_l0_sine"], _ = fit_curve(far[s], lc, y, fixed=(commonF["l0"], commonF["amp"], commonF["phi"]))
        e["near"] = {}
        for D in Ds:
            af = naf[(0.08, D, zms[s])]; cn = commonN[D]
            r_ref, curve = fit_curve(af, lc, y); r_com, _ = fit_curve(af, lc, y, fixed=(cn["l0"], cn["amp"], cn["phi"]))
            e["near"][str(D)] = dict(refit=r_ref, common=r_com, model_amp2_nosine=M.harm(np.interp(M.LC, TH1, af / af.mean(), period=360))["amp2"])
        e["far_model_amp2_nosine"] = M.harm(np.interp(M.LC, TH1, far[s] / far[s].mean(), period=360))["amp2"]
        e["dipole_R2"], _ = fit_curve(np.ones(360), lc, y)
        out[s] = e
        print(s, "n", e["n"], "B %.2f" % B, "amp1 %.3f+-%.3f (p %.3f) amp2 %.3f+-%.3f (p %.3f, iso95 %.3f) amp3 %.3f+-%.3f ph2 %.0f coh2 %.3f | Bfull amp2 %.3f" %
              (h["amp1"], he["amp1"], e["p_iso_amp1"], h["amp2"], he["amp2"], e["p_iso_amp2"], e["iso_amp2_p95"], h["amp3"], he["amp3"], h["phase2_deg"], e["coherent2"], hB["amp2"]), flush=True)
        print("   far refit %.2f common %.2f modelamp2 %.3f" % (e["far_refit"]["R2"], e["far_common_l0_sine"]["R2"], e["far_model_amp2_nosine"]),
              " ".join("| near D%.1f refit %.2f common %.2f mamp2 %.4f" % (D, e["near"][str(D)]["refit"]["R2"], e["near"][str(D)]["common"]["R2"], e["near"][str(D)]["model_amp2_nosine"]) for D in Ds),
              "| dipole %.2f" % e["dipole_R2"]["R2"], flush=True)
    out["common"] = dict(far=commonF, near={str(D): commonN[D] for D in Ds})
    json.dump(out, open(os.path.join(HERE, "slices2.json"), "w"), indent=1, default=float)
