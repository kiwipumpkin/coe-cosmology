"""Per-slice term-by-term solve (Chris's geometry, Newtonian integrator, Omega_s = 0.049 shell with the headline rule
R_in = D - 1.028 R, R_out = D + 1.475 R, R = 14.3 Gpc). Sampling radius Rs = slice median cz/71.6 (and a per-cluster distance-distribution
variant). Unknowns: Mc/Mcrit(D), l0, sine Amp, phi (D scanned: only weakly constrained at Rs << D). Equations: the slice's dipole and
quadrupole Fourier terms (C1, S1, C2, S2) from migkas_zslices_dipquad.json. Model = a_f(l - l0) [1 - Amp sin(l - l0 + phi)], mean-normalised,
compared linearly (point model vs the data's fitted terms, as the v3 term-by-term sheet). Variant 'cone': model terms multiplied by the
90-deg top-hat factors (1l 0.900, 2l 0.637) so they are compared like-for-like with the cone-binned profile."""
import json, sys, numpy as np, multiprocessing as mp
from scipy.optimize import least_squares
import tcore as T
sys.path.insert(0, "../cycle3"); import migkas_A as M
DQ = json.load(open("../migkas_zslices_dipquad.json")); S2 = json.load(open("../cycle3/slices2.json"))
z = M.C["z"]; rcl = 299792.458 * z / 71.6 / 1e3
SL = {"0.03-0.07": (0.03, 0.07), "z>0.07": (0.07, 1.0), "all": (0, 1.0)}
HM = 71.6; F1, F2 = np.sin(np.pi / 4) / (np.pi / 4), np.sin(np.pi / 2) / (np.pi / 2)
def af_curve(D, mcc, Rs_list):
    afs = []; Ht = []
    for Rs in Rs_list:
        r = T.run(D, mcc, Rs); afs.append(r["af"]); Ht.append(r["Ht"])
    return np.mean(afs, axis=0), float(np.mean(Ht))
def model(D, mcc, l0, amp, phi, Rs_list):
    af, Ht = af_curve(D, mcc, Rs_list); full = T.A_from(T.shifted(af, l0), T.TH1 - l0, amp, phi); return full, Ht
def solve(s, D, mode="median", cone=False, l0_branch=0):
    d = DQ[s]; tg = T.target(d["a1"], d["p1"], d["a2"], d["p2"]); lo, hi = SL[s]; sel = (z >= lo) & (z < hi)
    Rs_list = [float(np.median(rcl[sel]))] if mode == "median" else list(np.quantile(rcl[sel], [0.07, 0.21, 0.36, 0.5, 0.64, 0.79, 0.93]))
    f1, f2 = (F1, F2) if cone else (1.0, 1.0)
    def res(x):
        lm, l0, amp, phi = x; A, _ = model(D, 10**lm, l0, amp, phi, Rs_list); A = A / A.mean()
        c1, s1 = T.harm(A, 1); c2, s2 = T.harm(A, 2)
        return np.array([f1 * c1 - tg[1][0], f1 * s1 - tg[1][1], f2 * c2 - tg[2][0], f2 * s2 - tg[2][1]])
    l0g = d["p2"] + 180 * l0_branch; best = None
    for pg in np.linspace(-np.pi, np.pi, 4, endpoint=False):
        x0 = [np.log10(d["a2"] / f2 / 43.0), l0g, -d["a1"] / f1, pg]
        r = least_squares(res, x0, x_scale=[0.1, 3, 0.02, 0.3], diff_step=1e-4)
        if best is None or r.cost < best.cost: best = r
    lm, l0, amp, phi = best.x; mcc = 10**lm
    if amp > 0: amp, phi = -amp, phi + np.pi
    A, Ht = model(D, mcc, l0, amp, phi, Rs_list); raw_mean = A.mean(); An = A / raw_mean
    k = HM / raw_mean; HE = Ht * k; Mc = mcc * 4 * np.pi / 3 * T.rhoc(HE) * (D * T.GPCm)**3; Mc70 = mcc * 4 * np.pi / 3 * T.rhoc(70.0) * (D * T.GPCm)**3
    # geometry-only terms and sine-supplied dipole
    af, _ = af_curve(D, mcc, Rs_list); g = T.dq(T.shifted(af, l0)); m = T.dq(An)
    prof = np.array(S2[s]["A"], float); ok = np.isfinite(prof); y = prof[ok] / np.nanmean(prof)
    mod_at = np.interp(M.LC[ok], T.TH1, An, period=360)
    mod_cone = np.convolve(np.r_[An[-45:], An, An[:45]], np.ones(91) / 91, mode="valid"); mod_cone_at = np.interp(M.LC[ok], T.TH1, mod_cone, period=360)
    r2 = lambda yy, mm: 1 - np.sum((yy - mm)**2) / np.sum((yy - yy.mean())**2)
    Rs_med = float(np.median(rcl[sel])); v_tilt = abs(amp) * HM * Rs_med * 1e3
    return dict(slice=s, D=D, mode=mode, cone=cone, l0_branch=l0_branch, Rs=Rs_list, resid=float(np.sqrt(2 * best.cost)), mcc=mcc, l0=float(l0 % 360), amp=float(amp), phi=float(phi % (2 * np.pi)),
                model_terms=m, geom_terms=g, R2_point=float(r2(y, mod_at)), R2_conebinned=float(r2(y, mod_cone_at)), R2_dipquad_data=d["r2"],
                Ht=Ht, k=k, H_E=HE, Mc_kg=Mc, Mc_kg_H70=Mc70, v_tilt_kms=v_tilt, sine_dipole_peak_l=float((l0 - np.degrees(phi) + 90) % 360), curve=An.tolist())
JOBS = [(s, D, mode, cone, br) for s in SL for D in (20.0, 30.0, 49.1746, 94.0, 200.0) for mode in ("median",) for cone in (False,) for br in (0, 1)]
JOBS += [(s, 49.1746, "dist", False, 0) for s in SL] + [(s, 49.1746, "median", True, 0) for s in SL]
def job(a): return solve(*a)
if __name__ == "__main__":
    with mp.Pool(7) as P: R = P.map(job, JOBS, chunksize=1)
    json.dump(R, open("solve.json", "w"), default=float)
    for r in R:
        print(r["slice"], "D %.1f %s cone=%s br=%d" % (r["D"], r["mode"], r["cone"], r["l0_branch"]), "resid %.1e" % r["resid"], "mcc %.5f l0 %.1f amp %.4f phi %.3f" % (r["mcc"], r["l0"], r["amp"], r["phi"]),
              "| R2 point %.3f cone %.3f (dipquad %.2f)" % (r["R2_point"], r["R2_conebinned"], r["R2_dipquad_data"]), "| Mc %.3g kg (H70 %.3g) HE %.2f" % (r["Mc_kg"], r["Mc_kg_H70"], r["H_E"]),
              "| geom a1 %.4f a2 %.4f a3 %.4f | v_tilt %.0f km/s, dipole peak l %.0f" % (r["geom_terms"]["a1"], r["geom_terms"]["a2"], r["geom_terms"]["a3"], r["v_tilt_kms"], r["sine_dipole_peak_l"]), flush=True)
