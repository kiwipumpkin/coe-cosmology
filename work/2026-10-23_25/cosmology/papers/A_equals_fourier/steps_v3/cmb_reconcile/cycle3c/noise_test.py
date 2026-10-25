"""Cycle 3c: noise-aware consistency of Chris's fixed full-sample fits (Newtonian, Omega_s = 0.049 and Omega_s free; his integrator,
tables from cycle3b/af_tables.json) evaluated at each slice's own cluster distances, vs the rebuilt Migkas 2020 A(l) profiles.
Data profile: Migkas-style sliding 90-deg longitude cones (5-deg steps), L-T normalisation per cone with slope B fixed at the slice fit.
1. Bootstrap (resample clusters, 1000x) -> covariance of the 72-point profile. chi2 in two ways:
   (F) on the 1l-3l Fourier coefficients of the profile (6 numbers; the 90-deg cones carry little above 3l), full bootstrap covariance;
   (G) on the profile itself using the leading eigenmodes of the bootstrap covariance (99% of variance), Hartlap-corrected.
   Also the effective number of independent points (tr C)^2 / tr(C^2). Model = per-cluster prediction binned in the same cones
   (each cluster gets the model A at its own longitude and distance), normalised to the slice mean. Flat (isotropic) line for comparison.
2. Mocks (1000 per slice and model): L_mock = L-T prediction x model A(cluster) x a randomly permuted real residual (keeps the real scatter,
   removes the real sky pattern); rebuild the profile with the same pipeline. Gives the noise distribution of the 2l phase, the max/min
   longitudes and chi2 (mock-calibrated p-values)."""
import sys, json, numpy as np, multiprocessing as mp
from scipy.stats import chi2 as CHI2, norm
sys.path.insert(0, "../cycle3"); sys.path.insert(0, "/workspace/cosmology/papers/spreadsheet_model/integrated")
import migkas_A as M
from integ_core import A_from
T = json.load(open("../cycle3b/af_tables.json")); PAR = T["params"]
TH1 = np.arange(0, 360, 1.0); LC = M.LC; C = M.C; z = C["z"]; r_cl = 299792.458 * z / 71.6 / 1e3
SL = {"z<0.03": (0, 0.03), "0.03-0.07": (0.03, 0.07), "z>0.07": (0.07, 1.0), "all": (0, 1.0)}
NB, NM = 1000, 1000
sh = lambda af, l0: np.interp((TH1 - l0) % 360, np.append(TH1, 360), np.append(af, af[0]))
def A_curve(name, Rs):
    m = PAR[name]; R = np.array(T[name]["Rs"]); AF = np.array(T[name]["af"]); o = np.argsort(R); R, AF = R[o], AF[o]
    lr = np.log(np.clip(Rs, R[0], R[-1])); af = np.array([np.interp(lr, np.log(R), AF[:, j]) for j in range(360)])
    return A_from(sh(af, m["l0"]), TH1 - m["l0"], m["amp"], m["phi"])
def model_cluster(name, idx):
    if name == "iso": return np.ones(len(idx))
    Ai = np.array([np.interp(C["GLON"][i] % 360, TH1, A_curve(name, r_cl[i]), period=360) for i in idx]); return Ai / Ai.mean()
def cone_mean(vals, gl):
    out = np.full(len(LC), np.nan)
    for j, l in enumerate(LC):
        w = np.abs((gl - l + 180) % 360 - 180) <= 45
        if w.sum() >= 5: out[j] = vals[w].mean()
    return out
def four(y):
    t = np.radians(LC); ok = np.isfinite(y); cols = [np.ones(ok.sum())]
    for m in (1, 2, 3): cols += [np.cos(m * t[ok]), np.sin(m * t[ok])]
    c = np.linalg.lstsq(np.c_[tuple(cols)], y[ok], rcond=None)[0]; return c
def feats(y):
    c = four(y); fit = c[0] + sum(c[2 * m - 1] * np.cos(m * np.radians(TH1)) + c[2 * m] * np.sin(m * np.radians(TH1)) for m in (1, 2, 3))
    return dict(c=c[1:] / c[0], ph2=float((np.degrees(np.arctan2(c[4], c[3])) / 2) % 180), amp2=float(np.hypot(c[3], c[4]) / c[0]),
                max=float(TH1[fit.argmax()]), min=float(TH1[fit.argmin()]))
cd = lambda a, b, P=360: np.abs((a - b + P / 2) % P - P / 2)
def boot_job(args):
    s, seed = args; lo, hi = SL[s]; sel = (z >= lo) & (z < hi); B = M.fitAB_w(sel, sel.astype(int))[1]
    rng = np.random.default_rng(seed); idx0 = np.where(sel)[0]; out = []
    for _ in range(NB // 10):
        cnt = np.bincount(rng.choice(idx0, len(idx0)), minlength=len(sel)); p = M.profile_w(cnt, B=B); out.append(p / np.nanmean(p))
    return out
def mock_job(args):
    s, name, seed = args; lo, hi = SL[s]; sel = (z >= lo) & (z < hi); idx = np.where(sel)[0]
    a, B = M.fitAB_w(sel, sel.astype(int)); la = np.log10(a); res = M.Y[idx] - la - B * M.X[idx]
    Ai = model_cluster(name, idx); rng = np.random.default_rng(seed); Y0 = M.Y.copy(); out = []
    for _ in range(NM // 10):
        Y = Y0.copy(); Y[idx] = la + B * M.X[idx] + np.log10(Ai) + rng.permutation(res); M.Y = Y
        p = M.profile_w(sel.astype(int), B=B); out.append(p / np.nanmean(p))
    M.Y = Y0; return out
if __name__ == "__main__":
    R = {}
    with mp.Pool(7) as P:
        boots = {s: sum(P.map(boot_job, [(s, 100 + k) for k in range(10)]), []) for s in SL}
        mocks = {(s, n): sum(P.map(mock_job, [(s, n, 500 + k) for k in range(10)]), []) for s in SL for n in ("omb", "free", "iso")}
    for s, (lo, hi) in SL.items():
        sel = (z >= lo) & (z < hi); idx = np.where(sel)[0]; B = M.fitAB_w(sel, sel.astype(int))[1]
        d = M.profile_w(sel.astype(int), B=B); d = d / np.nanmean(d); ok = np.isfinite(d)
        Bt = np.array(boots[s]); Bt = Bt[np.all(np.isfinite(Bt[:, ok]), axis=1)]
        Cg = np.cov(Bt[:, ok].T); err = np.sqrt(np.diag(Cg)); lam, V = np.linalg.eigh(Cg); lam, V = lam[::-1], V[:, ::-1]
        neff = float(lam.sum() ** 2 / np.sum(lam ** 2)); k = int(np.searchsorted(np.cumsum(lam) / lam.sum(), 0.99) + 1)
        hart = (len(Bt) - k - 2) / (len(Bt) - 1)
        Fb = np.array([feats(b)["c"] for b in Bt]); CF = np.cov(Fb.T); ICF = np.linalg.inv(CF) * (len(Bt) - 6 - 2) / (len(Bt) - 1)
        fd = feats(d)
        def chis(mod):
            r = (d - mod)[ok]; rg = V[:, :k].T @ r; xg = float(np.sum(rg ** 2 / lam[:k]) * hart)
            rf = fd["c"] - feats(mod)["c"]; xf = float(rf @ ICF @ rf); xd = float(np.sum(r ** 2 / err ** 2))
            return dict(chi2_F=xf, dof_F=6, p_F=float(CHI2.sf(xf, 6)), chi2_G=xg, dof_G=k, p_G=float(CHI2.sf(xg, k)), chi2_diag=xd, npts=int(ok.sum()),
                        chi2_diag_per_neff=xd / ok.sum())
        flat = np.where(ok, 1.0, np.nan); e = dict(n=int(sel.sum()), B=float(B), neff=neff, k_modes=k, data=d.tolist(), err=np.where(ok, 0, np.nan).tolist(),
                                                   data_feats={kk: (v.tolist() if hasattr(v, "tolist") else v) for kk, v in fd.items()}, flat=chis(flat))
        ee = np.full(len(LC), np.nan); ee[ok] = err; e["err"] = ee.tolist()
        e["boot_band"] = [np.nanpercentile(Bt, 16, axis=0).tolist(), np.nanpercentile(Bt, 84, axis=0).tolist()]
        e["models"] = {}
        for n in ("omb", "free", "iso"):
            Ai = model_cluster(n, idx); mod = cone_mean(Ai, C["GLON"][idx]); mod = mod / np.nanmean(mod)
            Mk = np.array(mocks[(s, n)]); Mk = Mk[np.all(np.isfinite(Mk[:, ok]), axis=1)]; mmean = np.nanmean(Mk, axis=0)
            fm = feats(mod); fmm = feats(mmean); x = chis(mod)
            if n != "iso":
                Pc = np.array([A_curve(n, r) for r in r_cl[idx]]).mean(0); Pc = np.interp(LC, TH1, Pc / Pc.mean()); xpc = chis(np.where(ok, Pc, np.nan))
            else: xpc = None
            # mock-calibrated: chi2 of each mock vs the mock mean, same covariance; and observed vs mock mean
            def xF(y, ref): r = feats(y)["c"] - feats(ref)["c"]; return float(r @ ICF @ r)
            xo = xF(d, mmean); xm = np.array([xF(y, mmean) for y in Mk])
            MF = [feats(y) for y in Mk]; dph = np.array([cd(f["ph2"], fmm["ph2"], 180) for f in MF]); dmax = np.array([cd(f["max"], fmm["max"]) for f in MF])
            dmin = np.array([cd(f["min"], fmm["min"]) for f in MF])
            o_ph, o_mx, o_mn = cd(fd["ph2"], fmm["ph2"], 180), cd(fd["max"], fmm["max"]), cd(fd["min"], fmm["min"])
            amp_m = np.array([f["amp2"] for f in MF]); p_amp = float(np.mean(amp_m >= fd["amp2"]))
            pz = lambda p: float(norm.isf(p / 2)) if p > 0 else np.inf
            e["models"][n] = dict(model=mod.tolist(), mock_mean=mmean.tolist(), mock_band=[np.nanpercentile(Mk, 16, axis=0).tolist(), np.nanpercentile(Mk, 84, axis=0).tolist()],
                                  pred=dict(ph2=fm["ph2"], max=fm["max"], min=fm["min"], amp2=fm["amp2"]), mock_mean_feats=dict(ph2=fmm["ph2"], max=fmm["max"], min=fmm["min"], amp2=fmm["amp2"]),
                                  chi=x, chi2_F_vs_mockmean=xo, p_F_mockcal=float(np.mean(xm >= xo)), sigma_F_mockcal=pz(max(np.mean(xm >= xo), 1 / len(xm))),
                                  obs_dph2=float(o_ph), p_ph2=float(np.mean(dph >= o_ph)), obs_dmax=float(o_mx), p_max=float(np.mean(dmax >= o_mx)),
                                  obs_dmin=float(o_mn), p_min=float(np.mean(dmin >= o_mn)), mock_ph2_std=float(np.sqrt(np.mean(dph ** 2))), nmock=len(Mk),
                                  p_amp2=p_amp, mock_amp2_median=float(np.median(amp_m)), mock_ph2=[f["ph2"] for f in MF], chi_pointcurve=xpc)
            m = e["models"][n]
            print(s, n, "n %d neff %.1f k %d | chi2_F %.1f/6 p %.3g (mockcal p %.3g, %.1f sig) | chi2_G %.1f/%d p %.3g | flat chi2_F %.1f p %.3g, chi2_G %.1f p %.3g" %
                  (e["n"], neff, k, x["chi2_F"], x["p_F"], m["p_F_mockcal"], m["sigma_F_mockcal"], x["chi2_G"], k, x["p_G"], e["flat"]["chi2_F"], e["flat"]["p_F"], e["flat"]["chi2_G"], e["flat"]["p_G"]), flush=True)
            print("    pred ph2 %.0f (mockmean %.0f, noise rms %.0f) obs %.0f -> d %.0f p %.3f | max pred %.0f obs %.0f d %.0f p %.3f | min pred %.0f obs %.0f d %.0f p %.3f | amp2 pred %.3f mockmed %.3f obs %.3f p %.3f | pointcurve chi2_F %s" %
                  (fm["ph2"], fmm["ph2"], m["mock_ph2_std"], fd["ph2"], o_ph, m["p_ph2"], fmm["max"], fd["max"], o_mx, m["p_max"], fmm["min"], fd["min"], o_mn, m["p_min"], fm["amp2"], m["mock_amp2_median"], fd["amp2"], p_amp, "%.1f p %.3g" % (xpc["chi2_F"], xpc["p_F"]) if xpc else "-"), flush=True)
        R[s] = e
    json.dump(R, open("noise_test.json", "w"), default=float)
