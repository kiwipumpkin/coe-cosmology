"""Correct distance-dependent test: full-sample fits (Newtonian free D=93.7, Omega_s=0.049 D=49.2) held fixed; model A(l) evaluated at
each slice's own cluster distances (r_i = c z_i / 71.6 km/s/Mpc, the model <H>), normalised to the slice mean, vs the rebuilt Migkas 2020
slice profiles (cycle3/slices2.json, L-T slope fixed per slice). No refitting. Three model constructions:
 (P) point curve averaged over the slice's distance distribution (Chris's convention: model point curve vs binned data);
 (PB) the same curve binned with the 90-deg longitude top-hat used for the data;
 (C) per-cluster: each cluster gets A(l_i; r_i); in each 90-deg bin the mean over that bin's clusters, over the slice mean (same sampling as data)."""
import sys, json, numpy as np
sys.path.insert(0, "/workspace/cosmology/papers/spreadsheet_model/integrated"); sys.path.insert(0, "../cycle3")
from integ_core import A_from, r2
import migkas_A as M
T = json.load(open("af_tables.json")); PAR = T["params"]; S2 = json.load(open("../cycle3/slices2.json"))
TH1 = np.arange(0, 360, 1.0); LC = M.LC; C = M.C; z = C["z"]; r_cl = 299792.458 * z / 71.6 / 1e3
SL = {"z<0.03": (0, 0.03), "0.03-0.07": (0.03, 0.07), "z>0.07": (0.07, 1.0), "all": (0, 1.0)}
sh = lambda af, l0: np.interp((TH1 - l0) % 360, np.append(TH1, 360), np.append(af, af[0]))
def A_curve(name, Rs):                      # raw model A(l) at sampling radius Rs (log-interpolated in Rs), sine fixed
    m = PAR[name]; R = np.array(T[name]["Rs"]); AF = np.array(T[name]["af"]); o = np.argsort(R); R, AF = R[o], AF[o]
    lr = np.log(np.clip(Rs, R[0], R[-1])); af = np.array([np.interp(lr, np.log(R), AF[:, j]) for j in range(360)])
    return A_from(sh(af, m["l0"]), TH1 - m["l0"], m["amp"], m["phi"])
def binned(A1):
    k = np.ones(91) / 91; ext = np.r_[A1[-45:], A1, A1[:45]]; return np.convolve(ext, k, mode="valid")
def harm(y, lc):
    t = np.radians(lc); ok = np.isfinite(y); cols = [np.ones(ok.sum())]
    for mm in (1, 2, 3): cols += [np.cos(mm * t[ok]), np.sin(mm * t[ok])]
    c = np.linalg.lstsq(np.c_[tuple(cols)], y[ok], rcond=None)[0]
    fit = c[0] + sum(c[2 * mm - 1] * np.cos(mm * np.radians(TH1)) + c[2 * mm] * np.sin(mm * np.radians(TH1)) for mm in (1, 2, 3))
    return dict(amp1=float(np.hypot(c[1], c[2])), amp2=float(np.hypot(c[3], c[4])), phase2=float((np.degrees(np.arctan2(c[4], c[3])) / 2) % 180),
                max3=float(TH1[fit.argmax()]), min3=float(TH1[fit.argmin()]))
def ext(y, lc):
    ok = np.isfinite(y); return float(lc[ok][np.argmax(y[ok])]), float(lc[ok][np.argmin(y[ok])])
out = {}
for name in ("free", "omb"):
    out[name] = {}
    for s, (lo, hi) in SL.items():
        sel = (z >= lo) & (z < hi); d = np.array(S2[s]["A"]); ok = np.isfinite(d); lc, y = LC[ok], d[ok]
        Ais = np.array([A_curve(name, r) for r in r_cl[sel]])                         # per-cluster full curves
        P = Ais.mean(0); P = P / P.mean(); PB = binned(P)
        # per-cluster sampling
        Ai = np.array([np.interp(C["GLON"][i] % 360, TH1, A_curve(name, r_cl[i]), period=360) for i in np.where(sel)[0]])
        gl = C["GLON"][sel]; Cb = np.full(len(LC), np.nan)
        for j, l in enumerate(LC):
            w = np.abs((gl - l + 180) % 360 - 180) <= 45
            if w.sum() >= 5: Cb[j] = Ai[w].mean()
        Cb = Cb / np.nanmean(Ai)
        mP = np.interp(lc, TH1, P, period=360); mPB = np.interp(lc, TH1, PB, period=360); mC = Cb[ok]
        e = dict(n=int(sel.sum()), r_median_Gpc=float(np.median(r_cl[sel])), r_range_Gpc=[float(r_cl[sel].min()), float(r_cl[sel].max())],
                 R2_point=float(r2(y, mP)), R2_point_binned=float(r2(y, mPB)), R2_percluster=float(r2(y, mC[np.isfinite(mC)]) if np.all(np.isfinite(mC)) else np.nan),
                 model_point=dict(zip(("max", "min"), ext(P, TH1)), **harm(P, TH1)), model_binned=dict(zip(("max", "min"), ext(PB, TH1)), **harm(PB, TH1)),
                 model_percluster=dict(zip(("max", "min"), ext(Cb, LC)), **harm(Cb, LC)), data=dict(zip(("max", "min"), ext(d, LC)), **harm(d, LC)),
                 curve_point=P.tolist(), curve_binned=PB.tolist(), curve_percluster=Cb.tolist())
        out[name][s] = e
        print(name, s, "n", e["n"], "r_med %.3f" % e["r_median_Gpc"], "R2 point %.2f binned %.2f per-cluster %.2f" % (e["R2_point"], e["R2_point_binned"], e["R2_percluster"]),
              "| model max/min %d/%d ph2 %.0f (binned %d/%d ph2 %.0f; per-cl %d/%d ph2 %.0f)" % (e["model_point"]["max"], e["model_point"]["min"], e["model_point"]["phase2"],
              e["model_binned"]["max"], e["model_binned"]["min"], e["model_binned"]["phase2"], e["model_percluster"]["max"], e["model_percluster"]["min"], e["model_percluster"]["phase2"]),
              "| data max/min %d/%d (3-harm %d/%d) ph2 %.0f" % (e["data"]["max"], e["data"]["min"], e["data"]["max3"], e["data"]["min3"], e["data"]["phase2"]), flush=True)
# continuous: positions vs sampling distance
RSc = np.geomspace(0.02, 14.3, 80); cont = {}
for name in ("free", "omb"):
    rows = []
    for Rs in RSc:
        A = A_curve(name, Rs); A = A / A.mean(); B = binned(A); h = harm(A, TH1); hb = harm(B, TH1)
        rows.append(dict(Rs=float(Rs), max=float(TH1[A.argmax()]), min=float(TH1[A.argmin()]), ph2=h["phase2"], bmax=float(TH1[B.argmax()]), bmin=float(TH1[B.argmin()]), bph2=hb["phase2"]))
    cont[name] = rows
out["continuous"] = cont
json.dump(out, open("compare.json", "w"), indent=1)
