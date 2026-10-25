"""Part 1: redshift-split Migkas test with the real per-cluster data (Migkas+2020, VizieR J/A+A/636/A15, 313 clusters).
Per slice: rebuild A(l)/A_all exactly as Fig. 6 (Delta l = 90 deg bins, all b, 5-deg steps), bootstrap errors, Fourier amplitudes.
Fit (a) Chris's far core (Newtonian integrator unchanged, Omega_s=0.049 headline D=49.17, mcc=1.742e-3, kin, kout) sampled at the
slice's median distance, with (a1) headline l0/Amp/phi fixed and (a2) l0 + sine refit; (b) near-core variant (exact LTB, growing mode,
core + uniform dust) for D = 0.2/0.3/0.5 Gpc, mcc = 0.05/0.08/0.15, sampled at the slice's median z, l0 + sine refit.
R2 as in Chris's convention (model point curve vs binned data) and also with the model binned the same way (90-deg top-hat)."""
import sys, os, json, numpy as np, multiprocessing as mp
import migkas_A as M
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/workspace/cosmology/papers/spreadsheet_model/integrated")
from integ_core import r2, A_from, best_sine
TH1 = np.arange(0, 360, 1.0)
def shifted(af, l0): return np.interp((TH1 - l0) % 360, np.append(TH1, 360), np.append(af, af[0]))
def binned(A1):   # 90-deg top-hat in longitude (the Fig. 6 bins), on the 1-deg grid (approximation: latitude-integrated)
    k = np.ones(91) / 91; ext = np.r_[A1[-45:], A1, A1[:45]]; return np.convolve(ext, k, mode="valid")
def fit_curve(af1, lc, y, sine=True, fixed=None, bin_model=False):
    """af1: model on 1-deg grid with psi=0 at index 0 meaning 'direction l0' convention of Chris (A(l) = af(l - l0))."""
    best = (-np.inf,)
    l0s = [fixed[0]] if fixed else np.arange(0, 360, 1.0)
    for l0 in l0s:
        A = shifted(af1, l0); A = A / A.mean()
        if fixed: am, ph = fixed[1], fixed[2]
        elif sine: am, ph = best_sine(np.interp(lc, TH1, A, period=360), lc - l0, y)
        else: am, ph = 0.0, 0.0
        full = A_from(A, TH1 - l0, am, ph)
        if bin_model: full = binned(full)
        m = np.interp(lc, TH1, full, period=360); v = r2(y, m)
        if v > best[0]: best = (v, l0, am, ph, full)
    return dict(R2=best[0], l0=best[1], amp=best[2], phi=best[3]), best[4]
# ---------- far core: Chris's Newtonian integrator ----------
BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"
def far_af(Rs):
    sys.path.insert(0, BUB); cwd = os.getcwd(); os.chdir(BUB)
    import newt_omb as N
    r = N.af_of(N.geo(49.1746, 0.0017418, 1.0528, 1.0815, Rs)); os.chdir(cwd); return r[0]
HEAD = (13.2242, -0.0621, 0.9307)          # omb headline l0, Amp, phi
# ---------- near core: exact LTB ----------
def near_af(args):
    mcc, D, z = args
    sys.path.insert(0, HERE); from ltb_common import make_p, setup, C
    from numpy.polynomial.legendre import legfit, legval
    Mc = mcc * 4 * np.pi / 3 * 9.2042e-27 * (D * 3.0857e25)**3
    p = make_p(Mc, D, 0.049, 0.03 * D, 5000.0, cavity=False); p[11] = 0.01 * D; t0 = setup(p, D)
    nb_ = C.nearby_fast(p, t0, D, (z,), 24); DL = nb_["DL"][:, 0]
    cf = legfit(nb_["mu"], DL, len(DL) - 1)
    psi = TH1 - 180.0                          # af(l - l0) with psi = l - l0 - 180 (CoE at l0 + 180, as rm_core.lfits)
    return legval(np.cos(np.radians(psi)), cf) ** -1.0
SL = {"z<0.03": (0, 0.03), "0.03-0.07": (0.03, 0.07), "z>0.07": (0.07, 1.0), "all": (0, 1.0)}
if __name__ == "__main__":
    z = M.C["z"]; out = {"source": "Migkas+2020 VizieR J/A+A/636/A15/tablec1 (313 clusters)", "slices": {}}
    near_grid = [(m, D) for D in (0.2, 0.3, 0.5) for m in (0.05, 0.08, 0.15)]
    jobs = []
    for s, (lo, hi) in SL.items():
        sel = (z >= lo) & (z < hi); zm = float(np.median(z[sel]))
        jobs += [(m, D, zm) for m, D in near_grid]
    with mp.Pool(6) as P: naf = dict(zip(jobs, P.map(near_af, jobs)))
    for s, (lo, hi) in SL.items():
        sel = (z >= lo) & (z < hi); zm = float(np.median(z[sel])); Rs = 299792.458 * zm / 71.6 / 1e3   # Gpc (low-z Hubble law)
        prof = M.profile(sel); err, herr = M.bootstrap(sel, n=200)
        ok = np.isfinite(prof); lc, y = M.LC[ok], prof[ok]
        H = M.harm(prof); H.update({f"{k}_err": v for k, v in herr.items()})
        e = dict(n=int(sel.sum()), z_median=zm, Rs_Gpc=Rs, A=prof.tolist(), A_err=err.tolist(), harmonics=H, n_bins_valid=int(ok.sum()))
        af_far = far_af(Rs); af_far14 = far_af(14.30048287260589)
        e["far_core_at_slice_depth_headline_fixed"], _ = fit_curve(af_far, lc, y, fixed=HEAD)
        e["far_core_at_slice_depth_refit"], cfar = fit_curve(af_far, lc, y)
        e["far_core_at_slice_depth_refit_binned"], _ = fit_curve(af_far, lc, y, bin_model=True)
        e["far_core_headline_R14.3_fixed"], _ = fit_curve(af_far14, lc, y, fixed=HEAD)
        e["far_core_headline_R14.3_refit"], _ = fit_curve(af_far14, lc, y)
        e["far_core_model_harm"] = M.harm(cfar[::5] if False else np.interp(M.LC, TH1, cfar, period=360))
        e["near"] = []
        for m, D in near_grid:
            af = naf[(m, D, zm)]; f1, cn = fit_curve(af, lc, y); f2, _ = fit_curve(af, lc, y, bin_model=True)
            e["near"].append(dict(mcc=m, D=D, R2=f1["R2"], R2_binned=f2["R2"], l0=f1["l0"], amp=f1["amp"], model_amp2=M.harm(np.interp(M.LC, TH1, cn, period=360))["amp2"]))
        # sine only (dipole) reference and free 3-harmonic reference
        e["dipole_only_R2"], _ = fit_curve(np.ones(360), lc, y)
        out["slices"][s] = e
        print(s, "n", e["n"], "zmed %.3f" % zm, "amp1 %.3f+-%.3f amp2 %.3f+-%.3f amp3 %.3f+-%.3f phase2 %.0f" % (H["amp1"], herr["amp1"], H["amp2"], herr["amp2"], H["amp3"], herr["amp3"], H["phase2_deg"]),
              "| far fixed %.2f refit %.2f (binned %.2f) R14 fixed %.2f" % (e["far_core_at_slice_depth_headline_fixed"]["R2"], e["far_core_at_slice_depth_refit"]["R2"], e["far_core_at_slice_depth_refit_binned"]["R2"], e["far_core_headline_R14.3_fixed"]["R2"]),
              "| near best %.2f" % max(n["R2"] for n in e["near"]), "| dipole %.2f" % e["dipole_only_R2"]["R2"], flush=True)
    json.dump(out, open(os.path.join(HERE, "slices.json"), "w"), indent=1, default=float)
