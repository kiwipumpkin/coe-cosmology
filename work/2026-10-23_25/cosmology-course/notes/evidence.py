"""Bayesian evidence (ln Z) for CoE models vs Migkas A/A_all (71 bins), dynesty nested sampling.

Usage:  python3 evidence.py run <model> <priorset>     # one nested-sampling run -> evidence_runs/<model>_<priorset>.json
        python3 evidence.py combine                     # null + dipole + table -> evidence-results.json, tbl/evidence.tex/png
models: pointmass, exp_free, exp_lock, soft_lock, soft_lock_rot
priorsets: approved  (mc 1e-4..1e-2 log-uniform, as approved by the user)
           widemc    (SENSITIVITY ONLY: mc 1e-4..30 log-uniform; everything else identical)

Likelihood: L = prod_i N(y_i | m_i, sigma_i), sigma_i=(lo+hi)/2  ->  ln L = -chi2/2 + C,
C = -sum ln(sqrt(2 pi) sigma_i) (identical for every model, so ln B is unaffected).
Invalid integrator output (non-finite a_f, vrmin<=0, lock solve failure, code guards) -> ln L = -inf (zero likelihood).
Model evaluation is direct (no grid/interpolation): one call costs 0.005-0.05 s after numba JIT.
"""
import sys, os, json, time
from pathlib import Path
import numpy as np
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
NOTES = Path("/workspace/cosmology-course/notes")
code = (NOTES / "catwise-density-levelset-migkas.py").read_text().replace('if __name__ == "__main__":\n    main()', 'pass')
ns = {}; exec(compile(code, "levelset", "exec"), ns)
from bubble_core import run as bc_run, shifted, IMIG, TH1   # noqa (sys.path set by levelset module)
Y, SIG, W, RV, D_LOCK = ns["Y"], ns["SIG"], ns["W"], ns["RV"], ns["D_LOCK"]
CHI2_NULL = float(np.sum(W * (1 - Y) ** 2))
LNC = float(-np.sum(np.log(np.sqrt(2 * np.pi) * SIG)))
LBIN = np.radians(TH1[IMIG])
LD_PT = 1e-4
BAD = 1e300

def chi2(p): return float(np.sum(W * (p - Y) ** 2))

def norm_af(af):
    if not np.all(np.isfinite(af)): return None
    return af / af.mean()

def pred_pointmass(D, mc, l0):
    rr = bc_run(D / RV, LD_PT, mc, 0.0); af = norm_af(rr["af"])
    if af is None or rr["vrmin"] <= 0: return None
    return shifted(af, l0)[IMIG]

def pred_exp_free(D, L, mc, l0):
    rr = bc_run(D / RV, L / D, mc, 0.0); af = norm_af(rr["af"])
    if af is None or rr["vrmin"] <= 0: return None
    return shifted(af, l0)[IMIG]

def pred_locked(kind, D, mc, l0, shape=None):
    o = ns["model_bins"](kind, D, mc, l0, D_LOCK, shape)
    return None if o is None else o[0][IMIG]

def rot(A, phi): return 1 - A * np.sin(LBIN + np.radians(phi))

def lu(u, lo, hi): return 10 ** (np.log10(lo) + u * (np.log10(hi) - np.log10(lo)))
def un(u, lo, hi): return lo + u * (hi - lo)

MC = dict(approved=(1e-4, 1e-2), widemc=(1e-4, 30.0))
MODELS = {
 "pointmass": dict(names=["D", "mc", "l0"], label="Point mass at CoE"),
 "exp_free": dict(names=["D", "L", "mc", "l0"], label="Exponential, L free"),
 "exp_lock": dict(names=["D", "mc", "l0"], label="Exponential, CatWISE lock"),
 "soft_lock": dict(names=["D", "mc", "l0", "alpha"], label="Soft core, CatWISE lock"),
 "soft_lock_rot": dict(names=["D", "mc", "l0", "alpha", "A", "phi"], label="Soft core lock + rotation"),
}

def make(model, ps):
    mlo, mhi = MC[ps]
    def ptf(u):
        x = [lu(u[0], 15, 150), lu(u[1 if model != "exp_free" else 2], mlo, mhi)]
        if model == "exp_free":
            return np.array([x[0], lu(u[1], 1, 100), x[1], un(u[3], 0, 360)])
        out = [x[0], x[1], un(u[2], 0, 360)]
        if model.startswith("soft"): out.append(un(u[3], 0.3, 3.0))
        if model == "soft_lock_rot": out += [un(u[4], 0, 0.1), un(u[5], 0, 360)]
        return np.array(out)
    def c2(t):
        if model == "pointmass": p = pred_pointmass(*t)
        elif model == "exp_free": p = pred_exp_free(*t)
        elif model == "exp_lock": p = pred_locked("exp", *t)
        else:
            p = pred_locked("soft", t[0], t[1], t[2], dict(alpha=float(t[3]), beta=1.0))
            if p is not None and model == "soft_lock_rot": p = p * rot(t[4], t[5])
        return BAD if p is None else chi2(p)
    def lnl(t):
        c = c2(t)
        return -1e300 if c >= BAD else -0.5 * c + LNC
    return ptf, c2, lnl

def inprior(model, ps, t):
    lo = dict(D=15, L=1, mc=MC[ps][0], l0=-1e9, alpha=0.3, A=0, phi=-1e9)
    hi = dict(D=150, L=100, mc=MC[ps][1], l0=1e9, alpha=3.0, A=0.1, phi=1e9)
    return all(lo[n] <= v <= hi[n] for n, v in zip(MODELS[model]["names"], t))

_LNL = _PTF = None
def g_lnl(t): return _LNL(t)
def g_ptf(u): return _PTF(u)

def do_run(model, ps, nlive=200, npool=1):
    import dynesty
    from dynesty import utils as dyfunc
    from scipy.optimize import minimize
    ptf, c2, lnl = make(model, ps); nd = len(MODELS[model]["names"])
    lnl(ptf(np.full(nd, 0.5)))  # JIT warm-up
    t0 = time.time()
    global _LNL, _PTF
    _LNL, _PTF = lnl, ptf
    import multiprocessing as mp
    pool = mp.get_context("fork").Pool(npool) if npool > 1 else None
    kw = dict(pool=pool, queue_size=npool) if pool else {}
    s = dynesty.NestedSampler(g_lnl, g_ptf, nd, nlive=nlive,
                              bound="multi", sample="rwalk", walks=20, rstate=np.random.default_rng(42), **kw)
    for it, res in enumerate(s.sample(dlogz=0.5)):
        if it % 500 == 0:
            print(f"it={it} ncall={s.ncall} logz={res[8]:.3f} t={time.time()-t0:.0f}s", flush=True)
    s.add_live_points()
    r = s.results
    if pool: pool.close()
    import pickle; pickle.dump(r, open(NOTES / "evidence_runs" / f"{model}_{ps}.pkl", "wb"))
    lnz = float(r.logz[-1])
    # ln Z uncertainty: std over 100 jittered (prior-volume-resampled) realisations of the run (dynesty utils.jitter_run)
    jz = [dyfunc.jitter_run(r, rstate=np.random.default_rng(100 + k)).logz[-1] for k in range(100)]
    lnzerr = float(np.std(jz)); lnzerr_dynesty = float(r.logzerr[-1])
    w = np.exp(r.logwt - r.logz[-1]); post = dyfunc.resample_equal(r.samples, w, rstate=np.random.default_rng(1))
    stats = {}
    for i, n in enumerate(MODELS[model]["names"]):
        q = np.percentile(post[:, i], [16, 50, 84]); stats[n] = dict(median=float(q[1]), lo68=float(q[0]), hi68=float(q[2]))
    ib = int(np.argmax(r.logl)); tb = r.samples[ib].copy(); cb = c2(tb)
    logv = [n in ("D", "L", "mc") for n in MODELS[model]["names"]]
    enc = lambda t: np.array([np.log10(v) if lg else v for v, lg in zip(t, logv)])
    dec = lambda u: np.array([10 ** v if lg else v for v, lg in zip(u, logv)])
    f = lambda u: c2(dec(u)) if inprior(model, ps, dec(u)) else BAD
    m = minimize(f, enc(tb), method="Nelder-Mead", options=dict(maxiter=3000, xatol=1e-6, fatol=1e-7, adaptive=True))
    if m.fun < cb: cb, tb = float(m.fun), dec(m.x)
    out = dict(model=model, label=MODELS[model]["label"], priorset=ps, mc_prior=MC[ps], params=MODELS[model]["names"],
               lnZ=lnz, lnZ_err=lnzerr, lnZ_err_dynesty_internal=lnzerr_dynesty, lnZ_minus_C=lnz - LNC, ncall=int(np.sum(r.ncall)), niter=int(r.niter), nlive=nlive,
               runtime_s=time.time() - t0, posterior=stats, maxL_chi2_in_prior=float(cb),
               maxL_params=dict(zip(MODELS[model]["names"], map(float, tb))),
               maxL_at_prior_edge={n: bool(abs(np.log10(v / MC[ps][1])) < 0.02) for n, v in zip(MODELS[model]["names"], tb) if n == "mc"},
               lnL_max_sampled=float(r.logl[ib]))
    p = NOTES / "evidence_runs" / f"{model}_{ps}.json"; p.write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1), flush=True)

def dipole_lnz():
    a = np.linspace(0, 0.3, 601); ld = np.linspace(0, 360, 721)[:-1]
    A, LDd = np.meshgrid(a, ld, indexing="ij")
    pr = 1 + A[..., None] * np.cos(LBIN[None, None, :] - np.radians(LDd)[..., None])
    c = np.sum(W * (pr - Y) ** 2, axis=-1); ll = -0.5 * c
    mx = ll.max(); Lm = np.exp(ll - mx)
    I = np.trapz(Lm.mean(axis=1), a) / 0.3   # uniform in l_d (periodic -> mean), uniform in a
    lnz = float(np.log(I) + mx + LNC)
    from scipy.optimize import minimize
    i, j = np.unravel_index(np.argmax(ll), ll.shape)
    fm = lambda x: float(np.sum(W * (1 + x[0] * np.cos(LBIN - np.radians(x[1])) - Y) ** 2))
    m = minimize(fm, [a[i], ld[j]], method="Nelder-Mead")
    pa = Lm.sum(axis=1); pa /= pa.sum(); cdf = np.cumsum(pa)
    pl = Lm.sum(axis=0); pl /= pl.sum(); cl = np.cumsum(pl)
    qa = [float(np.interp(q, cdf, a)) for q in (.16, .5, .84)]; ql = [float(np.interp(q, cl, ld)) for q in (.16, .5, .84)]
    return dict(model="dipole", label="Pure dipole 1+a cos(l-l_d)", lnZ=lnz, lnZ_err=0.0, lnZ_minus_C=lnz - LNC,
                maxL_chi2_in_prior=float(m.fun), maxL_params=dict(a=float(m.x[0]), l_d=float(m.x[1] % 360)),
                posterior=dict(a=dict(median=qa[1], lo68=qa[0], hi68=qa[2]), l_d=dict(median=ql[1], lo68=ql[0], hi68=ql[2])),
                method="2D grid quadrature (601x720), exact to <1e-3 in ln Z")

def jeffreys(lnb):
    s = "favours model" if lnb > 0 else "favours comparison"; x = abs(lnb)
    g = "inconclusive" if x < 1 else "weak" if x < 2.5 else "moderate" if x < 5 else "strong"
    return f"{g} ({s})" if x >= 1 else "inconclusive"

def combine():
    null = dict(model="flat", label="Flat A/A_all = 1", lnZ=-0.5 * CHI2_NULL + LNC, lnZ_err=0.0, lnZ_minus_C=-0.5 * CHI2_NULL,
                maxL_chi2_in_prior=CHI2_NULL)
    dip = dipole_lnz()
    runs = {}
    for f in sorted((NOTES / "evidence_runs").glob("*.json")):
        d = json.loads(f.read_text()); runs[f"{d['model']}_{d['priorset']}"] = d
    for d in runs.values():
        for ref, key in ((null, "flat"), (dip, "dipole")):
            b = d["lnZ"] - ref["lnZ"]
            d[f"lnB_vs_{key}"] = b; d[f"lnB_vs_{key}_err"] = d["lnZ_err"]; d[f"jeffreys_vs_{key}"] = jeffreys(b)
    dip["lnB_vs_flat"] = dip["lnZ"] - null["lnZ"]; dip["jeffreys_vs_flat"] = jeffreys(dip["lnB_vs_flat"])
    res = dict(date="2026-10-07 (America/Chicago)", data="Migkas Fig.6 A/A_all, N=71, sigma=(lo+hi)/2, WLS chi2", chi2_null=CHI2_NULL,
               R_Gpc=RV, lnL_const_C=LNC, sampler="dynesty NestedSampler nlive=300 bound=multi sample=unif dlogz=0.1 seed=42",
               priors=dict(D="15-150 Gpc log-uniform", mc_approved="1e-4..1e-2 log-uniform", mc_widemc_SENSITIVITY="1e-4..30 log-uniform",
                           L="1-100 Gpc log-uniform (exp_free only)", alpha="0.3-3 uniform", l0="0-360 uniform", A="0-0.1 uniform",
                           phi="0-360 uniform", dipole="a 0-0.3 uniform, l_d 0-360 uniform"),
               jeffreys_scale="|lnB|<1 inconclusive, 1-2.5 weak, 2.5-5 moderate, >5 strong",
               approximation="none: model evaluated directly per likelihood call (no grid/interpolation)",
               flat=null, dipole=dip, runs=runs,
               prior_chi2_crosscheck=dict(pointmass=20.84, exp_free=8.48, exp_lock=10.57, soft_lock=10.43, soft_lock_rot=5.21))
    (NOTES / "evidence-results.json").write_text(json.dumps(res, indent=1))
    make_table(res)

def fmt_par(d):
    s = []
    for n, v in d["posterior"].items():
        m, lo, hi = v["median"], v["lo68"], v["hi68"]
        if n == "mc": s.append(f"$m_c={m:.2g}^{{+{hi-m:.1g}}}_{{-{m-lo:.1g}}}$")
        else:
            sym = dict(D="D", L="L", l0=r"\ell_0", alpha=r"\alpha", A="A", phi=r"\varphi", a="a", l_d=r"\ell_d")[n]
            dg = 3 if n in ("A", "a") else 2 if n == "alpha" else 0
            s.append(f"${sym}={m:.{dg}f}^{{+{hi-m:.{dg}f}}}_{{-{m-lo:.{dg}f}}}$")
    return ", ".join(s)

def make_table(res):
    import subprocess
    order = ["pointmass", "exp_free", "exp_lock", "soft_lock", "soft_lock_rot"]
    rows = []
    def row(lab, d, k):
        return (f"{lab} & {k} & {d['lnZ_minus_C']:.2f}$\\pm${d.get('lnZ_err',0):.2f} & {d.get('lnB_vs_flat',0):+.2f} & "
                f"{d.get('lnB_vs_dipole', d['lnZ']-res['dipole']['lnZ']):+.2f} & {d['maxL_chi2_in_prior']:.2f}\\\\")
    rows.append(row("Flat $A/A_{\\rm all}=1$", res["flat"], 0))
    rows.append(row("Pure dipole", res["dipole"], 2))
    blocks = {}
    for ps, title in (("approved", "Approved priors ($m_c$ $10^{-4}$--$10^{-2}$)"), ("widemc", "Sensitivity: $m_c$ $10^{-4}$--$30$")):
        rr = [res["runs"][f"{m}_{ps}"] for m in order if f"{m}_{ps}" in res["runs"]]
        if not rr: continue
        blocks[ps] = (title, [row(d["label"], d, len(d["params"])) for d in rr], rr)
    tex = [r"\documentclass[border=12pt,varwidth=24cm]{standalone}", r"\usepackage{booktabs,amsmath,array}",
           r"\renewcommand{\arraystretch}{1.35}", r"\begin{document}\Large",
           r"\textbf{Bayesian evidence, Migkas $A/A_{\mathrm{all}}$ ($N=71$, $\chi^2_{\rm null}=45.14$)}\\[6pt]",
           r"\begin{tabular}{p{8.5cm}cccccc}".replace("cccccc", "ccccc"), r"\toprule",
           r"Model & $k$ & $\ln Z - C$ & $\ln B_{\rm flat}$ & $\ln B_{\rm dip}$ & $\chi^2_{\min}$\\", r"\midrule"] + rows
    for ps, (title, rws, rr) in blocks.items():
        tex += [r"\midrule", rf"\multicolumn{{6}}{{l}}{{\textit{{{title}}}}}\\", r"\midrule"] + rws
    tex += [r"\bottomrule", r"\end{tabular}\\[8pt]", r"\normalsize"]
    tex.append(r"\begin{minipage}{23cm}\textbf{Posterior medians, 68\%:}\\")
    for ps, (title, rws, rr) in blocks.items():
        tex.append(rf"\textit{{{title}}}\\")
        for d in rr: tex.append(rf"{d['label']}: {fmt_par(d)}; {d['jeffreys_vs_flat']} vs flat\\")
    tex.append(rf"Dipole: {fmt_par(res['dipole'])}\\[4pt]")
    tex.append(r"$C=-\sum\ln(\sqrt{2\pi}\sigma_i)$ common to all models. Jeffreys: $|\ln B|<1$ inconclusive, 1--2.5 weak, 2.5--5 moderate, $>5$ strong. "
               r"Priors: $D$ 15--150 Gpc, $L$ 1--100 Gpc, $m_c$ log-uniform; $\alpha$ 0.3--3, $\ell_0,\varphi$ 0--360$^\circ$, $A$ 0--0.1, $a$ 0--0.3 uniform. "
               r"$\chi^2_{\min}$ is the best fit inside the prior. dynesty, 300 live points.\end{minipage}")
    tex.append(r"\end{document}")
    tdir = NOTES / "tbl"; (tdir / "evidence.tex").write_text("\n".join(tex))
    subprocess.run(["pdflatex", "-interaction=nonstopmode", "evidence.tex"], cwd=tdir, check=True, capture_output=True)
    subprocess.run(["pdftoppm", "-png", "-r", "150", "-singlefile", "evidence.pdf", "evidence"], cwd=tdir, check=True)

if __name__ == "__main__":
    import multiprocessing as _mp; _mp.set_start_method("fork", force=True)
    if sys.argv[1] == "run": do_run(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 200, int(sys.argv[5]) if len(sys.argv) > 5 else 1)
    else: combine()
