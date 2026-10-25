"""A: far core, exact GR (LTB dust + central point-like core, Lambda=0), CMB l>=2 residual for an off-centre observer at D.
Chris's LTB code unchanged (ltb_fill + cmb_redshift rays/analyze). Construction = the rel-fit convention: single bang t_B=0,
E(r) from R(r,t0)=r (areal label today), t0 from local mean H=70 at Earth. Core-only: the cavity is filled (uniform dust + core).
The dust is extended to label 5000 Gpc so the last-scattering surface lies in dust (in open Omega_s=0.049 GR the LSS sits at
areal-today label ~180 Gpc, proper radial ~19 Gpc - outside the fitted R_out).
Outputs per case: exact CMB multipoles (fixed proper-time emission 'a' with dipole boosted away, and fixed local density 'b'),
and the l=2 of ln(1+z) for emission at a ladder of earlier times (tests Phi ~ 1/a vs Phi ~ const)."""
import sys, json, time, numpy as np, multiprocessing as mp
from ltb_common import *
TL = (0.3, 0.1, 0.03, 0.01, 3e-3, 1e-3, 3e-4, 1e-4)
def run(case):
    name, f, cav, mscale = case
    p = make_p(f["Mc"] * mscale, f["D"], f["Os"], f["Rin"] if cav else f["Rin"], 5000.0, cavity=cav)
    t0 = setup(p, f["D"]); t = time.time()
    if t0 is None: return dict(name=name, fail="t0")
    try:
        c = analyze.cmb(p, t0, f["D"], N=24)
    except Exception as e:
        return dict(name=name, fail=str(e)[:200])
    mu, w = leggauss(24)
    sols = [rays.trace(p, t0, f["D"], ps, t0 * 1e-4 * 0.99) for ps in np.arccos(mu)]
    ladder = {}
    for tf in TL:
        z = np.array([rays.z_at_t(s, t0, t0 * tf) for s in sols]); L = rays.legendre_coeffs(mu, w, -np.log1p(z) + np.mean(np.log1p(z)), 4)
        ladder[str(tf)] = dict(zbar=float(0.5 * np.sum(w * z)), a1=L[1], a2=L[2], a3=L[3])
    Rt, Rr = LT.local_rates(f["D"], t0, p)
    out = dict(name=name, cavity=cav, mscale=mscale, mcc=f["mcc"] * mscale, D=f["D"], t0=t0, H_perp=Rt * C.K, H_r=Rr * C.K, null_err=c["null_err"],
               t_rec_Myr=c["t_rec_Myr"], c_a=c["c_a"].tolist(), cb_a=c["cb_a"].tolist(), beta_a=c["beta_a"],
               c_b=np.asarray(c["c_b"]).tolist(), cb_b=np.asarray(c["cb_b"]).tolist(), ladder=ladder, sec=time.time() - t)
    print(name, "a2(boosted, fixed t) =", out["cb_a"][2], "a2(fixed rho) =", out["cb_b"][2], "sec", out["sec"], flush=True)
    return out
if __name__ == "__main__":
    F = fitparams(); cases = []
    for k in ("omb", "free", "rel"):
        cases.append((f"{k}|core-only", F[k], False, 1.0))
    cases.append(("omb|as-fitted cavity", F["omb"], True, 1.0))
    cases.append(("rel|as-fitted cavity", F["rel"], True, 1.0))
    for ms in (1e-1, 1e-2, 1e-3):
        cases.append((f"omb|core-only|Mc x{ms:g}", F["omb"], False, ms))
    cases.append(("null|Mc x1e-9", F["omb"], False, 1e-9))
    with mp.Pool(5) as P: res = P.map(run, cases, chunksize=1)
    json.dump(res, open("A_ltb_cmb.json", "w"), indent=1, default=float)
