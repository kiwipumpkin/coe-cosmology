"""C (VARIANT, not a replacement): Chris's geometry with the core near, D ~ 1-2 Gpc, sampled at Migkas depths.
Newtonian: Chris's integrator is scale-free (A depends only on D/Rs, Mc/Mcrit(D), Omega_s, shell edges in units of D), so sampling at
Rs = 0.316 / 0.5 / 1.19 Gpc with the shell rule written in units of Rs (R_in = D - kin Rs) reproduces the far-core R2 exactly at D = (D/Rs) Rs.
Exact GR: LTB (single bang = pure growing mode, i.e. the core 'grows with the structure'), core + uniform Omega_s=0.049 dust everywhere
(no cavity), A = linear directional rate at fixed z (D_L^-1, Chris's LTB A) with l0 + sine, scanned over Mc/Mcrit and D/Rs; then exact CMB.
Writes C_near_core.json."""
import sys, json, time, numpy as np, multiprocessing as mp
from ltb_common import *
from bubble_core import MA
RHOC70 = 9.2042e-27
def mc_kg(mcc, D): return mcc * 4 * np.pi / 3 * RHOC70 * (D * 3.0857e25)**3
def Afit(args):
    mcc, D, z = args
    p = make_p(mc_kg(mcc, D), D, 0.049, 0.03 * D, 5000.0, cavity=False); p[11] = 0.01 * D
    t0 = setup(p, D)
    if t0 is None: return dict(mcc=mcc, D=D, z=z, fail="t0")
    try:
        nb_ = C.nearby_fast(p, t0, D, (z,), 24)
        b, Ad = C.lfits(nb_["DL"][:, 0], nb_["mu"])
    except Exception as e: return dict(mcc=mcc, D=D, z=z, fail=str(e)[:100])
    Hp, Hr = LT.local_rates(D, t0, p)
    return dict(mcc=mcc, D=D, z=z, R2_s=float(b["s"][0]), l0=float(b["s"][1]), amp=float(b["s"][2]), phi=float(b["s"][3]), R2_ns=float(b["ns"][0]), Hr_over_Hp=Rr_over(Hp, Hr))
def Rr_over(Hp, Hr): return float(Hr / Hp)
def cmb_case(args):
    mcc, D = args
    p = make_p(mc_kg(mcc, D), D, 0.049, 0.03 * D, 5000.0, cavity=False); p[11] = 0.01 * D
    t0 = setup(p, D)
    try: c = analyze.cmb(p, t0, D, N=24)
    except Exception as e: return dict(mcc=mcc, D=D, fail=str(e)[:150])
    return dict(mcc=mcc, D=D, t_rec_Myr=c["t_rec_Myr"], null_err=c["null_err"], c_a=c["c_a"].tolist(), cb_a=c["cb_a"].tolist(), beta_a=c["beta_a"],
                c_b=np.asarray(c["c_b"]).tolist(), cb_b=np.asarray(c["cb_b"]).tolist())
if __name__ == "__main__":
    z = 0.12                                     # ~0.5 Gpc, Migkas bulk-flow depth
    grid = [(m, D, z) for D in (1.0, 1.6, 2.5) for m in (0.02, 0.05, 0.08, 0.15, 0.3)]
    with mp.Pool(3) as P: A = P.map(Afit, grid, chunksize=1)
    for a in A: print(a, flush=True)
    ok = [a for a in A if "R2_s" in a]; best = max(ok, key=lambda a: a["R2_s"])
    extra = [(best["mcc"], best["D"], zz) for zz in (0.075, 0.3)]
    with mp.Pool(2) as P: A2 = P.map(Afit, extra, chunksize=1)
    cm = [(best["mcc"], best["D"]), (best["mcc"] / 10, best["D"]), (best["mcc"], 1.0 if best["D"] != 1.0 else 2.5)]
    with mp.Pool(3) as P: Cm = P.map(cmb_case, cm, chunksize=1)
    for c in Cm: print(c, flush=True)
    json.dump(dict(Afits=A, best=best, best_other_z=A2, cmb=Cm), open("C_near_core.json", "w"), indent=1, default=float)
