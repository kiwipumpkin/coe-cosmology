"""Part 2: far-core two-mode mix in exact LTB (growing mode = core, single bang; early-pull/decaying mode = bang-time bump
t_B(r) = Ab exp(-(r/w)^2) centred on the CoE). Base = relativistic headline geometry (D = 32.94 Gpc, Omega_s = 0.049), core only.
Every point: full-domain regularity (validity.full), A(l) R2 at z = 0.3 (l0 + sine, rm_core.lfits vs digitised Migkas Fig. 6),
exact CMB multipoles at fixed local age since the local bang (cmbfix.cmb_age), per-ray null-condition error.
Usage: python3 mix3.py <stage>  (reads mix3_points_<stage>.json: list of [ms, Ab_frac, w_Gpc, N]) -> mix3_<stage>.json"""
import sys, json, time, numpy as np, multiprocessing as mp
from ltb_common import *
import cmbfix, validity
F = fitparams()["rel"]; D = F["D"]
def build(ms, Ab_frac, w):
    p = make_p(F["Mc"] * ms, D, 0.049, F["Rin"], 5000.0, cavity=False)
    if setup(p, D) is None: return None
    for _ in range(3):
        p[6] = Ab_frac * p[5]; p[7] = w
        if setup(p, D) is None: return None
    return p
def evalp(args):
    ms, Ab_frac, w, N = args; T = time.time(); out = dict(ms=ms, mcc=F["mcc"] * ms, Ab_frac=Ab_frac, w=w, N=N)
    p = build(ms, Ab_frac, w)
    if p is None: out["fail"] = "t0"; return out
    t0 = p[5]; out["t0"] = t0; out.update(validity.full(p, t0))
    nb_ = C.nearby_fast(p, t0, D, (0.3,), 24); b, Ad = C.lfits(nb_["DL"][:, 0], nb_["mu"])
    out.update(R2_s=float(b["s"][0]), R2_ns=float(b["ns"][0]), l0=float(b["s"][1]))
    try:
        c = cmbfix.cmb_age(p, t0, D, N=N)
        out.update(psi=c["psi"].tolist(), z=c["z"].tolist(), nerr_max=float(np.max(c["nerr"])), tau_rec=c["tau_rec_over_t0"], stop=c["stop"])
        if "cb" in c: out.update(cb=np.asarray(c["cb"]).tolist(), beta=float(c["beta"]), res234=float(np.max(np.abs(c["cb"][2:5]))))
    except Exception as e: out["cmb_fail"] = str(e)[:150]
    out["sec"] = time.time() - T
    print({k: (float("%.4g" % v) if isinstance(v, float) else v) for k, v in out.items() if k not in ("psi", "z", "stop", "cb")},
          "a2-4", [float("%.3g" % x) for x in out.get("cb", [np.nan] * 5)[2:5]], flush=True)
    return out
if __name__ == "__main__":
    stage = sys.argv[1]; pts = [tuple(x) for x in json.load(open(f"mix3_points_{stage}.json"))]
    with mp.Pool(6) as P: res = P.map(evalp, pts, chunksize=1)
    json.dump(res, open(f"mix3_{stage}.json", "w"), indent=1, default=float)
