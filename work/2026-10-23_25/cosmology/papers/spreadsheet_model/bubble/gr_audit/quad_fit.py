"""Best achievable fit of an axisymmetric local-shear pattern A(l) = [1 + s cos^2(l-l0)]^p / mean  (p=2: A~H^2; p=1: A~H) to Migkas Fig.6,
with and without the (1 - amp sin(th+phi)) factor (|amp|<=0.15, th = l - l0). s>0 means H larger ALONG the CoE axis."""
import sys, json, numpy as np
sys.path.insert(0, ".."); sys.path.insert(0, "../../integrated")
from bubble_core import MA, IMIG, TH1, r2, A_from, best_sine
def pattern(s, l0, p=2):
    A = (1 + s * np.cos(np.radians(TH1 - l0))**2)**p; return A / A.mean()
def fit_s_l0(s, l0, p=2):
    A = pattern(s, l0, p); x = A[IMIG]; th = (TH1 - l0)[IMIG]; a, ph = best_sine(x, th, MA)
    return r2(MA, x), r2(MA, A_from(x, th, a, ph)), a, ph
if __name__ == "__main__":
    out = {}
    S = np.r_[-np.geomspace(0.8, 0.005, 30), np.geomspace(0.005, 3, 40)]; L = np.arange(0, 180, 3.0)
    for p in (2, 1):
        best = {}
        for mode, k in (("nosine", 0), ("sine", 1)):
            for sign, Sg in (("s>0 (along>across, Migkas-like)", S[S > 0]), ("s<0 (across>along, LTB GR-A-like)", S[S < 0])):
                b = max(((fit_s_l0(s, l0, p)[k], s, l0) for s in Sg for l0 in L))
                best[mode + " | " + sign] = dict(R2=b[0], s=b[1], l0=b[2], Hr_over_Hperp=1 + b[1])
        out["A~H^%d" % p] = best
    # the other worker's GR-A shear (s = 30.29/89.86 - 1) at its l0 = 12.82
    s_A = 30.287 / 89.856 - 1
    out["GR-A_local_s"] = dict(s=s_A, R2_nosine_sine=fit_s_l0(s_A, 12.8215)[:2])
    json.dump(out, open("quad_fit.json", "w"), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
