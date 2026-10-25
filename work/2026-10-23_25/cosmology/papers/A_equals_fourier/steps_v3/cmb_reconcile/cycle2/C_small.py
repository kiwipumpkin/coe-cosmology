"""C extension: shrink the near-core variant (D and sampling depth together, D/Rs ~ 3) to find where the exact-LTB CMB residual reaches 1e-5."""
import json, multiprocessing as mp
import C_near_core as CN
if __name__ == "__main__":
    A = [(0.08, 0.5, 0.04), (0.08, 0.3, 0.025), (0.08, 0.2, 0.016), (0.05, 0.3, 0.025)]
    with mp.Pool(4) as P: fits = P.map(CN.Afit, A); cm = P.map(CN.cmb_case, [(a[0], a[1]) for a in A])
    for f, c in zip(fits, cm): print(f, "\n   CMB a2..a4", c.get("cb_a", c)[2:5] if "cb_a" in c else c, flush=True)
    json.dump(dict(fits=fits, cmb=cm), open("C_small.json", "w"), indent=1, default=float)
