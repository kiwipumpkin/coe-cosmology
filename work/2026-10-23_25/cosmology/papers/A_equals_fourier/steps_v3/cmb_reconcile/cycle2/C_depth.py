"""C extension 2: the passing near-core variants (D = 0.2-0.3 Gpc) sampled at deeper Migkas depths (z = 0.05, 0.075, 0.12)."""
import json, multiprocessing as mp
import C_near_core as CN
if __name__ == "__main__":
    A = [(0.08, 0.2, z) for z in (0.03, 0.05, 0.075, 0.12)] + [(0.08, 0.3, z) for z in (0.05, 0.075, 0.12)]
    with mp.Pool(7) as P: fits = P.map(CN.Afit, A)
    for f in fits: print({k: f.get(k) for k in ("mcc", "D", "z", "R2_s", "R2_ns", "fail")}, flush=True)
    json.dump(fits, open("C_depth.json", "w"), indent=1, default=float)
