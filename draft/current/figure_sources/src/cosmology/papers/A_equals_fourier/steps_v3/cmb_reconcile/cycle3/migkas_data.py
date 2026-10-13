"""Load the per-cluster Migkas et al. 2020 (A&A 636, A15) eeHIFLUGCS-based sample from VizieR J/A+A/636/A15/tablec1 (313 clusters)."""
import numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))
def load():
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "data", "migkas2020_tablec1.tsv")) if l.strip() and not l.startswith("#")]
    hdr = rows[0]; body = [r for r in rows[3:] if len(r) == len(hdr)]
    d = {h.strip(): np.array([r[i].strip() for r in body]) for i, h in enumerate(hdr)}
    out = {k: d[k].astype(float) for k in ("z", "GLON", "GLAT", "T", "E_T", "e_T", "LX", "e_LX", "FX", "NHtot")}
    out["name"] = d["Cluster"]; out["note"] = d["Note"]; return out
if __name__ == "__main__":
    c = load(); z = c["z"]; print(len(z), "clusters; z range", z.min(), z.max(), "median", np.median(z))
    for lo, hi in ((0, 0.03), (0.03, 0.07), (0.07, 1)): print(lo, hi, ((z >= lo) & (z < hi)).sum())
    print("quantiles", np.quantile(z, [0.25, 0.5, 0.75, 0.9]))
