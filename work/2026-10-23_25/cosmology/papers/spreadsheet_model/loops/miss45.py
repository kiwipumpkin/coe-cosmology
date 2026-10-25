import json, sys, numpy as np
sys.path.insert(0, "/workspace/cosmology/figures")
from coe_half_loop_geometry import lb_to_n, ring_basis
F = json.load(open("loop_fit.json")); l0 = F["l0"]; b0 = F["fits"]["b0_only"]["b"]
ang = lambda u, v: np.degrees(np.arccos(np.clip(np.dot(u, v), -1, 1)))
D = lb_to_n(l0, b0); out = {}
for roll in (F["fits"]["b0_only"]["roll"], 45.0):
    Y, _ = ring_basis(D, roll); row = []
    for p in F["pairs"]:
        t = np.radians(p["theta"]); row.append(ang(np.cos(t) * D + np.sin(t) * Y, lb_to_n(*p["migkas"])))
    out[roll] = row; print(roll, [p["kind"] + str(p["migkas"]) for p in F["pairs"]], np.round(row, 1), "RMS", round(float(np.sqrt(np.mean(np.square(row)))), 1))
json.dump({"b0": b0, "rolls": {str(k): v for k, v in out.items()}, "pairs": [p["migkas"] for p in F["pairs"]]}, open("miss45.json", "w"), indent=1)
