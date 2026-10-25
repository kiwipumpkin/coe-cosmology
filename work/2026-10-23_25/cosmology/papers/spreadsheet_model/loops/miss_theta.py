import json, sys, numpy as np
sys.path.insert(0, "/workspace/cosmology/figures")
from coe_half_loop_geometry import lb_to_n, ring_basis
F = json.load(open("loop_fit.json")); l0 = F["l0"]; b0 = F["fits"]["b0_only"]["b"]
ang = lambda u, v: np.degrees(np.arccos(np.clip(np.dot(u, v), -1, 1)))
D = lb_to_n(l0, b0); rows = {}
for tag, roll in (("roll45", 45.0), ("roll80p8", F["fits"]["b0_only"]["roll"])):
    Y, _ = ring_basis(D, roll); ring, ray = [], []
    for p in F["pairs"]:
        th = p["l_model"]; P = lb_to_n(*p["migkas"]); t = np.radians(th)
        ring.append(ang(D, P) - np.degrees(np.arccos(np.cos(t)))); ray.append(ang(np.cos(t) * D + np.sin(t) * Y, P))
    rows[tag] = dict(ring=ring, ray=ray, rms_ring=float(np.sqrt(np.mean(np.square(ring)))), rms_ray=float(np.sqrt(np.mean(np.square(ray)))))
    print(tag, np.round(ring, 1), round(rows[tag]["rms_ring"], 1), np.round(ray, 1), round(rows[tag]["rms_ray"], 1))
json.dump({"b0": b0, "theta": [p["l_model"] for p in F["pairs"]], "migkas": [p["migkas"] for p in F["pairs"]], **rows}, open("miss_theta.json", "w"), indent=1)
