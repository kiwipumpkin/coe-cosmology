import time, rm_core as C, json
for Mc, D, R in ((2.472258e54, 93.72, 14.30), (6e56, 93.72, 14.30), (2e57, 93.72, 14.30), (1e58, 93.72, 14.3), (1e57, 20, 5), (1e55, 5, 2)):
    t = time.time(); o = C.evaluate(Mc, D, R)
    print(round(time.time() - t, 2), {k: (round(v, 4) if isinstance(v, float) else v) for k, v in o.items() if k != "fits"})
    if "fits" in o: print("   ", {z: (round(f["R2_ns"], 3), round(f["R2_s"], 3), round(f["m2"], 3), round(f["corr"], 3)) for z, f in o["fits"].items()})
