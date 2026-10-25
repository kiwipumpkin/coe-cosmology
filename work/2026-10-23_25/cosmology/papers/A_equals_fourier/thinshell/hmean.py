from tm_core import *
out = {}
for key in ("free", "omb"):
    a, d = newt_params(key); D = a["DR"] * RV; res = {}
    for tag, Rs in (("R=14.3", RV), ("z=0.3 (1.19 Gpc)", 1.19), ("z=0.075 (0.316 Gpc)", 0.316), ("local (0.01 Gpc)", 0.01)):
        b = dict(a, DR=D / Rs); rr, m1 = newt_run(b); af = rr["af"]          # code units: 1/t0 = k
        As = A_from(shifted(af, d["l0"]), TH1 - d["l0"], d["amp"], d["phi"])
        res[tag] = dict(geom=float(af.mean()), sine=float(As.mean()), min=float(af.min()), max=float(af.max()))
    res["H_E/k (Earth radial rate, code Ht)"] = rr["Ht"]; res["rE"] = rr["rE"]; res["vE"] = rr["vE"]
    out[key] = res; print(key, json.dumps(res, indent=0))
json.dump(out, open("/workspace/cosmology/papers/A_equals_fourier/thinshell/hmean.json", "w"), indent=1)
