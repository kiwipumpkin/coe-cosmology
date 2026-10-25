import json, B_hybrid_fit as B
out = {"no_cavity_free_mcc": B.fit(None, cavity=False)}; print(out, flush=True)
out["no_cavity_ladder"] = []
for mcc in (1.7e-3, 5e-4, 1e-4, 2e-6):
    d = B.fit(mcc, cavity=False); out["no_cavity_ladder"].append(d); print("nocav", d, flush=True)
json.dump(out, open(B.OUTF.replace(".json", "_nocav.json"), "w"), indent=1, default=float)
