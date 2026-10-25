"""Chris's Newtonian CoE integrator (ptshell + prof_core.run, exactly as coe_steps_v3.ipynb Step 2), full-sample fit parameters
held fixed (D, Mc, shell, l0, sine); only the sampling radius Rs changes (DR = D/Rs). Writes af(theta; Rs) tables."""
import sys, os, json, numpy as np, multiprocessing as mp
BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"
sys.path.insert(0, BUB); sys.path.insert(0, BUB + "/gr_ltb/refit_mass"); sys.path.insert(0, "/workspace/cosmology/papers/spreadsheet_model/integrated")
RV = 14.30048287260589
def models():
    hug = json.load(open(BUB + "/thin_b_free_hug.json")); omb = json.load(open(BUB + "/universe_newt_omb.json"))["Rs14.3|sine"]
    return {"free": dict(D=hug["DR"] * RV, Os=hug["Os"], mcc=hug["mcc"], rin=hug["rin"], rout=hug["rout"], l0=hug["l0"], amp=hug["amp"], phi=hug["phi"]),
            "omb": dict(D=omb["D"], Os=0.049, mcc=omb["mcc"], rin=omb["R_in"] / omb["D"], rout=omb["R_out"] / omb["D"], l0=omb["l0"], amp=omb["amp"], phi=omb["phi"])}
def af_at(args):
    name, Rs = args; cwd = os.getcwd(); os.chdir(BUB)
    import ptshell, prof_core as P
    m = models()[name]; a = dict(DR=m["D"] / Rs, Os=m["Os"], mcc=m["mcc"], rin=m["rin"], rout=m["rout"], Lc=0.005, mc=0.0, Ob=0.0)
    rr = P.run("ptshell", a, 0.0, 0.0); os.chdir(cwd)
    return name, Rs, np.asarray(rr["af"]).tolist(), float(rr["vrmin"])
if __name__ == "__main__":
    RS = sorted(set(np.round(np.geomspace(0.02, 2.0, 41), 5).tolist() + [0.05, 0.1, 0.2, 0.316, 0.5, 1.19, 3.0, 6.0, 10.0, RV]))
    jobs = [(n, r) for n in ("free", "omb") for r in RS]
    with mp.Pool(7) as Pl: res = Pl.map(af_at, jobs, chunksize=1)
    out = {n: {"Rs": [], "af": [], "vrmin": []} for n in ("free", "omb")}
    for n, r, af, v in res: out[n]["Rs"].append(r); out[n]["af"].append(af); out[n]["vrmin"].append(v)
    out["params"] = models(); json.dump(out, open("af_tables.json", "w"))
    print("done", len(res))
