from dscan import *
def os_sens(DR):
    D = DR * RV; out = []
    base = [r for r in json.load(open("/workspace/cosmology/papers/A_equals_fourier/steps_v3/dscan/dscan_newt.json")) if r["fit"] == "free" and abs(r["DR"] - DR) < 1e-6][0]
    for Os in (0.01, 0.049, 0.3, 1.0, 3.0, 10.0):
        f = lambda u: 1 - newt_eval(D, 10**u[0], Os, 1, 1, u[1])
        r = minimize(f, [np.log10(base["mcc"]), base["l0"]], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-6, maxiter=300))
        mcc = 10**r.x[0]; m = masses(D, mcc, Os, 1, 1)
        out.append(dict(DR=DR, Os=Os, R2=1 - r.fun, mcc=mcc, share=m["share"], Mc=m["Mc"], Mshell=m["Mshell"]))
    print(DR, [(o["Os"], round(o["R2"], 4), round(o["share"], 4)) for o in out], flush=True); return out
if __name__ == "__main__":
    with mp.Pool(7) as P:
        a = P.map_async(fit_newt, [("omb", 30), ("omb", 150), ("free", 30), ("free", 150)], chunksize=1)
        b = P.map_async(os_sens, [T["free"]["DR"], 10, 20, 50, 100, 200], chunksize=1)
        json.dump(dict(extra=a.get(), os_sens=b.get()), open("/workspace/cosmology/papers/A_equals_fourier/steps_v3/dscan/dscan_newt_extra.json", "w"), indent=1, default=float)
