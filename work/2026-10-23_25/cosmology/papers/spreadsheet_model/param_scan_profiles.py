import sys, json, numpy as np, time
sys.path.insert(0, ".")
exec(open("param_scan.py").read().split("BASE = [")[0])
BASE = ["DR", "LD", "M", "amp", "phi"]
res = json.load(open("param_scan_results.json"))
gb = res["global"]["mig"]["params"]
GRID = dict(DR=np.linspace(1.05, 10, 40), LD=np.geomspace(0.02, 3, 40), M=np.geomspace(0.1*M0, 1000*M0, 50),
            amp=np.linspace(-0.15, 0.15, 41), phi=np.linspace(0, 2*np.pi, 41), l0=np.linspace(0, 360, 73))
prof = {}
t0 = time.time()
for n, g in GRID.items():
    fixedc, profd, pars = [], [], []
    prev = dict(gb)
    for v in g:
        p = dict(CUR); p[n] = v; fixedc.append(TARGETS["mig"](p))
        free = [m for m in BASE if m != n]
        fx = {n: v} if n != "l0" else {"l0": v}
        if n != "l0": fx["l0"] = 0.0
        s1 = {m: prev[m] for m in free}; s2 = {m: gb[m] for m in free}
        val, pp = fit("mig", free, fx, starts=4, seeds=[s1, s2])
        profd.append(val); pars.append(pp); prev = pp
    # backward sweep to fix any path-dependence
    for i in range(len(g) - 2, -1, -1):
        free = [m for m in BASE if m != n]; fx = {n: g[i]}
        if n != "l0": fx["l0"] = 0.0
        val, pp = fit("mig", free, fx, starts=2, seeds=[{m: pars[i+1][m] for m in free}])
        if val > profd[i]: profd[i], pars[i] = val, pp
    prof[n] = dict(grid=list(map(float, g)), fixed=list(map(float, fixedc)), profile=list(map(float, profd)), pars=pars)
    print(n, "max prof", max(profd), "at", g[int(np.argmax(profd))], "max fixed", max(fixedc), "at", g[int(np.argmax(fixedc))], round(time.time()-t0))
json.dump(prof, open("param_scan_profiles.json", "w"))
# 2D DR x LD map, (M, amp, phi) profiled, l0 = 0
DRg = np.linspace(1.1, 10, 30); LDg = np.geomspace(0.03, 3, 30)
Z = np.zeros((len(LDg), len(DRg))); prevrow = None
for j, ld in enumerate(LDg):
    prev = dict(gb)
    for i, dr in enumerate(DRg):
        seeds = [{m: prev[m] for m in ["M", "amp", "phi"]}, {m: gb[m] for m in ["M", "amp", "phi"]}]
        if prevrow is not None: seeds.append({m: prevrow[i][m] for m in ["M", "amp", "phi"]})
        val, pp = fit("mig", ["M", "amp", "phi"], {"DR": dr, "LD": ld, "l0": 0.0}, starts=2, seeds=seeds)
        Z[j, i] = val; prev = pp
        if i == 0: row = []
        row.append(pp)
    prevrow = row
print("2D max", Z.max(), round(time.time()-t0))
np.savez("param_scan_map.npz", DR=DRg, LD=LDg, Z=Z)
