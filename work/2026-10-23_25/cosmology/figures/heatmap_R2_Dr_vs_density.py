#!/usr/bin/env python3
"""R^2 heat maps: D/r (x) vs a density parameter (y), with the OTHER parameters
(kappa, l0, b0 within ±30 deg, remaining density params) re-optimised at every cell
(Nelder-Mead, 1-2 starts: best solutions from R2_vs_Dr_density_families_raw.json at nearest D/r
with the y-parameter overwritten, (and for exp maps also the pure-exp best), maxfev=1200). Same model / direct R^2 as
r2_vs_Dr_density_families.py (kappa free; local-density shape convention of the code)."""
import json, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
import r2_vs_Dr_density_families as M
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT = M.OUT
raw = json.load(open(f"{OUT}/R2_vs_Dr_density_families_raw.json"))["results"]
DRG = np.round(np.linspace(1.2, 5.0, 12), 3)
MAPS = {
    "DL": dict(fam="a_exp_floor", idx=3, ys=np.round(np.geomspace(0.5, 30, 15), 3), log=True, lab=r"$D/L$ (exp core)"),
    "fu": dict(fam="a_exp_floor", idx=4, ys=np.round(np.linspace(0.0, 0.98, 15), 3), log=False, lab=r"uniform floor fraction $f_u$ (exp core)"),
    "n":  dict(fam="b_powerlaw_floor", idx=3, ys=np.round(np.linspace(0.0, 10.0, 15), 3), log=False, lab=r"power-law index $n$, $\rho\propto(1+R/R_c)^{-n}$"),
}

def starts(fam, Dr):
    S = []
    Drs = sorted(set(r["D_over_r"] for r in raw), key=lambda d: abs(d - Dr))[:1]
    for r in raw:
        if r["strict_local_density"] or r["D_over_r"] not in Drs: continue
        x = list(r["x"])
        if fam == "a_exp_floor" and r["family"] == "d_pure_exp": x = x + [0.0]
        if r["family"] == fam or (fam == "a_exp_floor" and r["family"] == "d_pure_exp"):
            S.append(np.array(x, float))
    return S

def cell(args):
    key, Dr, y, seed = args
    m = MAPS[key]; fam, idx = m["fam"], m["idx"]
    yv = np.log(y) if m["log"] else y
    eps = 1.0 / Dr
    def f(z):
        v = np.insert(z, idx, yv) if idx < 3 + len(M.FAM[fam]["names"]) else None
        return M.objective(v, fam, eps)
    b = M.bounds(fam); rng = np.random.default_rng(seed)
    S = [np.delete(s, idx) for s in starts(fam, Dr)]
    pass  # random start dropped for speed
    best = (10, None)
    for z0 in S:
        r = minimize(f, z0, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-8, maxfev=1200, adaptive=True))
        if r.fun < best[0]: best = (r.fun, r.x)
    v = np.insert(best[1], idx, yv)
    kap, l0, b0, dens = M.unpack(fam, v)
    return key, Dr, y, -best[0], dict(kappa=float(kap), l0=float(l0 % 360), b0=float(b0), density=dens.tolist())

if __name__ == "__main__":
    jobs = [(k, Dr, y, 7 * i + j) for k in MAPS for i, Dr in enumerate(DRG) for j, y in enumerate(MAPS[k]["ys"])]
    with Pool(8) as p: out = p.map(cell, jobs, chunksize=2)
    res = {}
    for k, m in MAPS.items():
        Z = np.full((len(m["ys"]), len(DRG)), np.nan); P = {}
        for (kk, Dr, y, R, par) in out:
            if kk != k: continue
            i, j = list(m["ys"]).index(y), list(DRG).index(Dr); Z[i, j] = R; P[f"{Dr}_{y}"] = par
        ys = m["ys"]
        def edges(v, lg):
            w = np.log(v) if lg else v; dv = np.diff(w).mean(); e = np.concatenate([w - dv / 2, [w[-1] + dv / 2]])
            return np.exp(e) if lg else e
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        pm = ax.pcolormesh(edges(DRG, False), edges(ys, m["log"]), Z, cmap="magma", vmin=max(0, np.nanmin(Z)), vmax=np.nanmax(Z), shading="flat")
        cs = ax.contour(DRG, ys, Z, levels=[0.6, 0.7, 0.72, 0.74, 0.75, 0.77], colors="c", linewidths=0.8)
        ax.clabel(cs, fontsize=7, fmt="%.2f")
        iy, ix = np.unravel_index(np.nanargmax(Z), Z.shape)
        ax.plot(DRG[ix], ys[iy], "*", ms=16, color="w", mec="k", label=f"best cell D/r={DRG[ix]:.2f}, y={ys[iy]:.3g}, $R^2$={Z[iy, ix]:.3f}")
        sub = Z[:, (DRG >= 3) & (DRG <= 4)]
        iy2, ix2 = np.unravel_index(np.nanargmax(sub), sub.shape); Dsub = DRG[(DRG >= 3) & (DRG <= 4)]
        ax.plot(Dsub[ix2], ys[iy2], "D", ms=8, color="#17becf", mec="k", label=f"best in 3≤D/r≤4: D/r={Dsub[ix2]:.2f}, y={ys[iy2]:.3g}, $R^2$={sub[iy2, ix2]:.3f}")
        ax.axvline(3, color="w", ls="--", lw=1.2); ax.axvline(4, color="w", ls="--", lw=1.2)
        if m["log"]: ax.set_yscale("log")
        ax.set_xlabel(r"$D/r$  ($d/D=r/D$; $D = 14.3\,D/r$ Gpc)"); ax.set_ylabel(m["lab"])
        ax.set_title(r"direct $R^2_A$ vs Migkas; $\kappa,\ell_0,b_0$ & other density params re-optimised per cell", fontsize=9)
        cb = fig.colorbar(pm, ax=ax); cb.set_label(r"$R^2_A$ (no rescale)")
        ax.legend(fontsize=7.5, loc="lower right", framealpha=0.9)
        fig.text(0.01, 0.01, "Grid 12×15; each cell = best of 1–2 warm-started Nelder–Mead runs (maxfev 1200). 71 sector points overlap (N_eff≈16).", fontsize=7, color="#444")
        fig.tight_layout(rect=(0, 0.03, 1, 1))
        for ext in ("png", "pdf"): fig.savefig(f"{OUT}/heatmap_R2_Dr_vs_{k}.{ext}", dpi=200)
        plt.close(fig)
        res[k] = dict(family=m["fam"], y=ys.tolist(), D_over_r=DRG.tolist(), Z=Z.tolist(), params=P,
                      best=dict(D_over_r=float(DRG[ix]), y=float(ys[iy]), R2=float(Z[iy, ix])),
                      best_in_3_4=dict(D_over_r=float(Dsub[ix2]), y=float(ys[iy2]), R2=float(sub[iy2, ix2])))
        print(k, res[k]["best"], res[k]["best_in_3_4"], flush=True)
    json.dump(dict(doc=__doc__, maps=res), open(f"{OUT}/heatmap_R2_Dr_vs_density.json", "w"), indent=1)
