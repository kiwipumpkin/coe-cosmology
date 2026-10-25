#!/usr/bin/env python3
"""R^2 heat maps, D/r in [1.05,10] (16 log-spaced columns) vs exp-core scale:
 (a) y = L [Gpc] (log 0.5-100), D/L = 14.3 (D/r) / L ; (b) y = D/L (log 0.5-40).
Every cell: kappa, l0, b0 (|b0|<=30), f_u re-optimised (Nelder-Mead, maxfev 1000) from 2 warm starts
(kappa-free and rho(D)=0.3rho_c best solutions of R2_vs_Dr_1to10.json at nearest D/r, D/L overwritten).
Then in every column the best cell is polished: 5 extra starts (the two warm starts, the best cell's own
solution, best solutions of the neighbouring cells in that column, 1 random) with maxfev 4000, twice.
Same model/direct R^2 as r2_vs_Dr_density_families.py (kappa free, code local-density convention).
u(r) grid [0.02,2]D, 80 pts; D/L bound widened to 700."""
import json, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
import r2_vs_Dr_density_families as M
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M.RG_LO, M.NRAD = 0.02, 80
M.FAM["a_exp_floor"]["hi"][0] = 700.0
OUT = M.OUT; RG = 14.3; FAM = "a_exp_floor"
scan = [r for r in json.load(open(f"{OUT}/R2_vs_Dr_1to10.json"))["results"] if r["family"] == FAM]
DRG = np.geomspace(1.05, 10, 16)
MAPS = {"L_Gpc": np.geomspace(0.5, 100, 15), "DL": np.geomspace(0.5, 40, 15)}

def DL_of(key, Dr, y): return RG * Dr / y if key == "L_Gpc" else y

def warm(Dr):
    near = min(set(r["D_over_r"] for r in scan), key=lambda d: abs(np.log(d / Dr)))
    return [np.delete(np.array(r["x"]), 3) for r in scan if r["D_over_r"] == near]

def opt(Dr, DL, starts, maxfev, passes=1):
    eps = 1 / Dr; lx = np.log(DL)
    f = lambda z: M.objective(np.insert(z, 3, lx), FAM, eps)
    best = (10, None)
    for z0 in starts:
        r = minimize(f, z0, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9, maxfev=maxfev, adaptive=True))
        for _ in range(passes - 1):
            r = minimize(f, r.x, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-10, maxfev=maxfev, adaptive=True))
        if r.fun < best[0]: best = (r.fun, r.x)
    return -best[0], best[1]

def cell(a):
    key, i, j = a; Dr, y = DRG[j], MAPS[key][i]
    R, z = opt(Dr, DL_of(key, Dr, y), warm(Dr), 1000)
    return key, i, j, R, z

def polish(a):
    key, j, i, zs, seed = a; Dr, y = DRG[j], MAPS[key][i]
    rng = np.random.default_rng(seed); b = M.bounds(FAM)
    rnd = np.delete(np.array([rng.uniform(lo, hi) for lo, hi in b]), 3)
    R, z = opt(Dr, DL_of(key, Dr, y), warm(Dr) + [np.array(s) for s in zs] + [rnd], 4000, passes=2)
    return key, i, j, R, z

def pars(key, i, j, z):
    Dr = DRG[j]; DL = DL_of(key, Dr, MAPS[key][i])
    return dict(D_over_r=float(Dr), D_over_L=float(DL), L_Gpc=float(RG * Dr / DL), kappa=float(np.exp(z[0])),
                l0=float(z[1] % 360), b0=float(z[2]), f_u=float(z[3]))

if __name__ == "__main__":
    jobs = [(k, i, j) for k in MAPS for i in range(15) for j in range(16)]
    with Pool(8) as p: out = p.map(cell, jobs, chunksize=4)
    Z = {k: np.full((15, 16), np.nan) for k in MAPS}; ZZ = {k: [[None] * 16 for _ in range(15)] for k in MAPS}
    for k, i, j, R, z in out: Z[k][i, j] = R; ZZ[k][i][j] = z
    pj = []
    for k in MAPS:
        for j in range(16):
            i = int(np.nanargmax(Z[k][:, j]))
            zs = [ZZ[k][ii][j] for ii in (i - 1, i, i + 1) if 0 <= ii < 15]
            pj.append((k, j, i, zs, 17 * j + (k == "DL")))
    with Pool(8) as p: pol = p.map(polish, pj, chunksize=1)
    for k, i, j, R, z in pol:
        if R > Z[k][i, j]: Z[k][i, j] = R; ZZ[k][i][j] = z
    res = {}
    for k, ys in MAPS.items():
        Zk = Z[k]; zmax = np.nanmax(Zk); iy, ix = np.unravel_index(np.nanargmax(Zk), Zk.shape)
        fig, ax = plt.subplots(figsize=(8.2, 5.6))
        def edges(v): w = np.log(v); d = np.diff(w).mean(); return np.exp(np.concatenate([w - d / 2, [w[-1] + d / 2]]))
        lo = max(np.nanpercentile(Zk, 5), zmax - 0.25)
        pm = ax.pcolormesh(edges(DRG), edges(ys), Zk, cmap="magma", vmin=lo, vmax=zmax, shading="flat")
        cs = ax.contour(DRG, ys, zmax - Zk, levels=[0.01, 0.03], colors=["c", "w"], linewidths=[1.6, 1.2], linestyles=["-", "--"])
        ax.clabel(cs, fmt={0.01: "ΔR²=0.01", 0.03: "ΔR²=0.03"}, fontsize=7)
        colbest = [int(np.nanargmax(Zk[:, j])) for j in range(16)]
        ax.plot(DRG, ys[colbest], "o", ms=4, color="#17becf", mec="k", mew=0.5, label="best cell per D/r column (polished)")
        dec_edges = [0.5, 1, 10, 100] if k == "L_Gpc" else [0.5, 1, 10, 40]
        mk = ["s", "D", "^"]; decs = []
        for q in range(3):
            sel = (ys >= dec_edges[q]) & (ys < dec_edges[q + 1] * 1.0001)
            if not sel.any(): continue
            sub = np.where(sel[:, None], Zk, -np.inf); a, bb = np.unravel_index(np.argmax(sub), sub.shape)
            decs.append(dict(y_decade=[dec_edges[q], dec_edges[q + 1]], **pars(k, a, bb, ZZ[k][a][bb]), R2=float(Zk[a, bb])))
            ax.plot(DRG[bb], ys[a], mk[q], ms=9, color="#ffd92f", mec="k",
                    label=f"best, y∈[{dec_edges[q]},{dec_edges[q+1]}]: D/r={DRG[bb]:.2f}, y={ys[a]:.3g}, R²={Zk[a, bb]:.3f}")
        ax.plot(DRG[ix], ys[iy], "*", ms=18, color="w", mec="k", label=f"best overall: D/r={DRG[ix]:.2f}, y={ys[iy]:.3g}, R²={zmax:.3f}")
        ax.axvspan(3, 4, color="w", alpha=0.12); ax.axvline(3, color="w", ls=":", lw=1); ax.axvline(4, color="w", ls=":", lw=1)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xticks([1.05, 1.5, 2, 3, 4, 5, 7, 10]); ax.set_xticklabels(["1.05", "1.5", "2", "3", "4", "5", "7", "10"])
        ax.set_xlabel(r"$D/r$  ($d=r$, $d/D=r/D$)")
        ax.set_ylabel(r"exp-core scale length $L$ [Gpc]" if k == "L_Gpc" else r"$D/L$ (exp core)")
        top = ax.secondary_xaxis("top", functions=(lambda x: x * RG, lambda x: x / RG)); top.set_xlabel(r"$D$ [Gpc]  ($r=14.3$ Gpc)")
        ax.set_title(r"direct $R^2_A$ vs Migkas; $\kappa,\ell_0,b_0,f_u$ re-optimised per cell", fontsize=9.5)
        cb = fig.colorbar(pm, ax=ax, extend="min"); cb.set_label(r"$R^2_A$ (no rescale)")
        ax.legend(fontsize=6.8, loc="lower left" if k == "DL" else "upper left", framealpha=0.9)
        fig.text(0.01, 0.01, "16×15 grid; warm-started Nelder–Mead per cell, column maxima polished. 71 overlapping sectors (N_eff≈16).", fontsize=7, color="#444")
        fig.tight_layout(rect=(0, 0.03, 1, 1))
        name = "heatmap_R2_Dr1to10_vs_L_Gpc" if k == "L_Gpc" else "heatmap_R2_Dr1to10_vs_DL"
        for ext in ("png", "pdf"): fig.savefig(f"{OUT}/{name}.{ext}", dpi=200)
        plt.close(fig)
        cols = []
        for j in range(16):
            i = colbest[j]; ok = Zk[:, j] > Zk[i, j] - 0.01
            Ls = [pars(k, ii, j, ZZ[k][ii][j])["L_Gpc"] for ii in range(15) if ok[ii]]
            cols.append(dict(**pars(k, i, j, ZZ[k][i][j]), R2=float(Zk[i, j]), L_Gpc_range_dR2_0p01_in_column=[min(Ls), max(Ls)]))
        res[k] = dict(y=ys.tolist(), D_over_r=DRG.tolist(), Z=Zk.tolist(), best=dict(**pars(k, iy, ix, ZZ[k][iy][ix]), R2=float(zmax)),
                      best_per_y_decade=decs, column_best=cols,
                      cell_params=[[pars(k, i, j, ZZ[k][i][j]) for j in range(16)] for i in range(15)])
        print(k, res[k]["best"], flush=True)
    json.dump(dict(doc=__doc__, maps=res), open(f"{OUT}/heatmap_R2_Dr1to10_vs_L.json", "w"), indent=1)
    print("done")
