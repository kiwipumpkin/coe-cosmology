#!/usr/bin/env python3
"""Composite of A6-A9 heatmaps (one figure, 2x2) in the style of Fig A16.
Reuses make_four_var_heatmaps.py grids/holds and _hm.r2_geometry_cell exactly."""
import sys, importlib.util
from pathlib import Path
import numpy as np, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import composite_layout as CL
MAKERS = Path("/workspace/figs_letter_exact_python/makers")
spec = importlib.util.spec_from_file_location("fv", MAKERS / "make_four_var_heatmaps.py")
fv = importlib.util.module_from_spec(spec); spec.loader.exec_module(fv)
from run_density_local_match import params_twoscale_local_match, load_migkas
BASE, edges, r2 = fv.BASE, fv.edges, fv.r2_geometry_cell
CACHE = Path(__file__).with_name("heat_2x2_cache.npz")
x_g = np.round(np.arange(3.0, 10.001, 0.25), 3)
eps_g = np.round(np.arange(0.02, 0.8001, 0.04), 3)
alpha_g = np.round(np.arange(0.2, 2.001, 0.1), 3)
kappa_g = np.round(np.arange(0.05, 1.0001, 0.05), 4)

def scan(mig, xv, yv, xk, yk, hold):
    ZA = np.full((len(yv), len(xv)), np.nan); ZB = ZA.copy()
    for j, y in enumerate(yv):
        for i, x in enumerate(xv):
            kw = dict(hold); kw[xk] = float(x); kw[yk] = float(y)
            rec = r2(mig, x=kw.get("x", BASE["x"]), f_u=kw.get("f_u", BASE["f_u"]),
                     alpha=kw.get("alpha", BASE["alpha"]), eps=kw.get("eps", BASE["eps"]),
                     kappa=kw.get("kappa"), af_mode=kw.get("af_mode", BASE["af_mode"]))
            ZA[j, i], ZB[j, i] = rec["R2_A"], rec["R2_B"]
    return ZA, ZB

def main():
    fv.style_four()
    K0 = float(params_twoscale_local_match(float(BASE["x"]), float(BASE["f_u"]))["kappa"])
    if CACHE.exists():
        d = np.load(CACHE); DdA, DdB, DaA, KaA = d["DdA"], d["DdB"], d["DaA"], d["KaA"]
    else:
        mig = load_migkas()
        com = {"f_u": BASE["f_u"], "af_mode": BASE["af_mode"]}
        DdA, DdB = scan(mig, x_g, eps_g, "x", "eps", {"alpha": BASE["alpha"], "kappa": K0, **com})
        DaA, _ = scan(mig, x_g, alpha_g, "x", "alpha", {"eps": BASE["eps"], "kappa": K0, **com})
        KaA, _ = scan(mig, kappa_g, alpha_g, "kappa", "alpha", {"x": BASE["x"], "eps": BASE["eps"], **com})
        np.savez(CACHE, DdA=DdA, DdB=DdB, DaA=DaA, KaA=KaA)
    panels = [
        ("A6  Distance × radius: $R^2_A$", x_g, eps_g, DdA, (BASE["x"], BASE["eps"]), None),
        ("A7  Distance × radius: $R^2_B$", x_g, eps_g, DdB, (BASE["x"], BASE["eps"]), None),
        (r"A8  Distance × velocity: $R^2_A$", x_g, alpha_g, DaA, (BASE["x"], BASE["alpha"]), None),
        (r"A9  Mass × velocity: $R^2_A$", kappa_g, alpha_g, KaA, (K0, BASE["alpha"]), 0.15),
    ]
    cmaps = ["magma", "plasma", "magma", "magma"]  # originals: A=magma, B=plasma
    fig = plt.figure(figsize=CL.FIGSIZE)
    axs = np.array([fig.add_axes(b) for b in CL.boxes()]).reshape(2, 2)
    best = {}
    for ax, cm, (ttl, xv, yv, Z, base, vl) in zip(axs.flat, cmaps, panels):
        pm = ax.pcolormesh(edges(xv), edges(yv), Z, cmap=cm, vmin=0.0, shading="flat",
                           edgecolors="face", linewidth=0, rasterized=True)
        iy, ix = np.unravel_index(np.nanargmax(Z), Z.shape)
        if vl is not None:
            ax.axvline(vl, color="#17becf", ls=":", lw=1.5, label=r"$\kappa=\Omega_m/2$")
        ax.plot(*base, marker="o", ms=8, color="#17becf", mec="k", mew=0.7, ls="none", label="baseline")
        ax.plot(xv[ix], yv[iy], marker="*", ms=14, color="w", mec="k", mew=0.8, ls="none",
                label=rf"best ({xv[ix]:.3g}, {yv[iy]:.3g}), $R^2$={Z[iy, ix]:.3f}")
        ax.set_title(ttl, fontsize=11)
        ax.legend(loc="lower right", fontsize=7.5, framealpha=0.9, handletextpad=0.3, borderpad=0.3)
        p = ax.get_position()
        cax = fig.add_axes([p.x1 + CL.CB_PAD, p.y0, CL.CB_W, p.height])
        cb = fig.colorbar(pm, cax=cax); cb.ax.tick_params(labelsize=8)
        best[ttl] = (float(xv[ix]), float(yv[iy]), float(Z[iy, ix]))
    # shared axes where meaningful
    axs[0, 1].sharey(axs[0, 0]); axs[0, 1].tick_params(labelleft=False)
    axs[1, 1].sharey(axs[1, 0]); axs[1, 1].tick_params(labelleft=False)
    axs[0, 0].set_ylabel(fv.LBL_RAD); axs[1, 0].set_ylabel(fv.LBL_VEL)
    axs[0, 0].set_xlabel(fv.LBL_DIST, labelpad=1); axs[0, 1].set_xlabel(fv.LBL_DIST, labelpad=1)
    axs[1, 0].set_xlabel(fv.LBL_DIST); axs[1, 1].set_xlabel(r"Mass $\kappa$")
    fig.suptitle("Geometry-pair $R^2$ heatmaps (AF=hratio, fixed Amp×$\\phi$; "
                 rf"held $\kappa$={K0:.4g} unless on axis)", fontsize=12.5, y=CL.SUP_Y, va="top")
    out = Path("/workspace/cosmo_preprint/figs/FigApp_heat_2x2.png")
    fig.savefig(out, dpi=200); print(out)
    for k, v in best.items(): print(k, v)

if __name__ == "__main__":
    main()
