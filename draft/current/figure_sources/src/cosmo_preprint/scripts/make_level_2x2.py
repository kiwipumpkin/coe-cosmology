#!/usr/bin/env python3
"""Composite of A10-A13 one-at-a-time level sets (2x2, Fig A16 style).
Reuses make_one_at_a_time_levelsets.py levels/holds/model_A_curve exactly
(incl. alpha=1.3 held for the d/D sweep, per the maker)."""
import importlib.util
from pathlib import Path
import numpy as np, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import composite_layout as CL
MAKERS = Path("/workspace/figs_letter_exact_python/makers")
spec = importlib.util.spec_from_file_location("ls", MAKERS / "make_one_at_a_time_levelsets.py")
L = importlib.util.module_from_spec(spec); spec.loader.exec_module(L)
from run_density_local_match import load_migkas

def main():
    L.style_ls()
    mig = load_migkas()
    figs = [("x", L.LEVELS_DL, 5.9, "A10"), ("alpha", L.LEVELS_ALPHA, 1.0, "A11"),
            ("kappa", L.LEVELS_KAPPA, 0.23, "A12"), ("eps", L.LEVELS_DD, 0.10, "A13")]
    fig = plt.figure(figsize=CL.FIGSIZE)
    bx = CL.boxes(); a0 = fig.add_axes(bx[0])
    axs = np.array([a0] + [fig.add_axes(b, sharex=a0, sharey=a0) for b in bx[1:]]).reshape(2, 2)
    for a in axs[0]: a.tick_params(labelbottom=False)
    for a in axs[:, 1]: a.tick_params(labelleft=False)
    for ax, (key, levels, hi, tag) in zip(axs.flat, figs):
        hold = L.baseline_hold()
        if key == "eps":
            hold["alpha"] = 1.3  # Chris: observer-radius levels need higher held velocity
        ax.errorbar(mig["l"], mig["A"], yerr=[mig["A_elo"], mig["A_ehi"]], fmt="o", ms=3.0,
                    color="k", ecolor="#999999", elinewidth=0.8, capsize=0, zorder=5, label="Migkas")
        ax.axhline(1.0, color="#888888", ls="--", lw=0.9)
        cols = [plt.cm.viridis(v) for v in np.linspace(0.08, 0.92, len(levels))]
        for val, c in zip(levels, cols):
            kw = dict(hold); kw[key] = float(val)
            lf, ma = L.model_A_curve(x=kw["x"], f_u=kw["f_u"], alpha=kw["alpha"],
                                     eps=kw["eps"], kappa=kw.get("kappa"))
            h = abs(float(val) - hi) < 1e-6
            ax.plot(lf, ma, color="#E69F00" if h else c, lw=2.6 if h else 1.3,
                    label=L.KEY_FMT[key](val) + (r" $\star$" if h else ""), zorder=4 if h else 3)
        held = ", ".join(L.KEY_FMT[k](hold[k] if k != "kappa" else round(hold[k], 4))
                         for k in ("x", "eps", "alpha", "kappa") if k != key)
        ax.set_title(f"{tag}  {L.VARY_TITLE[key]}\n" + r"{\small held} " + held if False else
                     f"{tag}  {L.VARY_TITLE[key]}   (held {held})", fontsize=8)
        ax.grid(True, alpha=0.35)
        ax.legend(loc="upper right", fontsize=7, ncol=2, framealpha=0.9,
                  handlelength=1.6, columnspacing=0.8, borderpad=0.3, labelspacing=0.25)
    ax = axs[0, 0]; ax.set_xlim(0, 360); ax.set_ylim(0.80, 1.30); ax.set_xticks(range(0, 361, 60))
    for a in axs[1]: a.set_xlabel(r"Galactic longitude $\ell$ [deg]")
    for a in axs[:, 0]: a.set_ylabel(r"$A/A_{\mathrm{all}}$")
    fig.suptitle("One-at-a-time level sets vs Migkas $A/A_{\\rm all}$ "
                 rf"(Amp$_A$={L.SHEET_AMP_A}, $\phi_A$={L.SHEET_PHI_RAD} rad; $\star$ = baseline)",
                 fontsize=12.5, y=CL.SUP_Y, va="top")
    out = Path("/workspace/cosmo_preprint/figs/FigApp_level_2x2.png")
    fig.savefig(out, dpi=200); print(out)

if __name__ == "__main__":
    main()
