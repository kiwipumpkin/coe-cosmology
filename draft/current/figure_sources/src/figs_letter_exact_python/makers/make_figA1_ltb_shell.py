#!/usr/bin/env python3
"""Fig A1 — LTB shell evolution: flat m+r; E=0; E>0; E<0."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))
from physics import save_fig
from paths import OUT
OUT_NAME = "FigA1_ltb_shell.png"


def evolve_dust(E: float, t_max: float = 8.0, a0: float = 1e-4):
    """Integrate ȧ² = 2E + 2GM/a with GM=1 units; a(0)=a0, ȧ>0."""
    GM = 1.0

    def f(t, y):
        a, = y
        rhs = 2 * E + 2 * GM / max(a, 1e-12)
        if rhs < 0:
            return [0.0]
        return [np.sqrt(rhs)]

    sol = solve_ivp(f, (0, t_max), [a0], dense_output=True, max_step=0.01, rtol=1e-6)
    t = np.linspace(0, min(t_max, sol.t[-1]), 600)
    a = sol.sol(t)[0]
    # stop if recollapse starts
    if E < 0:
        # analytic cycloid max a = GM/|E|
        amax = GM / abs(E)
        # continue through recollapse with signed root
        def f2(t, y):
            a, = y
            rhs = 2 * E + 2 * GM / max(a, 1e-12)
            # use phase from energy: after turnaround ȧ < 0
            return [-np.sqrt(max(rhs, 0.0))]

        # find turnaround
        i_turn = int(np.argmax(a))
        if a[i_turn] > 0.5 * amax:
            t_turn = t[i_turn]
            sol2 = solve_ivp(f2, (t_turn, t_max), [a[i_turn]], dense_output=True, max_step=0.01)
            t2 = np.linspace(t_turn, min(t_max, sol2.t[-1]), 400)
            a2 = sol2.sol(t2)[0]
            t = np.concatenate([t[: i_turn + 1], t2[1:]])
            a = np.concatenate([a[: i_turn + 1], a2[1:]])
    return t, a


def evolve_flat_mr(t_max=8.0):
    """Schematic flat matter+radiation: early a~t^{1/2}, late a~t^{2/3}."""
    t = np.linspace(1e-3, t_max, 800)
    # match at teq=1: continuous a and prefer radiation then matter
    teq = 1.0
    a_eq = teq ** 0.5
    a = np.where(t < teq, t ** 0.5, a_eq * (t / teq) ** (2 / 3))
    return t, a


def main():
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.5), dpi=160)
    panels = [
        (axes[0, 0], "Flat matter + radiation", evolve_flat_mr()),
        (axes[0, 1], r"Dust $E=0$ (marginally bound)", evolve_dust(0.0)),
        (axes[1, 0], r"Dust $E>0$ (expand forever)", evolve_dust(+0.15)),
        (axes[1, 1], r"Dust $E<0$ (recollapse)", evolve_dust(-0.25, t_max=10.0)),
    ]
    for ax, title, (t, a) in panels:
        ax.plot(t, a, color="C0", lw=2.0)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r"$t$ (code units)")
        ax.set_ylabel(r"$a(t)$ / areal factor")
        ax.grid(True, alpha=0.3)
        if "radiation" in title.lower():
            ax.axvline(1.0, color="0.5", ls=":", lw=1)
            ax.text(1.05, ax.get_ylim()[1] * 0.1 if False else 0.3, r"$t_{\mathrm{eq}}$", fontsize=8, color="0.4")

    # annotate asymptotics on first panel
    axes[0, 0].text(0.35, 0.55, r"$a\sim t^{1/2}$", fontsize=9, color="0.35", transform=axes[0, 0].transAxes)
    axes[0, 0].text(0.62, 0.75, r"$a\sim t^{2/3}$", fontsize=9, color="0.35", transform=axes[0, 0].transAxes)

    fig.suptitle("LTB shell evolution (Einstein/Friedmann shell constraint)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    png = save_fig(fig, "_figA1_ltb_shell", also_pdf=False)
    dest = OUT / OUT_NAME
    import shutil
    shutil.move(str(png), str(dest))
    p = dest
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()

