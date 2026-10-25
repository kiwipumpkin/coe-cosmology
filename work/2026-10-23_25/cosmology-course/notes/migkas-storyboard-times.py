"""A/A_all vs Galactic longitude at each universe-storyboard epoch. No refit, no sine.

D and R are taken from universe-storyboard.json. The crossing has R=0, so the
integrator is not called (a_f divides by R/D). Later epochs keep mcc, Omega_s
and ell_0. R_in is the 0.50 Gpc floor. R_out is D+R, except today, which uses
the frozen fit R_out.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"
sys.path.insert(0, BUB)
os.chdir(BUB)
import ptshell  # noqa: E402,F401
import prof_core as P  # noqa: E402
from bubble_core import MA, MHI, MLO, TH1, shifted  # noqa: E402

# Frozen angular parameters. Not refit. Match migkas-earlier-times.py.
MCC = 0.0017121308626185587
L0 = 15.1629523788697
OS = 0.049
RIN = 0.50
ROUT_TODAY = 72.83442650308822

STORY = Path("/workspace/cosmology-course/notes/universe-storyboard.json")
OUT_PNG = Path("/workspace/cosmology-course/notes/migkas-storyboard-times.png")
OUT_JSON = Path("/workspace/cosmology-course/notes/migkas-storyboard-times.json")
LDEG = np.loadtxt("/workspace/cosmology/migkas/fig6_digitized.csv", delimiter=",", skiprows=1, usecols=0)


def sky(D, R, Rin, Rout):
    a = dict(
        DR=D / R, Os=OS, mcc=float(MCC),
        rin=float(Rin / D), rout=float(Rout / D),
        Lc=0.005, mc=0.0, Ob=0.0,
    )
    rr = P.run("ptshell", a, 0.0, 0.0)
    af = rr["af"]
    if (not np.all(np.isfinite(af))) or rr["vrmin"] <= 0:
        raise RuntimeError(
            f"flow not usable D={D:.4f} R={R:.4f} vrmin={rr['vrmin']}"
        )
    return af / af.mean(), rr


def main():
    frames = json.loads(STORY.read_text())["frames"]
    rows = []
    curves = []
    skipped = None
    for fr in frames:
        t = float(fr["t_after_Gyr"])
        D = float(fr["D_Gpc"])
        R = float(fr["R_Gpc"])
        if R <= 0:
            skipped = dict(
                t_after_Gyr=t,
                D_Gpc=D,
                R_Gpc=R,
                note="R=0, no profile. Integrator divides by R/D, so the crossing is not plotted.",
            )
            print(f"skip t={t:g} D={D:.4f} R={R:.4f} (R=0, no profile)", flush=True)
            continue
        today = abs(t - 61.2) < 1e-6
        Rout = ROUT_TODAY if today else float(fr["D_plus_R_Gpc"])
        below = bool((D - R) < RIN)
        af, rr = sky(D, R, RIN, Rout)
        full = shifted(af, L0)
        row = dict(
            t_after_Gyr=t,
            D_Gpc=D,
            R_Gpc=R,
            D_minus_R_Gpc=float(D - R),
            D_plus_R_Gpc=float(D + R),
            R_in_Gpc=float(RIN),
            R_out_Gpc=float(Rout),
            R_out_is_frozen_fit=bool(today),
            D_minus_R_below_Rin=below,
            A_pp=float(full.max() - full.min()),
            A_min=float(full.min()),
            A_max=float(full.max()),
            today=bool(today),
        )
        rows.append(row)
        curves.append(full)
        print(
            f"t={t:g} D={D:.4f} R={R:.4f} DmR={D-R:.4f} "
            f"Rout={Rout:.4f} pp={row['A_pp']:.4f} below_Rin={below}",
            flush=True,
        )

    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    ax.errorbar(
        LDEG, MA, yerr=np.vstack([MLO, MHI]), fmt="o", ms=2.6, color="0.55",
        ecolor="0.78", elinewidth=0.55, capsize=0, alpha=0.45, zorder=2,
        label="Migkas, present-day only",
    )
    # Early epochs first, today last and black.
    cmap = plt.cm.viridis
    n_early = sum(1 for r in rows if not r["today"])
    k = 0
    for full, row in zip(curves, rows):
        if row["today"]:
            continue
        col = cmap(0.15 + 0.75 * k / max(n_early - 1, 1))
        k += 1
        ax.plot(
            TH1, full, color=col, lw=1.25, zorder=3,
            label=rf"$t={row['t_after_Gyr']:g}$ Gyr, $D={row['D_Gpc']:.2f}$, $R={row['R_Gpc']:.2f}$",
        )
    for full, row in zip(curves, rows):
        if not row["today"]:
            continue
        ax.plot(
            TH1, full, color="black", lw=2.0, zorder=4,
            label=rf"$t={row['t_after_Gyr']:g}$ Gyr (today), $D={row['D_Gpc']:.2f}$, $R={row['R_Gpc']:.2f}$",
        )
    note = "Crossing $t=0$: $R=0$, no profile"
    ax.text(0.99, 0.97, note, transform=ax.transAxes, ha="right", va="top", fontsize=8, color="0.25")
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$")
    ax.set_xticks(np.arange(0, 361, 45))
    ax.legend(fontsize=6.6, loc="lower left", frameon=True)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)

    payload = dict(
        formula="A/A_all = a_f / <a_f>, no sine, no extra multipoles",
        frozen=dict(ell0_deg=L0, mcc=MCC, Omega_s=OS, R_in_Gpc=RIN, R_out_today_Gpc=ROUT_TODAY),
        D_and_R_source=str(STORY),
        skipped_crossing=skipped,
        D_minus_R_below_Rin_any=any(r["D_minus_R_below_Rin"] for r in rows),
        epochs=[
            dict(
                t_after_Gyr=r["t_after_Gyr"],
                D_Gpc=r["D_Gpc"],
                R_Gpc=r["R_Gpc"],
                A_pp=r["A_pp"],
                D_minus_R_Gpc=r["D_minus_R_Gpc"],
                R_out_Gpc=r["R_out_Gpc"],
                D_minus_R_below_Rin=r["D_minus_R_below_Rin"],
            )
            for r in rows
        ],
    )
    OUT_JSON.write_text(json.dumps(payload, indent=1))
    print("png", OUT_PNG.stat().st_size, "json", OUT_JSON.stat().st_size, flush=True)


if __name__ == "__main__":
    main()
