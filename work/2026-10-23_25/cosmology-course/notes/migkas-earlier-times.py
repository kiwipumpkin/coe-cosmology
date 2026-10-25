"""A/A_all vs Galactic longitude at five epochs of the frozen shell. No refit, no sine.

Today is the WLS point in migkas-no-kdipole.py. Earlier epochs keep mcc, Omega_s
and ell_0. D follows the E=0 chart, anchored so t=61.2 Gyr is the frozen D.
R_vision scales with time after the crossing. R_in stays on the integrator
floor. R_out is the frozen value today and D+R at earlier epochs.
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
from bubble_core import GPC, IMIG, MA, MHI, MLO, TH1, shifted  # noqa: E402

# Frozen snapshot. Not refit.
D0 = 58.19699837268826
RV0 = 14.30048287260589
MCC = 0.0017121308626185587
KIN = 4.034634647911553
KOUT = 1.0235618098210888
L0 = 15.1629523788697
OS = 0.049
RIN0 = D0 - KIN * RV0
ROUT0 = D0 + KOUT * RV0

T_AFTER0 = 61.2
T_CROSS = 0.550
T_SING0 = T_AFTER0 + T_CROSS
# integ_core.GPC is 1 Gpc in km. c in km/s. Julian year, matching the storyboard's 0.3066 Gpc/Gyr.
C_KMS = 2.99792458e5
GYR_S = 365.25 * 86400.0 * 1e9
C_GPC_PER_GYR = C_KMS * GYR_S / GPC

OUT_PNG = Path("/workspace/cosmology-course/notes/migkas-earlier-times.png")
OUT_JSON = Path("/workspace/cosmology-course/notes/migkas-earlier-times.json")
LDEG = np.loadtxt("/workspace/cosmology/migkas/fig6_digitized.csv", delimiter=",", skiprows=1, usecols=0)

FRACTIONS = (1.0, 0.75, 0.5, 0.25, 0.125)
LABELS = ("t", "3t/4", "t/2", "t/4", "t/8")


def epoch(frac):
    t_after = frac * T_AFTER0
    t_sing = t_after + T_CROSS
    D = D0 * (t_sing / T_SING0) ** (2.0 / 3.0)
    R = RV0 * frac
    # E=0 worldline: D = [(3/2) t_sing sqrt(2GM)]^{2/3}, so dD/dt = (2/3) D/t_sing.
    # The a_f integrator has no separate shell-centre speed; this is that chart derivative.
    v_over_c = (2.0 / 3.0) * D / t_sing / C_GPC_PER_GYR
    return dict(
        frac=frac,
        t_after_Gyr=t_after,
        t_sing_Gyr=t_sing,
        D_Gpc=float(D),
        R_Gpc=float(R),
        D_minus_R_Gpc=float(D - R),
        D_plus_R_Gpc=float(D + R),
        v_over_c=float(v_over_c),
    )


def sky(D, R, Rin, Rout):
    a = dict(
        DR=D / R, Os=OS, mcc=float(MCC),
        rin=float(Rin / D), rout=float(Rout / D),
        Lc=0.005, mc=0.0, Ob=0.0,
    )
    rr = P.run("ptshell", a, 0.0, 0.0)
    af = rr["af"]
    if (not np.all(np.isfinite(af))) or rr["vrmin"] <= 0:
        raise RuntimeError(f"flow not usable D={D:.4f} R={R:.4f} vrmin={rr['vrmin']}")
    return af / af.mean(), rr


def main():
    rows = []
    curves = []
    for frac in FRACTIONS:
        ep = epoch(frac)
        if frac == 1.0:
            D_use, R_use, Rin, Rout = D0, RV0, RIN0, ROUT0
        else:
            D_use, R_use = ep["D_Gpc"], ep["R_Gpc"]
            Rin, Rout = RIN0, ep["D_plus_R_Gpc"]
        af, rr = sky(D_use, R_use, Rin, Rout)
        full = shifted(af, L0)
        below = bool(ep["D_minus_R_Gpc"] < Rin)
        row = dict(ep)
        row.update(
            D_used_Gpc=float(D_use),
            R_used_Gpc=float(R_use),
            R_in_Gpc=float(Rin),
            R_out_Gpc=float(Rout),
            D_minus_R_below_Rin=below,
            DR=float(D_use / R_use),
            rin=float(Rin / D_use),
            rout=float(Rout / D_use),
            vrmin=float(rr["vrmin"]),
            integrator_vE=float(rr["vE"]),
            integrator_Ht=float(rr["Ht"]),
            A_min=float(full.min()),
            A_max=float(full.max()),
            A_pp=float(full.max() - full.min()),
            ell0_deg=L0,
            mcc=MCC,
            Omega_s=OS,
        )
        rows.append(row)
        curves.append(full)
        print(
            f"{frac:g} D={ep['D_Gpc']:.4f} R={ep['R_Gpc']:.4f} "
            f"DmR={ep['D_minus_R_Gpc']:.4f} DpR={ep['D_plus_R_Gpc']:.4f} "
            f"v/c={ep['v_over_c']:.4f} Rin={Rin:.4f} Rout={Rout:.4f} "
            f"pp={row['A_pp']:.4f} below={below}",
            flush=True,
        )

    fig, ax = plt.subplots(figsize=(8.0, 3.25))
    ax.errorbar(
        LDEG, MA, yerr=np.vstack([MLO, MHI]), fmt="o", ms=2.6, color="0.55",
        ecolor="0.78", elinewidth=0.55, capsize=0, alpha=0.55, zorder=2,
        label="Migkas, present-day measurement",
    )
    colors = ("#000000", "#1d4e89", "#1a7a6d", "#e07a1f", "#c0392b")
    for full, row, col, lab in zip(curves, rows, colors, LABELS):
        ax.plot(
            TH1, full, color=col, lw=1.85 if row["frac"] == 1.0 else 1.35, zorder=4 if row["frac"] == 1.0 else 3,
            label=rf"{lab}   $D={row['D_Gpc']:.1f}$ Gpc, $R={row['R_Gpc']:.2f}$ Gpc",
        )
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$")
    ax.set_xticks(np.arange(0, 361, 45))
    ax.legend(fontsize=7.1, loc="lower left", frameon=True)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=160)

    payload = dict(
        scaling=(
            "R_vision scales as t_after/61.2; D scales as t_sing^{2/3} anchored to the "
            "frozen D(61.2 Gyr); R_in is held at the 0.50 Gpc integrator floor, not D-R; "
            "R_out is the frozen D+kout*R today and D+R at earlier epochs; "
            "mcc, Omega_s and ell_0 are unchanged. No sine."
        ),
        v_over_c=(
            "dD/dt of the E=0 worldline that defines these epochs, "
            "v=(2/3)*D/t_sing, divided by c using integ_core.GPC. "
            "The a_f integrator does not return a shell-centre speed."
        ),
        c_Gpc_per_Gyr=C_GPC_PER_GYR,
        formula="A/A_all(ell)=a_f(ell-l0)/<a_f>, a_f=|V_source-V_Earth|/(R/D)",
        frozen=dict(D=D0, R=RV0, mcc=MCC, kin=KIN, kout=KOUT, l0=L0, Rin=RIN0, Rout=ROUT0, Omega_s=OS),
        D_minus_R_drops_below_Rin=any(r["D_minus_R_below_Rin"] for r in rows),
        epochs=rows,
        longitude_deg=TH1.tolist(),
        A_over_Aall={LABELS[i]: curves[i].tolist() for i in range(5)},
    )
    OUT_JSON.write_text(json.dumps(payload, indent=1))
    print("c_Gpc_per_Gyr", C_GPC_PER_GYR, flush=True)
    print("png", OUT_PNG.stat().st_size, "json", OUT_JSON.stat().st_size, flush=True)


if __name__ == "__main__":
    main()
