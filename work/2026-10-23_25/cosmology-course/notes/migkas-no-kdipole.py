"""Same WLS geometry curve as migkas-geometry-fit.py, with the model's own m=1 removed.

The Newtonian integrator has no separate CMB kinematic boost. a_f is
|V_source - V_Earth| / (R/D), mean-normalised. The dipole in that curve is the
m=1 Fourier component of the shell geometry. This script subtracts that
component (cos psi and sin psi, psi = ell - l0) and does not refit D, l0, or
the density profile. No phase factor, sine template, or higher multipoles.
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
from bubble_core import IMIG, MA, MHI, MLO, TH1, shifted  # noqa: E402

RV = 14.30048287260589
OS = 0.049
# Exact WLS point behind the blue curve (rounded in the legend as D=58.2, l0=15.2).
D = 58.19699837268826
MCC = 0.0017121308626185587
KIN = 4.034634647911553
KOUT = 1.0235618098210888
L0 = 15.1629523788697

OUT_PNG = Path("/workspace/cosmology-course/notes/migkas-no-kdipole.png")
OUT_JSON = Path("/workspace/cosmology-course/notes/migkas-no-kdipole.json")

Y = MA
SIG = 0.5 * (MLO + MHI)
W = 1.0 / SIG ** 2
LDEG = np.loadtxt("/workspace/cosmology/migkas/fig6_digitized.csv", delimiter=",", skiprows=1, usecols=0)


def r2_u(y, m):
    return float(1.0 - np.sum((y - m) ** 2) / np.sum((y - y.mean()) ** 2))


def r2_w(y, m, w):
    yb = np.sum(w * y) / np.sum(w)
    return float(1.0 - np.sum(w * (y - m) ** 2) / np.sum(w * (y - yb) ** 2))


def sky():
    rin = 1.0 - KIN * RV / D
    rout = 1.0 + KOUT * RV / D
    a = dict(DR=D / RV, Os=OS, mcc=float(MCC), rin=float(rin), rout=float(rout),
             Lc=0.005, mc=0.0, Ob=0.0)
    rr = P.run("ptshell", a, 0.0, 0.0)
    af = rr["af"]
    if (not np.all(np.isfinite(af))) or rr["vrmin"] <= 0:
        raise RuntimeError("flow not usable")
    return af / af.mean()


def remove_m1(full, l0):
    """Subtract the m=1 Fourier component in psi = ell - l0. Mean is unchanged."""
    psi = np.radians((TH1 - l0) % 360.0)
    # Uniform 1-degree grid: a_c = (2/N) sum A cos, likewise sin. DC is separate.
    n = full.size
    ac = (2.0 / n) * np.sum(full * np.cos(psi))
    as_ = (2.0 / n) * np.sum(full * np.sin(psi))
    nodip = full - ac * np.cos(psi) - as_ * np.sin(psi)
    return nodip, float(ac), float(as_)


def metrics(m):
    return dict(
        R2_unweighted=r2_u(Y, m),
        R2_weighted=r2_w(Y, m, W),
        chi2=float(np.sum(W * (m - Y) ** 2)),
    )


def main():
    af = sky()
    full = shifted(af, L0)
    nodip, ac, as_ = remove_m1(full, L0)
    # Also quote the cosine series on the unshifted af (psi measured from anti-CoE).
    th = np.radians(TH1)
    coeffs = {}
    for m in range(0, 6):
        if m == 0:
            coeffs["a0"] = float(np.mean(af))
        else:
            coeffs[f"a{m}"] = float((2.0 / af.size) * np.sum(af * np.cos(m * th)))
            coeffs[f"b{m}"] = float((2.0 / af.size) * np.sum(af * np.sin(m * th)))
    base = metrics(full[IMIG])
    cut = metrics(nodip[IMIG])
    print("m1 on sky", ac, as_, flush=True)
    print("coeffs", coeffs, flush=True)
    print("original", base, flush=True)
    print("no m1", cut, flush=True)
    print("delta peak anti-CoE vs CoE", float(full[int(round(L0)) % 360]), float(full[int(round(L0 + 180)) % 360]),
          float(nodip[int(round(L0)) % 360]), float(nodip[int(round(L0 + 180)) % 360]), flush=True)

    fig, ax = plt.subplots(figsize=(9.2, 4.7))
    ax.errorbar(LDEG, Y, yerr=np.vstack([MLO, MHI]), fmt="o", ms=3.6, color="k",
                ecolor="0.45", elinewidth=0.7, capsize=0, zorder=3,
                label=r"Migkas Fig. 6  $A/A_{\mathrm{all}}$")
    ax.plot(TH1, full, color="0.55", lw=1.35, ls="--",
            label=(rf"geometry $a_f$ (with its dipole): "
                   rf"$R^2={base['R2_unweighted']:.3f}$"))
    ax.plot(TH1, nodip, color="C0", lw=1.8,
            label=(rf"same curve, model $m=1$ removed "
                   rf"($a_1={ac:+.3f}$): $R^2={cut['R2_unweighted']:.3f}$"))
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$")
    ax.set_title(r"Shell geometry, model's own dipole removed (not a CMB boost)")
    ax.legend(fontsize=8.0, loc="lower left")
    ax.set_xticks(np.arange(0, 361, 45))
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=160)
    payload = dict(
        what_removed=("No separate kinematic-dipole term exists in the integrator. "
                      "Subtracted the m=1 Fourier component of the geometry curve "
                      "a_f/<a_f> in psi=ell-l0: A -> A - a_c cos psi - a_s sin psi. "
                      "This is the model's own dipole, not a CMB velocity boost."),
        a_c=ac, a_s=as_,
        fourier_about_antiCoE=coeffs,
        D=D, l0=L0, mcc=MCC, kin=KIN, kout=KOUT,
        Rin=float(D - KIN * RV), Rout=float(D + KOUT * RV),
        original=base, m1_removed=cut,
    )
    OUT_JSON.write_text(json.dumps(payload, indent=1))
    print("wrote", OUT_PNG, flush=True)


if __name__ == "__main__":
    main()
