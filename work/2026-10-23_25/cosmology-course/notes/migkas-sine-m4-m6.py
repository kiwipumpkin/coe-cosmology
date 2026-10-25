"""m=4 and m=6 on the frozen geometry times the fixed paper dipole.

Geometry is the WLS Newtonian curve (D, l0, mcc, kin, kout frozen).
The dipole 1 - Amp sin(ell + phi), Amp=-0.06, phi=11.5 deg, is not refit.
Multipoles multiply that product: 1 + a cos(m ell) + b sin(m ell),
ell galactic longitude in radians. Weighted least squares on the 71 Migkas bins,
sigma = (A_err_low + A_err_high)/2. No m=8.
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
D = 58.19699837268826
MCC = 0.0017121308626185587
KIN = 4.034634647911553
KOUT = 1.0235618098210888
L0 = 15.1629523788697
AMP = -0.06
PHI_DEG = 11.5

OUT_PNG = Path("/workspace/cosmology-course/notes/migkas-sine-m4-m6.png")
OUT_JSON = Path("/workspace/cosmology-course/notes/migkas-sine-m4-m6.json")

Y = MA
SIG = 0.5 * (MLO + MHI)
W = 1.0 / SIG ** 2


def r2_u(y, m):
    return float(1.0 - np.sum((y - m) ** 2) / np.sum((y - y.mean()) ** 2))


def r2_w(y, m, w):
    yb = np.sum(w * y) / np.sum(w)
    return float(1.0 - np.sum(w * (y - m) ** 2) / np.sum(w * (y - yb) ** 2))


def metrics(m):
    return dict(
        R2_unweighted=r2_u(Y, m),
        R2_weighted=r2_w(Y, m, W),
        chi2=float(np.sum(W * (m - Y) ** 2)),
        n_bins=int(Y.size),
    )


def sky():
    rin = 1.0 - KIN * RV / D
    rout = 1.0 + KOUT * RV / D
    a = dict(DR=D / RV, Os=OS, mcc=float(MCC), rin=float(rin), rout=float(rout),
             Lc=0.005, mc=0.0, Ob=0.0)
    rr = P.run("ptshell", a, 0.0, 0.0)
    af = rr["af"]
    if (not np.all(np.isfinite(af))) or rr["vrmin"] <= 0:
        raise RuntimeError("flow not usable")
    # af is on the shell angle. shifted maps it onto Galactic longitude
    # the same way migkas-geometry-fit.py does: value at ell is af(ell - l0).
    norm = af / af.mean()
    return shifted(norm, L0)


def dipole_factor(ell_deg, amp=AMP, phi_deg=PHI_DEG):
    """1 - Amp sin(Omega * ell + phi), ell in degrees, Omega = pi/180."""
    return 1.0 - amp * np.sin(np.radians(ell_deg) + np.radians(phi_deg))


def coeff_block(a, b):
    c = float(np.hypot(a, b))
    psi = float(np.degrees(np.arctan2(b, a)))
    phi = float(np.degrees(np.arctan2(a, b)) % 360.0)
    return dict(a=float(a), b=float(b), C=c, psi_deg=psi, Amp=-c, phi_deg=phi)


def fit_multipoles(base_bins, ell_deg, orders):
    """Weighted LS of base * (1 + sum a cos(m ell) + b sin(m ell)) onto Y."""
    cols = []
    for m in orders:
        mr = np.radians(m * ell_deg)
        cols.append(base_bins * np.cos(mr))
        cols.append(base_bins * np.sin(mr))
    x = np.column_stack(cols)
    sw = np.sqrt(W)
    beta, *_ = np.linalg.lstsq(x * sw[:, None], (Y - base_bins) * sw, rcond=None)
    coefs = {}
    for i, m in enumerate(orders):
        coefs[str(m)] = coeff_block(beta[2 * i], beta[2 * i + 1])
    return coefs


def apply_multipoles(base, ell_deg, coefs):
    fac = np.ones_like(base, dtype=float)
    for key, c in coefs.items():
        m = int(key)
        mr = np.radians(m * ell_deg)
        fac = fac + c["a"] * np.cos(mr) + c["b"] * np.sin(mr)
    return base * fac


def case(name, nfree, pred_bins, coefs=None, note=None):
    out = dict(free_parameters=nfree, **metrics(pred_bins))
    if note:
        out["note"] = note
    if coefs is not None:
        out["coefficients"] = coefs
    return out


def main():
    geom = sky()
    ell = TH1
    bins = ell[IMIG]
    gbin = geom[IMIG]
    dip = dipole_factor(ell)
    base = geom * dip
    bbin = base[IMIG]

    c4 = fit_multipoles(bbin, bins, (4,))
    c6 = fit_multipoles(bbin, bins, (6,))
    c46 = fit_multipoles(bbin, bins, (4, 6))

    full4 = apply_multipoles(base, ell, c4)
    full6 = apply_multipoles(base, ell, c6)
    full46 = apply_multipoles(base, ell, c46)

    gmet = metrics(gbin)
    dmet = metrics(bbin)
    print("geometry", gmet, flush=True)
    print("dipole", dmet, flush=True)
    print("m4", metrics(full4[IMIG]), c4, flush=True)
    print("m6", metrics(full6[IMIG]), c6, flush=True)
    print("m4m6", metrics(full46[IMIG]), c46, flush=True)

    fig, ax = plt.subplots(figsize=(7.0, 3.675))
    ax.errorbar(bins, Y, yerr=np.vstack([MLO, MHI]), fmt="o", ms=3.2, color="k",
                ecolor="0.45", elinewidth=0.7, capsize=0, zorder=3,
                label="Migkas Fig. 6")
    ax.plot(ell, geom, color="0.45", lw=1.15, label="geometry", zorder=2)
    ax.plot(ell, base, color="C0", lw=1.6,
            label=r"geometry $\times[1-(-0.06)\sin]$", zorder=3)
    ax.plot(ell, full46, color="C1", lw=1.6,
            label=r"geometry $\times$ dipole $\times(m=4+m=6)$", zorder=4)
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$")
    ax.set_xticks(np.arange(0, 361, 45))
    ax.legend(fontsize=7.2, loc="upper right", framealpha=0.92)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=160)

    payload = dict(
        geometry=dict(
            D_Gpc=D, l0_deg=L0, mcc=MCC, Rin_Gpc=float(D - KIN * RV),
            Rout_Gpc=float(D + KOUT * RV), kin=KIN, kout=KOUT,
            Omega_s=OS, R_Gpc=RV,
            source="exact wls_fit in migkas-geometry-fit.json; not refit",
        ),
        weights="sigma=(A_err_low+A_err_high)/2, w=1/sigma^2, chi2=sum((model-data)/sigma)^2",
        R2="unweighted: 1-SS_res/SS_tot about the unweighted mean; weighted: same with w, about the weighted mean",
        modulation=("1 + a_m cos(m*ell) + b_m sin(m*ell), ell galactic longitude in radians. "
                    "Equivalent 1 + C cos(m*ell - psi), C=hypot(a,b), psi=atan2(b,a). "
                    "Paper sine 1 - Amp sin(m*ell + phi), Amp=-C, phi=atan2(a,b) wrapped to [0, 360). "
                    "Angles in degrees. Multipoles multiply geometry x the fixed dipole, not geometry alone."),
        fixed_dipole=dict(Amp=AMP, phi_deg=PHI_DEG,
                          convention="1 - Amp * sin(ell_rad + phi_rad)", refit=False),
        n_bins=71,
        cases=dict(
            geometry=case("geometry", 0, gbin, note="frozen WLS geometry"),
            geometry_x_fixed_dipole=case(
                "dipole", 0, bbin,
                note="Amp=-0.06 and phi=11.5 deg held fixed; not refit"),
            geometry_x_fixed_dipole_x_m4=case("m4", 2, full4[IMIG], c4),
            geometry_x_fixed_dipole_x_m6=case("m6", 2, full6[IMIG], c6),
            geometry_x_fixed_dipole_x_m4_m6=case("m46", 4, full46[IMIG], c46),
        ),
        diagnostics_not_proofs=True,
    )
    # stamp the fixed dipole numbers onto that case for the reader
    payload["cases"]["geometry_x_fixed_dipole"]["Amp"] = AMP
    payload["cases"]["geometry_x_fixed_dipole"]["phi_deg"] = PHI_DEG
    OUT_JSON.write_text(json.dumps(payload, indent=1))
    print("wrote", OUT_PNG, OUT_JSON, flush=True)


if __name__ == "__main__":
    main()
