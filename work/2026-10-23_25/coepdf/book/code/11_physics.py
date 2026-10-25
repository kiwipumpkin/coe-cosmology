"""Shared physics / data helpers for CoE letter figures.

Imports digitised Migkas CSVs and may call Newtonian helpers under the
cosmology bubble tree. Figure scripts should import from here, not copy
legacy plotting code.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

from paths import (
    BUBBLE,
    CANCEL_JSON,
    COE_ANALEMMA_JSON,
    COE_CARDINAL_JSON,
    DTT_JSON,
    LTB_CURVE_NPY,
    MIGKAS_CSV,
    MIGKAS_FIG7_CSV,
    NEWTON_CURVE_NPY,
    THIN_FREE,
    TWO_EMITTERS_JSON,
    UNIVERSE_NEWT,
    UNIVERSE_REL,
    VS_KILL_JSON,
)

# --- headline constants (letter / CoE_calculations.md) ---
R_VISION = 14.30048287260589  # Gpc
K_IN = 1.0276268339399446
K_OUT = 1.4754794866889613
H_MEAN = 71.6  # km/s/Mpc (cycle convention)
OS_OMB = 0.049
BAND_CENTRES = (29.5, 122.6, 187.5, 264.4)  # deg on model circle
BAND_HALF = 9.0  # deg
BAND_COLORS = ("#0097b2", "#6fa800", "#e0107a", "#2166d4")
BAND_KINDS = ("MAX", "MIN", "MAX", "MIN")
TILT_DEG = 15.0
ROLL_DEG = 80.8

# CoE sky: +15° tilt → Galactic (ℓ,b)=(180°, −15°)
COE_GAL_L, COE_GAL_B = 180.0, -15.0
PLANO_LAT, PLANO_LON, PLANO_H = 33.02, -96.70, 200.0

# A_F four-harmonic coefficients (bubble_core convention)
AF_COEFFS = {
    0: (0.9965, 0.0),
    1: (0.0296, 0.0434),
    2: (0.0636, 0.0366),
    3: (-0.0299, -0.0114),
}

TH1 = np.arange(0.0, 360.0, 1.0)


def A_F(ell: np.ndarray) -> np.ndarray:
    """Migkas Fourier reconstruction A_F(ℓ)."""
    ell = np.asarray(ell, dtype=float)
    rad = np.radians(ell)
    y = np.full_like(rad, AF_COEFFS[0][0], dtype=float)
    for m, (c, s) in AF_COEFFS.items():
        if m == 0:
            continue
        y = y + c * np.cos(m * rad) + s * np.sin(m * rad)
    return y


def load_migkas_fig6():
    """Digitised Migkas 2020 Fig. 6 A/A_all vs Galactic longitude.

    Uses fig6_digitized.csv (letter gold / bubble_core): ~5° sampling with
    asymmetric A_err_low / A_err_high. Returns (ell, A, err_low, err_high).
    """
    data = np.genfromtxt(MIGKAS_CSV, delimiter=",", names=True)
    return (
        np.asarray(data["l_deg"], float),
        np.asarray(data["A_over_Aall"], float),
        np.asarray(data["A_err_low"], float),
        np.asarray(data["A_err_high"], float),
    )


def load_migkas_fig7():
    data = np.genfromtxt(MIGKAS_FIG7_CSV, delimiter=",", names=True)
    return (
        np.asarray(data["galactic_longitude_deg"], float),
        np.asarray(data["A_over_Aall"], float),
        np.asarray(data["A_over_Aall_err_1sigma"], float),
    )


def load_fit_curves():
    """Return stored Newtonian / free / relativistic A(ℓ) on the 1° grid.

    Prefer recomputing the Ωs=0.049 Newtonian curve via newt_omb when available;
    fall back to the stored A array in universe_newt_omb.json.
    """
    newt_j = json.loads(UNIVERSE_NEWT.read_text())
    sine = newt_j["Rs14.3|sine"]
    A_newt = np.asarray(sine["A"], float)

    # free Ωs Newtonian (previous fit) from newton_curve.npy row 1
    A_free = np.load(NEWTON_CURVE_NPY)[1].astype(float)
    free_meta = json.loads(THIN_FREE.read_text())

    A_rel = np.load(LTB_CURVE_NPY).astype(float)
    rel_j = json.loads(UNIVERSE_REL.read_text())

    meta = {
        "newt": {
            "D": sine["D"],
            "R_in": sine["R_in"],
            "R_out": sine["R_out"],
            "Os": sine["Os"],
            "l0": sine["l0"],
            "amp": sine["amp"],
            "phi": sine["phi"],
            "R2": sine["R2"],
            "mcc": sine["mcc"],
        },
        "free": {
            "D": free_meta["DR"] * R_VISION,
            "Os": free_meta["Os"],
            "l0": free_meta["l0"],
            "amp": free_meta.get("amp"),
            "phi": free_meta.get("phi"),
            "R2": free_meta["R2"],
            "rin": free_meta["rin"],
            "rout": free_meta["rout"],
        },
        "rel": {
            "D": rel_j["D"],
            "R_in": rel_j["R_in"],
            "R_out": rel_j["R_out"],
            "Os": rel_j["Os"],
            "l0": rel_j["l0_s"],
            "amp": rel_j["amp"],
            "phi": rel_j["phi"],
            "R2": rel_j["R2_s"],
            "mcc": rel_j["mcc"],
            "z": rel_j.get("z", 0.3),
        },
    }
    return TH1.copy(), A_newt, A_free, A_rel, meta


def recompute_newt_omb_A(force: bool = False):
    """Optionally recompute Newtonian Ωs=0.049 A(ℓ) via bubble newt_omb."""
    cwd = os.getcwd()
    try:
        sys.path.insert(0, str(BUBBLE))
        os.chdir(BUBBLE)
        import newt_omb as N
        from bubble_core import A_from, shifted

        meta = json.loads(UNIVERSE_NEWT.read_text())["Rs14.3|sine"]
        a = N.geo(meta["D"], meta["mcc"], meta["kin"], meta["kout"], R_VISION)
        af, _rr = N.af_of(a)
        A = A_from(shifted(af, meta["l0"]), TH1 - meta["l0"], meta["amp"], meta["phi"])
        return TH1.copy(), np.asarray(A, float)
    finally:
        os.chdir(cwd)


def r2_score(y, yhat):
    y = np.asarray(y, float)
    yhat = np.asarray(yhat, float)
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan


def load_dtt():
    return json.loads(DTT_JSON.read_text())


def kill_from_doppler(dtt_row: dict) -> tuple[np.ndarray, np.ndarray]:
    """kill(θ) = −(Doppler − ⟨Doppler⟩); θ toward CoE."""
    th = np.asarray(dtt_row["theta_toward_deg"], float)
    dop = np.asarray(dtt_row["dtt_doppler"], float)
    kill = -(dop - dop.mean())
    return th, kill


def extend_0_360(th_0_180: np.ndarray, f_0_180: np.ndarray):
    """Mirror an axisymmetric 0–180° pattern onto a full 0–360° great circle."""
    th = np.asarray(th_0_180, float)
    f = np.asarray(f_0_180, float)
    # assume th goes 0..180
    th2 = np.concatenate([th, 360.0 - th[-2:0:-1]])
    f2 = np.concatenate([f, f[-2:0:-1]])
    return th2, f2


def unit_norm(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, float)
    m = np.max(np.abs(y))
    return y / m if m > 0 else y


def legendre_coeffs(f_theta: np.ndarray, theta_deg: np.ndarray, lmax: int = 4):
    """Legendre a_ℓ for axisymmetric f(θ), θ from CoE (deg)."""
    from numpy.polynomial.legendre import legval

    th = np.radians(theta_deg)
    mu = np.cos(th)
    o = np.argsort(mu)
    mu, f = mu[o], np.asarray(f_theta, float)[o]
    out = []
    for ell in range(lmax + 1):
        P = legval(mu, [0] * ell + [1])
        out.append(float((2 * ell + 1) / 2 * np.trapezoid(f * P, mu)))
    return out


def partner_cancel(doppler: np.ndarray, theta: np.ndarray, mode: str = "dip_p2"):
    """Cancel Doppler with a Legendre partner (dipole+P2 or ℓ≤4)."""
    a = legendre_coeffs(doppler, theta, lmax=4)
    from numpy.polynomial.legendre import legval

    mu = np.cos(np.radians(theta))
    if mode == "dip_p2":
        partner = -(a[1] * legval(mu, [0, 1]) + a[2] * legval(mu, [0, 0, 1]))
    elif mode == "l_le_4":
        partner = np.zeros_like(mu)
        for ell in range(1, 5):
            partner = partner - a[ell] * legval(mu, [0] * ell + [1])
    else:
        raise ValueError(mode)
    residual = doppler + partner - (doppler + partner).mean()
    return partner, residual, a


def load_coe_sky_meta():
    card = json.loads(COE_CARDINAL_JSON.read_text())
    return {
        "gal_lb": tuple(card["coe_galactic_lb_deg"]),
        "ra_deg": card["coe_icrs_ra_deg"],
        "dec_deg": card["coe_icrs_dec_deg"],
        "site_lat": card["site_lat_deg"],
        "site_lon": card["site_lon_deg"],
        "site_h": card["site_height_m"],
        "site_name": card["site_name"],
    }


def save_fig(fig, stem: str, also_pdf: bool = True, dpi: int = 200):
    from paths import OUT

    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / f"{stem}.png"
    fig.savefig(png, dpi=dpi, bbox_inches="tight", facecolor="white")
    if also_pdf:
        try:
            fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
        except Exception:
            pass
    return png


def shell_radii(D: float, kin: float = K_IN, kout: float = K_OUT, R: float = R_VISION):
    return D - kin * R, D + kout * R


# --- sky orientation geometry ---
def band_circle_xyz(theta_c: float, half: float = BAND_HALF, n: int = 360):
    """Unit-sphere band at polar angle theta_c from +x (anti-CoE), width ±half."""
    # In model frame: anti-CoE = +x, CoE = −x; observer at origin.
    # Band is a latitude belt at polar angle θ from +x.
    ph = np.linspace(0, 2 * np.pi, n, endpoint=False)
    # sample mid, inner, outer edges
    out = {}
    for name, th in (("mid", theta_c), ("lo", theta_c - half), ("hi", theta_c + half)):
        t = np.radians(th)
        # x = cos θ, cylindrical radius = sin θ
        x = np.full_like(ph, np.cos(t))
        y = np.sin(t) * np.cos(ph)
        z = np.sin(t) * np.sin(ph)
        out[name] = np.column_stack([x, y, z])
    return out


def rotate_xyz(xyz: np.ndarray, tilt_deg: float = 0.0, roll_deg: float = 0.0):
    """Apply tilt about y then roll about the CoE axis (−x after tilt).

    Untilted: anti-CoE at Galactic (ℓ,b)=(0,0) → +x in model.
    +tilt: rotate about y so anti-CoE goes to b=+tilt (and CoE to b=−tilt).
    +roll: rotate about the anti-CoE axis.
    """
    t = np.radians(tilt_deg)
    r = np.radians(roll_deg)
    # tilt about +y: (x,z) -> (x cos t - z sin t, z cos t + x sin t)
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    xt = x * np.cos(t) - z * np.sin(t)
    yt = y
    zt = x * np.sin(t) + z * np.cos(t)
    # roll about anti-CoE axis (tilted +x'): rotate (y,z) in the plane ⊥ to x'
    # After tilt, anti-CoE direction is (cos t, 0, sin t).
    # Use Rodrigues about that axis.
    ax = np.array([np.cos(t), 0.0, np.sin(t)])
    v = np.column_stack([xt, yt, zt])
    if abs(roll_deg) < 1e-12:
        return v
    k = ax
    c, s = np.cos(r), np.sin(r)
    # v' = v cos + (k×v) sin + k (k·v) (1-cos)
    kxv = np.cross(k, v)
    kdotv = v @ k
    return v * c + kxv * s + np.outer(kdotv, k) * (1 - c)


def xyz_to_galactic(xyz: np.ndarray):
    """Map model-frame unit vectors to Galactic (ℓ,b) deg.

    Convention after tilt/roll: anti-CoE direction → Galactic (ℓ,b).
    Untilted anti-CoE at (0,0); +15° tilt → anti-CoE at (0,+15°), CoE at (180,−15°).
    """
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    # Galactic: ℓ along equator, b latitude. Map x→ toward ℓ=0,b=0; y→ℓ=90; z→b=+90
    b = np.degrees(np.arcsin(np.clip(z, -1, 1)))
    ell = np.degrees(np.arctan2(y, x)) % 360.0
    return ell, b


def mollweide_xy(lon_deg: np.ndarray, lat_deg: np.ndarray, center_lon: float = 0.0):
    """Mollweide projection; lon increases to the left (astro convention)."""
    lon = np.radians(((lon_deg - center_lon + 180) % 360) - 180)
    lat = np.radians(lat_deg)
    # solve 2γ + sin 2γ = π sin φ
    gamma = lat.copy()
    for _ in range(12):
        gamma = gamma - (2 * gamma + np.sin(2 * gamma) - np.pi * np.sin(lat)) / (
            4 * np.cos(gamma) ** 2
        )
    x = (2 * np.sqrt(2) / np.pi) * lon * np.cos(gamma)
    # flip x so longitude increases to the left
    x = -x
    y = np.sqrt(2) * np.sin(gamma)
    return x, y


def load_cancel_json():
    return json.loads(CANCEL_JSON.read_text())


def load_vs_kill():
    return json.loads(VS_KILL_JSON.read_text())


def load_two_emitters():
    return json.loads(TWO_EMITTERS_JSON.read_text())

