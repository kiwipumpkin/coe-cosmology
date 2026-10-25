#!/usr/bin/env python3
"""CoE analemma from Plano, TX — altitude vs azimuth at fixed local mean time."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time
import astropy.units as u

OUT_DIR = Path(__file__).resolve().parent
PNG = OUT_DIR / "coe_analemma_texas.png"
JSON = OUT_DIR / "coe_analemma_texas.json"

COE_GAL_L, COE_GAL_B = 180.0, -15.0
coe_gal = SkyCoord(l=COE_GAL_L * u.deg, b=COE_GAL_B * u.deg, frame="galactic")
coe_icrs = coe_gal.icrs
COE_RA = float(coe_icrs.ra.deg)
COE_DEC = float(coe_icrs.dec.deg)

SITE_LAT, SITE_LON = 33.02, -96.70
SITE_HEIGHT = 200.0
site = EarthLocation(lat=SITE_LAT * u.deg, lon=SITE_LON * u.deg, height=SITE_HEIGHT * u.m)

YEAR = 2026
LMT_HOUR = 12.0
LON_HOURS = SITE_LON / 15.0


def lmt_to_time(year: int, doy: int, lmt_hour: float) -> Time:
    t0 = Time(f"{year}-01-01 00:00:00", scale="utc")
    lmt_day = (doy - 1) + lmt_hour / 24.0
    utc_day = lmt_day - LON_HOURS / 24.0
    return t0 + utc_day * u.day


def altaz_of(coord: SkyCoord, times: Time):
    aa = coord.transform_to(AltAz(obstime=times, location=site))
    return aa.az.deg.copy(), aa.alt.deg.copy()


def main():
    doys = np.arange(1, 366)
    times = Time([lmt_to_time(YEAR, int(d), LMT_HOUR) for d in doys])

    coe = SkyCoord(ra=COE_RA * u.deg, dec=COE_DEC * u.deg, frame="icrs")
    coe_az, coe_alt = altaz_of(coe, times)

    sun = get_sun(times)
    sun_az, sun_alt = altaz_of(sun, times)

    lst = times.sidereal_time("mean", longitude=site.lon)
    ha_coe = (lst - coe.ra).wrap_at(180 * u.deg).deg

    # One-day diurnal arc (Jan 15), every 2 min
    day_hours = np.arange(0, 24, 2 / 60)
    day_times = Time([lmt_to_time(YEAR, 15, h) for h in day_hours])
    day_az, day_alt = altaz_of(coe, day_times)

    month_doys = {
        "Jan": 15, "Feb": 46, "Mar": 74, "Apr": 105,
        "May": 135, "Jun": 166, "Jul": 196, "Aug": 227,
        "Sep": 258, "Oct": 288, "Nov": 319, "Dec": 349,
    }
    month_idx = {m: d - 1 for m, d in month_doys.items()}

    shape_note = (
        "closed diurnal oval / ellipse-like loop (full hour-angle sweep), "
        "not a tight Sun-like figure-8"
    )

    # ----- figure -----
    fig = plt.figure(figsize=(13.0, 7.6), dpi=170)
    fig.suptitle("CoE analemma from Plano, Texas", fontsize=16, fontweight="bold", y=0.985)
    fig.text(
        0.5, 0.942,
        (
            f"Sampled every day at {LMT_HOUR:.0f}:00 local mean time ({YEAR})  ·  "
            f"Plano {SITE_LAT:.2f}°N, {abs(SITE_LON):.2f}°W\n"
            f"CoE Galactic (ℓ,b)=({COE_GAL_L:.0f}°, {COE_GAL_B:.0f}°)  →  "
            f"ICRS RA {COE_RA:.2f}°, Dec {COE_DEC:+.2f}°"
        ),
        ha="center", va="top", fontsize=9, color="#333",
    )

    ax = fig.add_axes([0.065, 0.14, 0.56, 0.70])

    ax.axhline(0.0, color="#666", lw=1.0, ls="--", zorder=1)
    ax.axhspan(-90, 0, color="#f0f0f0", zorder=0)

    for az in (0, 90, 180, 270, 360):
        ax.axvline(az, color="#e6e6e6", lw=0.8, zorder=1)

    # Full 1-day track (locus context) — above/below horizon
    above = day_alt >= 0
    # Draw continuous segments without huge az jumps
    def plot_az_segments(az, alt, **kw):
        az = np.asarray(az)
        alt = np.asarray(alt)
        if len(az) < 2:
            return
        d = np.abs(np.diff(az))
        breaks = np.where(d > 180)[0]
        starts = np.r_[0, breaks + 1]
        ends = np.r_[breaks, len(az) - 1]
        for s, e in zip(starts, ends):
            if e > s:
                ax.plot(az[s : e + 1], alt[s : e + 1], **kw)

    plot_az_segments(
        day_az[above], day_alt[above],
        color="#c8c8c8", lw=7.0, solid_capstyle="round", zorder=2,
        label="1-day CoE track (Jan 15, above horizon)",
    )
    plot_az_segments(
        day_az[~above], day_alt[~above],
        color="#dcdcdc", lw=2.5, ls=":", zorder=2,
    )

    cmap = plt.cm.turbo
    norm = Normalize(vmin=1, vmax=365)

    # Annual CoE: scatter + HA-sorted smooth outline (avoids 0/360 DOY jump artifact)
    order = np.argsort(ha_coe)
    az_o = coe_az[order]
    alt_o = coe_alt[order]
    # split at az wrap for outline
    plot_az_segments(az_o, alt_o, color="k", lw=0.7, alpha=0.45, zorder=3)

    ax.scatter(
        coe_az, coe_alt, c=doys, cmap=cmap, norm=norm,
        s=16, zorder=4, linewidths=0, label="CoE @ 12:00 LMT each day",
    )

    # Month markers — offset labels away from loop center
    # Rough center of above-horizon arc
    for m, i in month_idx.items():
        ax.scatter(
            [coe_az[i]], [coe_alt[i]], s=55,
            facecolors="white", edgecolors="k", linewidths=1.1, zorder=6,
        )
        # Nudge label toward outside of loop
        ha_i = ha_coe[i]
        # outward in az-alt: push away from az=180, alt=20 centroid-ish
        dx = 8 if coe_az[i] >= 180 else -8
        dy = 8 if coe_alt[i] >= 20 else -10
        # special cases near edges
        if coe_az[i] < 40:
            dx, dy = 10, 4
        elif coe_az[i] > 320:
            dx, dy = -14, 4
        if abs(coe_alt[i]) < 5:
            dy = 10 if coe_alt[i] >= 0 else -12
        ax.annotate(
            m, (coe_az[i], coe_alt[i]),
            textcoords="offset points", xytext=(dx, dy),
            fontsize=8, fontweight="bold", zorder=7,
            arrowprops=dict(arrowstyle="-", color="#666", lw=0.5)
            if abs(dx) > 6 or abs(dy) > 6 else None,
        )

    # Sun analemma (tiny figure-8 near south)
    ax.plot(sun_az, sun_alt, color="#c44e00", lw=1.8, zorder=5, label="Sun analemma (same LMT)")
    ax.scatter(sun_az, sun_alt, c=doys, cmap=cmap, norm=norm, s=7, zorder=5, linewidths=0)
    for m in ("Jun", "Dec"):
        i = month_idx[m]
        ax.annotate(
            f"Sun {m}", (sun_az[i], sun_alt[i]),
            textcoords="offset points",
            xytext=(6, 8 if m == "Jun" else -12),
            fontsize=7.5, color="#a04000", fontweight="bold", zorder=7,
        )

    for az, lab in [(0, "N"), (90, "E"), (180, "S"), (270, "W"), (360, "N")]:
        ax.text(az, -55, lab, ha="center", va="top", fontsize=11, fontweight="bold", color="#333")

    ax.set_xlim(-8, 368)
    ax.set_ylim(-58, 95)
    ax.set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
    ax.set_xlabel("Azimuth (° from North, clockwise)", fontsize=10)
    ax.set_ylabel("Altitude (°)", fontsize=10)
    ax.set_title(
        "CoE annual locus at 12:00 local mean time\n"
        "closed loop = diurnal path of fixed CoE, sampled once per day over the year",
        fontsize=10, pad=8,
    )
    ax.grid(True, alpha=0.22)
    ax.legend(loc="upper right", fontsize=7.8, framealpha=0.95)

    cax = fig.add_axes([0.065, 0.055, 0.56, 0.018])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_label("Day of year", fontsize=9)

    # Sun zoom
    ax2 = fig.add_axes([0.70, 0.56, 0.27, 0.32])
    ax2.plot(sun_az, sun_alt, color="#c44e00", lw=1.8, zorder=3)
    ax2.scatter(sun_az, sun_alt, c=doys, cmap=cmap, norm=norm, s=11, zorder=4, linewidths=0)
    for m in ("Jan", "Mar", "Jun", "Sep", "Dec"):
        i = month_idx[m]
        ax2.scatter([sun_az[i]], [sun_alt[i]], s=32, facecolors="white",
                    edgecolors="#c44e00", zorder=5, linewidths=1.1)
        ax2.annotate(m, (sun_az[i], sun_alt[i]), textcoords="offset points",
                     xytext=(4, 4), fontsize=7.5)
    ax2.axvline(180, color="#bbb", lw=0.8, ls=":")
    ax2.set_xlabel("Azimuth (°)", fontsize=8)
    ax2.set_ylabel("Altitude (°)", fontsize=8)
    ax2.set_title("Sun analemma @ 12:00 LMT\n(classic figure-8)", fontsize=9.5)
    ax2.grid(True, alpha=0.3)
    pad_az = 1.5
    pad_alt = 2.5
    ax2.set_xlim(np.nanmin(sun_az) - pad_az, np.nanmax(sun_az) + pad_az)
    ax2.set_ylim(np.nanmin(sun_alt) - pad_alt, np.nanmax(sun_alt) + pad_alt)

    # CoE alt vs HA
    ax3 = fig.add_axes([0.70, 0.14, 0.27, 0.32])
    ax3.scatter(ha_coe, coe_alt, c=doys, cmap=cmap, norm=norm, s=11, linewidths=0, zorder=3)
    ax3.plot(ha_coe[order], coe_alt[order], color="k", lw=0.6, alpha=0.4, zorder=2)
    ax3.axhline(0, color="#666", lw=0.8, ls="--")
    for m in ("Jan", "Apr", "Jul", "Oct"):
        i = month_idx[m]
        ax3.annotate(m, (ha_coe[i], coe_alt[i]), textcoords="offset points",
                     xytext=(3, 3), fontsize=7.5)
    ax3.set_xlabel("CoE hour angle (°; − = east)", fontsize=8)
    ax3.set_ylabel("Altitude (°)", fontsize=8)
    ax3.set_title("CoE: altitude vs hour angle\n(same daily LMT samples)", fontsize=9.5)
    ax3.set_xlim(-180, 180)
    ax3.set_ylim(-50, 90)
    ax3.grid(True, alpha=0.3)

    fig.text(
        0.70, 0.085,
        (
            f"Shape: {shape_note}.\n"
            f"CoE alt ≈ [{np.nanmin(coe_alt):.1f}°, {np.nanmax(coe_alt):.1f}°]; "
            f"azimuth spans the full circle.\n"
            f"At fixed mean solar time a fixed ICRS direction walks through all HAs once/year,\n"
            f"so the “analemma” is the diurnal oval — Sun’s figure-8 is tiny by comparison."
        ),
        ha="left", va="top", fontsize=7.6, color="#333",
        transform=fig.transFigure,
    )

    fig.savefig(PNG, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {PNG}")

    payload = {
        "title": "CoE analemma from Plano, Texas",
        "png": str(PNG),
        "site": {
            "name": "Plano, Texas",
            "lat_deg": SITE_LAT,
            "lon_deg": SITE_LON,
            "height_m": SITE_HEIGHT,
        },
        "coe_galactic_lb_deg": [COE_GAL_L, COE_GAL_B],
        "coe_icrs_ra_deg": COE_RA,
        "coe_icrs_dec_deg": COE_DEC,
        "year": YEAR,
        "sample_clock": {
            "type": "local_mean_time",
            "lmt_hour": LMT_HOUR,
            "note": (
                "Each day at fixed Local Mean Time "
                f"(UTC = LMT − lon/15h; lon_hours={LON_HOURS:.4f}). "
                "Classic analemma sampling."
            ),
        },
        "shape": {
            "classification": "closed_diurnal_oval",
            "description": shape_note,
            "coe_alt_min_deg": float(np.nanmin(coe_alt)),
            "coe_alt_max_deg": float(np.nanmax(coe_alt)),
            "coe_az_min_deg": float(np.nanmin(coe_az)),
            "coe_az_max_deg": float(np.nanmax(coe_az)),
            "sun_alt_min_deg": float(np.nanmin(sun_alt)),
            "sun_alt_max_deg": float(np.nanmax(sun_alt)),
            "sun_az_min_deg": float(np.nanmin(sun_az)),
            "sun_az_max_deg": float(np.nanmax(sun_az)),
        },
        "month_markers_lmt": {
            m: {
                "doy": month_doys[m],
                "coe_az_deg": float(coe_az[i]),
                "coe_alt_deg": float(coe_alt[i]),
                "sun_az_deg": float(sun_az[i]),
                "sun_alt_deg": float(sun_alt[i]),
                "coe_ha_deg": float(ha_coe[i]),
            }
            for m, i in month_idx.items()
        },
        "daily": {
            "doy": doys.tolist(),
            "coe_az_deg": np.round(coe_az, 3).tolist(),
            "coe_alt_deg": np.round(coe_alt, 3).tolist(),
            "coe_ha_deg": np.round(ha_coe, 3).tolist(),
            "sun_az_deg": np.round(sun_az, 3).tolist(),
            "sun_alt_deg": np.round(sun_alt, 3).tolist(),
        },
    }
    JSON.write_text(json.dumps(payload, indent=2))
    print(f"Wrote {JSON}")
    print(
        f"CoE alt [{payload['shape']['coe_alt_min_deg']:.2f}, {payload['shape']['coe_alt_max_deg']:.2f}]  "
        f"Sun az span {payload['shape']['sun_az_max_deg']-payload['shape']['sun_az_min_deg']:.2f}°"
    )


if __name__ == "__main__":
    main()
