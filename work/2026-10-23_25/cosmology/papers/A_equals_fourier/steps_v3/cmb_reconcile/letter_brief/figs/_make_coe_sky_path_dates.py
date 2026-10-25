#!/usr/bin/env python3
"""Make a clean annual CoE sky-path chart for Plano, Texas."""
from __future__ import annotations

from pathlib import Path
import calendar

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.utils import iers

iers.conf.auto_download = False

i = Path(__file__).resolve().parent
OUT = i / "coe_sky_path_dates_texas.png"

# Inputs supplied in the request.
YEAR = 2026
LAT, LON, HEIGHT_M = 33.02, -96.70, 200.0
LMT_HOUR = 12.0
GAL_L, GAL_B = 180.0, -15.0

site = EarthLocation(lat=LAT * u.deg, lon=LON * u.deg, height=HEIGHT_M * u.m)
coe_gal = SkyCoord(l=GAL_L * u.deg, b=GAL_B * u.deg, frame="galactic")
coe = coe_gal.icrs

# Local mean solar time: UT = LMT - longitude/15 hours.
doys = np.arange(1, 366)
t0 = Time(f"{YEAR}-01-01 00:00:00", scale="utc")
times = t0 + (doys - 1 + LMT_HOUR / 24.0 - LON / 15.0 / 24.0) * u.day
altaz = coe.transform_to(AltAz(obstime=times, location=site))
az = np.asarray(altaz.az.deg)
alt = np.asarray(altaz.alt.deg)

# Month labels at the 15th, an uncluttered representative date for each month.
month_names = [calendar.month_abbr[m] for m in range(1, 13)]
month_doys = np.array([sum(calendar.monthrange(YEAR, m)[1] for m in range(1, mm)) + 15 for mm in range(1, 13)])
month_idx = month_doys - 1

# Build a daily colored path; omit the single segment that crosses the 0/360 seam.
points = np.column_stack([az, alt])
segments = np.stack([points[:-1], points[1:]], axis=1)
seam = np.abs(np.diff(az)) > 180.0
segments = segments[~seam]
seg_day = (doys[:-1] + doys[1:]) / 2.0
seg_day = seg_day[~seam]

fig, ax = plt.subplots(figsize=(12.4, 8.4), dpi=200)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Below-horizon region and horizon.
ax.axhspan(-45, 0, color="#f1f3f5", zorder=0)
ax.axhline(0, color="#30343b", lw=1.45, zorder=1)
ax.text(356.5, 1.7, "horizon", ha="right", va="bottom", fontsize=8.5, color="#50555c")

# Cardinal guides: restrained enough to keep the path dominant.
for x in (0, 90, 180, 270, 360):
    ax.axvline(x, color="#dfe3e7", lw=0.8, zorder=0)

norm = Normalize(vmin=1, vmax=365)
cmap = plt.get_cmap("turbo")
lc = LineCollection(segments, cmap=cmap, norm=norm, linewidths=2.5, alpha=0.96,
                    capstyle="round", joinstyle="round", zorder=3)
lc.set_array(seg_day)
ax.add_collection(lc)

# Daily sample dots make the once-per-day sampling explicit.
ax.scatter(az, alt, c=doys, cmap=cmap, norm=norm, s=11, linewidths=0,
           alpha=0.98, zorder=4)

# Month markers and labels. Offsets are chosen to keep all labels clear of the curve.
offsets = {
    "Jan": (-28, -18), "Feb": (-28, 10), "Mar": (-25, 10), "Apr": (-28, 9),
    "May": (-28, 10), "Jun": (10, 10), "Jul": (10, 9), "Aug": (11, 5),
    "Sep": (11, -1), "Oct": (10, -15), "Nov": (-35, -4), "Dec": (10, -7),
}
for name, idx in zip(month_names, month_idx):
    x, y = az[idx], alt[idx]
    ax.scatter([x], [y], s=68, facecolor="white", edgecolor="#20242a",
               linewidth=1.15, zorder=6)
    dx, dy = offsets[name]
    ax.annotate(
        f"{name} 15",
        (x, y), xytext=(dx, dy), textcoords="offset points",
        fontsize=9.2, fontweight="bold", color="#20242a", zorder=7,
        bbox=dict(boxstyle="round,pad=0.20", fc="white", ec="#b8bec6", lw=0.65, alpha=0.94),
        arrowprops=dict(arrowstyle="-", color="#747b84", lw=0.75, shrinkA=4, shrinkB=4),
    )

# Axes, grid, and cardinal labels.
ax.set_xlim(0, 360)
ax.set_ylim(-45, 85)
ax.set_xticks([0, 90, 180, 270, 360])
ax.set_xticklabels(["N\n0°", "E\n90°", "S\n180°", "W\n270°", "N\n360°"],
                   fontsize=10.5, fontweight="bold", color="#252a30")
ax.set_xticks(np.arange(45, 360, 45), minor=True)
ax.set_yticks(np.arange(-40, 81, 20))
ax.set_xlabel("Azimuth (clockwise from North)", fontsize=11, labelpad=13)
ax.set_ylabel("Altitude above horizon (°)", fontsize=11, labelpad=10)
ax.grid(which="major", color="#cdd2d8", lw=0.75, alpha=0.62)
ax.grid(which="minor", axis="x", color="#e5e8eb", lw=0.55, alpha=0.8)
ax.tick_params(axis="x", which="minor", length=3, color="#8d949c")
ax.tick_params(axis="both", which="major", labelsize=9.5)
for spine in ax.spines.values():
    spine.set_color("#858c94")
    spine.set_linewidth(0.8)

fig.suptitle("CoE annual sky path from Plano, Texas", fontsize=18, fontweight="bold",
             color="#171a1e", y=0.965)
fig.text(
    0.5, 0.925,
    "Daily position at 12:00 local mean solar time in 2026  •  fixed celestial direction",
    ha="center", va="center", fontsize=11.5, color="#3e444b",
)
fig.text(
    0.5, 0.895,
    f"Plano: {LAT:.2f}°N, {abs(LON):.2f}°W   •   Galactic (ℓ,b)=({GAL_L:.0f}°, {GAL_B:.0f}°)   →   "
    f"ICRS (RA, Dec)=({coe.ra.deg:.2f}°, {coe.dec.deg:+.2f}°)",
    ha="center", va="center", fontsize=9.4, color="#555c64",
)

# Color scale is a compact legend for the date progression.
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.045, pad=0.13,
                    aspect=55)
cbar.set_label("Date in 2026 (day of year)", fontsize=10, labelpad=5)
cbar.set_ticks(month_doys)
cbar.set_ticklabels(month_names)
cbar.ax.tick_params(labelsize=8.5, length=3)

fig.subplots_adjust(left=0.095, right=0.985, top=0.84, bottom=0.18)
fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Wrote {OUT}")
print(f"CoE ICRS RA={coe.ra.deg:.8f} deg, Dec={coe.dec.deg:+.8f} deg")
print(f"Altitude range: {alt.min():.3f} to {alt.max():.3f} deg")
