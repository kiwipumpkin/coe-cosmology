#!/usr/bin/env python3
"""Make a north-celestial-hemisphere map for the fixed CoE direction.

NCP-centered polar projection of the northern celestial hemisphere
(celestial equator at the rim; no southern sky).  Stick figures come from
Stellarium modern_st constellationship.fab (Hipparcos IDs) with positions
from VizieR I/239/hip_main (_RAJ2000/_DEJ2000).  CoE is the fixed ICRS
direction of Galactic (l,b)=(180°,-15°), fixed ICRS direction (no clock time — it does not move).
"""
from __future__ import annotations

from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

ROOT = Path("/workspace/figs_letter_exact_python")
OUT = ROOT / "coe_ncp_hemisphere_path.png"
DATA = Path(__file__).resolve().parent / "data"

# ---------------------------------------------------------------------------
# CoE: fixed Galactic direction → ICRS
# ---------------------------------------------------------------------------
GAL_L, GAL_B = 180.0, -15.0
coe = SkyCoord(l=GAL_L * u.deg, b=GAL_B * u.deg, frame="galactic").icrs
RA_DEG = float(coe.ra.deg)
DEC_DEG = float(coe.dec.deg)
coe_r = 90.0 - DEC_DEG

# Sidereal rotation rate relative to mean solar day
SIDEREAL_REV_PER_SOLAR_DAY = 1.00273790935

TZ = ZoneInfo("America/Chicago")
NOW = datetime.now(TZ)
NOW_LABEL = NOW.strftime("%Y-%m-%d %H:%M:%S %Z (UTC%z)").replace("UTC-0500", "UTC-5")

# Main northern constellations (reference list)
MAIN_CONST = [
    ("UMi", "Ursa Minor"),
    ("Dra", "Draco"),
    ("Cep", "Cepheus"),
    ("Cas", "Cassiopeia"),
    ("Cam", "Camelopardalis"),
    ("UMa", "Ursa Major"),
    ("Cyg", "Cygnus"),
    ("Lyr", "Lyra"),
    ("Her", "Hercules"),
    ("Boo", "Boötes"),
    ("Per", "Perseus"),
    ("Aur", "Auriga"),
    ("Lyn", "Lynx"),
    ("Tau", "Taurus"),
    ("Gem", "Gemini"),
    ("LMi", "Leo Minor"),
    ("Leo", "Leo"),
]

# Named bright stars to label (HIP → common name)
NAMED_STARS = {
    11767: r"Polaris ($\alpha$)",
    91262: r"Vega ($\alpha$)",
    102098: r"Deneb ($\alpha$)",
    69673: r"Arcturus ($\alpha$)",
    54061: r"Dubhe ($\alpha$)",
    72607: r"Kochab ($\beta$)",
    105199: r"Alderamin ($\alpha$)",
    3179: r"Schedar ($\alpha$)",
    68756: r"Thuban ($\alpha$)",
}


def load_hip_coords(path: Path) -> dict[int, tuple[float, float, float | None]]:
    out: dict[int, tuple[float, float, float | None]] = {}
    with path.open() as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            hip = int(parts[0])
            ra = float(parts[1])
            de = float(parts[2])
            vmag = float(parts[3]) if len(parts) > 3 and parts[3].strip() else None
            out[hip] = (ra, de, vmag)
    return out


def load_bright_field(path: Path) -> list[tuple[float, float, float, int]]:
    """Return list of (ra_deg, dec_deg, vmag, hip) for Dec >= 0."""
    rows: list[tuple[float, float, float, int]] = []
    with path.open() as f:
        for line in f:
            if line.startswith("#") or line.startswith("-") or not line.strip():
                continue
            if "_RAJ" in line or line.strip().startswith("deg"):
                continue
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) < 4:
                continue
            try:
                ra = float(parts[0])
                de = float(parts[1])
                hip = int(parts[2])
                vmag = float(parts[3])
            except ValueError:
                continue
            if de >= 0.0:
                rows.append((ra, de, vmag, hip))
    return rows


def load_stellarium_lines(path: Path, keep: set[str]) -> dict[str, list[tuple[int, int]]]:
    lines: dict[str, list[tuple[int, int]]] = {}
    with path.open() as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            abbr = parts[0]
            if abbr not in keep:
                continue
            nseg = int(parts[1])
            ids = [int(x) for x in parts[2:]]
            if len(ids) != 2 * nseg:
                raise ValueError(f"Bad constellationship line for {abbr}")
            lines[abbr] = [(ids[i], ids[i + 1]) for i in range(0, len(ids), 2)]
    return lines


def ra_dec_to_polar(ra_deg: float, dec_deg: float) -> tuple[float, float]:
    """Map ICRS RA/Dec to polar (theta, r) for NCP-centered looking-up view.

    Looking toward NCP from Earth: RA increases clockwise.
    r = 90° − Dec (equator at rim).
    """
    return np.deg2rad(ra_deg % 360.0), 90.0 - dec_deg


def circular_mean_deg(angles_deg: np.ndarray) -> float:
    a = np.deg2rad(angles_deg)
    return float(np.rad2deg(np.arctan2(np.mean(np.sin(a)), np.mean(np.cos(a)))) % 360.0)


def main() -> None:
    hip = load_hip_coords(DATA / "hip_coords.tsv")
    field = load_bright_field(DATA / "hip_bright_north.tsv")
    stick = load_stellarium_lines(
        DATA / "constellationship.fab", {a for a, _ in MAIN_CONST}
    )

    fig = plt.figure(figsize=(11.0, 11.6), dpi=220, facecolor="white")
    ax = fig.add_axes([0.07, 0.14, 0.86, 0.72], projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)  # RA increases clockwise (looking toward NCP)
    ax.set_rlim(0.0, 90.0)
    ax.set_rorigin(0.0)
    ax.set_facecolor("white")

    # Declination rings
    dec_ticks = [60, 30, 0]
    r_ticks = [90 - d for d in dec_ticks]
    ax.set_rticks(r_ticks)
    ax.set_rlabel_position(135)
    ax.set_yticklabels(
        [f"Dec +{d}°" if d else "Dec 0°" for d in dec_ticks],
        color="#475569",
        fontsize=8.5,
    )
    ax.grid(True, which="major", color="#94a3b8", lw=0.7, alpha=0.85)
    ax.spines["polar"].set_color("#64748b")
    ax.spines["polar"].set_linewidth(1.6)

    ra_hours = np.arange(0, 24, 2)
    ax.set_xticks(np.deg2rad(ra_hours * 15.0))
    ax.set_xticklabels(
        [f"{h}h" for h in ra_hours], fontsize=10.0, fontweight="bold", color="#0f172a"
    )
    ax.tick_params(axis="x", pad=8, colors="#0f172a")

    # Background bright stars (northern only)
    for ra, de, vmag, _hip in field:
        th, r = ra_dec_to_polar(ra, de)
        if r > 90.0:
            continue
        size = max(4.0, 55.0 * 10 ** (-0.35 * (vmag - 0.0)))
        ax.scatter(
            [th],
            [r],
            s=size,
            c="#0f172a",
            alpha=0.55,
            linewidths=0,
            zorder=2,
            rasterized=True,
        )

    # Stick figures + constellation member stars
    line_color = "#6d28d9"
    label_ras: dict[str, list[float]] = {a: [] for a, _ in MAIN_CONST}
    label_decs: dict[str, list[float]] = {a: [] for a, _ in MAIN_CONST}
    plotted_hips: set[int] = set()

    for abbr, _name in MAIN_CONST:
        segs = stick[abbr]
        for hip_a, hip_b in segs:
            if hip_a not in hip or hip_b not in hip:
                continue
            ra1, de1, v1 = hip[hip_a]
            ra2, de2, v2 = hip[hip_b]
            # Skip any segment that dips south of the equator
            if de1 < 0.0 or de2 < 0.0:
                continue
            th1, r1 = ra_dec_to_polar(ra1, de1)
            th2, r2 = ra_dec_to_polar(ra2, de2)
            # Avoid wrapping across the plot when RA span is large: draw short arc in theta
            dth = ((th2 - th1 + np.pi) % (2 * np.pi)) - np.pi
            th2_draw = th1 + dth
            ax.plot(
                [th1, th2_draw],
                [r1, r2],
                color=line_color,
                lw=1.35,
                alpha=0.95,
                solid_capstyle="round",
                zorder=4,
            )
            for h, ra, de, v in ((hip_a, ra1, de1, v1), (hip_b, ra2, de2, v2)):
                if h in plotted_hips:
                    continue
                plotted_hips.add(h)
                th, r = ra_dec_to_polar(ra, de)
                vmag = 4.5 if v is None else v
                size = max(8.0, 70.0 * 10 ** (-0.32 * vmag))
                ax.scatter(
                    [th],
                    [r],
                    s=size,
                    c="#0f172a",
                    edgecolors="#4c1d95",
                    linewidths=0.4,
                    zorder=5,
                )
                label_ras[abbr].append(ra)
                label_decs[abbr].append(de)

    # Constellation name labels near stick-figure centroids
    for abbr, name in MAIN_CONST:
        if not label_ras[abbr]:
            continue
        ra_c = circular_mean_deg(np.asarray(label_ras[abbr], dtype=float))
        de_c = float(np.mean(label_decs[abbr]))
        if de_c < 2.0:
            continue
        th, r = ra_dec_to_polar(ra_c, de_c)
        # Slight outward nudge for readability near pole / edge
        r_lab = min(86.0, max(4.0, r + (2.5 if r < 25 else 0.0)))
        ax.text(
            th,
            r_lab,
            name.upper(),
            ha="center",
            va="center",
            fontsize=8.2,
            color="#1e3a5f",
            fontweight="bold",
            alpha=0.95,
            zorder=6,
        )

    # Named bright-star labels
    for hip_id, label in NAMED_STARS.items():
        if hip_id not in hip:
            continue
        ra, de, _v = hip[hip_id]
        if de < 0.0:
            continue
        th, r = ra_dec_to_polar(ra, de)
        ax.annotate(
            label,
            xy=(th, r),
            xytext=(8, 6),
            textcoords="offset points",
            fontsize=7.2,
            color="#5b4636",
            zorder=7,
        )

    # CoE constant-Dec diurnal parallel
    theta = np.linspace(0.0, 2.0 * np.pi, 721)
    path_color = "#0284c7"
    ax.plot(
        theta,
        np.full_like(theta, coe_r),
        color=path_color,
        lw=2.6,
        alpha=0.95,
        solid_capstyle="round",
        zorder=8,
    )


    # Fixed CoE point
    coe_theta = np.deg2rad(RA_DEG)
    ax.scatter(
        [coe_theta],
        [coe_r],
        s=220,
        marker="o",
        color="#ff5a36",
        edgecolor="#0f172a",
        linewidth=1.8,
        zorder=11,
    )
    ax.scatter(
        [coe_theta], [coe_r], s=48, marker="o", color="#1a0d0a", linewidth=0, zorder=12
    )

    # NCP marker
    ax.scatter([0], [0], s=36, c="#b45309", marker="x", linewidths=1.6, zorder=10)
    ax.text(
        0.0,
        3.5,
        "NCP +90°",
        ha="center",
        va="bottom",
        fontsize=9.5,
        fontweight="bold",
        color="#b45309",
        zorder=10,
    )

    # Legend / notes
    legend_handles = [
        Line2D([0], [0], color=line_color, lw=1.6, label="Constellation stick figures"),
        Line2D(
            [0],
            [0],
            color=path_color,
            lw=2.4,
            label=f"CoE Dec parallel · Dec {DEC_DEG:+.2f}°",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#ff5a36",
            markeredgecolor="#0f172a",
            markeredgewidth=1.1,
            markersize=9,
            label=f"CoE fixed ICRS · RA {RA_DEG:.2f}° ({RA_DEG / 15:.2f}h)",
        ),
    ]
    leg = ax.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        frameon=False,
        fontsize=8.6,
        ncol=1,
        handlelength=2.4,
        labelcolor="#0f172a",
    )

    fig.savefig(OUT, dpi=220, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)

    from PIL import Image

    with Image.open(OUT) as im:
        print(f"Wrote {OUT}")
        print(f"Pixel size: {im.width} x {im.height}")
    print(f"Timestamp: {NOW_LABEL}")
    print(f"CoE ICRS RA={RA_DEG:.8f} deg ({RA_DEG / 15:.8f} h), Dec={DEC_DEG:+.8f} deg")
    print(f"CoE diurnal parallel radius from NCP={coe_r:.8f} deg")
    print(
        f"Constellations: {', '.join(n for _, n in MAIN_CONST)} "
        f"({sum(len(stick[a]) for a,_ in MAIN_CONST)} segments, {len(plotted_hips)} HIP stars)"
    )
    print(
        "Sources: Stellarium modern_st constellationship.fab @ eb47095; "
        "VizieR I/239/hip_main; astropy SkyCoord galactic→icrs"
    )


if __name__ == "__main__":
    main()
