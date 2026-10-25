#!/usr/bin/env python3
"""
Clean 2D Galactic Mollweide maps of Migkas et al. 2020 Table 2 A extrema.

Deliverables (this directory):
  migkas_A_extrema_table2_mollweide.png / .pdf
  migkas_fig8_theta75_A_extrema_only.png
  migkas_A_extrema_table2.json

Convention (match Migkas Fig. 8 / existing overlay):
  lon_disp = -wrap(ℓ)  →  ℓ=0 at center, ℓ increases LEFT
  (ℓ=90 left of center, ℓ=270 right of center)

No CoE rings, no model curves — published Table 2 directions only.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, FancyBboxPatch
from PIL import Image

OUT_DIR = Path(__file__).resolve().parent
DPI = 220

# ---------------------------------------------------------------------------
# Migkas Table 2 — Galactic (ℓ, b) degrees. Authoritative from user brief.
# ---------------------------------------------------------------------------

PANEL_A = [
    # Our sample All θ cones (primary paper maps)
    {
        "id": "all_90",
        "label": "90°",
        "sample": "Our sample All 90°",
        "panel": "A",
        "low": {"ell": 272.0, "b": -8.0},
        "high": {"ell": 47.0, "b": 22.0},
        "dipole": {"ell": 230.0, "b": -20.0},
    },
    {
        "id": "all_75",
        "label": "75°",
        "sample": "Our sample All 75° (Fig. 8 primary)",
        "panel": "A",
        "low": {"ell": 274.0, "b": -22.0, "A": 0.999, "A_err": 0.050},
        "high": {"ell": 17.0, "b": -9.0, "A": 1.288, "A_err": 0.061},
        "dipole": {"ell": 263.0, "b": -21.0},
        "primary": True,
    },
    {
        "id": "all_60",
        "label": "60°",
        "sample": "Our sample All 60°",
        "panel": "A",
        "low": {"ell": 281.0, "b": -16.0},
        "high": {"ell": 34.0, "b": 4.0},
        "dipole": {"ell": 260.0, "b": -36.0},
    },
    {
        "id": "all_45",
        "label": "45°",
        "sample": "Our sample All 45°",
        "panel": "A",
        "low": {"ell": 280.0, "b": 1.0},
        "high": {"ell": 32.0, "b": 14.0},
        "dipole": {"ell": 255.0, "b": -53.0},
    },
]

PANEL_B = [
    {
        "id": "all_75_freeB_c1",
        "label": "75° free B C1",
        "sample": "All 75° free B Case 1",
        "panel": "B",
        "low": {"ell": 272.0, "b": -21.0},
        "high": {"ell": 26.0, "b": -13.0},
    },
    {
        "id": "all_75_freeB_c2",
        "label": "75° free B C2",
        "sample": "All 75° free B Case 2",
        "panel": "B",
        "low": {"ell": 269.0, "b": -17.0},
        "high": {"ell": 23.0, "b": -10.0},
    },
    {
        "id": "T2p5_z003_75",
        "label": "T>2.5 z>0.03",
        "sample": "T>2.5 keV z>0.03 (75°)",
        "panel": "B",
        "low": {"ell": 288.0, "b": -35.0},
        "high": {"ell": 10.0, "b": 16.0},
    },
    {
        "id": "T3_z005_75",
        "label": "T>3 z>0.05",
        "sample": "T>3 keV z>0.05 (75°)",
        "panel": "B",
        "low": {"ell": 286.0, "b": -36.0},
        "high": {"ell": 9.0, "b": 15.0},
    },
    {
        "id": "Zcore_75",
        "label": "Zcore≤0.59",
        "sample": "Zcore≤0.59 (75°)",
        "panel": "B",
        "low": {"ell": 264.0, "b": -18.0},
        "high": {"ell": 30.0, "b": 23.0},
    },
    {
        "id": "Zout_75",
        "label": "Zout≤0.426",
        "sample": "Zout≤0.426 (75°)",
        "panel": "B",
        "low": {"ell": 270.0, "b": -14.0},
        "high": {"ell": 24.0, "b": 15.0},
    },
    {
        "id": "NH_75",
        "label": "NH≤5.16e20",
        "sample": "NH≤5.16e20 (75°)",
        "panel": "B",
        "low": {"ell": 242.0, "b": -27.0},
        "high": {"ell": 35.0, "b": -15.0},
    },
    {
        "id": "ACC_75",
        "label": "ACC 75°",
        "sample": "ACC 75°",
        "panel": "B",
        "low": {"ell": 314.0, "b": -17.0},
        "high": {"ell": 77.0, "b": 15.0},
    },
    {
        "id": "XCS_75",
        "label": "XCS-DR1 75°",
        "sample": "XCS-DR1 75°",
        "panel": "B",
        "low": {"ell": 281.0, "b": 24.0},
        "high": {"ell": 31.0, "b": 25.0},
    },
    {
        "id": "joint_ACC_XCS_75",
        "label": "Joint+ACC+XCS",
        "sample": "Joint+ACC+XCS 75°",
        "panel": "B",
        "low": {"ell": 303.0, "b": -27.0},
        "high": {"ell": 34.0, "b": 26.0},
        "dipole": {"ell": 265.0, "b": -20.0},
    },
    {
        "id": "joint_freeB_c1",
        "label": "Joint free B C1",
        "sample": "Joint free B Case 1",
        "panel": "B",
        "low": {"ell": 312.0, "b": -21.0},
        "high": {"ell": 45.0, "b": 21.0},
    },
]

ALL_ROWS = PANEL_A + PANEL_B

# Colors
C_HIGH = "#c0392b"   # warm red — highest A
C_LOW = "#2471a3"    # cool blue — lowest A
C_DIP = "#b7950b"    # gold — dipole
C_ARC = "#7f8c8d"    # gray arcs

# Fig. 8 registration (from make_migkas_fig8_coe_rings_overlay.py)
GOOD_REG = {"cx": 863.0, "cy": 763.25, "Rx": 693.9927, "Ry": 392.542425}
SOURCE_CANDIDATES = [
    Path("/workspace/cosmology/map_compare/FULL_migkas_theta75.png"),
    Path("/workspace/cosmology/map_compare/migkas_fig8_theta75_crop.png"),
    Path("/workspace/cosmology/migkas/pub_maps/Fig8_full_hires.png"),
]


# ---------------------------------------------------------------------------
# Coordinate helpers (same as Fig. 8 overlay)
# ---------------------------------------------------------------------------

def wrap_lon_deg(lon):
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def lon_disp_from_ell(ell_deg):
    """Display longitude: ℓ=0 center, ℓ increases LEFT (Migkas Fig. 8)."""
    return -wrap_lon_deg(ell_deg)


def mollweide_xy(lon_rad, lat_rad, n_iter=64):
    lon = np.asarray(lon_rad, dtype=float)
    lat = np.asarray(lat_rad, dtype=float)
    target = np.pi * np.sin(np.clip(lat, -np.pi / 2, np.pi / 2))
    lo = np.full_like(lat, -np.pi / 2, dtype=float)
    hi = np.full_like(lat, np.pi / 2, dtype=float)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        g = 2.0 * mid + np.sin(2.0 * mid)
        lo = np.where(g < target, mid, lo)
        hi = np.where(g >= target, mid, hi)
    theta = 0.5 * (lo + hi)
    theta = np.where(lat >= np.pi / 2 - 1e-15, np.pi / 2, theta)
    theta = np.where(lat <= -np.pi / 2 + 1e-15, -np.pi / 2, theta)
    x = (2.0 * np.sqrt(2.0) / np.pi) * lon * np.cos(theta)
    y = np.sqrt(2.0) * np.sin(theta)
    return x, y


def lb_to_px(ell_deg, b_deg, cx, cy, Rx, Ry):
    lon = np.deg2rad(lon_disp_from_ell(np.asarray(ell_deg, dtype=float)))
    lat = np.deg2rad(np.asarray(b_deg, dtype=float))
    x, y = mollweide_xy(lon, lat)
    Xn = x / (2.0 * np.sqrt(2.0))
    Yn = y / np.sqrt(2.0)
    return cx + Xn * Rx, cy - Yn * Ry


def lb_to_n(l_deg, b_deg):
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    return np.array(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], dtype=float
    )


def n_to_lb(n):
    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    b = float(np.rad2deg(np.arcsin(np.clip(n[2], -1.0, 1.0))))
    l = float(np.rad2deg(np.arctan2(n[1], n[0])) % 360.0)
    return l, b


def great_circle_lb(ell1, b1, ell2, b2, n=181):
    """Short great-circle arc from (ℓ1,b1) to (ℓ2,b2)."""
    n1 = lb_to_n(ell1, b1)
    n2 = lb_to_n(ell2, b2)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    cos_ang = float(np.clip(np.dot(n1, n2), -1.0, 1.0))
    ang = np.arccos(cos_ang)
    if ang < 1e-8:
        return np.array([ell1]), np.array([b1])
    t = np.linspace(0.0, 1.0, n)
    # slerp
    pts = (
        (np.sin((1 - t) * ang)[:, None] * n1 + np.sin(t * ang)[:, None] * n2)
        / np.sin(ang)
    )
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    return ls, bs


# ---------------------------------------------------------------------------
# Matplotlib Mollweide helpers
# ---------------------------------------------------------------------------

def to_moll_rad(ell_deg, b_deg):
    """Return (lon_disp_rad, lat_rad) for ax.plot on projection='mollweide'."""
    lon = np.deg2rad(lon_disp_from_ell(ell_deg))
    lat = np.deg2rad(b_deg)
    return lon, lat


def plot_gc_arc(ax, ell1, b1, ell2, b2, **kw):
    ls, bs = great_circle_lb(ell1, b1, ell2, b2, n=241)
    lon, lat = to_moll_rad(ls, bs)
    # Split on antimeridian jumps in display longitude
    dlon = np.diff(lon)
    cuts = np.where(np.abs(dlon) > np.pi)[0] + 1
    idx = np.concatenate([[0], cuts, [len(lon)]])
    for i0, i1 in zip(idx[:-1], idx[1:]):
        if i1 - i0 >= 2:
            ax.plot(lon[i0:i1], lat[i0:i1], **kw)


def setup_moll_ax(ax, title=None):
    ax.set_longitude_grid(45)
    ax.set_latitude_grid(30)
    ax.set_longitude_grid_ends(90)
    ax.grid(True, linestyle=":", color="#bbbbbb", linewidth=0.6, alpha=0.9)
    # Tick labels: show Galactic ℓ at display positions
    # matplotlib mollweide shows display lon; we remapped so labels need care.
    # Set custom tick labels for longitude every 45°.
    tick_ells = [180, 135, 90, 45, 0, 315, 270, 225]  # left → right on map
    tick_disp = [lon_disp_from_ell(e) for e in tick_ells]
    ax.set_xticks(np.deg2rad(tick_disp))
    ax.set_xticklabels([rf"${e}^\circ$" for e in tick_ells], fontsize=8)
    ax.set_yticks(np.deg2rad([-60, -30, 0, 30, 60]))
    ax.set_yticklabels(
        [r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$+30^\circ$", r"$+60^\circ$"],
        fontsize=8,
    )
    if title:
        ax.set_title(title, fontsize=11, pad=10)


def annotate_point(ax, ell, b, text, color, offset=(6, 6), fontsize=7.5, weight="normal"):
    lon, lat = to_moll_rad(ell, b)
    ax.annotate(
        text,
        xy=(lon, lat),
        xytext=offset,
        textcoords="offset points",
        fontsize=fontsize,
        color=color,
        fontweight=weight,
        ha="left",
        va="bottom",
        bbox=dict(
            boxstyle="round,pad=0.15",
            facecolor="white",
            edgecolor="none",
            alpha=0.75,
        ),
        zorder=10,
    )


# ---------------------------------------------------------------------------
# Figure 1: two-panel Table 2 Mollweide
# ---------------------------------------------------------------------------

def draw_panel_a(ax):
    setup_moll_ax(
        ax,
        title=(
            r"Panel A — Our sample All $\theta$ cones "
            r"($90^\circ$, $75^\circ$, $60^\circ$, $45^\circ$)"
        ),
    )
    # Distinct sizes / slight offset labels per cone
    offsets = {
        "90°": {"high": (8, 8), "low": (8, -12)},
        "75°": {"high": (8, -14), "low": (-55, -14)},
        "60°": {"high": (8, 10), "low": (8, 8)},
        "45°": {"high": (-50, 8), "low": (8, 10)},
    }
    for row in PANEL_A:
        lo, hi = row["low"], row["high"]
        lab = row["label"]
        # pair arc
        plot_gc_arc(
            ax, lo["ell"], lo["b"], hi["ell"], hi["b"],
            color=C_ARC, lw=0.9, ls="--", alpha=0.65, zorder=2,
        )
        # low
        lon_l, lat_l = to_moll_rad(lo["ell"], lo["b"])
        ax.plot(
            lon_l, lat_l, marker="v", ms=9, color=C_LOW,
            markeredgecolor="k", markeredgewidth=0.6,
            linestyle="None", zorder=5,
        )
        # high
        lon_h, lat_h = to_moll_rad(hi["ell"], hi["b"])
        ax.plot(
            lon_h, lat_h, marker="*", ms=13, color=C_HIGH,
            markeredgecolor="k", markeredgewidth=0.55,
            linestyle="None", zorder=6,
        )
        # dipole (smaller)
        if "dipole" in row:
            d = row["dipole"]
            lon_d, lat_d = to_moll_rad(d["ell"], d["b"])
            ax.plot(
                lon_d, lat_d, marker="o", ms=4.5, color=C_DIP,
                markeredgecolor="k", markeredgewidth=0.4,
                linestyle="None", zorder=4, alpha=0.9,
            )
        off = offsets.get(lab, {"high": (6, 6), "low": (6, -10)})
        weight = "bold" if row.get("primary") else "normal"
        annotate_point(
            ax, hi["ell"], hi["b"], f"{lab} max", C_HIGH,
            offset=off["high"], fontsize=7.5 if not row.get("primary") else 8.5,
            weight=weight,
        )
        annotate_point(
            ax, lo["ell"], lo["b"], f"{lab} min", C_LOW,
            offset=off["low"], fontsize=7.5 if not row.get("primary") else 8.5,
            weight=weight,
        )


def draw_panel_b(ax):
    setup_moll_ax(
        ax,
        title=(
            "Panel B — Other Table 2 rows "
            "(subsamples + ACC / XCS / Joint)"
        ),
    )
    # Compact markers; legend identifies samples rather than crowding labels
    for i, row in enumerate(PANEL_B):
        lo, hi = row["low"], row["high"]
        lon_l, lat_l = to_moll_rad(lo["ell"], lo["b"])
        lon_h, lat_h = to_moll_rad(hi["ell"], hi["b"])
        ax.plot(
            lon_l, lat_l, marker="v", ms=7, color=C_LOW,
            markeredgecolor="k", markeredgewidth=0.45,
            linestyle="None", zorder=5, alpha=0.95,
        )
        ax.plot(
            lon_h, lat_h, marker="*", ms=10, color=C_HIGH,
            markeredgecolor="k", markeredgewidth=0.45,
            linestyle="None", zorder=6, alpha=0.95,
        )
        if "dipole" in row:
            d = row["dipole"]
            lon_d, lat_d = to_moll_rad(d["ell"], d["b"])
            ax.plot(
                lon_d, lat_d, marker="o", ms=3.5, color=C_DIP,
                markeredgecolor="k", markeredgewidth=0.35,
                linestyle="None", zorder=4, alpha=0.85,
            )
        # Tiny numeric index near high-A for legend cross-ref
        ax.annotate(
            str(i + 1),
            xy=(lon_h, lat_h),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=6.5,
            color="#333333",
            ha="left",
            va="bottom",
            zorder=10,
            bbox=dict(boxstyle="round,pad=0.1", facecolor="white", edgecolor="none", alpha=0.7),
        )


def make_table2_mollweide():
    fig = plt.figure(figsize=(11.5, 12.0), dpi=DPI)
    ax_a = fig.add_subplot(2, 1, 1, projection="mollweide")
    ax_b = fig.add_subplot(2, 1, 2, projection="mollweide")
    draw_panel_a(ax_a)
    draw_panel_b(ax_b)

    # Shared legend (below panel B)
    handles = [
        Line2D(
            [0], [0], marker="*", color="w", markerfacecolor=C_HIGH,
            markeredgecolor="k", linestyle="None", ms=12,
            label="Highest A (max)",
        ),
        Line2D(
            [0], [0], marker="v", color="w", markerfacecolor=C_LOW,
            markeredgecolor="k", linestyle="None", ms=9,
            label="Lowest A (min)",
        ),
        Line2D(
            [0], [0], marker="o", color="w", markerfacecolor=C_DIP,
            markeredgecolor="k", linestyle="None", ms=6,
            label="Dipole direction (where listed)",
        ),
        Line2D(
            [0], [0], color=C_ARC, lw=1.0, ls="--",
            label=r"Great-circle arc linking All-$\theta$ min–max pair (Panel A)",
        ),
    ]
    leg1 = fig.legend(
        handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.055),
        ncol=4, fontsize=8.5, frameon=True, fancybox=False,
        edgecolor="#333333", framealpha=1.0, facecolor="white",
    )
    leg1.get_frame().set_linewidth(0.8)

    # Panel B numbered legend
    b_lines = []
    for i, row in enumerate(PANEL_B):
        b_lines.append(f"{i + 1}. {row['sample']}")
    # Two-column text block
    col1 = "\n".join(b_lines[:6])
    col2 = "\n".join(b_lines[6:])
    fig.text(
        0.08, 0.012, col1, fontsize=7.0, family="sans-serif",
        va="bottom", ha="left", color="#222222",
        linespacing=1.35,
    )
    fig.text(
        0.52, 0.012, col2, fontsize=7.0, family="sans-serif",
        va="bottom", ha="left", color="#222222",
        linespacing=1.35,
    )
    fig.text(
        0.08, 0.048,
        "Panel B index (number next to each high-A marker):",
        fontsize=7.5, fontweight="bold", color="#222222",
    )

    fig.suptitle(
        "Migkas et al. 2020 — published Table 2 $A$ extrema on the Galactic sky\n"
        r"(Galactic Mollweide; $\ell\!=\!0$ center, $\ell$ increases LEFT — Fig. 8 convention; "
        "not a CoE / model map)",
        fontsize=12.5,
        y=0.985,
    )
    fig.subplots_adjust(left=0.04, right=0.98, top=0.91, bottom=0.12, hspace=0.28)

    png = OUT_DIR / "migkas_A_extrema_table2_mollweide.png"
    pdf = OUT_DIR / "migkas_A_extrema_table2_mollweide.pdf"
    fig.savefig(png, dpi=DPI, facecolor="white", bbox_inches="tight")
    fig.savefig(pdf, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {png}")
    print(f"wrote {pdf}")
    return png, pdf


# ---------------------------------------------------------------------------
# Figure 2: Fig. 8 θ=75° overlay — published max/min/dipole only
# ---------------------------------------------------------------------------

def choose_source():
    for p in SOURCE_CANDIDATES:
        if p.is_file():
            return p
    raise FileNotFoundError("No Migkas Fig.8 theta=75 source image found")


def make_fig8_overlay():
    source = choose_source()
    rgb = np.asarray(Image.open(source).convert("RGB"))
    h, w, _ = rgb.shape
    cx, cy, Rx, Ry = GOOD_REG["cx"], GOOD_REG["cy"], GOOD_REG["Rx"], GOOD_REG["Ry"]

    # Spot-check: 75° max should land in warmer pixels than 75° min
    def sample(ell, b, r=5):
        px, py = lb_to_px(ell, b, cx, cy, Rx, Ry)
        xi, yi = int(round(float(px))), int(round(float(py)))
        patch = rgb[max(0, yi - r): yi + r + 1, max(0, xi - r): xi + r + 1]
        return patch.mean(axis=(0, 1)), (xi, yi)

    cmax, pmax = sample(17.0, -9.0)
    cmin, pmin = sample(274.0, -22.0)
    cdip, pdip = sample(263.0, -21.0)
    warm_ok = (cmax[0] + cmax[1]) > (cmin[0] + cmin[1]) + 40
    print(
        f"Fig8 source={source.name} size={w}x{h}; "
        f"max px={pmax} RGB={cmax.round(1)}; "
        f"min px={pmin} RGB={cmin.round(1)}; "
        f"dipole px={pdip}; warm_check={'OK' if warm_ok else 'WARN'}"
    )

    dpi = 140
    fig, ax = plt.subplots(figsize=(w / dpi, (h / dpi) + 0.95), dpi=dpi)
    ax.imshow(rgb, origin="upper", extent=(0, w, h, 0), interpolation="nearest")
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    clip_oval = Ellipse(
        (cx, cy), width=2.0 * Rx, height=2.0 * Ry,
        facecolor="none", edgecolor="none", transform=ax.transData,
    )
    ax.add_patch(clip_oval)

    def mark(ell, b, marker, color, ms, label, offset):
        px, py = lb_to_px(ell, b, cx, cy, Rx, Ry)
        px, py = float(px), float(py)
        ax.plot(
            px, py, marker=marker, ms=ms + 2.2, color="white",
            linestyle="None", zorder=6, clip_path=clip_oval, clip_on=True,
        )
        ax.plot(
            px, py, marker=marker, ms=ms, color=color,
            markeredgecolor="k", markeredgewidth=0.7,
            linestyle="None", zorder=7, clip_path=clip_oval, clip_on=True,
        )
        ax.annotate(
            label,
            xy=(px, py),
            xytext=offset,
            textcoords="offset points",
            fontsize=10,
            fontweight="bold",
            color=color,
            ha="left",
            va="center",
            bbox=dict(
                boxstyle="round,pad=0.25",
                facecolor="white",
                edgecolor=color,
                alpha=0.92,
                linewidth=1.0,
            ),
            zorder=8,
            clip_on=False,
        )

    # Optional thin dashed great-circle between min and max (in image pixels)
    ls, bs = great_circle_lb(274.0, -22.0, 17.0, -9.0, n=301)
    px_arc, py_arc = lb_to_px(ls, bs, cx, cy, Rx, Ry)
    px_arc = np.asarray(px_arc, dtype=float)
    py_arc = np.asarray(py_arc, dtype=float)
    # Split on large pixel jumps (antimeridian)
    segs_x, segs_y = [], []
    cx_list, cy_list = [], []
    for x, y in zip(px_arc, py_arc):
        # keep inside oval
        Xn = (x - cx) / Rx
        Yn = (y - cy) / Ry
        if Xn * Xn + Yn * Yn > 0.98:
            if len(cx_list) >= 2:
                segs_x.append(cx_list)
                segs_y.append(cy_list)
            cx_list, cy_list = [], []
            continue
        if cx_list and np.hypot(x - cx_list[-1], y - cy_list[-1]) > 40:
            if len(cx_list) >= 2:
                segs_x.append(cx_list)
                segs_y.append(cy_list)
            cx_list, cy_list = [x], [y]
        else:
            cx_list.append(x)
            cy_list.append(y)
    if len(cx_list) >= 2:
        segs_x.append(cx_list)
        segs_y.append(cy_list)
    for sx, sy in zip(segs_x, segs_y):
        ax.plot(
            sx, sy, color="white", lw=2.0, ls="--", alpha=0.85, zorder=4,
            clip_path=clip_oval, clip_on=True,
        )
        ax.plot(
            sx, sy, color="#555555", lw=1.1, ls="--", alpha=0.9, zorder=5,
            clip_path=clip_oval, clip_on=True,
        )

    mark(17.0, -9.0, "*", C_HIGH, 16, r"$75^\circ$ max $(17^\circ,-9^\circ)$", (14, 10))
    mark(274.0, -22.0, "v", C_LOW, 12, r"$75^\circ$ min $(274^\circ,-22^\circ)$", (14, -8))
    mark(263.0, -21.0, "o", C_DIP, 8, r"dipole $(263^\circ,-21^\circ)$", (14, 18))

    handles = [
        Line2D(
            [0], [0], marker="*", color="w", markerfacecolor=C_HIGH,
            markeredgecolor="k", linestyle="None", ms=14,
            label=r"Highest $A$: $(17^\circ,-9^\circ)$, $A=1.288\pm0.061$",
        ),
        Line2D(
            [0], [0], marker="v", color="w", markerfacecolor=C_LOW,
            markeredgecolor="k", linestyle="None", ms=10,
            label=r"Lowest $A$: $(274^\circ,-22^\circ)$, $A=0.999\pm0.050$",
        ),
        Line2D(
            [0], [0], marker="o", color="w", markerfacecolor=C_DIP,
            markeredgecolor="k", linestyle="None", ms=7,
            label=r"Dipole: $(263^\circ,-21^\circ)$",
        ),
    ]
    leg = ax.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.01),
        ncol=3, fontsize=9, frameon=True, fancybox=False,
        edgecolor="#222222", framealpha=1.0, facecolor="white",
    )
    leg.get_frame().set_linewidth(1.0)

    fig.suptitle(
        r"Migkas et al. 2020 Fig. 8 — $\theta=75^\circ$ $A/A_{\rm all}$ map "
        r"with published Table 2 extrema only"
        "\n"
        r"(no CoE rings / model; registration: lon$_{\rm disp}=-\mathrm{wrap}(\ell)$, "
        r"GOOD_REG from working Fig. 8 overlay)",
        fontsize=11,
        y=0.99,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.91, bottom=0.10)
    out = OUT_DIR / "migkas_fig8_theta75_A_extrema_only.png"
    fig.savefig(out, dpi=dpi, facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")
    return out, {
        "source_image": str(source),
        "registration": dict(GOOD_REG),
        "pixel_checks": {
            "max_75": {"ell": 17.0, "b": -9.0, "px": list(pmax), "rgb_mean": [float(x) for x in cmax]},
            "min_75": {"ell": 274.0, "b": -22.0, "px": list(pmin), "rgb_mean": [float(x) for x in cmin]},
            "dipole_75": {"ell": 263.0, "b": -21.0, "px": list(pdip)},
            "max_warmer_than_min": bool(warm_ok),
        },
    }


# ---------------------------------------------------------------------------
# JSON
# ---------------------------------------------------------------------------

def write_json(fig8_meta):
    points = []
    for row in ALL_ROWS:
        entry = {
            "id": row["id"],
            "sample": row["sample"],
            "panel": row["panel"],
            "label": row["label"],
            "low_A": {
                "ell_deg": row["low"]["ell"],
                "b_deg": row["low"]["b"],
            },
            "high_A": {
                "ell_deg": row["high"]["ell"],
                "b_deg": row["high"]["b"],
            },
        }
        if "A" in row["low"]:
            entry["low_A"]["A"] = row["low"]["A"]
            entry["low_A"]["A_err"] = row["low"]["A_err"]
        if "A" in row["high"]:
            entry["high_A"]["A"] = row["high"]["A"]
            entry["high_A"]["A_err"] = row["high"]["A_err"]
        if "dipole" in row:
            entry["dipole"] = {
                "ell_deg": row["dipole"]["ell"],
                "b_deg": row["dipole"]["b"],
            }
        if row.get("primary"):
            entry["primary"] = True
        points.append(entry)

    meta = {
        "source": "Migkas et al. 2020, Table 2 (published A extrema directions)",
        "coordinate_system": "Galactic (ℓ, b) degrees",
        "mollweide_convention": (
            "lon_disp = -wrap(ℓ); ℓ=0 at center; ℓ increases LEFT "
            "(matches Migkas Fig. 8)"
        ),
        "note": (
            "These are published Migkas 2D map extrema from Table 2 — "
            "not CoE model rings or Amp×sin curves."
        ),
        "style": {
            "high_A": "red star",
            "low_A": "blue inverted triangle",
            "dipole": "gold circle (where listed)",
        },
        "points": points,
        "outputs": {
            "mollweide_png": str(OUT_DIR / "migkas_A_extrema_table2_mollweide.png"),
            "mollweide_pdf": str(OUT_DIR / "migkas_A_extrema_table2_mollweide.pdf"),
            "fig8_overlay_png": str(OUT_DIR / "migkas_fig8_theta75_A_extrema_only.png"),
        },
        "fig8_overlay": fig8_meta,
        "script": str(Path(__file__).resolve()),
    }
    path = OUT_DIR / "migkas_A_extrema_table2.json"
    path.write_text(json.dumps(meta, indent=2))
    print(f"wrote {path}")
    return path


def main():
    png, pdf = make_table2_mollweide()
    fig8_png, fig8_meta = make_fig8_overlay()
    jpath = write_json(fig8_meta)
    print("--- done ---")
    for p in (png, pdf, fig8_png, jpath):
        print(p)


if __name__ == "__main__":
    main()
