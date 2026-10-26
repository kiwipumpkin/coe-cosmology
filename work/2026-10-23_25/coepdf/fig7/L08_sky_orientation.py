# PROVENANCE COPY of /workspace/figs_from_makers/makers/make_fig05_sky_orientation.py (original generating script, rev6 figure; run in its original directory).
# Output figure copied here as L08_sky_orientation.jpg. Not re-run for rev7: see Appendix B.
#!/usr/bin/env python3
"""Fig05 — Four-row sky orientation (untilted → tilt → gradient → roll).

Regenerates analytic band circles; writes JPEG to match letter delivery format.
Gold target ~1370×1420.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths import OUT  # noqa: E402
from physics import (  # noqa: E402
    BAND_CENTRES,
    BAND_COLORS,
    BAND_HALF,
    BAND_KINDS,
    ROLL_DEG,
    TILT_DEG,
    band_circle_xyz,
    mollweide_xy,
    rotate_xyz,
    xyz_to_galactic,
)

OUT_NAME = "Fig05_sky_orientation_four_row.jpg"


def _draw_model_frame(ax, tilt, roll, gradient=False):
    ax.set_aspect("equal")
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.axis("off")
    ax.add_patch(Circle((0, 0), 1.0, fill=False, ec="0.4", lw=1.4))
    ax.plot(0, 0, "ko", ms=5)
    ax.annotate("", xy=(1.15, 0), xytext=(-1.15, 0),
                arrowprops=dict(arrowstyle="<->", color="#6c3483", lw=1.5))
    ax.text(-1.28, 0.12, "CoE (−x)", color="#6c3483", fontsize=7, ha="left")
    ax.text(1.28, 0.12, "anti-CoE (+x)", color="#6c3483", fontsize=7, ha="right")
    for th, col, kind in zip(BAND_CENTRES, BAND_COLORS, BAND_KINDS):
        t = np.radians(th)
        alpha = 0.35
        if gradient:
            dist = min(abs(th - 180), abs(th - 180 + 360), abs(th - 180 - 360))
            alpha = 0.15 + 0.55 * (1 - dist / 180.0)
        for sgn in (+1, -1):
            x = np.cos(t)
            z = sgn * np.sin(t)
            tt = np.radians(tilt)
            xr = x * np.cos(tt) - z * np.sin(tt)
            zr = x * np.sin(tt) + z * np.cos(tt)
            rr = np.radians(roll)
            zr2 = zr * np.cos(rr)
            ax.plot([0.82 * xr, 1.05 * xr], [0.82 * zr2, 1.05 * zr2],
                    color=col, lw=4.5, alpha=alpha, solid_capstyle="round")
        xl = 1.18 * np.cos(t)
        zl = 1.18 * np.sin(t)
        ax.text(xl, zl, rf"$\theta$={th:g}°", color=col, fontsize=6.5,
                ha="center", va="center", alpha=min(1, alpha + 0.3))


def _draw_mollweide(ax, tilt, roll, gradient=False):
    t = np.linspace(0, 2 * np.pi, 400)
    ox = 2 * np.sqrt(2) * np.cos(t)
    oy = np.sqrt(2) * np.sin(t)
    ax.plot(-ox, oy, "k-", lw=0.8)
    lon = np.linspace(0, 360, 400)
    x, y = mollweide_xy(lon, np.zeros_like(lon))
    ax.plot(x, y, color="#1a7f37", lw=0.9, ls="--", alpha=0.8)
    for th, col, kind in zip(BAND_CENTRES, BAND_COLORS, BAND_KINDS):
        band = band_circle_xyz(th, BAND_HALF, n=720)
        for edge_name, lw in (("mid", 2.2), ("lo", 0.9), ("hi", 0.9)):
            xyz = rotate_xyz(band[edge_name], tilt_deg=tilt, roll_deg=roll)
            ell, b = xyz_to_galactic(xyz)
            x, y = mollweide_xy(ell, b)
            d = np.hypot(np.diff(x), np.diff(y))
            breaks = np.where(d > 0.35)[0]
            idx = np.r_[0, breaks + 1, len(x)]
            for i0, i1 in zip(idx[:-1], idx[1:]):
                if i1 - i0 < 2:
                    continue
                xs, ys = x[i0:i1], y[i0:i1]
                if gradient:
                    mid = ell[i0:i1].mean()
                    d0 = min(mid % 360, 360 - (mid % 360))
                    alpha = 0.2 + 0.75 * (1 - d0 / 180.0)
                else:
                    alpha = 0.9 if edge_name == "mid" else 0.45
                ax.plot(xs, ys, color=col, lw=lw if edge_name == "mid" else 0.8, alpha=alpha)
    ac = rotate_xyz(np.array([[1.0, 0.0, 0.0]]), tilt_deg=tilt, roll_deg=roll)
    ell0, b0 = xyz_to_galactic(ac)
    x0, y0 = mollweide_xy(ell0, b0)
    ax.plot(x0, y0, "o", color="#7b241c", ms=6, zorder=5)
    ax.text(x0[0], y0[0] + 0.12, "anti-CoE", color="#7b241c", fontsize=7, ha="center")
    ax.set_aspect("equal")
    ax.set_xlim(-3.0, 3.0)
    ax.set_ylim(-1.7, 1.7)
    ax.axis("off")
    for L, xpos in ((180, -2.0), (270, -1.0), (0, 0.0), (90, 1.0), (180, 2.0)):
        ax.text(xpos, -1.55, "180" if L == 180 else str(L), fontsize=7, ha="center", color="0.35")
    ax.text(0, -1.72, r"Galactic $\ell$ (increases ←)", fontsize=7, ha="center", color="0.35")


def _draw_sphere3d(ax, tilt, roll, gradient=False):
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 20)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, color="0.75", lw=0.3, alpha=0.5)
    ax.plot([-1.2, 1.2], [0, 0], [0, 0], color="#6c3483", lw=1.5)
    for th, col in zip(BAND_CENTRES, BAND_COLORS):
        band = band_circle_xyz(th, 0.0, n=200)["mid"]
        xyz = rotate_xyz(band, tilt_deg=tilt, roll_deg=roll)
        alpha = 0.9
        if gradient:
            alpha = 0.95 if th <= 90 else 0.35
        ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=col, lw=2.0, alpha=alpha)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(-1.2, 1.2)
    ax.set_axis_off()
    ax.view_init(elev=18, azim=35)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [
        (0.0, 0.0, False, "Untilted, no roll"),
        (TILT_DEG, 0.0, False, rf"+{TILT_DEG:g}° tilt"),
        (TILT_DEG, 0.0, True, rf"+{TILT_DEG:g}° tilt, gradient"),
        (TILT_DEG, ROLL_DEG, True, rf"+{TILT_DEG:g}° tilt, gradient, +{ROLL_DEG:g}° roll"),
    ]
    # Aim near gold 1370×1420
    fig = plt.figure(figsize=(9.13, 9.47), dpi=150)
    for i, (tilt, roll, grad, title) in enumerate(rows):
        ax0 = fig.add_subplot(4, 3, 3 * i + 1)
        _draw_model_frame(ax0, tilt, roll, gradient=grad)
        ax0.set_title("Model frame" if i == 0 else "", fontsize=10)
        ax0.text(0.0, 1.22, title, ha="center", fontsize=9, transform=ax0.transData, clip_on=False)
        ax1 = fig.add_subplot(4, 3, 3 * i + 2)
        _draw_mollweide(ax1, tilt, roll, gradient=grad)
        ax1.set_title("Mollweide (Galactic)" if i == 0 else "", fontsize=10)
        ax2 = fig.add_subplot(4, 3, 3 * i + 3, projection="3d")
        _draw_sphere3d(ax2, tilt, roll, gradient=grad)
        ax2.set_title("Sphere 3D" if i == 0 else "", fontsize=10)
    fig.suptitle(
        "Sky orientation sequence: untilted → +15° tilt → gradient → +80.8° roll",
        fontsize=12, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    png = OUT / "_Fig05_sky_orientation_four_row.png"
    fig.savefig(png, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    # JPEG delivery (letter uses jpg)
    dest = OUT / OUT_NAME
    im = Image.open(png).convert("RGB")
    # Resize width to gold 1370 keeping aspect
    w = 1370
    h = round(im.height * w / im.width)
    im.resize((w, h), Image.LANCZOS).save(dest, quality=92)
    print("wrote", dest, dest.stat().st_size, im.size, "->", (w, h))
    return dest


if __name__ == "__main__":
    main()

