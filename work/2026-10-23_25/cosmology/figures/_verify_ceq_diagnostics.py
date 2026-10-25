#!/usr/bin/env python3
"""Save before/after Dec=0 vs published CEQ diagnostics."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from make_migkas_fig8_coe_rings_overlay import (
    choose_source,
    detect_oval_bbox,
    celestial_equator_lb,
    lon_disp_from_ell,
    mollweide_xy,
    published_ceq_mask,
    skeletonize_stroke,
    chamfer_distance_to_dark,
    refine_registration,
)

OUT = Path(__file__).resolve().parent


def main():
    source = choose_source()
    rgb = np.asarray(Image.open(source).convert("RGB"))
    seed = detect_oval_bbox(rgb)
    scx, scy, sRx, sRy = (float(seed[k]) for k in ("cx", "cy", "Rx", "Ry"))
    h, w, _ = rgb.shape

    ceq = published_ceq_mask(rgb, seed)
    skel, _ = skeletonize_stroke(ceq)
    dist = chamfer_distance_to_dark(skel)
    ys, xs = np.where(ceq)

    ell, b = celestial_equator_lb(2800)
    mx, my = mollweide_xy(np.deg2rad(lon_disp_from_ell(ell)), np.deg2rad(b))
    Xn = mx / (2 * np.sqrt(2))
    Yn = my / np.sqrt(2)
    lon = lon_disp_from_ell(ell)

    def curve(cx, cy, Rx, Ry):
        return cx + Xn * Rx, cy - Yn * Ry

    def stats(cx, cy, Rx, Ry, rmax=0.95):
        px, py = curve(cx, cy, Rx, Ry)
        m = ((px - scx) / sRx) ** 2 + ((py - scy) / sRy) ** 2 <= rmax ** 2
        m &= (px >= 1) & (px < w - 1) & (py >= 1) & (py < h - 1)
        xi = np.clip(np.round(px[m]).astype(int), 0, w - 1)
        yi = np.clip(np.round(py[m]).astype(int), 0, h - 1)
        d = dist[yi, xi]
        return {
            "median": float(np.median(d)),
            "p75": float(np.percentile(d, 75)),
            "p90": float(np.percentile(d, 90)),
            "mean": float(np.mean(d)),
            "n": int(d.size),
            "lon": lon[m],
            "d": d,
            "px": px,
            "py": py,
            "m": m,
        }

    # BEFORE = previous (edge-locking) registration from prior JSON
    old = dict(cx=863.0, cy=763.25, Rx=693.9927, Ry=392.542425)
    before = stats(**old)
    print("BEFORE skel residual", {k: before[k] for k in ("median", "p75", "p90", "mean", "n")})

    # AFTER = new refine
    reg = refine_registration(rgb, seed)
    after = stats(reg["cx"], reg["cy"], reg["Rx"], reg["Ry"])
    print("AFTER skel residual", {k: after[k] for k in ("median", "p75", "p90", "mean", "n")})
    print("reg", {k: reg[k] for k in ("cx", "cy", "Rx", "Ry", "ceq_median_px")})

    def make_fig(tag, st, other=None, title=""):
        fig, axes = plt.subplots(2, 2, figsize=(13, 10))
        crops = [
            (scx - 200, scx + 200, scy + 140, scy - 140, f"{tag} center"),
            (scx - 90, scx + 90, scy - 40, scy - 210, f"{tag} upper arc"),
            (220, 520, 920, 560, f"{tag} left"),
        ]
        for ax, (x0, x1, y0, y1, ttl) in zip(axes.ravel()[:3], crops):
            ax.imshow(rgb)
            ax.scatter(xs[::2], ys[::2], s=1.5, c="yellow", alpha=0.45, zorder=2)
            m = st["m"]
            ax.plot(st["px"][m], st["py"][m], color="lime", lw=1.5, label="Dec=0 geom", zorder=3)
            if other is not None:
                mo = other["m"]
                ax.plot(
                    other["px"][mo],
                    other["py"][mo],
                    color="red",
                    lw=1.0,
                    ls="--",
                    alpha=0.7,
                    label="other",
                    zorder=3,
                )
            if tag == "BEFORE":
                # show what inset 0.97 did
                pxi = old["cx"] + (st["px"] - old["cx"]) * 0.97
                pyi = old["cy"] + (st["py"] - old["cy"]) * 0.97
                ax.plot(pxi[m], pyi[m], "c--", lw=1.0, label="inset 0.97", zorder=3)
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            ax.set_title(ttl)
            ax.legend(fontsize=7)
        ax = axes[1, 1]
        order = np.argsort(st["lon"])
        ax.plot(st["lon"][order], st["d"][order], "g.", ms=1.2, label="this")
        if other is not None:
            oo = np.argsort(other["lon"])
            ax.plot(other["lon"][oo], other["d"][oo], "r.", ms=1.0, alpha=0.4, label="other")
        ax.set_ylim(0, 20)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        ax.set_xlabel("lon_disp deg")
        ax.set_ylabel("px to pub CEQ skeleton")
        ax.set_title(
            f'{tag} med={st["median"]:.2f} p90={st["p90"]:.2f}'
        )
        fig.suptitle(title, fontsize=11)
        fig.tight_layout()
        out = OUT / f"_verify_ceq_{tag.lower()}.png"
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print("wrote", out)

    make_fig(
        "BEFORE",
        before,
        title="Dec=0 BEFORE retune (yellow=pub CEQ stroke; lime=old geom; cyan=inset 0.97)",
    )
    make_fig(
        "AFTER",
        after,
        other=before,
        title="Dec=0 AFTER retune (lime=new geom inset=1; red=old geom)",
    )

    meta = {
        "before_old_registration": old,
        "before_skel_residual": {k: before[k] for k in ("median", "p75", "p90", "mean", "n")},
        "after_registration": {
            "cx": reg["cx"],
            "cy": reg["cy"],
            "Rx": reg["Rx"],
            "Ry": reg["Ry"],
        },
        "after_skel_residual": {k: after[k] for k in ("median", "p75", "p90", "mean", "n")},
        "dofs": ["cx", "cy", "Rx", "Ry"],
        "target": "DT-ridge medial axis of achromatic published CEQ stroke",
        "note": (
            "Previous gray<55 chamfer reported ~1 px by locking onto thick stroke EDGE; "
            "honest centerline median was ~5.2 px. Drawing inset=0.97 added ~8 px visual offset."
        ),
    }
    (OUT / "_verify_ceq_stats.json").write_text(json.dumps(meta, indent=2))
    print("wrote", OUT / "_verify_ceq_stats.json")


if __name__ == "__main__":
    main()
