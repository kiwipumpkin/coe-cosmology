#!/usr/bin/env python3
"""Tune Mollweide registration so Dec=0 tracks published CEQ centerline."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_migkas_fig8_coe_rings_overlay import (
    choose_source,
    detect_oval_bbox,
    celestial_equator_lb,
    lon_disp_from_ell,
    mollweide_xy,
    chamfer_distance_to_dark,
)

OUT = Path(__file__).resolve().parent


def published_ceq_mask(rgb, seed):
    """Achromatic black stroke = published celestial equator (not dark magma)."""
    rgb = np.asarray(rgb, dtype=np.float64)
    h, w, _ = rgb.shape
    achromatic = (rgb.max(2) < 30) & ((rgb.max(2) - rgb.min(2)) < 25)
    scx, scy, sRx, sRy = [float(seed[k]) for k in ("cx", "cy", "Rx", "Ry")]
    yy, xx = np.mgrid[0:h, 0:w]
    oval = ((xx - scx) / sRx) ** 2 + ((yy - scy) / sRy) ** 2 <= 0.995**2
    return achromatic & oval


def mollweide_norm_ceq(n=4000):
    ell, b = celestial_equator_lb(n)
    mx, my = mollweide_xy(np.deg2rad(lon_disp_from_ell(ell)), np.deg2rad(b))
    Xn = mx / (2.0 * np.sqrt(2.0))
    Yn = my / np.sqrt(2.0)
    return ell, b, Xn, Yn


def extract_centerline(Xn, Yn, cx, cy, Rx, Ry, mask, seed, search=30):
    """Pair each Dec=0 sample with published-stroke center along local normal."""
    h, w = mask.shape
    scx, scy, sRx, sRy = [float(seed[k]) for k in ("cx", "cy", "Rx", "Ry")]
    px = cx + Xn * Rx
    py = cy - Yn * Ry
    t = np.gradient(np.stack([px, py], 1), axis=0)
    tn = t / (np.linalg.norm(t, axis=1, keepdims=True) + 1e-9)
    nm = np.stack([-tn[:, 1], tn[:, 0]], 1)
    tx, ty, ix, off = [], [], [], []
    for i in range(len(px)):
        x, y = float(px[i]), float(py[i])
        if ((x - scx) / sRx) ** 2 + ((y - scy) / sRy) ** 2 > 0.955**2:
            continue
        if not (2 <= x < w - 2 and 2 <= y < h - 2):
            continue
        hits = []
        for s in range(-search, search + 1):
            xi = int(round(x + s * nm[i, 0]))
            yi = int(round(y + s * nm[i, 1]))
            if 0 <= xi < w and 0 <= yi < h and mask[yi, xi]:
                hits.append(s)
        if not hits:
            continue
        hits = np.asarray(hits, dtype=float)
        # contiguous runs; pick run closest to s=0
        runs, a, prev = [], hits[0], hits[0]
        for v in hits[1:]:
            if v - prev > 1.5:
                runs.append((a, prev))
                a = v
            prev = v
        runs.append((a, prev))
        a, b = min(runs, key=lambda ab: abs(0.5 * (ab[0] + ab[1])))
        s_c = 0.5 * (a + b)
        tx.append(x + s_c * nm[i, 0])
        ty.append(y + s_c * nm[i, 1])
        ix.append(i)
        off.append(s_c)
    return (
        np.asarray(tx),
        np.asarray(ty),
        np.asarray(ix, dtype=int),
        np.asarray(off),
    )


def fit_diagonal(Xn, Yn, tx, ty, cx0, cy0, Rx0, Ry0):
    """Median residual grid search over (cx,cy,Rx,Ry)."""

    def med(cx, cy, Rx, Ry):
        return float(np.median(np.hypot(cx + Rx * Xn - tx, cy - Ry * Yn - ty)))

    best = dict(cx=cx0, cy=cy0, Rx=Rx0, Ry=Ry0, sc=med(cx0, cy0, Rx0, Ry0))
    for dcx in np.linspace(-16, 16, 33):
        for dcy in np.linspace(-16, 16, 33):
            for sRx in np.linspace(0.95, 1.05, 21):
                for sRy in np.linspace(0.95, 1.05, 21):
                    cx, cy = cx0 + dcx, cy0 + dcy
                    Rx, Ry = Rx0 * sRx, Ry0 * sRy
                    sc = med(cx, cy, Rx, Ry)
                    if sc < best["sc"]:
                        best = dict(cx=cx, cy=cy, Rx=Rx, Ry=Ry, sc=sc)
    b = best
    for dcx in np.linspace(-2.5, 2.5, 26):
        for dcy in np.linspace(-2.5, 2.5, 26):
            for sRx in np.linspace(0.99, 1.01, 21):
                for sRy in np.linspace(0.99, 1.01, 21):
                    cx = b["cx"] + dcx
                    cy = b["cy"] + dcy
                    Rx = b["Rx"] * sRx
                    Ry = b["Ry"] * sRy
                    sc = med(cx, cy, Rx, Ry)
                    if sc < best["sc"]:
                        best = dict(cx=cx, cy=cy, Rx=Rx, Ry=Ry, sc=sc)
    b = best
    for dcx in np.linspace(-0.6, 0.6, 13):
        for dcy in np.linspace(-0.6, 0.6, 13):
            for sRx in np.linspace(0.998, 1.002, 9):
                for sRy in np.linspace(0.998, 1.002, 9):
                    cx = b["cx"] + dcx
                    cy = b["cy"] + dcy
                    Rx = b["Rx"] * sRx
                    Ry = b["Ry"] * sRy
                    sc = med(cx, cy, Rx, Ry)
                    if sc < best["sc"]:
                        best = dict(cx=cx, cy=cy, Rx=Rx, Ry=Ry, sc=sc)
    return best


def fit_affine(Xn, Yn, tx, ty):
    M = np.column_stack([Xn, Yn, np.ones(len(Xn))])
    bx, *_ = np.linalg.lstsq(M, tx, rcond=None)
    by, *_ = np.linalg.lstsq(M, ty, rcond=None)
    px = bx[2] + bx[0] * Xn + bx[1] * Yn
    py = by[2] + by[0] * Xn + by[1] * Yn
    err = np.hypot(px - tx, py - ty)
    return {
        "Axx": float(bx[0]),
        "Axy": float(bx[1]),
        "cx": float(bx[2]),
        "Ayx": float(by[0]),
        "Ayy": float(by[1]),
        "cy": float(by[2]),
        "median": float(np.median(err)),
        "p90": float(np.percentile(err, 90)),
        "mean": float(np.mean(err)),
    }


def residual_stats(Xn, Yn, tx, ty, cx, cy, Rx, Ry):
    err = np.hypot(cx + Rx * Xn - tx, cy - Ry * Yn - ty)
    return {
        "median": float(np.median(err)),
        "p75": float(np.percentile(err, 75)),
        "p90": float(np.percentile(err, 90)),
        "mean": float(np.mean(err)),
        "n": int(len(err)),
    }


def main():
    source = choose_source()
    rgb = np.asarray(Image.open(source).convert("RGB"))
    seed = detect_oval_bbox(rgb)
    mask = published_ceq_mask(rgb, seed)
    print("source", source, "ceq_dark", int(mask.sum()))
    print("seed", {k: seed[k] for k in ("cx", "cy", "Rx", "Ry")})

    ell, b, Xn0, Yn0 = mollweide_norm_ceq(4000)

    # --- BEFORE: previous refine_registration params from JSON ---
    old = dict(cx=863.0, cy=763.25, Rx=693.9927, Ry=392.542425)
    tx_o, ty_o, ix_o, off_o = extract_centerline(
        Xn0, Yn0, old["cx"], old["cy"], old["Rx"], old["Ry"], mask, seed, search=24
    )
    before = residual_stats(
        Xn0[ix_o], Yn0[ix_o], tx_o, ty_o, old["cx"], old["cy"], old["Rx"], old["Ry"]
    )
    before["off_abs_med"] = float(np.median(np.abs(off_o)))
    print("BEFORE (old reg vs centerline)", before)

    # Freeze targets from seed prior (wide search)
    tx, ty, ix, off = extract_centerline(
        Xn0, Yn0, seed["cx"], seed["cy"], seed["Rx"], seed["Ry"], mask, seed, search=32
    )
    print(
        "seed extract n",
        len(tx),
        "off_med",
        float(np.median(off)),
        "abs_med",
        float(np.median(np.abs(off))),
    )
    Xn, Yn = Xn0[ix], Yn0[ix]

    diag0 = fit_diagonal(Xn, Yn, tx, ty, seed["cx"], seed["cy"], seed["Rx"], seed["Ry"])
    print("diag fit from seed targets", diag0)

    # Re-extract with refined prior (tighter), refit once
    tx2, ty2, ix2, off2 = extract_centerline(
        Xn0,
        Yn0,
        diag0["cx"],
        diag0["cy"],
        diag0["Rx"],
        diag0["Ry"],
        mask,
        seed,
        search=16,
    )
    print(
        "re-extract n",
        len(tx2),
        "abs_med",
        float(np.median(np.abs(off2))),
        "abs_p90",
        float(np.percentile(np.abs(off2), 90)),
    )
    Xn2, Yn2 = Xn0[ix2], Yn0[ix2]
    diag = fit_diagonal(
        Xn2, Yn2, tx2, ty2, diag0["cx"], diag0["cy"], diag0["Rx"], diag0["Ry"]
    )
    print("diag refit", diag)

    after = residual_stats(
        Xn2, Yn2, tx2, ty2, diag["cx"], diag["cy"], diag["Rx"], diag["Ry"]
    )
    after["off_abs_med"] = float(np.median(np.abs(off2)))
    print("AFTER diagonal", after)

    aff = fit_affine(Xn2, Yn2, tx2, ty2)
    print("affine LS", {k: aff[k] for k in aff})
    # How much shear vs pure diagonal?
    print(
        "affine vs diag: Axy",
        aff["Axy"],
        "Ayx",
        aff["Ayx"],
        "Axx~Rx",
        aff["Axx"],
        " -Ayy~Ry",
        -aff["Ayy"],
    )

    # Optional: one more extract using affine-projected prior for scoring only
    # Decide: if affine median clearly better (>20% relative) keep affine DOFs
    use_affine = aff["median"] < 0.85 * after["median"]
    print("use_affine", use_affine)

    # --- Diagnostics figures ---
    ys, xs = np.where(mask)
    scx, scy, sRx, sRy = [float(seed[k]) for k in ("cx", "cy", "Rx", "Ry")]

    def draw_curve(ax, cx, cy, Rx, Ry, **kw):
        px = cx + Xn0 * Rx
        py = cy - Yn0 * Ry
        m = ((px - scx) / sRx) ** 2 + ((py - scy) / sRy) ** 2 <= 0.96**2
        ax.plot(px[m], py[m], **kw)

    def draw_affine(ax, aff, **kw):
        px = aff["cx"] + aff["Axx"] * Xn0 + aff["Axy"] * Yn0
        py = aff["cy"] + aff["Ayx"] * Xn0 + aff["Ayy"] * Yn0
        m = ((px - scx) / sRx) ** 2 + ((py - scy) / sRy) ** 2 <= 0.96**2
        ax.plot(px[m], py[m], **kw)

    # BEFORE figure
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    for ax, (x0, x1, y0, y1, title) in zip(
        axes.ravel()[:3],
        [
            (scx - 200, scx + 200, scy + 140, scy - 140, "BEFORE center"),
            (scx - 90, scx + 90, scy - 40, scy - 210, "BEFORE upper arc"),
            (220, 520, 920, 560, "BEFORE left"),
        ],
    ):
        ax.imshow(rgb)
        ax.scatter(xs[::2], ys[::2], s=1.5, c="yellow", alpha=0.45, zorder=2)
        draw_curve(
            ax,
            old["cx"],
            old["cy"],
            old["Rx"],
            old["Ry"],
            color="lime",
            lw=1.4,
            label="old Dec=0 (geom)",
            zorder=3,
        )
        # inset cosmetic as previously drawn
        px = old["cx"] + Xn0 * old["Rx"]
        py = old["cy"] - Yn0 * old["Ry"]
        pxi = old["cx"] + (px - old["cx"]) * 0.97
        pyi = old["cy"] + (py - old["cy"]) * 0.97
        m = ((px - scx) / sRx) ** 2 + ((py - scy) / sRy) ** 2 <= 0.96**2
        ax.plot(
            pxi[m],
            pyi[m],
            "c--",
            lw=1.0,
            label="old inset 0.97",
            zorder=3,
        )
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_title(title)
        ax.legend(fontsize=7, loc="best")
    ax = axes[1, 1]
    lon_o = lon_disp_from_ell(ell)[ix_o]
    err_o = np.hypot(
        old["cx"] + old["Rx"] * Xn0[ix_o] - tx_o,
        old["cy"] - old["Ry"] * Yn0[ix_o] - ty_o,
    )
    order = np.argsort(lon_o)
    ax.plot(lon_o[order], err_o[order], "g.", ms=1.2)
    ax.set_ylim(0, 20)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("lon_disp deg")
    ax.set_ylabel("px to pub centerline")
    ax.set_title(
        f'BEFORE residual med={before["median"]:.2f} p90={before["p90"]:.2f}'
    )
    fig.suptitle(
        "Dec=0 BEFORE retune (yellow = published CEQ black stroke)", fontsize=11
    )
    fig.tight_layout()
    fig.savefig(OUT / "_verify_ceq_before.png", dpi=140)
    plt.close(fig)
    print("wrote", OUT / "_verify_ceq_before.png")

    # AFTER figure
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    for ax, (x0, x1, y0, y1, title) in zip(
        axes.ravel()[:3],
        [
            (scx - 200, scx + 200, scy + 140, scy - 140, "AFTER center"),
            (scx - 90, scx + 90, scy - 40, scy - 210, "AFTER upper arc"),
            (220, 520, 920, 560, "AFTER left"),
        ],
    ):
        ax.imshow(rgb)
        ax.scatter(xs[::2], ys[::2], s=1.5, c="yellow", alpha=0.45, zorder=2)
        draw_curve(
            ax,
            diag["cx"],
            diag["cy"],
            diag["Rx"],
            diag["Ry"],
            color="lime",
            lw=1.5,
            label="new Dec=0 (cx,cy,Rx,Ry)",
            zorder=3,
        )
        draw_curve(
            ax,
            old["cx"],
            old["cy"],
            old["Rx"],
            old["Ry"],
            color="red",
            lw=1.0,
            ls="--",
            alpha=0.7,
            label="old",
            zorder=3,
        )
        if use_affine:
            draw_affine(ax, aff, color="cyan", lw=1.0, ls=":", label="affine", zorder=3)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_title(title)
        ax.legend(fontsize=7, loc="best")
    ax = axes[1, 1]
    lon2 = lon_disp_from_ell(ell)[ix2]
    err2 = np.hypot(
        diag["cx"] + diag["Rx"] * Xn2 - tx2,
        diag["cy"] - diag["Ry"] * Yn2 - ty2,
    )
    order = np.argsort(lon2)
    ax.plot(lon2[order], err2[order], "g.", ms=1.2, label="new")
    ax.plot(lon_o[np.argsort(lon_o)], err_o[np.argsort(lon_o)], "r.", ms=1.0, alpha=0.4, label="old")
    ax.set_ylim(0, 20)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_xlabel("lon_disp deg")
    ax.set_title(
        f'AFTER med={after["median"]:.2f} p90={after["p90"]:.2f} (was {before["median"]:.2f})'
    )
    fig.suptitle(
        "Dec=0 AFTER retune to published CEQ centerline", fontsize=11
    )
    fig.tight_layout()
    fig.savefig(OUT / "_verify_ceq_after.png", dpi=140)
    plt.close(fig)
    print("wrote", OUT / "_verify_ceq_after.png")

    result = {
        "before": before,
        "after_diagonal": after,
        "affine": aff,
        "use_affine": use_affine,
        "registration_diagonal": {
            "cx": float(diag["cx"]),
            "cy": float(diag["cy"]),
            "Rx": float(diag["Rx"]),
            "Ry": float(diag["Ry"]),
            "median_obj": float(diag["sc"]),
        },
        "old_registration": old,
        "seed": {k: float(seed[k]) if k in ("cx", "cy", "Rx", "Ry") else seed[k]
                 for k in ("cx", "cy", "Rx", "Ry", "x0", "x1", "y0", "y1")},
    }
    (OUT / "_ceq_tune_result.json").write_text(json.dumps(result, indent=2))
    print("wrote", OUT / "_ceq_tune_result.json")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
