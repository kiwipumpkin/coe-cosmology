#!/usr/bin/env python3
"""A14+A15 as one 1x2 figure whose plot spines align with Fig A16's columns.
The maker (make_amp_phi_and_geometry_heatmaps.py part1) needs /tmp/chris_sheet_series.csv,
which is gone, so the plot interiors (heatmap + baked-in markers/legend) are cropped from the
existing PNGs and re-placed on fresh axes with the original data extents; ticks, labels and
colorbars are redrawn (A: magma 0..0.737, B: plasma 0..~0.795 measured from the old colorbar)."""
import numpy as np, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mi
F = "/workspace/cosmo_preprint/figs/"
A16_W = 1280.0; COLS = [(72, 629), (697, 1254)]        # A16 spine x (px of 1280)
W, H = 10.0, 5.0
PB, PT = 0.10, 0.75                                     # panel bottom/top (fig frac)
panels = [
  dict(f="FigApp_heat_Amp_phi_A.png", crop=(277, 1597, 126, 1189), ext=(0, 180, 0.002, -0.152),
       cmap="magma", vmax=0.7371, xl=r"rotation $\phi$ [deg]", yl="Amp$_A$",
       t=r"A14  $R^2_A$ over Amp$_A\times\phi$", cl=r"$R^2_A$"),
  dict(f="FigApp_heat_Amp_phi_B.png", crop=(278, 1692, 192, 1249), ext=(0, 350, -0.125, 0.125),
       cmap="plasma", vmax=0.795, xl=r"$\phi$ [deg]", yl="Amp$_B$",
       t=r"A15  $R^2_B$ over Amp$_B\times\phi$", cl=r"$R^2_B$"),
]
mpl.rcParams.update({"font.size": 11, "axes.labelsize": 11, "xtick.labelsize": 9.5, "ytick.labelsize": 9})
fig = plt.figure(figsize=(W, H))
for p, (x0, x1) in zip(panels, COLS):
    ax = fig.add_axes([x0 / A16_W, PB, (x1 - x0) / A16_W, PT - PB])
    img = mi.imread(F + p["f"]); c = p["crop"]
    l, r, b, t = p["ext"]
    ax.imshow(img[c[2]:c[3], c[0]:c[1]], extent=(l, r, b, t), aspect="auto", interpolation="lanczos")
    ax.set_xlim(l, r); ax.set_ylim(b, t)
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    ax.set_xlabel(p["xl"]); ax.set_ylabel(p["yl"], labelpad=0); ax.tick_params(axis="y", pad=2)
    pos = ax.get_position()
    cax = fig.add_axes([pos.x0, PT + 0.03, pos.width, 0.035])
    cb = fig.colorbar(mpl.cm.ScalarMappable(mpl.colors.Normalize(0, p["vmax"]), p["cmap"]),
                      cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("top"); cb.ax.tick_params(labelsize=9, pad=1)
    cax.set_title(p["t"] + "   (" + p["cl"] + ")", fontsize=12, pad=14)
fig.suptitle(r"Amp $\times\ \phi$ $R^2$ heatmaps (sheet AF/AN)", fontsize=14, y=0.995, va="top")
fig.savefig(F + "FigApp_heat_Amp_phi_1x2.png", dpi=256)   # 2560 px wide = 2x A16
