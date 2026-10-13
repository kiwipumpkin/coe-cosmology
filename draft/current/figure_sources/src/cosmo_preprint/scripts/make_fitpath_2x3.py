#!/usr/bin/env python3
"""Fit-path composite (old A1-A5) as one 2x3 grid. No makers/data exist for these panels,
so the existing JPGs are cropped and uniformly scaled (no aspect change) so that every
panel's plot spines share the same width and left x-position within its column."""
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import matplotlib
F = "/workspace/cosmo_preprint/figs/"
font = ImageFont.truetype(matplotlib.get_data_path() + "/fonts/ttf/DejaVuSans-Bold.ttf", 44)
# (file, x-crop, spine_left, spine_right) spine x measured on the source
cells = [
    ("FigApp_AF_sine_rotation.jpg", (0, 1280), 114, 1267, "(a)"),
    ("FigApp_2ell_correction.jpg", (0, 1280), 121, 1266, "(b)"),
    ("FigApp_bestA_assocB.jpg", (0, 1280), 103, 1266, "(c)"),
    ("FigApp_bestB_assocA.jpg", (0, 1280), 103, 1266, "(d)"),
    ("FigApp_sheet_scatter_AB.jpg", (0, 660), 81, 624, "(e)"),
    ("FigApp_sheet_scatter_AB.jpg", (662, 1280), 725, 1268, "(f)"),
]
SW = 1160                     # common spine width (px)
def prep(f, xc, sl, sr):
    im = Image.open(F + f).convert("RGB").crop((xc[0], 0, xc[1], Image.open(F + f).size[1]))
    if f.startswith("FigApp_sheet"):   # drop the shared suptitle band
        im = im.crop((0, 45, im.size[0], im.size[1]))
    im = im.crop((0, 0, im.size[0], im.size[1] - 8))   # drop source bottom border line
    a = np.asarray(im.convert("L")); nz = np.where((a < 235).any(1))[0]
    im = im.crop((0, nz.min(), im.size[0], nz.max() + 1))
    s = SW / (sr - sl)
    im = im.resize((round(im.size[0] * s), round(im.size[1] * s)), Image.LANCZOS)
    return im, (sl - xc[0]) * s, (xc[1] - sr) * s
P = [prep(c[0], c[1], c[2], c[3]) for c in cells]
LM = max(p[1] for p in P); RM = max(p[2] for p in P)
COLW = int(LM + SW + RM); GAP = 10; VGAP = 40
rows = [max(P[2*r][0].size[1], P[2*r+1][0].size[1]) for r in range(3)]
canvas = Image.new("RGB", (2 * COLW + GAP, sum(rows) + 2 * VGAP), "white")
d = ImageDraw.Draw(canvas); y = 0
for r in range(3):
    for c in range(2):
        im, lm, _ = P[2*r+c]; x0 = c * (COLW + GAP) + int(LM - lm)
        canvas.paste(im, (x0, y))
        d.text((c * (COLW + GAP) + 4, y + 2), cells[2*r+c][4], fill="black", font=font)
    y += rows[r] + VGAP
canvas.save(F + "FigApp_fitpath_2x3.png", dpi=(300, 300)); print(canvas.size, LM, RM)
