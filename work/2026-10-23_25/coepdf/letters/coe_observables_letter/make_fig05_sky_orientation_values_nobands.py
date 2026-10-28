#!/usr/bin/env python3
"""Fig 4 — VALUES, NO BANDS version (Chris: gradient band lines confuse the colour scale).
MAX/MIN band ribbons + mid-strokes removed in panels A/B/C (draw_panelA_ribbon,
fill_mollweide_ribbon, fill_sphere_ribbon -> no-op); colour alpha raised.

Based on VALUES version:

Same four rows (untilted / +15 tilt / +15 tilt+gradient / +80.8 roll) built by
Chris's four ring scripts in /workspace/cosmology/figures, loaded unchanged;
their panel functions are wrapped to ADD, underneath the existing MAX/MIN bands,
the actual Migkas A/A_all value at every scan angle theta:

  data : /workspace/cosmology/migkas/fig6_digitized.csv (71 pts, l = 5..355 deg)
  curve: connect-the-dots, np.interp(theta, l, A, period=360)  (no splines)
  theta = Galactic-longitude scan angle of the Migkas curve, measured from +x
          (anti-CoE), exactly as the bands' theta in the original figure.
  Sphere: a direction n has x = n.D = cos(theta_abs); the half of each cone on
          the +y side of the ring basis gets A(theta_abs), the -y half gets
          A(360 - theta_abs)  (panel A: top of circle = theta, bottom = 360-theta).
          y is the row's own ring_basis (incl. +80.8 roll in row 4); rows without
          ring_basis use the same unrolled construction.
Panel A: circle rim coloured by A(theta) at angle theta. Panel B: Mollweide sky
coloured (alpha .55). Panel C: sphere surface coloured (alpha .4).
One shared colourbar (same norm, data min..max) appended under the legend.
Original PNGs/scripts are not touched; outputs go to *_values files here.
"""
import importlib.util, sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from PIL import Image

HERE = Path(__file__).resolve().parent
FIG_DIR = Path("/workspace/cosmology/figures")
sys.path.insert(0, str(HERE))
import types; sys.modules.setdefault('paths', types.SimpleNamespace(OUT=HERE))
import make_fig05_sky_orientation as base  # reuse compose_four_row unchanged

d = pd.read_csv("/workspace/cosmology/migkas/fig6_digitized.csv")
L, AV = d.l_deg.values.astype(float), d.A_over_Aall.values.astype(float)
def A_of(theta):
    return np.interp(np.mod(theta, 360.0), L, AV, period=360)
NORM = matplotlib.colors.Normalize(vmin=AV.min(), vmax=AV.max())
CMAP = plt.get_cmap("coolwarm")

def basis_y(mod, D):
    a = np.asarray(D, float); a /= np.linalg.norm(a)
    if hasattr(mod, "ring_basis"):
        return mod.ring_basis(a)[0]
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(a @ tmp) > 0.9: tmp = np.array([0.0, 1.0, 0.0])
    y = tmp - (tmp @ a) * a
    return y / np.linalg.norm(y)

def A_on_dirs(mod, D, n):  # n (...,3)
    a = np.asarray(D, float); a /= np.linalg.norm(a)
    y = basis_y(mod, a)
    th = np.rad2deg(np.arccos(np.clip(n @ a, -1, 1)))
    th = np.where(n @ y >= 0, th, 360.0 - th)
    return A_of(th)

def wrap(mod, tag):
    mod.draw_panelA_ribbon = lambda *a, **k: None
    mod.fill_mollweide_ribbon = lambda *a, **k: None
    mod.fill_sphere_ribbon = lambda *a, **k: None
    oA, oB, oC = mod.panel_A_model_frame, mod.panel_B_mollweide, mod.panel_C_sphere
    def pA(ax, geo):
        oA(ax, geo)
        t = np.linspace(0, 360, 721)
        x, y = np.cos(np.deg2rad(t)), np.sin(np.deg2rad(t))
        seg = np.stack([np.column_stack([x[:-1], y[:-1]]), np.column_stack([x[1:], y[1:]])], 1)
        ax.add_collection(LineCollection(seg, colors=CMAP(NORM(A_of(0.5*(t[:-1]+t[1:])))),
                                         linewidths=6, zorder=1.5, capstyle="butt"))
    def pB(ax, geo):
        oB(ax, geo)
        lon = np.linspace(-np.pi, np.pi, 361); lat = np.linspace(-np.pi/2, np.pi/2, 181)
        LO, LA = np.meshgrid(lon, lat)
        ell = np.mod(-np.rad2deg(LO), 360.0)  # module: lon_disp = -l  (l increases LEFT)
        # verify against module convention
        n = np.stack([np.cos(LA)*np.cos(np.deg2rad(ell)), np.cos(LA)*np.sin(np.deg2rad(ell)), np.sin(LA)], -1)
        ax.pcolormesh(LO, LA, A_on_dirs(mod, geo["D"], n), cmap=CMAP, norm=NORM,
                      shading="gouraud", alpha=0.85, zorder=0.5, rasterized=True)
    def pC(ax, geo):
        u, v = np.meshgrid(np.linspace(0, 2*np.pi, 121), np.linspace(0, np.pi, 61))
        n = np.stack([np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)], -1)
        fc = CMAP(NORM(A_on_dirs(mod, geo["D"], n))); fc[..., 3] = 0.6
        ax.plot_surface(n[..., 0]*0.995, n[..., 1]*0.995, n[..., 2]*0.995, facecolors=fc,
                        rstride=1, cstride=1, linewidth=0, shade=False, zorder=1)
        oC(ax, geo)
    mod.panel_A_model_frame, mod.panel_B_mollweide, mod.panel_C_sphere = pA, pB, pC

def load(script, i):
    if str(FIG_DIR) not in sys.path: sys.path.insert(0, str(FIG_DIR))
    spec = importlib.util.spec_from_file_location(f"vals_{i}", FIG_DIR / script)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

def main():
    pngs = []
    for i, s in enumerate(base.RING_SCRIPTS):
        m = load(s, i)
        # display-lon convention check
        print(s, "lon_disp(90)=", m.lon_disp_from_ell(np.array([90.0])))
        m.PNG_PATH = HERE / f"fig8_row{i}_values_nobands.png"
        m.JSON_PATH = HERE / f"fig8_row{i}_values_nobands.json"
        wrap(m, i); m.main(); pngs.append(m.PNG_PATH)
    comp = base.compose_four_row(pngs)
    fig = plt.figure(figsize=(comp.width/160, 0.75), dpi=160)
    cax = fig.add_axes([0.2, 0.62, 0.6, 0.22])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(NORM, CMAP), cax=cax, orientation="horizontal")
    cb.set_label(r"Migkas $A/A_{\rm all}$ at scan angle $\theta$ (fig6 digitized, linear interp.)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    tmp = HERE / "_cbar_values_nobands_tmp.png"; fig.savefig(tmp, dpi=160, facecolor="white", bbox_inches="tight", pad_inches=0.05); plt.close(fig)
    cbi = Image.open(tmp).convert("RGB"); cbi.thumbnail((int(comp.width*0.7), 10**4))
    out = Image.new("RGB", (comp.width, comp.height + cbi.height), "white")
    out.paste(comp, (0, 0)); out.paste(cbi, ((comp.width-cbi.width)//2, comp.height)); tmp.unlink()
    out.save(HERE / "fig8_orientation_seq_values_nobands.png"); print("wrote", HERE / "fig8_orientation_seq_values_nobands.png", out.size)

if __name__ == "__main__":
    main()
