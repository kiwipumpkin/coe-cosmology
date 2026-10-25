#!/usr/bin/env python3
"""Full-size plots: T_obs vs angle for plates & shells. Labels in open space, no leader lines."""
from __future__ import annotations
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(HERE, "two_emitters.json")))

plt.rcParams.update({"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13})


def tag(ax, x, y, text, color="k", ha="left", va="center", fs=11.5):
    ax.text(
        x, y, text, color=color, fontsize=fs, va=va, ha=ha, zorder=6,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="none", alpha=0.9),
    )


# =====================================================================
# Figure 1 — Parallel plates: T_obs along a full great circle
# =====================================================================
P = DATA["plates"]
th = np.asarray(P["theta_deg"])
fig, axes = plt.subplots(3, 1, figsize=(15, 16.5), sharex=True)

curves = [
    (axes[0], P["mid_sameT"], P["offset_sameT"],
     "same rest-frame T on both plates",
     [("mid (z=0)", "C0", P["mid_sameT"]),
      ("offset z = +d/2", "C1", P["offset_sameT"])],
     "PERFECTLY ISOTROPIC\n(mid and offset identical)",
     55, 1.045),
    (axes[1], P["mid_diffT"], None,
     r"different T: $T_{+}=1.1\,T_0$, $T_{-}=0.9\,T_0$  (mid)",
     [("T_obs", "C3", P["mid_diffT"])],
     "step: upper vs lower\nhemisphere",
     40, 1.12),
    (axes[2], P["mid_sameT_doppler"], None,
     r"same T + observer motion $\beta=v/c=0.01$ along $+z$",
     [("T_obs (Doppler)", "C2", P["mid_sameT_doppler"])],
     "classic dipole\non flat T sky",
     55, 1.018),
]

for ax, _, _, title, series, note, nx, ny in curves:
    for lab, col, y in series:
        ax.plot(th, y, color=col, lw=2.8, label=lab)
    ax.axhline(1.0, color="0.5", lw=0.8, ls=":")
    ax.axvline(90, color="0.75", lw=0.8, ls="--")
    ax.axvline(270, color="0.75", lw=0.8, ls="--")
    ax.set_ylabel(r"$T_{\rm obs}/T_0$")
    ax.set_title(title, loc="left", pad=8)
    ax.grid(alpha=0.25)
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
    if len(series) > 1:
        ax.legend(loc="upper right", fontsize=11, framealpha=0.92)
    tag(ax, nx, ny, note, color="0.25", ha="left")

# hemisphere labels on top panel
tag(axes[0], 25, 0.955, "+plate\nhemisphere", "C0", ha="center", fs=10.5)
tag(axes[0], 200, 0.955, "−plate\nhemisphere", "C0", ha="center", fs=10.5)

axes[-1].set_xlabel(
    r"great-circle angle $\theta$ (degrees)      "
    r"$0^\circ$/$360^\circ$ = toward $+z$ (upper plate)      ·      "
    r"$180^\circ$ = toward $-z$ (lower plate)"
)
fig.suptitle(
    "Infinite parallel plates at $z=\\pm d$, blackbody emitters at rest-frame $T$\n"
    "Observer between plates · every ray hits exactly one plate · "
    r"$I_\nu=B_\nu(T_{\rm plate})$ $\Rightarrow$ $T_{\rm obs}=T_{\rm plate}$"
    "\n"
    "Same $T$ $\\Rightarrow$ isotropic sky (offset does nothing).  "
    "Anisotropy only from unequal $T$ or motion.",
    fontsize=13.5,
)
fig.tight_layout(rect=[0, 0, 1, 0.91])
p1 = os.path.join(HERE, "plates_T_vs_theta.png")
fig.savefig(p1, dpi=160)
plt.close(fig)
print("wrote", p1)


# =====================================================================
# Figure 2 — Concentric shells: T_obs along great circle
# =====================================================================
S = DATA["shells"]
th = np.asarray(S["theta_deg"])
alpha = S["alpha_deg"]

fig, axes = plt.subplots(3, 1, figsize=(15, 16.5), sharex=True)

# panel 0: same T
ax = axes[0]
ax.plot(th, S["sameT"], color="C0", lw=2.8)
ax.axhline(1.0, color="0.5", lw=0.8, ls=":")
for a in (alpha, 360 - alpha, 180 - alpha, 180 + alpha):
    ax.axvline(a, color="0.75", lw=0.9, ls="--")
ax.set_ylabel(r"$T_{\rm obs}/T_0$")
ax.set_title(
    r"same rest-frame $T$ on both shells"
    f"     (inner angular radius $\\alpha=\\arcsin(R_{{\\rm in}}/D)={alpha:.1f}^\\circ$)",
    loc="left", pad=8,
)
ax.set_ylim(0.94, 1.08)
ax.grid(alpha=0.25)
tag(ax, 55, 1.055,
    "PERFECTLY ISOTROPIC\n(inner disk + outer annulus\nboth at same T)",
    "0.25", ha="left")
tag(ax, 8, 0.965, "inner\nshell", "C0", ha="center", fs=10.5)
tag(ax, 100, 0.965, "outer shell", "C0", ha="center", fs=10.5)

# panel 1: different T
ax = axes[1]
ax.plot(th, S["diffT"], color="C3", lw=2.8)
ax.axhline(1.0, color="0.5", lw=0.8, ls=":")
for a in (alpha, 360 - alpha, 180 - alpha, 180 + alpha):
    ax.axvline(a, color="0.75", lw=0.9, ls="--")
ax.set_ylabel(r"$T_{\rm obs}/T_0$")
ax.set_title(
    r"different $T$: $T_{\rm in}=1.1\,T_0$, $T_{\rm out}=0.9\,T_0$",
    loc="left", pad=8,
)
ax.grid(alpha=0.25)
tag(ax, 8, 1.13, "inner disk\n" + r"$T=1.1\,T_0$", "C3", ha="center", fs=10.5)
tag(ax, 100, 0.93, "outer annulus\n" + r"$T=0.9\,T_0$", "C3", ha="center", fs=10.5)
tag(ax, 200, 1.08,
    "anisotropy = step at\nlimb of inner shell",
    "0.25", ha="left")

# panel 2: Doppler
ax = axes[2]
ax.plot(th, S["doppler_comoving_obs"], color="C2", lw=2.8, label="comoving observer (Hubble flow)")
ax.plot(th, S["doppler_coe_rest_obs"], color="C4", lw=2.2, ls="--",
        label="observer at rest vs CoE")
ax.axhline(1.0, color="0.5", lw=0.8, ls=":")
for a in (alpha, 360 - alpha, 180 - alpha, 180 + alpha):
    ax.axvline(a, color="0.75", lw=0.9, ls="--")
ax.set_ylabel(r"$T_{\rm obs}/T_0$")
ax.set_title(
    r"same rest-frame $T$ + recession Doppler  "
    r"($\beta_{\rm rec}=0.02$ at $R_{\rm out}$; $v \propto r$ from CoE)",
    loc="left", pad=8,
)
ax.grid(alpha=0.25)
ax.legend(loc="upper right", fontsize=11, framealpha=0.92)
tag(ax, 40, 1.012,
    "smooth anisotropy from\nline-of-sight velocity",
    "0.25", ha="left")

for ax in axes:
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])

axes[-1].set_xlabel(
    r"great-circle angle $\theta$ (degrees)      "
    r"$0^\circ$/$360^\circ$ = toward CoE      ·      "
    r"$180^\circ$ = away from CoE"
)
fig.suptitle(
    "Two concentric spherical shells about CoE  ·  "
    r"$R_{\rm in}=1$, $R_{\rm out}=3$, observer at $D=1.8$"
    "\n"
    r"Each surface emits uniform rest-frame $T$ toward the observer · "
    r"ray hits inner disk ($|\theta|<\alpha$) or outer shell"
    "\n"
    r"Equal $T$ $\Rightarrow$ isotropic.  "
    r"Anisotropy from unequal $T$ (hard step at $\alpha$) or Doppler $V_{\rm rec}$.",
    fontsize=13.5,
)
fig.tight_layout(rect=[0, 0, 1, 0.90])
p2 = os.path.join(HERE, "shells_T_vs_theta.png")
fig.savefig(p2, dpi=160)
plt.close(fig)
print("wrote", p2)


# =====================================================================
# Figure 3 — Mollweide sky maps for shells
# =====================================================================
# rebuild maps from the compute module logic (avoid huge JSON)
import importlib.util
spec = importlib.util.spec_from_file_location(
    "compute_two_emitters",
    os.path.join(HERE, "compute_two_emitters.py"),
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

cases = [
    ("same rest-frame T", dict(T_in=1.0, T_out=1.0, beta_rec=0.0, beta_obs_radial=0.0),
     "isotropic — flat map"),
    (r"$T_{{\rm in}}=1.1$, $T_{{\rm out}}=0.9$",
     dict(T_in=1.1, T_out=0.9, beta_rec=0.0, beta_obs_radial=0.0),
     "hot inner disk"),
    (r"Doppler $V_{{\rm rec}}$ (comoving obs.)",
     dict(T_in=1.0, T_out=1.0, beta_rec=0.02,
          beta_obs_radial=0.02 * (mod.D_SHELL / mod.R_OUT)),
     "smooth dipole-like"),
]

fig, axes = plt.subplots(1, 3, figsize=(18, 6.2),
                         subplot_kw=dict(projection="mollweide"))
for ax, (title, kw, note) in zip(axes, cases):
    Lon, Lat, Tmap, geo = mod.sky_map_shells(nside_lon=720, nside_lat=361, **kw)
    # mollweide wants lon in radians in (−π, π), lat in (−π/2, π/2)
    Lon_r = np.deg2rad(Lon)
    Lat_r = np.deg2rad(Lat)
    vmin, vmax = Tmap.min(), Tmap.max()
    if vmax - vmin < 1e-9:
        vmin, vmax = 0.99, 1.01  # show flat
    pcm = ax.pcolormesh(Lon_r, Lat_r, Tmap, cmap="RdBu_r",
                        norm=Normalize(vmin=vmin, vmax=vmax),
                        shading="auto")
    ax.set_title(title, pad=12, fontsize=13)
    # limb of inner shell: θ=α ⇒ lat = 90−α (circle about north pole)
    # draw as parallel at lat = 90−α
    a = geo["alpha_deg"]
    lat_limb = np.deg2rad(90.0 - a)
    lon_c = np.linspace(-np.pi, np.pi, 360)
    ax.plot(lon_c, np.full_like(lon_c, lat_limb), color="k", lw=1.0, ls="--", alpha=0.7)
    # south: the great-circle map has only one cone toward CoE = north
    tag_y = 0.55 if note.startswith("iso") else (-0.15 if "disk" in note else 0.15)
    ax.text(0.0, tag_y, note, transform=ax.transAxes, ha="center", va="center",
            fontsize=11.5,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="none", alpha=0.9))
    divider = make_axes_locatable(ax)
    # colorbar via fig — simpler:
    cbar = fig.colorbar(pcm, ax=ax, orientation="horizontal", pad=0.08, fraction=0.06)
    cbar.set_label(r"$T_{\rm obs}/T_0$")
    ax.grid(True, alpha=0.35)

fig.suptitle(
    "Shell sky maps (Mollweide)  ·  toward CoE = north pole  ·  "
    r"dashed = inner-shell limb $\alpha=\arcsin(R_{\rm in}/D)$"
    "\n"
    "Equal rest-frame $T$: blank isotropic sky.  "
    "Unequal $T$: sharp disk.  "
    "Recession Doppler: smooth anisotropy.",
    fontsize=13.5,
)
fig.tight_layout(rect=[0, 0, 1, 0.88])
p3 = os.path.join(HERE, "shells_mollweide.png")
fig.savefig(p3, dpi=160)
plt.close(fig)
print("wrote", p3)


# =====================================================================
# Figure 4 — Plates Mollweide (hemispheres)
# =====================================================================
fig, axes = plt.subplots(1, 3, figsize=(18, 6.2),
                         subplot_kw=dict(projection="mollweide"))
lon = np.linspace(-180, 180, 720, endpoint=False)
lat = np.linspace(-90, 90, 361)
Lon, Lat = np.meshgrid(lon, lat)
# +z = north pole ⇒ upper plate = northern hemisphere (lat>0)
mu = np.sin(np.deg2rad(Lat))  # cos(θ_from_+z) = sin(lat)

maps = [
    ("same T (any offset)", np.ones_like(mu), (0.99, 1.01), "isotropic — flat map"),
    (r"$T_{+}=1.1$, $T_{-}=0.9$",
     np.where(mu >= 0, 1.1, 0.9), (0.85, 1.15), "hemisphere step"),
    (r"same $T$ + $\beta=0.01$ along $+z$",
     mod.plates_T_obs(90.0 - Lat, z_obs=0.0, T_plus=1.0, T_minus=1.0, beta_obs=0.01),
     None, "classic dipole"),
]

for ax, (title, Tmap, vr, note) in zip(axes, maps):
    if vr is None:
        vmin, vmax = Tmap.min(), Tmap.max()
    else:
        vmin, vmax = vr
    pcm = ax.pcolormesh(np.deg2rad(Lon), np.deg2rad(Lat), Tmap,
                        cmap="RdBu_r", norm=Normalize(vmin=vmin, vmax=vmax),
                        shading="auto")
    ax.set_title(title, pad=12, fontsize=13)
    ax.grid(True, alpha=0.35)
    ax.text(0.0, 0.5, note, transform=ax.transAxes, ha="center", va="center",
            fontsize=11.5,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="none", alpha=0.9))
    cbar = fig.colorbar(pcm, ax=ax, orientation="horizontal", pad=0.08, fraction=0.06)
    cbar.set_label(r"$T_{\rm obs}/T_0$")

fig.suptitle(
    "Parallel-plate sky maps (Mollweide)  ·  $+z$ (upper plate) = north pole"
    "\n"
    "Infinite plates always fill exact hemispheres — offset never changes the map.  "
    "Equal $T$: isotropic.  Unequal $T$ / motion: step or dipole.",
    fontsize=13.5,
)
fig.tight_layout(rect=[0, 0, 1, 0.88])
p4 = os.path.join(HERE, "plates_mollweide.png")
fig.savefig(p4, dpi=160)
plt.close(fig)
print("wrote", p4)
print("done")
