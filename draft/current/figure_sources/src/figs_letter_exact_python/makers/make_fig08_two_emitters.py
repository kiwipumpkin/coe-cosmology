#!/usr/bin/env python3
"""Fig08 — Two-emitter baseline regenerated from cycle6 JSON + schematic.

Primary: cycle6_two_emitters/make_two_emitters_comoving.py
Fallback: self-contained schematic + T_obs from two_emitters.json.
"""
from __future__ import annotations

import json
import os
import runpy
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

OUT = Path("/workspace/figs_letter_exact_python")
CYCLE6 = Path(
    "/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle6_two_emitters"
)
SCRIPT = CYCLE6 / "make_two_emitters_comoving.py"
JSON_PATH = CYCLE6 / "two_emitters.json"
OUT_NAME = "Fig08_two_emitters.png"


def _fallback():
    data = json.loads(JSON_PATH.read_text()) if JSON_PATH.exists() else {}
    fig = plt.figure(figsize=(7.4, 11.9), dpi=160)  # portrait ~letter 1181×1908
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.25)
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-2.2, 2.4)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Shared bulk velocity (schematic)", loc="left", fontsize=12)
    ax.add_patch(Rectangle((-2.8, 1.2), 5.6, 0.35, color="#5dade2", alpha=0.85, ec="#1a5276"))
    ax.add_patch(Rectangle((-2.8, -1.55), 5.6, 0.35, color="#5dade2", alpha=0.85, ec="#1a5276"))
    ax.text(0, 1.55, "+ plate / shell", ha="center", fontsize=9, color="#1a5276")
    ax.text(0, -1.7, "− plate / shell", ha="center", fontsize=9, color="#1a5276")
    ax.plot(0, 0, "o", color="#2471a3", ms=10)
    ax.text(0.15, -0.35, "observer", color="#2471a3", fontsize=9)
    for y in (1.38, 0.0, -1.38):
        ax.annotate("", xy=(1.4, y), xytext=(-0.2, y),
                    arrowprops=dict(arrowstyle="->", color="#c0392b", lw=2.2))
    ax.text(1.55, 0.15, r"shared $\vec{V}$", color="#c0392b", fontsize=10)
    P = data.get("plates", {})
    th = np.asarray(P.get("theta_deg", np.linspace(0, 360, 361)), float)
    mid = np.asarray(P.get("mid_sameT", np.ones_like(th)), float)
    ax2.plot(th, mid, color="C0", lw=2.4, label=r"same $T$, shared $V$ (mid)")
    if "mid_sameT_doppler" in P:
        ax2.plot(th, P["mid_sameT_doppler"], color="C2", lw=1.6, ls="--",
                 label=r"observer boost $\beta=0.01$ (dipole)")
    ax2.axhline(1.0, color="0.5", lw=0.8, ls=":")
    ax2.set_xlim(0, 360)
    ax2.set_xlabel(r"great-circle angle $\theta$ (deg)")
    ax2.set_ylabel(r"$T_{\mathrm{obs}}/T_0$")
    ax2.set_title("Kinematic sky before gravity", loc="left", fontsize=12)
    ax2.legend(fontsize=8, frameon=False)
    ax2.grid(True, alpha=0.25)
    fig.suptitle("Two-emitter baseline", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    dest = OUT / OUT_NAME
    fig.savefig(dest, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return dest


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if SCRIPT.exists():
        cwd = os.getcwd()
        try:
            os.chdir(CYCLE6)
            runpy.run_path(str(SCRIPT), run_name="__main__")
        finally:
            os.chdir(cwd)
        cands = sorted(CYCLE6.glob("*.png"), key=lambda p: p.stat().st_mtime, reverse=True)
        prefer = [p for p in cands if "emitter" in p.name.lower() or "two" in p.name.lower()]
        src = (prefer or cands)[0] if (prefer or cands) else None
        if src is not None:
            dest = OUT / OUT_NAME
            shutil.copy2(src, dest)
            print("wrote", dest, "from", src.name)
            return dest
    return _fallback()


if __name__ == "__main__":
    main()

