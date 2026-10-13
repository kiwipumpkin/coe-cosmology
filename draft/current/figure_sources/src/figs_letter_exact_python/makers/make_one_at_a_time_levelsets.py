#!/usr/bin/env python3
"""
One-at-a-time A/A_all level-set curves vs Galactic longitude.

Each figure varies ONE geometry variable across several levels while the other
three stay fixed at baseline (NOT at heatmap-best partner values).

Baseline:
  D/L x = 5.9, d/D eps = 0.1, alpha = 1.0,
  kappa = 0.23058828928336284 (rescaled via scale_kappa when not varying Mass),
  f_u = 0.94, af_mode = "hratio",
  Amp_A = -0.06, phi_A = 0.35 rad
  Model A = (1 - Amp*sin(2πℓ/360+φ)) * AF

Four figures (A only, Migkas A points+errors):
  1. Vary Distance D/L
  2. Vary Radius of observer d/D
  3. Vary Velocity α
  4. Vary Mass κ
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path("/workspace/cosmology")
OUT = Path("/workspace/figs_letter_exact_python")
MAKERS = Path("/workspace/figs_letter_exact_python/makers")
sys.path.insert(0, str(ROOT / "density_local_match"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(MAKERS))

from run_density_local_match import (  # noqa: E402
    params_twoscale_local_match,
    load_migkas,
)

import importlib.util

_spec = importlib.util.spec_from_file_location(
    "amp_phi_geo", MAKERS / "make_amp_phi_and_geometry_heatmaps.py"
)
_hm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_hm)

scale_kappa = _hm.scale_kappa
af_an_alpha = _hm.af_an_alpha
model_AB_fixed = _hm.model_AB_fixed
style = _hm.style
BASE = dict(_hm.BASE)
SHEET_AMP_A = _hm.SHEET_AMP_A
SHEET_PHI_RAD = _hm.SHEET_PHI_RAD
SHEET_PHI_DEG = _hm.SHEET_PHI_DEG

KAPPA0 = 0.23058828928336284

# Level grids (include baseline as thicker highlighted orange curve)
LEVELS_DL = [3.0, 5.0, 5.9, 8.0, 10.0]
LEVELS_DD = [0.06, 0.10, 0.26, 0.50, 0.74]
LEVELS_ALPHA = [0.2, 0.5, 1.0, 1.3, 2.0]
LEVELS_KAPPA = [0.15, 0.23, 0.30, 0.45, 0.70]

KEY_LABEL = {
    "x": r"$D/L$",
    "eps": r"$d/D$",
    "alpha": r"$\alpha$",
    "kappa": r"$\kappa$",
}

KEY_FMT = {
    "x": lambda v: rf"$D/L={v:g}$",
    "eps": lambda v: rf"$d/D={v:g}$",
    "alpha": lambda v: rf"$\alpha={v:g}$",
    "kappa": lambda v: rf"$\kappa={v:g}$",
}

VARY_TITLE = {
    "x": "Vary Distance $D/L$",
    "eps": "Vary Radius of observer $d/D$",
    "alpha": r"Vary Velocity $\alpha$",
    "kappa": r"Vary Mass $\kappa$",
}


def style_ls():
    style()
    mpl.rcParams.update({"savefig.dpi": 220})


def baseline_hold() -> dict:
    """Full baseline for the three variables not being varied."""
    return {
        "x": float(BASE["x"]),
        "eps": float(BASE["eps"]),
        "alpha": float(BASE["alpha"]),
        "kappa": float(KAPPA0),
        "f_u": float(BASE["f_u"]),
        "af_mode": BASE["af_mode"],
    }


def model_A_curve(x, f_u, alpha, eps, kappa):
    """Return (l_deg, model_A) on fine longitude grid."""
    pmeta = params_twoscale_local_match(float(x), float(f_u))
    if kappa is not None:
        pmeta = scale_kappa(pmeta, float(kappa))
    lf, AF, AN = af_an_alpha(
        pmeta, eps=float(eps), alpha=float(alpha), af_mode=BASE["af_mode"]
    )
    ma, _mb = model_AB_fixed(AF, AN, lf)
    return lf, ma


def levelset_A_only(
    mig,
    vary_key: str,
    vary_vals: list[float],
    hold: dict,
    highlight: float,
    out_stem: str,
):
    """
    Single-panel A/A_all vs ℓ with Migkas points + level curves.
    Thick orange highlight for the baseline level of the varied axis.
    """
    l = mig["l"]
    cmap = plt.cm.viridis
    colors = [cmap(x) for x in np.linspace(0.08, 0.92, len(vary_vals))]

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ax.errorbar(
        l, mig["A"], yerr=[mig["A_elo"], mig["A_ehi"]],
        fmt="o", ms=3.5, color="k",
        ecolor="#999999", elinewidth=0.8, capsize=0, zorder=5,
        label="Migkas",
    )
    ax.axhline(1.0, color="#888888", ls="--", lw=0.9)

    held_parts = []
    for k in ("x", "eps", "alpha", "kappa"):
        if k == vary_key:
            continue
        if k in hold and hold[k] is not None:
            held_parts.append(KEY_FMT[k](hold[k]))
    held_str = ", ".join(held_parts)

    for val, c in zip(vary_vals, colors):
        kw = dict(hold)
        kw[vary_key] = float(val)
        lf, ma = model_A_curve(
            x=kw.get("x", BASE["x"]),
            f_u=kw.get("f_u", BASE["f_u"]),
            alpha=kw.get("alpha", BASE["alpha"]),
            eps=kw.get("eps", BASE["eps"]),
            kappa=kw.get("kappa", None),
        )
        is_hi = abs(float(val) - float(highlight)) < 1e-6
        lw = 2.8 if is_hi else 1.4
        col = "#E69F00" if is_hi else c
        lab = KEY_FMT[vary_key](val) + (" (baseline)" if is_hi else "")
        ax.plot(lf, ma, color=col, lw=lw, label=lab, zorder=4 if is_hi else 3)

    ax.set_xlabel(r"Galactic longitude $\ell$ [deg]")
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$")
    ax.set_xlim(0, 360)
    ax.set_ylim(0.80, 1.30)
    ax.grid(True, alpha=0.35)
    ax.legend(loc="upper right", fontsize=7.5, ncol=2, framealpha=0.92)
    ax.set_title(
        rf"{VARY_TITLE[vary_key]} (others at baseline)" "\n"
        rf"held: {held_str}; Amp$_A$={SHEET_AMP_A}, "
        rf"$\phi_A$={SHEET_PHI_RAD} rad ({SHEET_PHI_DEG:.1f}$^\circ$)",
        fontsize=11,
    )
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{out_stem}.{ext}", dpi=300)
    plt.close(fig)
    print(f"Wrote {out_stem}.png/.pdf  vary={vary_key} highlight={highlight}")


def main():
    style_ls()
    OUT.mkdir(parents=True, exist_ok=True)
    mig = load_migkas()

    # Confirm baseline kappa from params matches KAPPA0
    p0 = params_twoscale_local_match(float(BASE["x"]), float(BASE["f_u"]))
    print(f"params κ(x={BASE['x']})={float(p0['kappa']):.6g}; "
          f"using KAPPA0={KAPPA0:.6g}")

    report = {
        "amp_phi": {"Amp_A": SHEET_AMP_A, "phi_A_rad": SHEET_PHI_RAD},
        "baseline": {
            "x": BASE["x"],
            "eps": BASE["eps"],
            "alpha": BASE["alpha"],
            "kappa": KAPPA0,
            "f_u": BASE["f_u"],
            "af_mode": BASE["af_mode"],
        },
        "figures": [],
    }

    figures = [
        ("x", LEVELS_DL, 5.9, "levelset_vary_Distance"),
        ("eps", LEVELS_DD, 0.10, "levelset_vary_RadiusObs"),
        ("alpha", LEVELS_ALPHA, 1.0, "levelset_vary_Velocity"),
        ("kappa", LEVELS_KAPPA, 0.23, "levelset_vary_Mass"),
    ]

    for vary_key, levels, highlight, stem in figures:
        hold = baseline_hold()
        # Chris: observer-radius levels need higher held velocity.
        if vary_key == "eps":
            hold["alpha"] = 1.3
        # When varying kappa, still pass kappa through levels; hold keeps
        # x/eps/alpha at baseline. When not varying kappa, hold["kappa"]=KAPPA0
        # so Mass stays fixed while Distance (x) changes.
        levelset_A_only(
            mig, vary_key, levels, hold, highlight=highlight, out_stem=stem,
        )
        report["figures"].append({
            "stem": stem,
            "vary": vary_key,
            "levels": list(levels),
            "highlight": highlight,
            "hold": {k: hold[k] for k in ("x", "eps", "alpha", "kappa")},
        })

    expected = [stem for _, _, _, stem in figures]
    verify = {}
    for stem in expected:
        png = OUT / f"{stem}.png"
        pdf = OUT / f"{stem}.pdf"
        verify[stem] = {
            "png_exists": png.exists(),
            "pdf_exists": pdf.exists(),
            "png_bytes": png.stat().st_size if png.exists() else 0,
            "pdf_bytes": pdf.stat().st_size if pdf.exists() else 0,
        }
        print(f"  {stem}: png={verify[stem]['png_bytes']} pdf={verify[stem]['pdf_bytes']}")

    report["verify"] = verify
    out_json = OUT / "levelset_one_at_a_time_report.json"

    def _ser(obj):
        if isinstance(obj, dict):
            return {k: _ser(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_ser(v) for v in obj]
        if isinstance(obj, (np.floating, float)):
            return float(obj)
        if isinstance(obj, (np.integer, int)):
            return int(obj)
        return obj

    out_json.write_text(json.dumps(_ser(report), indent=2))
    print(f"Wrote {out_json}")
    return report


if __name__ == "__main__":
    main()
