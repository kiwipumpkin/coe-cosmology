#!/usr/bin/env python3
"""
Four-variable R² heatmaps vs Migkas Fig.6 (Chris CoE / Migkas fit).

Exactly four physical variables; each heatmap = exactly TWO on the axes.
All C(4,2)=6 pairs × A/B = 12 heatmaps.

Variables:
  1. Distance — D/L (model param x). Label: Distance $D/L$
  2. Radius of observer — d/D (eps). Label: Radius of observer $d/D$
  3. Velocity — α. Label: Velocity scale $\\alpha$
  4. Mass — κ = GM(<D)/(D³H²). Label: Mass $\\kappa = GM(\\lt D)/(D^3 H^2)$

When Mass is held: fix κ to baseline κ(x=5.9) so Distance and Mass stay independent.
"""
from __future__ import annotations

import json
import sys
import time
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

# Reuse helpers from existing Amp×φ / geometry maker
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "amp_phi_geo", MAKERS / "make_amp_phi_and_geometry_heatmaps.py"
)
_hm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_hm)

scale_kappa = _hm.scale_kappa
af_an_alpha = _hm.af_an_alpha
model_AB_fixed = _hm.model_AB_fixed
r2_geometry_cell = _hm.r2_geometry_cell
heat_ab = _hm.heat_ab
edges = _hm.edges
style = _hm.style
BASE = dict(_hm.BASE)
SHEET_AMP_A = _hm.SHEET_AMP_A
SHEET_PHI_DEG = _hm.SHEET_PHI_DEG
BEST_AMP_B = _hm.BEST_AMP_B
BEST_PHI_B_DEG = _hm.BEST_PHI_B_DEG

# Axis labels (Chris names)
LBL_DIST = r"Distance $D/L$"
LBL_RAD = r"Radius of observer $d/D$"
LBL_VEL = r"Velocity scale $\alpha$"
LBL_MASS = r"Mass $\kappa = GM(<D)/(D^3 H^2)$"


def style_four():
    style()
    mpl.rcParams.update({"savefig.dpi": 220})


def hold_title(held: dict) -> str:
    parts = []
    if "x" in held:
        parts.append(rf"$D/L={held['x']}$")
    if "eps" in held:
        parts.append(rf"$d/D={held['eps']}$")
    if "alpha" in held:
        parts.append(rf"$\alpha={held['alpha']}$")
    if "kappa" in held:
        parts.append(rf"$\kappa={held['kappa']:.4g}$ (held)")
    parts.append(r"AF=hratio, fixed Amp×φ")
    return ", ".join(parts)


def scan_pair(mig, x_vals, y_vals, x_key, y_key, hold: dict, name, xlab, ylab, title,
              marks=(), hlines=(), vlines=()):
    """
    Scan a 2D grid. keys in {x, eps, alpha, kappa}.
    hold: fixed values for the other two (may include kappa baseline).
    """
    ZA = np.full((len(y_vals), len(x_vals)), np.nan)
    ZB = np.full_like(ZA, np.nan)
    t0 = time.time()
    n = len(y_vals) * len(x_vals)
    print(f"Scanning {name}: {len(x_vals)}×{len(y_vals)}={n} ...", flush=True)
    for j, yv in enumerate(y_vals):
        for i, xv in enumerate(x_vals):
            kw = dict(hold)
            kw[x_key] = float(xv)
            kw[y_key] = float(yv)
            rec = r2_geometry_cell(
                mig,
                x=kw.get("x", BASE["x"]),
                f_u=kw.get("f_u", BASE["f_u"]),
                alpha=kw.get("alpha", BASE["alpha"]),
                eps=kw.get("eps", BASE["eps"]),
                kappa=kw.get("kappa", None),
                af_mode=kw.get("af_mode", BASE["af_mode"]),
            )
            ZA[j, i] = rec["R2_A"]
            ZB[j, i] = rec["R2_B"]
        if (j + 1) % max(1, len(y_vals) // 4) == 0 or j == len(y_vals) - 1:
            print(
                f"  {name} row {j+1}/{len(y_vals)}  ({time.time()-t0:.1f}s)",
                flush=True,
            )
    meta = heat_ab(
        x_vals, y_vals, ZA, ZB, xlab, ylab, title, name,
        marks=marks, hlines=hlines, vlines=vlines,
        out_prefix="heatmap",
    )
    meta["ZA_max"] = float(np.nanmax(ZA))
    meta["ZB_max"] = float(np.nanmax(ZB))
    meta["grid"] = {
        "x_key": x_key, "y_key": y_key,
        "nx": len(x_vals), "ny": len(y_vals),
        "hold": {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                 for k, v in hold.items()},
    }
    meta["runtime_s"] = time.time() - t0
    return meta


def main():
    style_four()
    OUT.mkdir(parents=True, exist_ok=True)
    mig = load_migkas()

    # Baseline κ at x=5.9 so Mass can be held independent of Distance
    p0 = params_twoscale_local_match(float(BASE["x"]), float(BASE["f_u"]))
    KAPPA0 = float(p0["kappa"])
    print(f"Baseline κ(x={BASE['x']}, f_u={BASE['f_u']}) = {KAPPA0:.6g}")
    print(f"Amp/φ fixed: A Amp={SHEET_AMP_A} φ={SHEET_PHI_DEG:.2f}°; "
          f"B Amp={BEST_AMP_B} φ={BEST_PHI_B_DEG}°")

    # Grids
    x_g = np.round(np.arange(3.0, 10.001, 0.25), 3)       # D/L
    eps_g = np.round(np.arange(0.02, 0.8001, 0.04), 3)    # d/D
    alpha_g = np.round(np.arange(0.2, 2.001, 0.1), 3)     # α
    kappa_g = np.round(np.arange(0.05, 1.0001, 0.05), 4)  # κ

    print(f"Grids: D/L={len(x_g)}, d/D={len(eps_g)}, α={len(alpha_g)}, κ={len(kappa_g)}")
    t_all = time.time()
    report = {
        "baseline": dict(BASE),
        "baseline_kappa": KAPPA0,
        "amp_phi": {
            "Amp_A": SHEET_AMP_A, "phi_A_deg": SHEET_PHI_DEG,
            "Amp_B": BEST_AMP_B, "phi_B_deg": BEST_PHI_B_DEG,
        },
        "grids": {
            "D/L": [float(x_g[0]), float(x_g[-1]), float(np.diff(x_g).mean()), len(x_g)],
            "d/D": [float(eps_g[0]), float(eps_g[-1]), float(np.diff(eps_g).mean()), len(eps_g)],
            "alpha": [float(alpha_g[0]), float(alpha_g[-1]), float(np.diff(alpha_g).mean()), len(alpha_g)],
            "kappa": [float(kappa_g[0]), float(kappa_g[-1]), float(np.diff(kappa_g).mean()), len(kappa_g)],
        },
        "pairs": {},
    }

    # 1. Distance × Mass  (D/L × κ); hold eps, alpha
    # When Mass is ON axis we vary κ via scale_kappa; Distance varies x.
    report["pairs"]["Distance_Mass"] = scan_pair(
        mig, x_g, kappa_g, "x", "kappa",
        hold={"eps": BASE["eps"], "alpha": BASE["alpha"], "f_u": BASE["f_u"],
              "af_mode": BASE["af_mode"]},
        name="Distance_Mass",
        xlab=LBL_DIST, ylab=LBL_MASS,
        title="Distance × Mass\n" + hold_title({"eps": BASE["eps"], "alpha": BASE["alpha"]}),
        marks=[(BASE["x"], KAPPA0, "o",
                rf"baseline $D/L={BASE['x']}$, $\kappa={KAPPA0:.3f}$")],
        hlines=[(0.15, r"$\kappa=\Omega_m/2=0.15$")],
    )

    # 2. Radius of observer × Mass  (d/D × κ); hold x (and thus native κ unless we set κ), alpha
    # Hold Distance at baseline x; Mass varies on axis so κ is free.
    report["pairs"]["RadiusObs_Mass"] = scan_pair(
        mig, eps_g, kappa_g, "eps", "kappa",
        hold={"x": BASE["x"], "alpha": BASE["alpha"], "f_u": BASE["f_u"],
              "af_mode": BASE["af_mode"]},
        name="RadiusObs_Mass",
        xlab=LBL_RAD, ylab=LBL_MASS,
        title="Radius of observer × Mass\n" + hold_title(
            {"x": BASE["x"], "alpha": BASE["alpha"]}),
        marks=[(BASE["eps"], KAPPA0, "o",
                rf"baseline $d/D={BASE['eps']}$, $\kappa={KAPPA0:.3f}$")],
        hlines=[(0.15, r"$\kappa=\Omega_m/2=0.15$")],
    )

    # 3. Velocity × Radius of observer  (α × d/D); hold Distance + Mass (κ=KAPPA0)
    report["pairs"]["Velocity_RadiusObs"] = scan_pair(
        mig, alpha_g, eps_g, "alpha", "eps",
        hold={"x": BASE["x"], "kappa": KAPPA0, "f_u": BASE["f_u"],
              "af_mode": BASE["af_mode"]},
        name="Velocity_RadiusObs",
        xlab=LBL_VEL, ylab=LBL_RAD,
        title="Velocity × Radius of observer\n" + hold_title(
            {"x": BASE["x"], "kappa": KAPPA0}),
        marks=[(BASE["alpha"], BASE["eps"], "o",
                rf"baseline $\alpha={BASE['alpha']}$, $d/D={BASE['eps']}$")],
    )

    # 4. Distance × Radius of observer  (D/L × d/D); hold Velocity + Mass (κ=KAPPA0)
    # IMPORTANT: rescale every D/L cell to KAPPA0 so Mass stays independent of Distance
    report["pairs"]["Distance_RadiusObs"] = scan_pair(
        mig, x_g, eps_g, "x", "eps",
        hold={"alpha": BASE["alpha"], "kappa": KAPPA0, "f_u": BASE["f_u"],
              "af_mode": BASE["af_mode"]},
        name="Distance_RadiusObs",
        xlab=LBL_DIST, ylab=LBL_RAD,
        title="Distance × Radius of observer\n" + hold_title(
            {"alpha": BASE["alpha"], "kappa": KAPPA0}),
        marks=[(BASE["x"], BASE["eps"], "o",
                rf"baseline $D/L={BASE['x']}$, $d/D={BASE['eps']}$")],
    )

    # 5. Distance × Velocity  (D/L × α); hold Radius + Mass (κ=KAPPA0)
    report["pairs"]["Distance_Velocity"] = scan_pair(
        mig, x_g, alpha_g, "x", "alpha",
        hold={"eps": BASE["eps"], "kappa": KAPPA0, "f_u": BASE["f_u"],
              "af_mode": BASE["af_mode"]},
        name="Distance_Velocity",
        xlab=LBL_DIST, ylab=LBL_VEL,
        title="Distance × Velocity\n" + hold_title(
            {"eps": BASE["eps"], "kappa": KAPPA0}),
        marks=[(BASE["x"], BASE["alpha"], "o",
                rf"baseline $D/L={BASE['x']}$, $\alpha={BASE['alpha']}$")],
    )

    # 6. Mass × Velocity  (κ × α); hold Distance + Radius
    report["pairs"]["Mass_Velocity"] = scan_pair(
        mig, kappa_g, alpha_g, "kappa", "alpha",
        hold={"x": BASE["x"], "eps": BASE["eps"], "f_u": BASE["f_u"],
              "af_mode": BASE["af_mode"]},
        name="Mass_Velocity",
        xlab=LBL_MASS, ylab=LBL_VEL,
        title="Mass × Velocity\n" + hold_title(
            {"x": BASE["x"], "eps": BASE["eps"]}),
        marks=[(KAPPA0, BASE["alpha"], "o",
                rf"baseline $\kappa={KAPPA0:.3f}$, $\alpha={BASE['alpha']}$")],
        vlines=[(0.15, r"$\kappa=\Omega_m/2=0.15$")],
    )

    report["runtime_s"] = time.time() - t_all

    # Expected PNG paths
    expected = []
    for stem in (
        "Distance_Mass", "RadiusObs_Mass", "Velocity_RadiusObs",
        "Distance_RadiusObs", "Distance_Velocity", "Mass_Velocity",
    ):
        for tag in ("A", "B"):
            expected.append(str(OUT / f"heatmap_{stem}_{tag}.png"))
    report["png_paths"] = expected
    report["png_exist"] = {p: Path(p).is_file() for p in expected}

    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()
                    if not isinstance(v, np.ndarray)}
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, Path):
            return str(o)
        return o

    report = clean(report)
    out_json = OUT / "heatmap_four_vars_report.json"
    json.dump(report, open(out_json, "w"), indent=2)
    print(f"\nWrote {out_json}")
    print(f"Total runtime {report['runtime_s']:.1f}s")
    missing = [p for p, ok in report["png_exist"].items() if not ok]
    if missing:
        print("MISSING:", missing)
    else:
        print(f"All {len(expected)} PNGs present.")
    # compact best summary
    print("\n=== BEST PER PANEL ===")
    for pname, meta in report["pairs"].items():
        for tag in ("A", "B"):
            m = meta[tag]
            print(f"  {pname}_{tag}: x={m['best_x']:.4g} y={m['best_y']:.4g} R2={m['R2']:.4f}")
    return report


if __name__ == "__main__":
    main()
