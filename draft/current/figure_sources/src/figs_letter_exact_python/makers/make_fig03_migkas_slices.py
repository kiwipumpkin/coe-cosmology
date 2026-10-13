#!/usr/bin/env python3
"""Fig03 — Migkas redshift slices + D–m_cc degeneracy.

PRIMARY (letter-era): cmb_reconcile/termslice/make_fig.py
Writes termslice_summary.png (2x2: three A/slice-mean panels + Degeneracy),
then copies to Figs_from_makers/Fig03_migkas_slices.png.

Do NOT use coe-figures/src/fig03_migkas_slices.py (wrong reconstruction).
"""
from __future__ import annotations

import os
import runpy
import shutil
from pathlib import Path

OUT = Path("/workspace/figs_letter_exact_python")
OUT_NAME = "Fig03_migkas_slices.png"
TERMSLICE = Path(
    "/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/termslice"
)
SCRIPT = TERMSLICE / "make_fig.py"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # make_fig.py uses relative imports (../cycle3, solve.json, …) — must run in termslice cwd
    prev = os.getcwd()
    try:
        os.chdir(TERMSLICE)
        runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        os.chdir(prev)
    src = TERMSLICE / "termslice_summary.png"
    if not src.is_file():
        raise SystemExit(f"missing output from letter-era maker: {src}")
    dest = OUT / OUT_NAME
    shutil.copy2(src, dest)
    print("wrote", dest, dest.stat().st_size)
    return dest


if __name__ == "__main__":
    main()
