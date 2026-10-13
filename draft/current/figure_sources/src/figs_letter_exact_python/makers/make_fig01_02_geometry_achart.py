#!/usr/bin/env python3
"""Fig01 geometry + Fig02 Achart — rerun spreadsheet_model/bubble/universe_thin.py.

Pixel-matches letter gold embeds (PNG encoding may differ; RGBA max abs == 0).
"""
from __future__ import annotations

import runpy
import shutil
from pathlib import Path

OUT = Path("/workspace/figs_letter_exact_python")
BUBBLE = Path("/workspace/cosmology/papers/spreadsheet_model/bubble")
SM = BUBBLE.parent
SCRIPT = BUBBLE / "universe_thin.py"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # universe_thin.py writes ../universe_thin_drawing.png and ../universe_thin_Achart.png
    # relative to bubble/; run in that directory via runpy with argv.
    import os, sys
    cwd = os.getcwd()
    sys.path.insert(0, str(BUBBLE))
    try:
        os.chdir(BUBBLE)
        runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        os.chdir(cwd)
    src1 = SM / "universe_thin_drawing.png"
    src2 = SM / "universe_thin_Achart.png"
    d1 = OUT / "Fig01_geometry.png"
    d2 = OUT / "Fig02_Achart.png"
    shutil.copy2(src1, d1)
    shutil.copy2(src2, d2)
    print("wrote", d1, d1.stat().st_size)
    print("wrote", d2, d2.stat().st_size)
    return d1, d2


if __name__ == "__main__":
    main()

