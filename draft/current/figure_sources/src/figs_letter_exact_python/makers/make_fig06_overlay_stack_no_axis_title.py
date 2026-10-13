#!/usr/bin/env python3
"""Fig06 overlay stack with baked-in 'Galactic latitude' axis title removed.

Method (a): white out the left title strip on a same-size copy of
FULL_migkas_theta75.png before overlay so GOOD_REG still registers.
Does NOT modify make_migkas_fig8_overlay_matched.py or CoE_calculations.md
or the letter-exact make_fig06_overlay_stack.py.

Ticks/numbers (+90°, …) are left intact; only the soft vertical title
string (cols ~70–81, AA ~55–100) is painted white.
"""
from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

OUT = Path("/workspace/figs_letter_exact_python")
MAKERS = Path(__file__).resolve().parent
SM = Path("/workspace/cosmology/papers/spreadsheet_model")
OVERLAY_PY = Path("/workspace/cosmology/figures/make_migkas_fig8_overlay_matched.py")
SOURCE_FULL = Path("/workspace/cosmology/map_compare/FULL_migkas_theta75.png")
SOURCE_NO_TITLE = MAKERS / "_FULL_migkas_theta75_no_axis_title.png"

# Title column in FULL_migkas_theta75.png (1975×1353):
# dark title glyphs ~x=70–81, y=622–885; AA soft pixels ~x=55–100, y=620–888.
# Tick numbers begin at x≈108 — keep them.
TITLE_X0, TITLE_X1 = 50, 100
TITLE_Y0, TITLE_Y1 = 600, 910

PANEL_NAMES = [
    "fig8_overlay_matched.png",
    "fig8_overlay_matched_grad.png",
    "fig8_overlay_matched_grad_rot81.png",
]
PANEL_LABELS = [
    "+15° tilt",
    "+15° tilt + gradient",
    "+15° tilt + gradient + 80.8° roll",
]
OUT_PANELS = [
    "Fig06_overlay_matched.png",
    "Fig06_overlay_matched_grad.png",
    "Fig06_overlay_matched_grad_rot81.png",
]
OUT_STACK = "Fig06_overlay_stack.png"
OUT_STACK_ALIAS = "Fig06_overlay_stack_NO_AXIS_TITLE.png"
QA_LEFT = "_qa_fig06_no_title_left.png"


def make_whitened_source() -> Path:
    """Paint white over the vertical 'Galactic latitude' strip; keep dims."""
    img = Image.open(SOURCE_FULL).convert("RGBA")
    arr = np.array(img)
    arr[TITLE_Y0:TITLE_Y1, TITLE_X0:TITLE_X1, :3] = 255
    arr[TITLE_Y0:TITLE_Y1, TITLE_X0:TITLE_X1, 3] = 255
    out = Image.fromarray(arr)
    assert out.size == img.size, (out.size, img.size)
    SOURCE_NO_TITLE.parent.mkdir(parents=True, exist_ok=True)
    out.save(SOURCE_NO_TITLE)
    gray = arr[:, :, :3].mean(axis=2)
    title_dark = int((gray[TITLE_Y0:TITLE_Y1, TITLE_X0:TITLE_X1] < 200).sum())
    print(
        f"whitened source {SOURCE_NO_TITLE} size={out.size} "
        f"title_strip_dark_px={title_dark} "
        f"rect=({TITLE_X0},{TITLE_Y0})-({TITLE_X1},{TITLE_Y1})"
    )
    if title_dark != 0:
        raise RuntimeError(f"title strip still has {title_dark} dark pixels")
    return SOURCE_NO_TITLE


def _run_overlay_on_whitened(whitened: Path):
    """Load overlay engine; point source at whitened PNG; SHOW_LEGEND=False."""
    spec = importlib.util.spec_from_file_location("migkas_overlay_matched_notitle", OVERLAY_PY)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    mod.SOURCE_CANDIDATES = [whitened] + [
        p for p in mod.SOURCE_CANDIDATES if Path(p) != whitened
    ]
    mod.SHOW_LEGEND = False
    mod.VERSIONS = [
        ("matched", False, 0.0),
        ("matched_grad", True, 0.0),
        ("matched_grad_rot81", True, 80.8),
    ]
    mod.main()
    return mod


def _stack_panels():
    OUT.mkdir(parents=True, exist_ok=True)
    parts = []
    dims = []
    for src_name, out_name in zip(PANEL_NAMES, OUT_PANELS):
        src = SM / src_name
        if not src.exists():
            raise FileNotFoundError(src)
        dest = OUT / out_name
        shutil.copy2(src, dest)
        im = Image.open(src).convert("RGB")
        parts.append(im)
        dims.append((out_name, im.size))
        print("panel", dest, im.size)

    try:
        font = ImageFont.truetype(
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 28
        )
    except Exception:
        font = ImageFont.load_default()

    label_h = 48
    gap = 12
    W = max(p.width for p in parts)
    H = sum(p.height + label_h + gap for p in parts) + 20
    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    y = 10
    for lab, p in zip(PANEL_LABELS, parts):
        bbox = draw.textbbox((0, 0), lab, font=font)
        tw = bbox[2] - bbox[0]
        draw.text(((W - tw) // 2, y + 8), lab, fill="black", font=font)
        y += label_h
        canvas.paste(p, ((W - p.width) // 2, y))
        y += p.height + gap

    stack = canvas.crop((0, 0, W, y))
    dest = OUT / OUT_STACK
    stack.save(dest)
    alias = OUT / OUT_STACK_ALIAS
    shutil.copy2(dest, alias)
    print("wrote stack", dest, stack.size, dest.stat().st_size)
    print("wrote alias", alias)
    return dest, dims, stack.size


def _qa_left_strip(stack_path: Path):
    """Crop left margin of the stacked figure to prove title string is gone."""
    stack = Image.open(stack_path).convert("RGB")
    # Wide enough to show ticks if present; title lived in far-left ~50–100 of source
    # which maps into the left margin of each matplotlib panel.
    qa = stack.crop((0, 0, min(280, stack.width), stack.height))
    dest = OUT / QA_LEFT
    qa.save(dest)
    arr = np.asarray(qa)
    gray = arr.mean(axis=2)
    # Heuristic: look for the soft vertical title column (no "Galactic" glyphs)
    # Scan left 120 px for tall vertical dark runs typical of rotated title.
    left = gray[:, :120]
    dark = left < 80
    # Count dark pixels that form a narrow vertical band (title signature)
    col_counts = dark.sum(axis=0)
    print(
        "QA left strip", dest, qa.size,
        "dark_in_left120=", int(dark.sum()),
        "max_col_dark=", int(col_counts.max()) if len(col_counts) else 0,
    )
    return dest


def main():
    print("=== Fig06 no-axis-title: whiten source title strip ===")
    whitened = make_whitened_source()
    print("=== regenerating Migkas Fig.8 overlays on whitened source ===")
    _run_overlay_on_whitened(whitened)
    print("=== stacking panels ===")
    stack_path, dims, stack_size = _stack_panels()
    print("=== QA left strip ===")
    qa = _qa_left_strip(stack_path)
    print("DONE")
    print("panels:", dims)
    print("stack:", stack_path, stack_size)
    print("qa:", qa)
    return stack_path


if __name__ == "__main__":
    main()
