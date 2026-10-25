# FIGURE_CODEBOOK_exact — SUPERSEDED

> **SUPERSEDED 2026-10-01:** Chris rejected copy-asset delivery.
> **Authoritative regenerating cookbook:** [`FIGURE_CODEBOOK_makers.md`](FIGURE_CODEBOOK_makers.md)
> (also `/workspace/figs_from_makers/FIGURE_CODEBOOK_makers.md`).
> Outputs: `/workspace/figs_from_makers/`. Run: `python3 /workspace/figs_from_makers/makers/run_all.py`.
> This file is retained only as historical inventory of the prior pixel-copy pass.

---


**Audience:** another bot / Retard Bot — code cookbook (not a paper doc).  
**Goal:** pixel-match letter PDF embeds. Gold extracts live in  
`/workspace/figs_exact_from_code/gold_extracts/` (from `pdfimages -all` on the Sep 30 2026 14-page letter; soft masks skipped).  
**Authoritative delivery:** `/workspace/figs_exact_from_code/delivery/` (see `INVENTORY.md` there).  
**Do not modify/delete:** `letter_brief/CoE_calculations.md`.

Honesty rule (from delivery INVENTORY): prefer re-run of original script when it pixel-matches; else **copy** the verified on-disk asset. Do **not** invent fake scripts that “should” match. Several letter-era one-offs are missing; current `make_*.py` / `coe-figures` recomputes often differ in size/layout.

---

## Table of contents

1. [Overview & verify recipe](#1-overview--verify-recipe)
2. [Per-figure cookbook](#2-per-figure-cookbook)
3. [Full sources — pixel-matching reruns](#3-full-sources--pixel-matching-reruns)
4. [Full sources — related (panels / non-exact)](#4-full-sources--related-panels--non-exact)
5. [Full sources — best-available approximates (do not pixel-match)](#5-full-sources--best-available-approximates-do-not-pixel-match)
6. [Dependencies](#6-dependencies)
7. [Gold embed map](#7-gold-embed-map)

---

## 1. Overview & verify recipe

| | |
|--|--|
| Delivery root | `/workspace/figs_exact_from_code/delivery/` |
| Gold extracts | `/workspace/figs_exact_from_code/gold_extracts/` |
| Letter refs | `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/refs_from_chris/` |
| Letter figs tree | `.../letter_brief/figs/` |

### Pixel / bitwise verification

```bash
python3 - <<'PY'
from pathlib import Path
import hashlib
import numpy as np
from PIL import Image

def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()

def max_abs(a, b):
    A = np.asarray(Image.open(a).convert("RGBA"), dtype=np.int16)
    B = np.asarray(Image.open(b).convert("RGBA"), dtype=np.int16)
    assert A.shape == B.shape, (A.shape, B.shape)
    d = np.abs(A - B)
    return int(d.max()), float(d.mean())

gold = Path("/workspace/figs_exact_from_code/gold_extracts")
deliv = Path("/workspace/figs_exact_from_code/delivery")
print("md5 gold", md5(gold / "gold-000.png"))
print("md5 deliv", md5(deliv / "Fig01_geometry.png"))
print("max/mean", max_abs(gold / "gold-000.png", deliv / "Fig01_geometry.png"))
# JPEG bitwise: md5 equality is sufficient when dims match
PY
```

- **PNG:** require `max abs diff == 0` after RGBA align (encoding may differ → md5 may differ).
- **JPEG (Fig05, Fig09, FigA2, FigA2b, comic):** require **md5 equality** with gold extract (bitwise).

---

## 2. Per-figure cookbook

### Fig01 — `Fig01_geometry.png`

- **Label:** Fig1 Geometry CoE/shell/vision
- **Gold extract:** `gold-000.png`
- **Expected dims:** 1460×488
- **Method (delivery):** `rerun:universe_thin.py`
- **Pixel=gold:** True · **Bitwise=gold:** False
- **Gold md5:** `ef40ce3d854574591d09d1edb96a38cc`
- **Delivered md5:** `d9fd6acd60b35a72f61515b87e18ee95`
- **Verified matching asset:** `/workspace/cosmology/papers/spreadsheet_model/universe_thin_drawing.png`
- **Notes:** Pixel-identical to PDF embed; bitwise-identical to spreadsheet_model/universe_thin_drawing.png. Script also writes Achart (Fig02).

- **Generating script (pixel-match on rerun):** `/workspace/cosmology/papers/spreadsheet_model/bubble/universe_thin.py` — full source in §3.

**Shell commands:**

```bash
cd /workspace/cosmology/papers/spreadsheet_model/bubble
python universe_thin.py
# outputs: ../universe_thin_drawing.png and ../universe_thin_Achart.png
cp ../universe_thin_drawing.png /workspace/figs_exact_from_code/delivery/Fig01_geometry.png
```

**Verify:**

```bash
python3 - <<'PY'
from PIL import Image; import numpy as np
A=np.asarray(Image.open('/workspace/figs_exact_from_code/gold_extracts/gold-000.png').convert('RGBA'),np.int16)
B=np.asarray(Image.open('/workspace/figs_exact_from_code/delivery/Fig01_geometry.png').convert('RGBA'),np.int16)
assert A.shape==B.shape, (A.shape,B.shape)
d=np.abs(A-B); print('max',int(d.max()),'mean',float(d.mean())); assert d.max()==0
PY
```

---

### Fig02 — `Fig02_Achart.png`

- **Label:** Fig2 Migkas A/A_all + model curves
- **Gold extract:** `gold-002.png`
- **Expected dims:** 1182×1287
- **Method (delivery):** `rerun:universe_thin.py`
- **Pixel=gold:** True · **Bitwise=gold:** False
- **Gold md5:** `5178d14e2707186ad6926d3a92ef20f7`
- **Delivered md5:** `6dac4950b3b907e5b421cd1a256df2f2`
- **Verified matching asset:** `/workspace/cosmology/papers/spreadsheet_model/universe_thin_Achart.png`
- **Notes:** Same universe_thin.py run as Fig01.

- **Generating script (pixel-match on rerun):** `/workspace/cosmology/papers/spreadsheet_model/bubble/universe_thin.py` — full source in §3.

**Shell commands:**

```bash
cd /workspace/cosmology/papers/spreadsheet_model/bubble
python universe_thin.py
cp ../universe_thin_Achart.png /workspace/figs_exact_from_code/delivery/Fig02_Achart.png
```

**Verify:**

```bash
python3 - <<'PY'
from PIL import Image; import numpy as np
A=np.asarray(Image.open('/workspace/figs_exact_from_code/gold_extracts/gold-002.png').convert('RGBA'),np.int16)
B=np.asarray(Image.open('/workspace/figs_exact_from_code/delivery/Fig02_Achart.png').convert('RGBA'),np.int16)
assert A.shape==B.shape, (A.shape,B.shape)
d=np.abs(A-B); print('max',int(d.max()),'mean',float(d.mean())); assert d.max()==0
PY
```

---

### Fig03 — `Fig03_migkas_slices.png`

- **Label:** Fig3 Migkas redshift slices + D-m_cc degeneracy
- **Gold extract:** `gold-004.png`
- **Expected dims:** 1123×1014
- **Method (delivery):** `copy:exact_on_disk_asset`
- **Pixel=gold:** True · **Bitwise=gold:** False
- **Gold md5:** `aee1382160adba01ab10ec5cc9c37323`
- **Delivered md5:** `6e97e4c31d507d3921aa8b702adefae6`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/migkas_slices_degeneracy_clean.png`
- **Notes:** Letter-era one-off plot script NOT FOUND. Exact delivery = copy of migkas_slices_degeneracy_clean.png. coe-figures/src/fig03_migkas_slices.py is a recompute that does NOT pixel-match (prior figs_regenerated_from_code Fig3 was WRONG vs PDF).

- **Letter-era one-off:** **missing** (or does not emit gold size).
- **(a) Verified asset path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/migkas_slices_degeneracy_clean.png`
- **(b) Best available current script:** `/workspace/coe-figures/src/fig03_migkas_slices.py` — full source in §4/§5; **does not guarantee pixel match**.

**Shell commands:**

```bash
# Exact match: COPY only (no matching generator found)
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/migkas_slices_degeneracy_clean.png \
   /workspace/figs_exact_from_code/delivery/Fig03_migkas_slices.png
# Optional non-matching recompute (do NOT expect pixel match):
cd /workspace/coe-figures/src && ../.venv/bin/python fig03_migkas_slices.py
```

**Verify:**

```bash
python3 - <<'PY'
from PIL import Image; import numpy as np
A=np.asarray(Image.open('/workspace/figs_exact_from_code/gold_extracts/gold-004.png').convert('RGBA'),np.int16)
B=np.asarray(Image.open('/workspace/figs_exact_from_code/delivery/Fig03_migkas_slices.png').convert('RGBA'),np.int16)
assert A.shape==B.shape, (A.shape,B.shape)
d=np.abs(A-B); print('max',int(d.max()),'mean',float(d.mean())); assert d.max()==0
PY
```

---

### Fig04 — `Fig04_amp_full.png`

- **Label:** Fig4 Full amp + residual
- **Gold extract:** `gold-006.png`
- **Expected dims:** 1800×1385
- **Method (delivery):** `copy:exact_on_disk_asset`
- **Pixel=gold:** True · **Bitwise=gold:** False
- **Gold md5:** `771d0afc75dd3aa8a128de519df58216`
- **Delivered md5:** `8ecabb7ddbd69a9c605cd45c459b9c6c`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/money_shot_amp_full.png`
- **Notes:** One-off script NOT FOUND. Use money_shot_amp_full.png (NOT money_shot_amp_full_SINE_BAD.png). coe-figures fig04_amp_full.py recomputes a different layout.

- **Letter-era one-off:** **missing** (or does not emit gold size).
- **(a) Verified asset path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/money_shot_amp_full.png`
- **(b) Best available current script:** `/workspace/coe-figures/src/fig04_amp_full.py` — full source in §4/§5; **does not guarantee pixel match**.

**Shell commands:**

```bash
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/money_shot_amp_full.png \
   /workspace/figs_exact_from_code/delivery/Fig04_amp_full.png
# Avoid: money_shot_amp_full_SINE_BAD.png
```

**Verify:**

```bash
python3 - <<'PY'
from PIL import Image; import numpy as np
A=np.asarray(Image.open('/workspace/figs_exact_from_code/gold_extracts/gold-006.png').convert('RGBA'),np.int16)
B=np.asarray(Image.open('/workspace/figs_exact_from_code/delivery/Fig04_amp_full.png').convert('RGBA'),np.int16)
assert A.shape==B.shape, (A.shape,B.shape)
d=np.abs(A-B); print('max',int(d.max()),'mean',float(d.mean())); assert d.max()==0
PY
```

---

### Fig05 — `Fig05_sky_orientation_four_row.jpg`

- **Label:** Fig5 Four-row sky orientation
- **Gold extract:** `gold-007.jpg`
- **Expected dims:** 1370×1420
- **Method (delivery):** `copy:bitwise_match_asset`
- **Pixel=gold:** True · **Bitwise=gold:** True
- **Gold md5:** `ed41b905e2563efb4df2c92b1e13e3b4`
- **Delivered md5:** `ed41b905e2563efb4df2c92b1e13e3b4`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/money_shot_four_row_only.jpg`
- **Notes:** Bitwise-identical to PDF JPEG. Letter-era four-row compose one-off not found as a single script. Row generators exist (make_coe_rings_*), and make_stack_clean.py is a related poster stacker (not the letter four_row). figs_regenerated recompute produced 1370x1441 (not 1370x1420) — close but not gold. Exact delivery = copy money_shot_four_row_only.jpg.

- **Letter-era one-off compose:** **missing** as a single preserved script.
- **(a) Verified asset path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/money_shot_four_row_only.jpg`
- **(b) Best available related sources:**
  - `/workspace/cosmology/figures/make_coe_rings_notilt.py`
  - `/workspace/cosmology/figures/make_coe_rings_migkas_pin_axis_matched.py`
  - `/workspace/cosmology/figures/make_coe_rings_migkas_pin_axis_matched_grad.py`
  - `/workspace/cosmology/figures/make_coe_rings_migkas_pin_axis_matched_grad_rot_fit.py`
  - `/workspace/cosmology/figures/make_stack_clean.py`
  - `/workspace/coe-figures/src/fig05_sky_orientation.py`

**Shell commands:**

```bash
# Exact match: COPY bitwise asset
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/money_shot_four_row_only.jpg \
   /workspace/figs_exact_from_code/delivery/Fig05_sky_orientation_four_row.jpg
# Best available row rebuild (outputs differ from gold size):
cd /workspace/cosmology/figures
python make_coe_rings_notilt.py
python make_coe_rings_migkas_pin_axis_matched.py
python make_coe_rings_migkas_pin_axis_matched_grad.py
python make_coe_rings_migkas_pin_axis_matched_grad_rot_fit.py
# Then crop/stack four loop_threepanel_*_web.jpg rows similar to make_stack_clean.py row_content()
# Regenerated compose historically wrote Fig05 at 1370x1441 != gold 1370x1420.
```

**Verify:**

```bash
md5sum /workspace/figs_exact_from_code/gold_extracts/gold-007.jpg /workspace/figs_exact_from_code/delivery/Fig05_sky_orientation_four_row.jpg
# both digests must match
```

---

### Fig06 — `Fig06_overlay_stack.png`

- **Label:** Fig6 Migkas Fig.8 overlay stack
- **Gold extract:** `gold-008.png`
- **Expected dims:** 2164×3222
- **Method (delivery):** `copy:bitwise_match_asset`
- **Pixel=gold:** True · **Bitwise=gold:** True
- **Gold md5:** `fd9f1c1c95194a3fb7864de08698704f`
- **Delivered md5:** `fd9f1c1c95194a3fb7864de08698704f`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/att14_overlay_stack.png`
- **Notes:** Bitwise-identical to PDF. Exact stack asset: att14_overlay_stack.png (2164x3222). Panel generator make_migkas_fig8_overlay_matched.py writes fig8_overlay_matched{,_grad,_grad_rot81}.png at 1975x1515 (with SHOW_LEGEND=True) into spreadsheet_model/. A naive vertical paste of those three panels yields ~1975x4689 (figs_regenerated), NOT the letter 2164x3222 crop/relabel stack. Letter-era stack compose one-liner not preserved as a standalone .py; use the verified stack PNG for exact match. Related: figs/fig8_overlays_three_stack.png is same dims but NOT bitwise identical to att14 (bak_before_restack).

- **Letter-era one-off:** **missing** (or does not emit gold size).
- **(a) Verified asset path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/att14_overlay_stack.png`
- **(b) Best available current script:** `/workspace/cosmology/figures/make_migkas_fig8_overlay_matched.py` — full source in §4/§5; **does not guarantee pixel match**.

**Shell commands:**

```bash
# Exact match: COPY
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/att14_overlay_stack.png \
   /workspace/figs_exact_from_code/delivery/Fig06_overlay_stack.png
# Rebuild panels (not the letter stack crop):
cd /workspace/cosmology/figures
# Temporarily set SHOW_LEGEND=True in make_migkas_fig8_overlay_matched.py for full legend panels,
# run, then restore SHOW_LEGEND=False (script default is False for letter single-panel use).
python make_migkas_fig8_overlay_matched.py
# Naive stack (NON-exact vs letter gold 2164x3222):
python3 -c "from PIL import Image; from pathlib import Path; P=Path('/workspace/cosmology/papers/spreadsheet_model'); parts=[Image.open(P/n).convert('RGB') for n in ('fig8_overlay_matched.png','fig8_overlay_matched_grad.png','fig8_overlay_matched_grad_rot81.png')]; W=max(p.width for p in parts); H=sum(p.height for p in parts); im=Image.new('RGB',(W,H),'white'); y=0
for p in parts:
 im.paste(p,((W-p.width)//2,y)); y+=p.height
im.save('/tmp/Fig06_naive_stack.png'); print(im.size)"
```

**Verify:**

```bash
md5sum /workspace/figs_exact_from_code/gold_extracts/gold-008.png /workspace/figs_exact_from_code/delivery/Fig06_overlay_stack.png
# both digests must match
```

---

### Fig07 — `Fig07_storyboard.png`

- **Label:** Fig7 Donut-universe storyboard A-D
- **Gold extract:** `gold-009.png`
- **Expected dims:** 1045×1795
- **Method (delivery):** `copy:exact_on_disk_asset`
- **Pixel=gold:** True · **Bitwise=gold:** False
- **Gold md5:** `05a980ff995d5a43f75de00c14fc750e`
- **Delivered md5:** `f2b0d04b62f80b9cb32c41a1e0d397ac`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/att14_storyboard.png`
- **Notes:** att14_storyboard.png pixel-identical. Current make_donut_universe_storyboard.py emits different size/layout (A-E, narrower) — not used for exact delivery.

- **Letter-era one-off:** **missing** (or does not emit gold size).
- **(a) Verified asset path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/att14_storyboard.png`
- **(b) Best available current script:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle7_recombination/make_donut_universe_storyboard.py` — full source in §4/§5; **does not guarantee pixel match**.

**Shell commands:**

```bash
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/att14_storyboard.png \
   /workspace/figs_exact_from_code/delivery/Fig07_storyboard.png
# Current script (different output):
cd /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle7_recombination
python make_donut_universe_storyboard.py
```

**Verify:**

```bash
python3 - <<'PY'
from PIL import Image; import numpy as np
A=np.asarray(Image.open('/workspace/figs_exact_from_code/gold_extracts/gold-009.png').convert('RGBA'),np.int16)
B=np.asarray(Image.open('/workspace/figs_exact_from_code/delivery/Fig07_storyboard.png').convert('RGBA'),np.int16)
assert A.shape==B.shape, (A.shape,B.shape)
d=np.abs(A-B); print('max',int(d.max()),'mean',float(d.mean())); assert d.max()==0
PY
```

---

### Fig08 — `Fig08_two_emitters.png`

- **Label:** Fig8 Two-emitter baseline
- **Gold extract:** `gold-011.png`
- **Expected dims:** 1181×1908
- **Method (delivery):** `copy:exact_on_disk_asset`
- **Pixel=gold:** True · **Bitwise=gold:** False
- **Gold md5:** `6953c5ae7902824f2fa24c854884cef2`
- **Delivered md5:** `41ddd33dbefce0695c7e225c72f26b58`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/att14_two_emitters.png`
- **Notes:** att14_two_emitters.png pixel-identical. Current make_two_emitters_comoving.py output differs — not used for exact delivery.

- **Letter-era one-off:** **missing** (or does not emit gold size).
- **(a) Verified asset path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/att14_two_emitters.png`
- **(b) Best available current script:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle6_two_emitters/make_two_emitters_comoving.py` — full source in §4/§5; **does not guarantee pixel match**.

**Shell commands:**

```bash
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/att14_two_emitters.png \
   /workspace/figs_exact_from_code/delivery/Fig08_two_emitters.png
cd /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle6_two_emitters
python make_two_emitters_comoving.py  # different size vs gold
```

**Verify:**

```bash
python3 - <<'PY'
from PIL import Image; import numpy as np
A=np.asarray(Image.open('/workspace/figs_exact_from_code/gold_extracts/gold-011.png').convert('RGBA'),np.int16)
B=np.asarray(Image.open('/workspace/figs_exact_from_code/delivery/Fig08_two_emitters.png').convert('RGBA'),np.int16)
assert A.shape==B.shape, (A.shape,B.shape)
d=np.abs(A-B); print('max',int(d.max()),'mean',float(d.mean())); assert d.max()==0
PY
```

---

### Fig09 — `Fig09_kill.jpg`

- **Label:** Fig9 Kill function vs plate/shell/shear
- **Gold extract:** `gold-013.jpg`
- **Expected dims:** 1206×1280
- **Method (delivery):** `copy:bitwise_match_asset`
- **Pixel=gold:** True · **Bitwise=gold:** True
- **Gold md5:** `15818532df3ff7b9bde79a200e779402`
- **Delivered md5:** `15818532df3ff7b9bde79a200e779402`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/kill_vs_cycle6_shapes.jpg`
- **Notes:** Letter JPEG one-off script NOT FOUND. kill_vs_cycle6_shapes.jpg bitwise=gold. Related cycle4 kill plots (plot_kill_term*.py) and coe-figures fig09_kill.py exist but layout differs.

- **Letter-era one-off:** **missing** (or does not emit gold size).
- **(a) Verified asset path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/kill_vs_cycle6_shapes.jpg`
- **(b) Best available current script:** `/workspace/coe-figures/src/fig09_kill.py` — full source in §4/§5; **does not guarantee pixel match**.

**Shell commands:**

```bash
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/kill_vs_cycle6_shapes.jpg \
   /workspace/figs_exact_from_code/delivery/Fig09_kill.jpg
# Related (non-matching): cycle4/dtt/plot_kill_term_0_360.py ; coe-figures/src/fig09_kill.py
```

**Verify:**

```bash
md5sum /workspace/figs_exact_from_code/gold_extracts/gold-013.jpg /workspace/figs_exact_from_code/delivery/Fig09_kill.jpg
# both digests must match
```

---

### Fig10 — `Fig10_cancel_partner.png`

- **Label:** Fig10 Cancel partner residuals
- **Gold extract:** `gold-014.png`
- **Expected dims:** 2084×2212
- **Method (delivery):** `copy:exact_on_disk_asset`
- **Pixel=gold:** True · **Bitwise=gold:** False
- **Gold md5:** `7f6701fcf63d81511470aea0c3626dc0`
- **Delivered md5:** `e8de308591d54d8c6de1b5481729b658`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/cancel_partner_residual.png`
- **Alt identical asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle6_two_emitters/simple_cmb_two_quad_cancel.png`
- **Notes:** cancel_partner_residual.png pixel-identical (= simple_cmb_two_quad_cancel.png). Letter one-off script not found; coe-figures fig10 is a recompute.

- **Letter-era one-off:** **missing** (or does not emit gold size).
- **(a) Verified asset path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/cancel_partner_residual.png`
- **(b) Best available current script:** `/workspace/coe-figures/src/fig10_cancel_partner.py` — full source in §4/§5; **does not guarantee pixel match**.

**Shell commands:**

```bash
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/cancel_partner_residual.png \
   /workspace/figs_exact_from_code/delivery/Fig10_cancel_partner.png
```

**Verify:**

```bash
python3 - <<'PY'
from PIL import Image; import numpy as np
A=np.asarray(Image.open('/workspace/figs_exact_from_code/gold_extracts/gold-014.png').convert('RGBA'),np.int16)
B=np.asarray(Image.open('/workspace/figs_exact_from_code/delivery/Fig10_cancel_partner.png').convert('RGBA'),np.int16)
assert A.shape==B.shape, (A.shape,B.shape)
d=np.abs(A-B); print('max',int(d.max()),'mean',float(d.mean())); assert d.max()==0
PY
```

---

### FigA1 — `FigA1_ltb_shell.png`

- **Label:** FigA1 LTB shell evolution
- **Gold extract:** `gold-016.png`
- **Expected dims:** 1444×969
- **Method (delivery):** `copy:exact_on_disk_asset`
- **Pixel=gold:** True · **Bitwise=gold:** False
- **Gold md5:** `717250da18b391623eb4618f233969fe`
- **Delivered md5:** `09462730cede0bb8dbb0052aef3eda1a`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/ltb_a_vs_t_all.png`
- **Notes:** ltb_a_vs_t_all.png pixel-identical. Letter one-off not found; coe-figures figA1_ltb_shell.py is best available recompute (may differ).

- **Letter-era one-off:** **missing** (or does not emit gold size).
- **(a) Verified asset path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/ltb_a_vs_t_all.png`
- **(b) Best available current script:** `/workspace/coe-figures/src/figA1_ltb_shell.py` — full source in §4/§5; **does not guarantee pixel match**.

**Shell commands:**

```bash
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/ltb_a_vs_t_all.png \
   /workspace/figs_exact_from_code/delivery/FigA1_ltb_shell.png
```

**Verify:**

```bash
python3 - <<'PY'
from PIL import Image; import numpy as np
A=np.asarray(Image.open('/workspace/figs_exact_from_code/gold_extracts/gold-016.png').convert('RGBA'),np.int16)
B=np.asarray(Image.open('/workspace/figs_exact_from_code/delivery/FigA1_ltb_shell.png').convert('RGBA'),np.int16)
assert A.shape==B.shape, (A.shape,B.shape)
d=np.abs(A-B); print('max',int(d.max()),'mean',float(d.mean())); assert d.max()==0
PY
```

---

### FigA2 — `FigA2_coe_sun.jpg`

- **Label:** FigA2 CoE/Sun timing 2026
- **Gold extract:** `gold-018.jpg`
- **Expected dims:** 1280×636
- **Method (delivery):** `copy:bitwise_match_asset`
- **Pixel=gold:** True · **Bitwise=gold:** True
- **Gold md5:** `9cb249b71dbb3ea478495deb47837921`
- **Delivered md5:** `9cb249b71dbb3ea478495deb47837921`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/appendix_coe_sun_timing_2026.jpg`
- **Notes:** Bitwise JPEG match. Letter-era plot script not found; coe-figures figA2_coe_sun.py is best available recompute (PNG, may differ).

- **Letter-era one-off:** **missing** (or does not emit gold size).
- **(a) Verified asset path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/appendix_coe_sun_timing_2026.jpg`
- **(b) Best available current script:** `/workspace/coe-figures/src/figA2_coe_sun.py` — full source in §4/§5; **does not guarantee pixel match**.

**Shell commands:**

```bash
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/appendix_coe_sun_timing_2026.jpg \
   /workspace/figs_exact_from_code/delivery/FigA2_coe_sun.jpg
```

**Verify:**

```bash
md5sum /workspace/figs_exact_from_code/gold_extracts/gold-018.jpg /workspace/figs_exact_from_code/delivery/FigA2_coe_sun.jpg
# both digests must match
```

---

### FigA2b — `FigA2b_sun_coe_crossings.jpg`

- **Label:** TableA1 Sun/CoE crossings 2027
- **Gold extract:** `gold-019.jpg`
- **Expected dims:** 1239×570
- **Method (delivery):** `copy:bitwise_match_asset`
- **Pixel=gold:** True · **Bitwise=gold:** True
- **Gold md5:** `6cf8f10855878c559c54e97ea06cb029`
- **Delivered md5:** `6cf8f10855878c559c54e97ea06cb029`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/appendix_sun_coe_crossings_2027.jpg`
- **Notes:** Bitwise JPEG match. Letter-era script not found; coe-figures figA_table_crossings.py is best available recompute.

- **Letter-era one-off:** **missing** (or does not emit gold size).
- **(a) Verified asset path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/appendix_sun_coe_crossings_2027.jpg`
- **(b) Best available current script:** `/workspace/coe-figures/src/figA_table_crossings.py` — full source in §4/§5; **does not guarantee pixel match**.

**Shell commands:**

```bash
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/appendix_sun_coe_crossings_2027.jpg \
   /workspace/figs_exact_from_code/delivery/FigA2b_sun_coe_crossings.jpg
```

**Verify:**

```bash
md5sum /workspace/figs_exact_from_code/gold_extracts/gold-019.jpg /workspace/figs_exact_from_code/delivery/FigA2b_sun_coe_crossings.jpg
# both digests must match
```

---

### FigA3 — `FigA3_coe_sky.png`

- **Label:** FigA3 Annual CoE sky path from Plano
- **Gold extract:** `gold-020.png`
- **Expected dims:** 2415×1432
- **Method (delivery):** `rerun:_make_coe_sky_path_dates.py`
- **Pixel=gold:** True · **Bitwise=gold:** False
- **Gold md5:** `1c0042622ee084daa22f2dfaf7f4a086`
- **Delivered md5:** `6d68e94d82d64d87d409f67620259ce4`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/coe_sky_path_dates_texas.png`
- **Notes:** Re-ran _make_coe_sky_path_dates.py; pixel-identical to gold. Needs astropy; iers auto_download disabled in script.

- **Generating script (pixel-match on rerun):** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/_make_coe_sky_path_dates.py` — full source in §3.

**Shell commands:**

```bash
cd /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs
python _make_coe_sky_path_dates.py
cp coe_sky_path_dates_texas.png /workspace/figs_exact_from_code/delivery/FigA3_coe_sky.png
```

**Verify:**

```bash
python3 - <<'PY'
from PIL import Image; import numpy as np
A=np.asarray(Image.open('/workspace/figs_exact_from_code/gold_extracts/gold-020.png').convert('RGBA'),np.int16)
B=np.asarray(Image.open('/workspace/figs_exact_from_code/delivery/FigA3_coe_sky.png').convert('RGBA'),np.int16)
assert A.shape==B.shape, (A.shape,B.shape)
d=np.abs(A-B); print('max',int(d.max()),'mean',float(d.mean())); assert d.max()==0
PY
```

---

### FigA_comic — `FigA_comic_HANDDRAWN_from_PDF.jpg`

- **Label:** p.14 hand-drawn comic
- **Gold extract:** `gold-022.jpg`
- **Expected dims:** 2080×944
- **Method (delivery):** `extract/copy:handdrawn_not_from_code`
- **Pixel=gold:** True · **Bitwise=gold:** True
- **Gold md5:** `262874d6412bd4b1d38b7905e237e72b`
- **Delivered md5:** `262874d6412bd4b1d38b7905e237e72b`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/appendix_coe_comic.jpg`
- **Notes:** HAND-DRAWN — no generating code. Extract from letter PDF (pdfimages) or copy appendix_coe_comic.jpg / gold-022.jpg.

- **Code:** none (hand-drawn / PDF extract only).

**Shell commands:**

```bash
# From gold extract (pdfimages):
cp /workspace/figs_exact_from_code/gold_extracts/gold-022.jpg \
   /workspace/figs_exact_from_code/delivery/FigA_comic_HANDDRAWN_from_PDF.jpg
# Or: pdfimages -all <letter.pdf> gold  -> use gold-022.jpg
```

**Verify:**

```bash
md5sum /workspace/figs_exact_from_code/gold_extracts/gold-022.jpg /workspace/figs_exact_from_code/delivery/FigA_comic_HANDDRAWN_from_PDF.jpg
# both digests must match
```

---

### eclipse — `eclipse_inner_shell.png`

- **Label:** Inner-shell eclipse geometry (supplemental; not a letter-page embed)
- **Gold extract:** `(n/a — not a letter embed)`
- **Expected dims:** 883×945
- **Method (delivery):** `rerun:attachment make_inner_shell_eclipse_sketch.py`
- **Pixel=gold:** None · **Bitwise=gold:** None
- **Delivered md5:** `f76d8526ed668a5df5855f9cb66314cc`
- **Verified matching asset:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/inner_shell_eclipse_sketch.png`
- **Notes:** Ran refs_from_chris/make_inner_shell_eclipse_sketch.py (also in cycle6_two_emitters/). Matches letter_brief/figs/inner_shell_eclipse_sketch.png.

- **Generating script (pixel-match on rerun):** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/refs_from_chris/make_inner_shell_eclipse_sketch.py` — full source in §3.

**Shell commands:**

```bash
cd /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/refs_from_chris
python make_inner_shell_eclipse_sketch.py
# writes inner_shell_eclipse_sketch.png next to the script (or edit OUT path)
cp /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/inner_shell_eclipse_sketch.png \
   /workspace/figs_exact_from_code/delivery/eclipse_inner_shell.png
```

**Verify:** compare to `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/inner_shell_eclipse_sketch.png` (md5 / max abs).

---

## 3. Full sources — pixel-matching reruns

These scripts were used for exact delivery (rerun → pixel match).

### 3.1 `universe_thin.py` (Fig01 + Fig02)

**Path:** `/workspace/cosmology/papers/spreadsheet_model/bubble/universe_thin.py`  
**Cwd:** that `bubble/` directory (imports `bubble_core`, `analysis_curve`; needs thin_*.json / universe_*.json next to it).  
**Deps:** matplotlib, numpy, scipy, pillow; local bubble package.

```python
# source: /workspace/cosmology/papers/spreadsheet_model/bubble/universe_thin.py
import json, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
from scipy.optimize import brentq
from PIL import Image
from bubble_core import AFfour, ML, MA, MLO, MHI, TH1, H0, GSI, RKM, GPC, IMIG, r2
import analysis_curve as AC
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
c = 2.99792458e8; Msun = 1.989e30; Gm = GPC * 1e3
def phys(k):
    d = json.load(open("thin_%s.json" % k.replace("|", "_")))
    Dm = d["D_Gpc"] * Gm; rhoc = 3 * H0**2 / (8 * np.pi * GSI); rhos = d["Os"] * rhoc; Mcrit = 4 * np.pi / 3 * rhoc * Dm**3
    Mc = d["mcc"] * Mcrit; rin = d["rin"] * Dm; rout = d["rout"] * Dm
    Msh = lambda r: 4 * np.pi / 3 * rhos * (np.clip(r, rin, rout)**3 - rin**3)
    def Phi(r):
        if r < rin: ps = -GSI * 2 * np.pi * rhos * (rout**2 - rin**2)
        elif r < rout: ps = -GSI * (4 * np.pi / 3 * rhos * (r**3 - rin**3) / r + 2 * np.pi * rhos * (rout**2 - r**2))
        else: ps = -GSI * Msh(rout) / r
        return -GSI * Mc / r + ps
    v2 = d["vesc2"] * 2 * abs(Phi(Dm)); E = v2 / 2 + Phi(Dm)
    rt = brentq(lambda r: Phi(r) - E, Dm, 1e6 * Dm) if E < 0 else np.inf
    rs = lambda M: 2 * GSI * M / c**2 / Gm
    Mt = Mc + Msh(min(rt, 1e300))
    return dict(d=d, D=d["D_Gpc"], R=RKM / GPC, DR=d["DR"], Os=d["Os"], rhos=rhos, Mc=Mc, Mcrit=Mcrit, mcc=d["mcc"], Msh=Msh(Dm), Mtot=Mc + Msh(Dm),
                Mhole=4 * np.pi / 3 * rhos * rin**3, Rin=rin / Gm, Rout=rout / Gm, rs_c=rs(Mc), rs_D=rs(Mc + Msh(Dm)), rt=rt / Gm, rs_t=rs(Mt),
                vE=np.sqrt(v2) / c, vesc2=d["vesc2"], age=d["age"], t_turn=d["t_turn"], amp=d["amp"], phi=d["phi"], l0=d["l0"], R2=d["R2_mig"], R2F=d["R2_four"],
                Lc=d["Lc"] * d["D_Gpc"], Dkm=d["DR"] * RKM, Lckm=d["Lc"] * d["DR"] * RKM)
P1, P2 = phys("b|free|hug"), phys("b|0.05|hug")
json.dump({k: {kk: vv for kk, vv in p.items() if kk != "d"} for k, p in (("free", P1), ("baryons", P2))}, open("universe_thin.json", "w"), indent=1, default=float)
for p in (P1, P2): print({k: (round(v, 4) if isinstance(v, float) and abs(v) < 1e6 else v) for k, v in p.items() if k != "d"})
# ---- relativistic (exact LTB) thick ordinary-matter shell, z = 0.3 case (gr_ltb/refit_mass) ----
def rel():
    o = json.load(open("universe_rel_mcfree.json"))   # all free: D, Mc, R_in, R_out, l0, Amp, phi (old Mc-constrained fit: best_rows_thick.json)
    D = o["D"]; Dm = D * Gm; rhoc = 3 * H0**2 / (8 * np.pi * GSI); rhos = o["Os"] * rhoc
    rs = lambda M: 2 * GSI * M / c**2 / Gm; Rv = RKM / GPC; Rin, Rout = o["R_in"], o["R_out"]
    cur = np.array(o["A"])
    return dict(D=D, R=Rv, DR=D / Rv, Os=o["Os"], rhos=rhos, Mc=o["Mc_kg"], Mcrit=4 * np.pi / 3 * rhoc * Dm**3, mcc=o["Mc_kg"] / (4 * np.pi / 3 * rhoc * Dm**3),
                Msh=o["M_D_kg"] - o["Mc_kg"], Mtot=o["M_D_kg"], Mshtot=4 * np.pi / 3 * rhos * ((Rout * Gm)**3 - (Rin * Gm)**3), Mhole=4 * np.pi / 3 * rhos * (Rin * Gm)**3,
                Mvis=4 * np.pi / 3 * rhos * (Rv * Gm)**3, Rin=Rin, Rout=Rout, rs_c=rs(o["Mc_kg"]), rs_D=rs(o["M_D_kg"]), vE=H0 * Dm / c, twoE=o["twoE"],
                age=o["age_Gyr"], amp=o["amp"], phi=o["phi"], l0=o["l0_s"], R2=o["R2_s"], R2ns=o["R2_ns"], R2F=r2(AFfour(TH1), cur), H_r=o["H_r"], H_perp=o["H_perp"],
                z=o["z"], Lc=0.005 * D, Dkm=D * GPC, Lckm=0.005 * D * GPC, cur=cur)
P3 = rel()
def newt_omb():   # Newtonian, Omega_s = 0.049 forced (newt_omb.py)
    o = json.load(open("universe_newt_omb.json"))["Rs14.3|sine"]; D = o["D"]; Dm = D * Gm; rhoc = 3 * H0**2 / (8 * np.pi * GSI); rhos = o["Os"] * rhoc
    rs = lambda M: 2 * GSI * M / c**2 / Gm; Rv = RKM / GPC; cur = np.array(o["A"])
    return dict(D=D, R=Rv, DR=D / Rv, Os=o["Os"], rhos=rhos, Mc=o["Mc_kg"], Mcrit=o["Mc_kg"] / o["mcc"], mcc=o["mcc"], Msh=o["shell_inside_D_kg"], Mtot=o["M_inside_D_kg"],
                Mshtot=o["shell_total_kg"], Mhole=4 * np.pi / 3 * rhos * (o["R_in"] * Gm)**3, Mvis=4 * np.pi / 3 * rhos * (Rv * Gm)**3, Rin=o["R_in"], Rout=o["R_out"],
                rs_c=rs(o["Mc_kg"]), rs_D=rs(o["M_inside_D_kg"]), vE=o["VO_over_c"], vesc2=o["vesc2_full_potential"], age=o["age_Gyr"], amp=o["amp"], phi=o["phi"], l0=o["l0"],
                R2=o["R2"], R2geo=json.load(open("universe_newt_omb.json"))["Rs14.3|geom"]["R2"], R2F=r2(AFfour(TH1), cur), Lc=0.005 * D, Dkm=D * GPC, Lckm=0.005 * D * GPC, cur=cur)
P4 = newt_omb()
for q in (P1, P2): q["Mshtot"] = 4 * np.pi / 3 * q["rhos"] * ((q["Rout"] * Gm)**3 - (q["Rin"] * Gm)**3); q["Mvis"] = 4 * np.pi / 3 * q["rhos"] * (q["R"] * Gm)**3
json.dump({k: {kk: vv for kk, vv in q.items() if kk not in ("d", "cur")} for k, q in (("free", P1), ("baryons", P2), ("rel", P3), ("newt_omb", P4))}, open("universe_thin.json", "w"), indent=1, default=float)
print("REL", {k: (round(v, 4) if isinstance(v, float) and abs(v) < 1e6 else v) for k, v in P3.items() if k != "cur"})
p = P3; Dg, Rg, l0 = p["D"], p["R"], p["l0"]
lpk = (90 - np.degrees(p["phi"]) + l0) % 360
Ev = np.array([Dg, 0.0])
def sky(l, r): t = np.radians(l - l0); return Ev + r * np.array([np.cos(t), np.sin(t)])
# ---------------- Figure 1: one figure, three axes, unified style ----------------
from matplotlib.gridspec import GridSpec
from itertools import combinations
FS_A, FS_L, FS_T = 7.5, 8.5, 9            # annotations, axis labels/ticks, titles
plt.rcParams.update({"mathtext.fontset": "cm", "axes.linewidth": 0.8, "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.labelsize": FS_L, "axes.titlesize": FS_T,
                     "xtick.major.width": 0.8, "ytick.major.width": 0.8})
CRED, CBLU, CDC, CGRN, CPUR, CSH = "#d62728", "#0055d4", "#1f77b4", "#2ca02c", "#9467bd", "0.82"
VEL = dict(arrowstyle="-|>", color=CBLU, lw=1.3, mutation_scale=9, shrinkA=0, shrinkB=0)
POS = lambda ls="-": dict(arrowstyle="-|>", color="0.2", lw=0.9, ls=ls, mutation_scale=8, shrinkA=0, shrinkB=0)
exA = (-0.6, 5.6, -2.6, 3.4); exB = (-7, 82.7, -34.5, 34.5); E2 = 230
ka = (exA[1] - exA[0]) / (exA[3] - exA[2]); kb = (exB[1] - exB[0]) / (exB[3] - exB[2])
FW = 7.4; L0, G1, G2, R0 = 0.04, 0.58, 0.54, 0.10
H = (FW - L0 - G1 - G2 - R0) / (ka + kb + 1); BOT, TOP = 0.42, 0.24; FH = H + BOT + TOP
fig = plt.figure(figsize=(FW, FH))
gs = GridSpec(1, 5, figure=fig, width_ratios=[ka, G1 / H, kb, G2 / H, 1], wspace=0, left=L0 / FW, right=1 - R0 / FW, bottom=BOT / FH, top=1 - TOP / FH)
vx, ax, bx = fig.add_subplot(gs[0]), fig.add_subplot(gs[2]), fig.add_subplot(gs[4])
def T(axx, *a, **k): k.setdefault("fontsize", FS_A); return axx.text(*a, **k)
# (a) Hubble velocity vectors (schematic)
C0_, O_, P_ = np.array([0, 0.9]), np.array([2.6, 0.9]), np.array([3.6, 2.3]); h_ = 0.4
VO, VP = h_ * (O_ - C0_), h_ * (P_ - C0_); VH = VP - VO
vx.annotate("", xy=O_, xytext=C0_, arrowprops=POS()); vx.annotate("", xy=P_, xytext=C0_, arrowprops=POS("--")); vx.annotate("", xy=P_, xytext=O_, arrowprops=POS(":"))
for base, vec in ((O_, VO), (P_, VP), (P_, VH)): vx.annotate("", xy=base + vec, xytext=base, arrowprops=VEL)
vx.plot(*C0_, "o", color=CRED, ms=5, zorder=6); vx.plot(*O_, "o", color=CBLU, mec="k", mew=0.6, ms=5, zorder=6); vx.plot(*P_, "o", color="0.3", ms=4.5, zorder=6)
T(vx, *(C0_ + [0, -0.18]), "CoE", color=CRED, ha="center", va="top")
T(vx, *(O_ + [0.02, -0.18]), "us ($O$)", ha="left", va="top")
T(vx, *(P_ + [0.1, -0.1]), "galaxy", ha="left", va="top")
T(vx, *((C0_ + O_) / 2 + [0, -0.12]), "$\\vec D$", ha="center", va="top")
T(vx, *((C0_ + P_) / 2 + [-0.12, 0.12]), "$\\vec D+\\vec r$", ha="right", va="bottom")
T(vx, *((O_ + P_) / 2 + [0.1, -0.05]), "$\\vec r$", ha="left", va="center")
T(vx, *(O_ + VO + [0.08, 0]), "$\\vec V_O$", color=CBLU, ha="left", va="bottom")
T(vx, *(P_ + VP + [0.04, -0.06]), "$\\vec V_P$", color=CBLU, ha="left", va="top")
T(vx, *(P_ + VH + [0.07, -0.02]), "$\\vec V_H$", color=CBLU, ha="left", va="bottom")
A_ = np.array([1.1, -1.45]); s_ = 1.6
for base, vec in ((A_, s_ * VO), (A_, s_ * VP), (A_ + s_ * VO, s_ * VH)): vx.annotate("", xy=base + vec, xytext=base, arrowprops=VEL)
T(vx, *(A_ + s_ * VO / 2 + [0, -0.1]), "$\\vec V_O$", color=CBLU, ha="center", va="top")
T(vx, *(A_ + s_ * VP / 2 + [-0.1, 0.1]), "$\\vec V_P$", color=CBLU, ha="right", va="bottom")
T(vx, *(A_ + s_ * VO + s_ * VH / 2 + [0.1, 0]), "$\\vec V_H$", color=CBLU, ha="left", va="center")
T(vx, 2.5, -2.55, "$\\vec V_H=\\vec V_P-\\vec V_O=H\\vec r$", ha="center", va="bottom")
vx.set_xlim(exA[0], exA[1]); vx.set_ylim(exA[2], exA[3]); vx.set_aspect("equal"); vx.set_xticks([]); vx.set_yticks([])
vx.set_title("(a) Hubble velocity vectors")
# (b) Hubble bubble scale
ax.add_patch(Circle((0, 0), p["Rout"], fc=CSH, ec="0.35", lw=0.8, zorder=0)); ax.add_patch(Circle((0, 0), p["Rin"], fc="white", ec="0.35", lw=0.8, zorder=0.5))
ax.add_patch(Circle((0, 0), Dg, fill=False, color=CDC, lw=1.0, ls="--", zorder=2))
ax.add_patch(Circle(Ev, Rg, fc=CGRN, alpha=0.18, ec=CGRN, lw=1.2, zorder=2))
ax.add_patch(Circle((0, 0), p["rs_c"], fill=False, color=CRED, lw=1.4, zorder=3))
ax.plot(0, 0, "o", color=CRED, ms=5, zorder=6); T(ax, -5.5, 2.4, "CoE core", color=CRED, ha="left", va="bottom")
ax.plot(*Ev, "o", color=CBLU, mec="k", mew=0.6, ms=5, zorder=7); T(ax, p["Rin"] - 1.4, -1.0, "Earth", ha="right", va="top", zorder=7)
T(ax, 5.5, 11.7, "empty\ninterior", color="0.4", ha="center", va="center", linespacing=1.1)
T(ax, 0.5 * (Dg + p["Rout"]) * np.cos(np.radians(42)), 0.5 * (Dg + p["Rout"]) * np.sin(np.radians(42)), "shell", color="0.25", ha="center", va="center", rotation=-48, rotation_mode="anchor")
tq = np.radians(125); q = Ev + Rg * np.array([np.cos(tq), np.sin(tq)])
ax.annotate("", xy=q, xytext=Ev, arrowprops=dict(arrowstyle="-|>", color=CGRN, lw=1.0, mutation_scale=8, shrinkA=0, shrinkB=0)); T(ax, (Ev[0] + q[0]) / 2 - 1.0, (Ev[1] + q[1]) / 2 - 1.0, "$R$", color=CGRN, ha="right", va="center")
td = np.radians(-38); qd = Dg * np.array([np.cos(td), np.sin(td)])
ax.annotate("", xy=qd, xytext=(0, 0), arrowprops=dict(arrowstyle="<|-|>", color=CDC, lw=0.9, mutation_scale=7, shrinkA=0, shrinkB=0)); T(ax, *(qd * 0.58 + [1.4, -1.4]), "$D$", color=CDC, ha="left", va="top")
T(ax, 43.5, 33, "core $r_s$ = %.0f Gpc\n(%.1f D), off scale" % (p["rs_c"], p["rs_c"] / Dg), color=CRED, ha="left", va="top", linespacing=1.1)
VX1 = 81; ax.annotate("", xy=(VX1, 0), xytext=Ev, arrowprops=VEL, zorder=8)
vlab = T(ax, VX1, -2.8, "$V_O\\!\\approx\\!%.1fc$" % p["vE"], color=CBLU, ha="right", va="top")
ax.set_xlim(exB[0], exB[1]); ax.set_ylim(exB[2], exB[3]); ax.set_aspect("equal")
ax.set_xlabel("$x$ [Gpc]"); ax.set_ylabel("$y$ [Gpc]"); ax.set_title("(b) Hubble bubble scale")
# (c) Cosmic shell scale
bx.add_patch(Circle((0, 0), p["Rout"], fc=CSH, ec="0.35", lw=0.6, zorder=1)); bx.add_patch(Circle((0, 0), p["Rin"], fc="white", ec="0.35", lw=0.5, zorder=1.5))
bx.add_patch(Circle((0, 0), p["rs_c"], fill=False, color=CRED, lw=1.4, zorder=3))
bx.add_patch(Circle(Ev, Rg, fc=CGRN, alpha=0.18, ec=CGRN, lw=1.2, zorder=2))   # vision sphere R, to scale (as in b)
T(bx, 0, 0.87 * p["rs_c"], "core $r_s$\n= %.0f Gpc" % p["rs_c"], color=CRED, ha="center", va="top", linespacing=1.1); T(bx, p["Rout"] + 6, 25, "shell", color="0.25", ha="left", va="center")
T(bx, 0, -0.996 * E2, "$r_s(<D)$ = %s Gpc\n(%.1f D), off scale" % ("{:,.0f}".format(p["rs_D"]), p["rs_D"] / Dg), color=CRED, ha="center", va="bottom", linespacing=1.0)
for ang in (45, 135, 225, 315):
    u = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))]); bx.annotate("", xy=u * 0.87 * p["rs_c"] * 0.75, xytext=u * p["Rout"] * 1.05, arrowprops=VEL, zorder=4)
bx.set_xlim(-E2, E2); bx.set_ylim(-E2, E2); bx.set_aspect("equal"); bx.set_xticks([-200, 0, 200]); bx.set_yticks([-200, 0, 200])
bx.set_xlabel("$x$ [Gpc]"); bx.set_title("(c) Cosmic shell scale")
# checks
fig.canvas.draw(); rr = fig.canvas.get_renderer()
for A_x in (vx, ax, bx):
    bbs = [(t.get_text()[:22], t.get_window_extent(rr)) for t in A_x.texts if t.get_text()]
    for (n1, b1), (n2, b2) in combinations(bbs, 2):
        if b1.overlaps(b2): print("OVERLAP", repr(n1), "|", repr(n2))
    W_ = A_x.get_window_extent(rr)
    for n1, b1 in bbs:
        if b1.x0 < W_.x0 or b1.x1 > W_.x1 or b1.y0 < W_.y0 or b1.y1 > W_.y1: print("OUTSIDE", repr(n1))
for A_x, circ, pad in ((ax, (p["Rin"], p["Rout"], p["rs_c"], Dg), 1.2), (bx, (p["Rin"], p["Rout"], p["rs_c"]), 17)):
    inv_ = A_x.transData.inverted()
    for t in A_x.texts:
        (x0, y0_), (x1, y1) = inv_.transform(t.get_window_extent(rr)); xs, ys = np.meshgrid(np.linspace(x0, x1, 60), np.linspace(y0_, y1, 20)); rr_ = np.hypot(xs, ys)
        for R_ in circ:
            if rr_.min() - pad < R_ < rr_.max() + pad: print("EDGE", repr(t.get_text()[:20]), round(R_, 1), np.round([x0, y0_, x1, y1], 1))
print("axes heights (in)", [round(a_.get_window_extent(rr).height / fig.dpi, 3) for a_ in (vx, ax, bx)], "H", round(H, 3), "tight W", round(fig.get_tightbbox(rr).width, 3))
fig.savefig("../universe_thin_drawing.png", dpi=200, bbox_inches="tight", pad_inches=0.02); plt.close(fig)
plt.rcParams.update({"mathtext.fontset": "dejavusans", "axes.linewidth": 0.8, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.labelsize": 9, "axes.titlesize": 12})
AH = float(os.environ.get("AH", "7.58"))
Af = AFfour(TH1); th2 = np.append(TH1, 360); Af2 = AFfour(th2)
fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, AH), sharex=True, gridspec_kw={"height_ratios": [2.0, 1], "hspace": 0.06})
a.errorbar(ML, MA, yerr=[MLO, MHI], fmt="o", ms=3, color="0.45", ecolor="0.75", lw=0.8, label="Migkas Fig. 6")
a.plot(th2, Af2, color="C0", lw=1.5, label="$A_F$")
for pp, col, ls, lab in [(P1, "C3", "-", "Newtonian thin shell + core, $\\Omega_s$=%.2f"), (P4, "C2", ":", "Newtonian shell + core, $\\Omega_s$=%.3f forced"), (P3, "C1", "--", "relativistic (LTB) shell + core, $\\Omega_s$=%.3f, $z$=0.3")]:
    A = AC.curve("ptshell", pp["d"]) if "d" in pp else pp["cur"]; assert abs(r2(MA, A[IMIG]) - pp["R2"]) < 2e-3, (r2(MA, A[IMIG]), pp["R2"])
    A2 = np.append(A, A[0])
    a.plot(th2, A2, color=col, ls=ls, lw=1.6 if ls != "--" else 1.3, label=(lab % pp["Os"]) + ", $\\ell_0$=%.1f° (fit $R^2$=%.3f; vs $A_F$ %.3f)" % (pp["l0"], pp["R2"], pp["R2F"]))
    b.plot(th2, 100 * (A2 - Af2) / Af2, color=col, ls=ls, lw=1.4 if ls == "-" else 1.2, label={"-": "Newtonian", ":": "Newt. $\\Omega_s$=0.049", "--": "relativistic"}[ls])
a.set_ylabel("$A/A_{\\rm all}$"); a.set_ylim(0.62, 1.22); a.legend(fontsize=7.2, loc="lower left", frameon=False)
b.axhline(0, color="k", lw=0.6); b.plot(ML, 100 * (MA - AFfour(ML)) / AFfour(ML), "o", ms=2.5, color="0.45", label="data")
b.set_ylabel("$100\\,(\\cdot-A_F)/A_F$ [%]", fontsize=8); b.set_xlabel("Galactic longitude $\\ell$ (deg)")
b.set_xlim(0, 360); b.set_xticks(range(0, 361, 45)); b.legend(fontsize=7, frameon=False, loc="upper right", ncol=4)
fig.savefig("../universe_thin_Achart.png", dpi=200, bbox_inches="tight")
for f in ["universe_thin_drawing", "universe_thin_Achart"]:
    im = Image.open("../" + f + ".png").convert("RGB"); im.thumbnail((1700, 1700)); im.save("../" + f + "_web.jpg", quality=88)
```

### 3.2 `_make_coe_sky_path_dates.py` (FigA3)

**Path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/_make_coe_sky_path_dates.py`

```python
# source: /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs/_make_coe_sky_path_dates.py
#!/usr/bin/env python3
"""Make a clean annual CoE sky-path chart for Plano, Texas."""
from __future__ import annotations

from pathlib import Path
import calendar

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.utils import iers

iers.conf.auto_download = False

i = Path(__file__).resolve().parent
OUT = i / "coe_sky_path_dates_texas.png"

# Inputs supplied in the request.
YEAR = 2026
LAT, LON, HEIGHT_M = 33.02, -96.70, 200.0
LMT_HOUR = 12.0
GAL_L, GAL_B = 180.0, -15.0

site = EarthLocation(lat=LAT * u.deg, lon=LON * u.deg, height=HEIGHT_M * u.m)
coe_gal = SkyCoord(l=GAL_L * u.deg, b=GAL_B * u.deg, frame="galactic")
coe = coe_gal.icrs

# Local mean solar time: UT = LMT - longitude/15 hours.
doys = np.arange(1, 366)
t0 = Time(f"{YEAR}-01-01 00:00:00", scale="utc")
times = t0 + (doys - 1 + LMT_HOUR / 24.0 - LON / 15.0 / 24.0) * u.day
altaz = coe.transform_to(AltAz(obstime=times, location=site))
az = np.asarray(altaz.az.deg)
alt = np.asarray(altaz.alt.deg)

# Month labels at the 15th, an uncluttered representative date for each month.
month_names = [calendar.month_abbr[m] for m in range(1, 13)]
month_doys = np.array([sum(calendar.monthrange(YEAR, m)[1] for m in range(1, mm)) + 15 for mm in range(1, 13)])
month_idx = month_doys - 1

# Build a daily colored path; omit the single segment that crosses the 0/360 seam.
points = np.column_stack([az, alt])
segments = np.stack([points[:-1], points[1:]], axis=1)
seam = np.abs(np.diff(az)) > 180.0
segments = segments[~seam]
seg_day = (doys[:-1] + doys[1:]) / 2.0
seg_day = seg_day[~seam]

fig, ax = plt.subplots(figsize=(12.4, 8.4), dpi=200)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Below-horizon region and horizon.
ax.axhspan(-45, 0, color="#f1f3f5", zorder=0)
ax.axhline(0, color="#30343b", lw=1.45, zorder=1)
ax.text(356.5, 1.7, "horizon", ha="right", va="bottom", fontsize=8.5, color="#50555c")

# Cardinal guides: restrained enough to keep the path dominant.
for x in (0, 90, 180, 270, 360):
    ax.axvline(x, color="#dfe3e7", lw=0.8, zorder=0)

norm = Normalize(vmin=1, vmax=365)
cmap = plt.get_cmap("turbo")
lc = LineCollection(segments, cmap=cmap, norm=norm, linewidths=2.5, alpha=0.96,
                    capstyle="round", joinstyle="round", zorder=3)
lc.set_array(seg_day)
ax.add_collection(lc)

# Daily sample dots make the once-per-day sampling explicit.
ax.scatter(az, alt, c=doys, cmap=cmap, norm=norm, s=11, linewidths=0,
           alpha=0.98, zorder=4)

# Month markers and labels. Offsets are chosen to keep all labels clear of the curve.
offsets = {
    "Jan": (-28, -18), "Feb": (-28, 10), "Mar": (-25, 10), "Apr": (-28, 9),
    "May": (-28, 10), "Jun": (10, 10), "Jul": (10, 9), "Aug": (11, 5),
    "Sep": (11, -1), "Oct": (10, -15), "Nov": (-35, -4), "Dec": (10, -7),
}
for name, idx in zip(month_names, month_idx):
    x, y = az[idx], alt[idx]
    ax.scatter([x], [y], s=68, facecolor="white", edgecolor="#20242a",
               linewidth=1.15, zorder=6)
    dx, dy = offsets[name]
    ax.annotate(
        f"{name} 15",
        (x, y), xytext=(dx, dy), textcoords="offset points",
        fontsize=9.2, fontweight="bold", color="#20242a", zorder=7,
        bbox=dict(boxstyle="round,pad=0.20", fc="white", ec="#b8bec6", lw=0.65, alpha=0.94),
        arrowprops=dict(arrowstyle="-", color="#747b84", lw=0.75, shrinkA=4, shrinkB=4),
    )

# Axes, grid, and cardinal labels.
ax.set_xlim(0, 360)
ax.set_ylim(-45, 85)
ax.set_xticks([0, 90, 180, 270, 360])
ax.set_xticklabels(["N\n0°", "E\n90°", "S\n180°", "W\n270°", "N\n360°"],
                   fontsize=10.5, fontweight="bold", color="#252a30")
ax.set_xticks(np.arange(45, 360, 45), minor=True)
ax.set_yticks(np.arange(-40, 81, 20))
ax.set_xlabel("Azimuth (clockwise from North)", fontsize=11, labelpad=13)
ax.set_ylabel("Altitude above horizon (°)", fontsize=11, labelpad=10)
ax.grid(which="major", color="#cdd2d8", lw=0.75, alpha=0.62)
ax.grid(which="minor", axis="x", color="#e5e8eb", lw=0.55, alpha=0.8)
ax.tick_params(axis="x", which="minor", length=3, color="#8d949c")
ax.tick_params(axis="both", which="major", labelsize=9.5)
for spine in ax.spines.values():
    spine.set_color("#858c94")
    spine.set_linewidth(0.8)

fig.suptitle("CoE annual sky path from Plano, Texas", fontsize=18, fontweight="bold",
             color="#171a1e", y=0.965)
fig.text(
    0.5, 0.925,
    "Daily position at 12:00 local mean solar time in 2026  •  fixed celestial direction",
    ha="center", va="center", fontsize=11.5, color="#3e444b",
)
fig.text(
    0.5, 0.895,
    f"Plano: {LAT:.2f}°N, {abs(LON):.2f}°W   •   Galactic (ℓ,b)=({GAL_L:.0f}°, {GAL_B:.0f}°)   →   "
    f"ICRS (RA, Dec)=({coe.ra.deg:.2f}°, {coe.dec.deg:+.2f}°)",
    ha="center", va="center", fontsize=9.4, color="#555c64",
)

# Color scale is a compact legend for the date progression.
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.045, pad=0.13,
                    aspect=55)
cbar.set_label("Date in 2026 (day of year)", fontsize=10, labelpad=5)
cbar.set_ticks(month_doys)
cbar.set_ticklabels(month_names)
cbar.ax.tick_params(labelsize=8.5, length=3)

fig.subplots_adjust(left=0.095, right=0.985, top=0.84, bottom=0.18)
fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Wrote {OUT}")
print(f"CoE ICRS RA={coe.ra.deg:.8f} deg, Dec={coe.dec.deg:+.8f} deg")
print(f"Altitude range: {alt.min():.3f} to {alt.max():.3f} deg")
```

### 3.3 `make_inner_shell_eclipse_sketch.py` (eclipse supplemental)

**Path:** `/workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/refs_from_chris/make_inner_shell_eclipse_sketch.py`  
(Also present under `cycle6_two_emitters/`.)

```python
# source: /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/refs_from_chris/make_inner_shell_eclipse_sketch.py
#!/usr/bin/env python3
"""Inner-shell eclipse geometry sketch (cycle-6 concentric shells).

R_in=1, R_out=3, D=1.8 → α = arcsin(R_in/D) ≈ 33.7°.
Toward CoE, lines of sight with |θ|<α hit the inner shell only
(eclipsing the outer); outside that disk they hit the outer shell.

Regenerates cycle6_two_emitters/inner_shell_eclipse_sketch.png
to match the approved reference used in letter drafts.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Circle, FancyArrowPatch, Polygon, Wedge

HERE = Path(__file__).resolve().parent
OUT = HERE / "inner_shell_eclipse_sketch.png"

R_IN, R_OUT, D = 1.0, 3.0, 1.8
ALPHA = float(np.arcsin(R_IN / D))  # rad
ALPHA_DEG = float(np.rad2deg(ALPHA))

# Colors matching the approved sketch
C_SHELL = "#8B1A1A"
C_SHELL_FILL = "#F5E6E0"
C_US = "#1E5AA8"
C_ECLIPSE = "#B8C4D0"
C_CAPTION = "#333333"


def main() -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 7.8), dpi=150)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-3.55, 3.55)
    ax.set_ylim(-3.85, 4.05)

    # Title + params
    ax.text(
        0.0,
        3.72,
        "Inner shell eclipses outer toward CoE",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
        color="0.05",
    )
    ax.text(
        0.0,
        3.35,
        rf"$R_{{\rm in}}=1,\ R_{{\rm out}}=3,\ D=1.8$"
        rf"$\qquad\alpha=\arcsin(R_{{\rm in}}/D)\approx{ALPHA_DEG:.1f}^\circ$",
        ha="center",
        va="center",
        fontsize=11,
        color="0.15",
    )

    # Outer fill + rings
    ax.add_patch(
        Circle((0, 0), R_OUT, facecolor=C_SHELL_FILL, edgecolor=C_SHELL, lw=2.2, zorder=1)
    )
    ax.add_patch(
        Circle((0, 0), R_IN, facecolor="white", edgecolor=C_SHELL, lw=2.2, zorder=2)
    )

    # Tangent geometry from US at (D, 0)
    # Tangent length and contact angle from CoE: sin(α)=R_in/D
    # Contact points on inner shell: angle ±(90°-α) from +x toward CoE side...
    # Observer at (D,0); tangents touch at angle φ where cos(φ_from_x_at_CoE related):
    # Contact: point = R_in * (cos ψ, sin ψ) with ψ = ±(π/2 + (π/2-α))? 
    # Standard: contact angles from +x (CoE→US): ψ = ±arccos(R_in/D) = ±(90°-α)... 
    # Actually: vector CoE→contact perpendicular to tangent; contact at angle ±β from US direction
    # where β = arcsin?  Contact angle from +x axis: ±(π - π/2 - α wait)
    # US at (D,0). Contact C on circle |C|=R_in, (C-US)·C = 0 ⇒ C·C = C·US ⇒ R_in^2 = D * C_x
    # ⇒ C_x = R_in^2/D, C_y = ±R_in*sqrt(1-(R_in/D)^2) = ±R_in*cos(α)
    cx = (R_IN**2) / D
    cy = R_IN * np.cos(ALPHA)
    c_top = np.array([cx, cy])
    c_bot = np.array([cx, -cy])
    us = np.array([D, 0.0])

    # Extend tangents past contact (toward CoE / beyond)
    def extend(p0, p1, scale=2.8):
        v = p1 - p0
        v = v / np.linalg.norm(v)
        return p0 + scale * v * R_OUT

    t_top_end = extend(us, c_top, scale=2.6)
    t_bot_end = extend(us, c_bot, scale=2.6)

    # Eclipse wedge (sector between tangents, toward CoE)
    # polygon: US → along top tangent far → arc or through contacts → bottom
    # Fill triangle-ish: US, far along top, far along bot (clipped visually by shells)
    eclipse_pts = np.vstack(
        [
            us,
            t_top_end,
            t_bot_end,
        ]
    )
    ax.add_patch(
        Polygon(
            eclipse_pts,
            closed=True,
            facecolor=C_ECLIPSE,
            edgecolor="none",
            alpha=0.55,
            zorder=0,
        )
    )

    # Redraw shells above eclipse fill for clean edges on near side
    ax.add_patch(
        Circle((0, 0), R_OUT, facecolor="none", edgecolor=C_SHELL, lw=2.2, zorder=3)
    )
    ax.add_patch(
        Circle((0, 0), R_IN, facecolor="white", edgecolor=C_SHELL, lw=2.2, zorder=4)
    )
    # re-fill eclipse only behind inner (clip): draw eclipse wedge again under labels
    # Keep simple: shade between tangents inside outer disk via wedge approx
    # Re-shade using Wedge from US is hard; use polygon clipped conceptually.

    # Tangents (dashed blue)
    for end in (t_top_end, t_bot_end):
        ax.plot(
            [us[0], end[0]],
            [us[1], end[1]],
            ls="--",
            color=C_US,
            lw=1.6,
            zorder=5,
        )

    # Direct LOS toward CoE (arrow to inner shell)
    ax.annotate(
        "",
        xy=(R_IN, 0),
        xytext=(D - 0.02, 0),
        arrowprops=dict(arrowstyle="-|>", color=C_US, lw=1.8, mutation_scale=12),
        zorder=6,
    )

    # Outer hit ray (solid red) — angle slightly outside α
    ang = ALPHA + np.deg2rad(18)
    # Ray from US in direction toward CoE-ish but above: direction = (-cos(ang), sin(ang))
    direction = np.array([-np.cos(ang), np.sin(ang)])
    # Intersect outer circle |US + t d|^2 = R_out^2
    # (us + t d)·(us + t d) = R_OUT^2
    a = np.dot(direction, direction)
    b = 2 * np.dot(us, direction)
    c = np.dot(us, us) - R_OUT**2
    disc = b**2 - 4 * a * c
    t1 = (-b - np.sqrt(disc)) / (2 * a)
    t2 = (-b + np.sqrt(disc)) / (2 * a)
    t_hit = t1 if t1 > 0 else t2
    hit = us + t_hit * direction
    ax.annotate(
        "",
        xy=hit,
        xytext=us,
        arrowprops=dict(arrowstyle="-|>", color=C_SHELL, lw=1.7, mutation_scale=11),
        zorder=6,
    )
    ax.plot(hit[0], hit[1], "o", color=C_SHELL, ms=5, zorder=7)
    ax.text(
        hit[0] + 0.12,
        hit[1] + 0.28,
        "hits outer",
        color=C_SHELL,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    # Contact dots
    ax.plot(c_top[0], c_top[1], "o", color=C_SHELL, ms=5, zorder=7)
    ax.plot(c_bot[0], c_bot[1], "o", color=C_SHELL, ms=5, zorder=7)

    # Alpha arc
    arc_r = 0.55
    ax.add_patch(
        Arc(
            us,
            2 * arc_r,
            2 * arc_r,
            angle=0,
            theta1=180 - ALPHA_DEG,
            theta2=180,
            color=C_US,
            lw=1.5,
            zorder=8,
        )
    )
    ax.text(
        us[0] - arc_r * 0.85,
        0.22,
        r"$\alpha$",
        color=C_US,
        fontsize=13,
        ha="center",
        va="bottom",
        fontweight="bold",
    )

    # Labels
    ax.plot(0, 0, "o", color="k", ms=7, zorder=8)
    ax.text(0, -0.28, "CoE", ha="center", va="top", fontsize=12, fontweight="bold")
    ax.plot(D, 0, "o", color=C_US, ms=11, zorder=8)
    ax.text(D, -0.32, "us", color=C_US, ha="center", va="top", fontsize=13, fontweight="bold")

    ax.text(
        0,
        R_OUT + 0.18,
        "outer shell",
        color=C_SHELL,
        ha="center",
        va="bottom",
        fontsize=11,
        fontweight="bold",
    )
    ax.text(
        0,
        R_IN + 0.12,
        "inner shell",
        color=C_SHELL,
        ha="center",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        -1.15,
        0.55,
        "hits inner\n(eclipses outer)",
        color=C_US,
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        zorder=9,
    )

    # Footer caption
    ax.text(
        0.0,
        -3.55,
        r"Sky disk $|\theta|<\alpha$ (and $360^\circ-\alpha$ twin): "
        r"lines of sight hit inner shell only",
        ha="center",
        va="center",
        fontsize=10.5,
        color=C_CAPTION,
        fontweight="bold",
    )

    fig.savefig(OUT, dpi=150, bbox_inches="tight", facecolor="white", pad_inches=0.15)
    plt.close(fig)
    print(f"wrote {OUT}  α={ALPHA_DEG:.2f}°")
    return OUT


if __name__ == "__main__":
    main()
```

## 4. Full sources — related (panels / non-exact)

### 4.1 Fig06 panel generator — `make_migkas_fig8_overlay_matched.py`

Produces `fig8_overlay_matched{,_grad,_grad_rot45,_grad_rot81}.png` (1975×1515).  
**Does not by itself emit letter Fig06 stack 2164×3222** — use `att14_overlay_stack.png` for exact.

**Stacking notes (Fig06):**
1. Run this script (optionally with `SHOW_LEGEND=True` temporarily) from `/workspace/cosmology/figures` with Migkas Fig.8 raster available (`/workspace/cosmology/map_compare/FULL_migkas_theta75.png` per run logs).
2. Outputs land under `/workspace/cosmology/papers/spreadsheet_model/` by default (`SM_DIR` in script).
3. Letter PDF embed is a **special stacked/cropped composition** `att14_overlay_stack.png` (2164×3222), bitwise=gold.
4. Naive `PIL` vertical concat of the three legend panels ≈ 1975×4689 (what `figs_regenerated_from_code` wrote) — **not** gold.
5. `figs/fig8_overlays_three_stack.png` is same dimensions as letter stack but **not** bitwise identical (matches `att14_overlay_stack.png.bak_before_restack`).

```python
# source: /workspace/cosmology/figures/make_migkas_fig8_overlay_matched.py
#!/usr/bin/env python3
"""
Superimpose CoE MAX/MIN highlight bands + Migkas markers onto the REAL Migkas
Fig. 8 theta=75 A/A_all Mollweide map.

Matches make_coe_rings_migkas_pin_axis.py (three-panel Panel B) geometry:
  Exactly FOUR features — each ONE transparent highlight band, 18° end-to-end
  (= ±9° = 5% of 360°) about the center θ + ONE solid mid-stroke:
    MAX 30°  → [21°, 39°]   color #00e5ff (bright cyan)
    MAX 185° → [176°, 194°] color #ff2d95 (magenta/hot pink)
    MIN 115° → [106°, 124°] color #b8f200 (lime)
    MIN 265° → [256°, 274°] color #4fc3f7 (sky blue)
  Drawing: clean PolyCollection fill + one mid-stroke (optional thin white
  underlay on mid only). NO black rim, NO thick white halo stacks, NO dual
  boundary-edge strokes. No 283° companion. No Dec=0 overlay.
  MAX triangles: (17°, −9°) and (185°, 0°)
  MIN squares:   (115°, 0°) and (274°, −22°)
  CoE axis tilted +30°: anti-CoE (0°,+30°), CoE (180°,−30°) (purple).

Registration (unchanged from working Fig. 8):
  GOOD_REG = (cx, cy, Rx, Ry) = (863.0, 763.25, 693.9927, 392.542425)
  Source: FULL_migkas_theta75.png
  CoE axis tilted +30°; clip/inset so overlays stay inside colored oval
  lon_disp = -wrap(ℓ)

Outputs:
  migkas_fig8_theta75_with_coe_rings.png
  migkas_fig8_theta75_with_coe_rings.json
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
from PIL import Image

OUT_DIR = Path(__file__).resolve().parent
SM_DIR = Path("/workspace/cosmology/papers/spreadsheet_model")
PNG_PATH = SM_DIR / "fig8_overlay_matched.png"
JSON_PATH = SM_DIR / "fig8_overlay_matched.json"

# Versions = panel B of each three-panel figure (grad: opacity gradient, roll: deg CCW about CoE axis)
VERSIONS = [
    ("matched", False, 0.0),
    ("matched_grad", True, 0.0),
    ("matched_grad_rot45", True, 45.0),
    ("matched_grad_rot81", True, 80.8),
    ("matched_rot81", False, 80.8),  # roll + all loops, no opacity gradient
]
GRAD = False
ROLL = 0.0
SHOW_LEGEND = False  # letter Fig.6: no legend box/strip

SOURCE_CANDIDATES = [
    Path("/workspace/cosmology/map_compare/FULL_migkas_theta75.png"),
    Path("/workspace/cosmology/map_compare/migkas_fig8_theta75_crop.png"),
    Path("/workspace/cosmology/migkas/pub_maps/Fig8_full_hires.png"),
]

# Strip A(ℓ) scan stops from +x (after-2ℓ A_after extrema) — four band centers.
# 5% of 360° = 18° end-to-end → half-width ±9° about each center (match three-panel).
BAND_HALF = 9.0
BAND_FULL = 2.0 * BAND_HALF  # 18°
THETA_MAX = [29.5, 187.5]
THETA_MIN = [122.6, 264.4]
MAX_BANDS = [(t - BAND_HALF, t + BAND_HALF) for t in THETA_MAX]  # [21,39], [176,194]
MIN_BANDS = [(t - BAND_HALF, t + BAND_HALF) for t in THETA_MIN]  # [106,124], [256,274]


def mu_pair_for_band(lo, hi):
    """Edge μ=cos(θ) for a scan-angle band (μ-space, not θ-space)."""
    return float(np.cos(np.deg2rad(lo))), float(np.cos(np.deg2rad(hi)))

# Visually-good registration (pre medial-axis retune). Hard-set as working default.
GOOD_REG = {
    "cx": 863.0,
    "cy": 763.25,
    "Rx": 693.9927,
    "Ry": 392.542425,
}

# Sky markers (ℓ, b): map extrema + strip stand-ins (same as 3-panel Panel B).
MIGKAS_MAXES = [(17.0, -9.0), (185.0, 0.0)]
MIGKAS_MINS = [(115.0, 0.0), (274.0, -22.0)]

# Galactic pose (Chris 2026-09-22): anti-CoE (+D/+x) at (0°, +30°); CoE at (180°, −30°)
ANTI_COE_LB = (0.0, 15.0)
COE_LB = (180.0, -15.0)
AXIS_TILT_B_DEG = 15.0


def lb_to_n(l_deg, b_deg):
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    return np.array(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], dtype=float
    )


D_HAT = lb_to_n(*ANTI_COE_LB)  # +D = anti-CoE = +x


def great_circle_through_poles(D, n=721):
    """Great circle through ±D (CoE axis meridian on the sky)."""
    a = D / np.linalg.norm(D)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([1.0, 0.0, 0.0])
    e2 = np.cross(a, tmp)
    e2 /= np.linalg.norm(e2)
    t = np.linspace(0.0, 2.0 * np.pi, n)
    pts = np.outer(np.cos(t), a) + np.outer(np.sin(t), e2)
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    return ls, bs

# High-contrast loop colors that pop on magma (avoid reds/oranges/corals)
C_MAX_30 = "#0097b2"    # cyan — MAX 29.5° (three-panel matched)
C_MAX_185 = "#e0107a"   # magenta — MAX 187.5°
C_MIN_115 = "#6fa800"   # green — MIN 122.6°
C_MIN_265 = "#2166d4"   # blue — MIN 264.4°
COLORS_MAX = [C_MAX_30, C_MAX_185]
COLORS_MIN = [C_MIN_115, C_MIN_265]
C_COE = "#6b3fa0"
C_ANTI = "#8e1b1b"  # dark red marker; distinct from loop reds (match three-panel)
C_AXIS = "#5b2c6f"
C_EQ = "#1e8449"


def wrap_lon_deg(lon):
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def lon_disp_from_ell(ell_deg):
    return -wrap_lon_deg(ell_deg)


def n_to_lb(n):
    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    b = float(np.rad2deg(np.arcsin(np.clip(n[2], -1.0, 1.0))))
    l = float(np.rad2deg(np.arctan2(n[1], n[0])) % 360.0)
    return l, b


def celestial_equator_lb(n=2500):
    ra_ngp = np.deg2rad(192.85948)
    dec_ngp = np.deg2rad(27.12825)
    l_ncp = np.deg2rad(122.93192)
    ra = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    dec = 0.0
    sin_b = np.clip(
        np.sin(dec) * np.sin(dec_ngp)
        + np.cos(dec) * np.cos(dec_ngp) * np.cos(ra - ra_ngp),
        -1.0,
        1.0,
    )
    b = np.arcsin(sin_b)
    y = np.cos(dec) * np.sin(ra - ra_ngp)
    x = (
        np.sin(dec) * np.cos(dec_ngp)
        - np.cos(dec) * np.sin(dec_ngp) * np.cos(ra - ra_ngp)
    )
    l = (l_ncp - np.arctan2(y, x)) % (2.0 * np.pi)
    return np.rad2deg(l), np.rad2deg(b)


def small_circle_about_axis(axis, mu, n=721):
    a = axis / np.linalg.norm(axis)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    tmp = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(a, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    t = np.linspace(0.0, 2.0 * np.pi, n)
    return (mu * a)[None, :] + rho * (
        np.outer(np.cos(t), e1) + np.outer(np.sin(t), e2)
    )


def mollweide_xy(lon_rad, lat_rad, n_iter=64):
    """Standard Mollweide; stable bracketed solve of 2θ+sin(2θ)=π sin(lat)."""
    lon = np.asarray(lon_rad, dtype=float)
    lat = np.asarray(lat_rad, dtype=float)
    target = np.pi * np.sin(np.clip(lat, -np.pi / 2, np.pi / 2))
    lo = np.full_like(lat, -np.pi / 2, dtype=float)
    hi = np.full_like(lat, np.pi / 2, dtype=float)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        g = 2.0 * mid + np.sin(2.0 * mid)
        lo = np.where(g < target, mid, lo)
        hi = np.where(g >= target, mid, hi)
    theta = 0.5 * (lo + hi)
    theta = np.where(lat >= np.pi / 2 - 1e-15, np.pi / 2, theta)
    theta = np.where(lat <= -np.pi / 2 + 1e-15, -np.pi / 2, theta)
    x = (2.0 * np.sqrt(2.0) / np.pi) * lon * np.cos(theta)
    y = np.sqrt(2.0) * np.sin(theta)
    return x, y


def lb_to_px(ell_deg, b_deg, cx, cy, Rx, Ry):
    lon = np.deg2rad(lon_disp_from_ell(np.asarray(ell_deg, dtype=float)))
    lat = np.deg2rad(np.asarray(b_deg, dtype=float))
    x, y = mollweide_xy(lon, lat)
    Xn = x / (2.0 * np.sqrt(2.0))
    Yn = y / np.sqrt(2.0)
    return cx + Xn * Rx, cy - Yn * Ry


def as_xy(px, py):
    return float(np.asarray(px).reshape(-1)[0]), float(np.asarray(py).reshape(-1)[0])


def chamfer_distance_to_dark(dark):
    h, w = dark.shape
    INF = 1e6
    D = np.where(dark, 0.0, INF)
    for y in range(h):
        for x in range(w):
            v = D[y, x]
            if x > 0:
                v = min(v, D[y, x - 1] + 1.0)
            if y > 0:
                v = min(v, D[y - 1, x] + 1.0)
            if x > 0 and y > 0:
                v = min(v, D[y - 1, x - 1] + 1.414)
            if x + 1 < w and y > 0:
                v = min(v, D[y - 1, x + 1] + 1.414)
            D[y, x] = v
    for y in range(h - 1, -1, -1):
        for x in range(w - 1, -1, -1):
            v = D[y, x]
            if x + 1 < w:
                v = min(v, D[y, x + 1] + 1.0)
            if y + 1 < h:
                v = min(v, D[y + 1, x] + 1.0)
            if x + 1 < w and y + 1 < h:
                v = min(v, D[y + 1, x + 1] + 1.414)
            if x > 0 and y + 1 < h:
                v = min(v, D[y + 1, x - 1] + 1.414)
            D[y, x] = v
    return D


def detect_oval_bbox(rgb):
    h, w, _ = rgb.shape
    white = (rgb[:, :, 0] > 250) & (rgb[:, :, 1] > 250) & (rgb[:, :, 2] > 250)
    oval_left = np.full(h, np.nan)
    oval_right = np.full(h, np.nan)
    for y in range(h):
        row = ~white[y]
        d = np.diff(np.r_[0, row.astype(int), 0])
        starts = np.where(d == 1)[0]
        ends = np.where(d == -1)[0] - 1
        if len(starts) == 0:
            continue
        lengths = ends - starts + 1
        candidates = [
            (int(lengths[i]), int(starts[i]), int(ends[i]))
            for i in range(len(starts))
            if starts[i] < 0.72 * w and lengths[i] > 40
        ]
        if not candidates:
            continue
        candidates.sort(reverse=True)
        _, L, R = candidates[0]
        oval_left[y] = L
        oval_right[y] = R
    widths = np.where(np.isfinite(oval_left), oval_right - oval_left, 0.0)
    yeq = int(np.argmax(widths))
    thresh = 0.08 * widths[yeq]
    y0 = yeq
    while y0 > 0 and widths[y0 - 1] >= thresh:
        y0 -= 1
    y1 = yeq
    while y1 < h - 1 and widths[y1 + 1] >= thresh:
        y1 += 1
    x0 = float(oval_left[yeq])
    x1 = float(oval_right[yeq])
    return {
        "x0": int(x0),
        "x1": int(x1),
        "y0": int(y0),
        "y1": int(y1),
        "yeq": int(yeq),
        "cx": 0.5 * (x0 + x1),
        "cy": 0.5 * (y0 + y1),
        "Rx": 0.5 * (x1 - x0),
        "Ry": 0.5 * (y1 - y0 + 1),
    }


def refine_registration(rgb, seed):
    """Gray/dark Dec=0 lock (pre medial-axis). Prefer GOOD_REG hard-set in main."""
    h, w, _ = rgb.shape
    dark = rgb.mean(axis=2) < 55
    dist = chamfer_distance_to_dark(dark)

    ell_eq, b_eq = celestial_equator_lb(2500)
    lon = np.deg2rad(lon_disp_from_ell(ell_eq))
    lat = np.deg2rad(b_eq)
    mx, my = mollweide_xy(lon, lat)
    Xn = mx / (2.0 * np.sqrt(2.0))
    Yn = my / np.sqrt(2.0)

    def median_dist(cx, cy, Rx, Ry):
        px = cx + Xn * Rx
        py = cy - Yn * Ry
        m = (px >= 1) & (px < w - 1) & (py >= 1) & (py < h - 1)
        xi = np.clip(np.round(px[m]).astype(int), 0, w - 1)
        yi = np.clip(np.round(py[m]).astype(int), 0, h - 1)
        return float(np.median(dist[yi, xi]))

    base = {k: float(seed[k]) for k in ("cx", "cy", "Rx", "Ry")}
    best_sc = median_dist(**base)
    best = dict(base)
    for dcx in np.linspace(-10, 10, 21):
        for dcy in np.linspace(-10, 10, 21):
            for sRx in (0.98, 0.99, 1.0, 1.01, 1.02):
                for sRy in (0.98, 0.99, 1.0, 1.01, 1.02):
                    p = dict(
                        cx=base["cx"] + dcx,
                        cy=base["cy"] + dcy,
                        Rx=base["Rx"] * sRx,
                        Ry=base["Ry"] * sRy,
                    )
                    sc = median_dist(**p)
                    if sc < best_sc:
                        best_sc = sc
                        best = p
    b = best
    for dcx in np.linspace(-2, 2, 17):
        for dcy in np.linspace(-2, 2, 17):
            for sRx in np.linspace(0.995, 1.005, 11):
                for sRy in np.linspace(0.995, 1.005, 11):
                    p = dict(
                        cx=b["cx"] + dcx,
                        cy=b["cy"] + dcy,
                        Rx=b["Rx"] * sRx,
                        Ry=b["Ry"] * sRy,
                    )
                    sc = median_dist(**p)
                    if sc < best_sc:
                        best_sc = sc
                        best = p
    best["ceq_median_px"] = float(best_sc)
    return best


def choose_source():
    for p in SOURCE_CANDIDATES:
        if p.is_file():
            return p
    raise FileNotFoundError("No Migkas Fig.8 theta=75 source image found")


def sample_rgb(rgb, ell, b, cx, cy, Rx, Ry, r=4):
    px, py = lb_to_px(ell, b, cx, cy, Rx, Ry)
    xi, yi = int(round(as_xy(px, py)[0])), int(round(as_xy(px, py)[1]))
    patch = rgb[max(0, yi - r) : yi + r + 1, max(0, xi - r) : xi + r + 1]
    return patch.mean(axis=(0, 1))


def validate_landmarks(rgb, reg):
    cx, cy, Rx, Ry = reg["cx"], reg["cy"], reg["Rx"], reg["Ry"]
    pts = {}
    landmarks = [
        ("center", 0.0, 0.0),
        ("ell90", 90.0, 0.0),
        ("ell270", 270.0, 0.0),
        ("north", 0.0, 90.0),
        ("south", 0.0, -90.0),
    ]
    for i, (ell, b) in enumerate(MIGKAS_MAXES):
        landmarks.append((f"migkas_max{i}", ell, b))
    for i, (ell, b) in enumerate(MIGKAS_MINS):
        landmarks.append((f"migkas_min{i}", ell, b))
    for name, ell, b in landmarks:
        x, y = as_xy(*lb_to_px(ell, b, cx, cy, Rx, Ry))
        pts[name] = {"ell": ell, "b": b, "px": x, "py": y}

    assert pts["ell90"]["px"] < pts["center"]["px"], "l=90 must be LEFT of center"
    assert pts["ell270"]["px"] > pts["center"]["px"], "l=270 must be RIGHT of center"
    assert abs(pts["center"]["px"] - cx) < 1.5
    assert abs(pts["center"]["py"] - cy) < 1.5
    assert pts["north"]["py"] < pts["center"]["py"]
    assert pts["south"]["py"] > pts["center"]["py"]

    # Map-extrema pair: (17,-9) should sit in a warmer patch than (274,-22).
    cmax = sample_rgb(rgb, *MIGKAS_MAXES[0], cx, cy, Rx, Ry)
    cmin = sample_rgb(rgb, *MIGKAS_MINS[1], cx, cy, Rx, Ry)
    assert (cmax[0] + cmax[1]) > (cmin[0] + cmin[1]) + 50
    return {
        "landmarks_px": pts,
        "migkas_max0_rgb_mean": [float(x) for x in cmax],
        "migkas_min1_rgb_mean": [float(x) for x in cmin],
    }


def inside_ellipse(px, py, ecx, ecy, eRx, eRy, r_lim=1.0):
    """Boolean mask: points inside axis-aligned ellipse (image y down)."""
    px = np.asarray(px, dtype=float)
    py = np.asarray(py, dtype=float)
    Xn = (px - ecx) / eRx
    Yn = (py - ecy) / eRy
    return Xn * Xn + Yn * Yn <= float(r_lim) * float(r_lim)


def plot_lb_curve(
    ax, ell, b, cx, cy, Rx, Ry, inset=1.0, jump_px=45.0,
    clip_r2=None, clip_path=None, color_oval=None, **kw
):
    """Draw curve in oval; inset<1 keeps strokes inside the published Mollweide edge."""
    px, py = lb_to_px(ell, b, cx, cy, Rx, Ry)
    px = np.asarray(px, dtype=float).reshape(-1)
    py = np.asarray(py, dtype=float).reshape(-1)
    # Pull points slightly toward oval center so thick strokes do not spill past the rim.
    px = cx + (px - cx) * inset
    py = cy + (py - cy) * inset
    Xn = (px - cx) / Rx
    Yn = (cy - py) / Ry
    r2_lim = (inset * inset) if clip_r2 is None else float(clip_r2)
    inside = Xn * Xn + Yn * Yn <= r2_lim
    # Also require published colored oval (seed) so registration oversize cannot bleed.
    if color_oval is not None:
        ocx, ocy, oRx, oRy, o_r = color_oval
        inside = inside & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)
    px = np.where(inside, px, np.nan)
    py = np.where(inside, py, np.nan)

    segs = []
    cur_x, cur_y = [], []
    for x, y in zip(px, py):
        if not (np.isfinite(x) and np.isfinite(y)):
            if len(cur_x) >= 2:
                segs.append((cur_x, cur_y))
            cur_x, cur_y = [], []
            continue
        if cur_x and np.hypot(x - cur_x[-1], y - cur_y[-1]) > jump_px:
            if len(cur_x) >= 2:
                segs.append((cur_x, cur_y))
            cur_x, cur_y = [float(x)], [float(y)]
        else:
            cur_x.append(float(x))
            cur_y.append(float(y))
    if len(cur_x) >= 2:
        segs.append((cur_x, cur_y))
    if clip_path is not None:
        kw["clip_path"] = clip_path
        kw["clip_on"] = True
    for sx, sy in segs:
        (line,) = ax.plot(sx, sy, **kw)
        if clip_path is not None:
            line.set_clip_path(clip_path)
            line.set_clip_on(True)


def plot_closed_ring(
    ax, mu, cx, cy, Rx, Ry, color, lw, inset=1.0, jump_px=80.0,
    clip_path=None, color_oval=None, ls="-", wash_alpha=None,
):
    """
    Exact projected small-circle about D_HAT=+x: n·x̂ = μ.

    Parameterize on the sphere, convert to Galactic (ℓ,b), apply
    lon_disp = -wrap(ℓ), Mollweide-project, and stroke continuous polylines.
    Split wherever display-longitude jumps across the antimeridian (or pixel
    distance jumps). Never invent interior chords, ellipse fits, or decorative
    oval-rim bridges.

    Notes on θ labels used elsewhere:
      μ=cos(30°)  → r=30° about +x (closed oval near map center)
      μ=cos(185°)=cos(175°) → r=175° about +x ≡ 5° about CoE
        → small arcs on BOTH left and right map edges (not a rim loop)
      μ=cos(115°) → r=115° about +x (large seam-crossing curve)
      μ=cos(265°)=cos(95°) → r=95° about +x (near-meridian, seam-crossing)
    """
    pts = small_circle_about_axis(D_HAT, mu, n=4001)
    # Keep sampling order around the axis (azimuth) — stable for splitting.
    phi = np.arctan2(pts[:, 2], pts[:, 1])
    order = np.argsort(phi)
    pts = pts[order]

    ell = np.empty(len(pts))
    b = np.empty(len(pts))
    for i, p in enumerate(pts):
        ell[i], b[i] = n_to_lb(p)

    lon = np.deg2rad(lon_disp_from_ell(ell))  # λ ∈ (-π, π]
    lat = np.deg2rad(b)

    # Split on antimeridian jumps in display longitude BEFORE pixel projection.
    dlon = np.diff(lon)
    cuts = np.where(np.abs(dlon) > np.pi)[0] + 1
    # Also close the periodic sample: check wrap between last and first.
    # Build open segments from azimuth-ordered points.
    idx = np.concatenate([[0], cuts, [len(lon)]])
    lon_segs = []
    for i0, i1 in zip(idx[:-1], idx[1:]):
        if i1 - i0 >= 2:
            lon_segs.append((lon[i0:i1].copy(), lat[i0:i1].copy()))

    # Circle sampling is cyclic. Merge the artificial last→first break when
    # that periodic join is not itself an antimeridian crossing. Otherwise a
    # seam-crossing circle falsely appears as three pieces with a tiny hole.
    if len(lon_segs) > 1 and abs(lon[0] - lon[-1]) <= np.pi:
        first_lo, first_la = lon_segs[0]
        last_lo, last_la = lon_segs[-1]
        lon_segs = [(
            np.concatenate([last_lo, first_lo]),
            np.concatenate([last_la, first_la]),
        )] + lon_segs[1:-1]

    # Only now pin TRUE seam endpoints to ±π. Never snap the arbitrary
    # parameterization start/end to the boundary.
    snapped = []
    for lo, la in lon_segs:
        if abs(abs(lo[0]) - np.pi) < np.deg2rad(10):
            lo[0] = np.sign(lo[0] if lo[0] != 0 else 1.0) * np.pi
        if abs(abs(lo[-1]) - np.pi) < np.deg2rad(10):
            lo[-1] = np.sign(lo[-1] if lo[-1] != 0 else 1.0) * np.pi
        snapped.append((lo, la))
    lon_segs = snapped

    def stroke_xy(xs, ys):
        if len(xs) < 2:
            return
        if wash_alpha is not None:
            # Soft μ-wash stroke: no white underlay; translucent color only.
            layers = ((color, lw, 2, float(wash_alpha)),)
        else:
            # Clean mid-stroke: optional thin white underlay only (lw+1.5), then color.
            # No black outer rim, no thick white halo stack.
            layers = (
                ("white", lw + 1.5, 3, 0.95),
                (color, lw, 4, 1.0),
            )
        for col, width, z, a in layers:
            (line,) = ax.plot(
                xs, ys, color=col, lw=width, linestyle=ls,
                solid_capstyle="round", solid_joinstyle="round",
                zorder=z, alpha=a,
                clip_path=clip_path, clip_on=True,
            )
            if clip_path is not None:
                line.set_clip_path(clip_path)
                line.set_clip_on(True)

    r2_lim = inset * inset + 1e-10
    drawn = 0
    for lo, la in lon_segs:
        # Mollweide in normalized oval coords, then pixels.
        x_m, y_m = mollweide_xy(lo, la)
        Xn = x_m / (2.0 * np.sqrt(2.0))
        Yn = y_m / np.sqrt(2.0)
        px = cx + Xn * Rx
        py = cy - Yn * Ry
        # Mild inset so thick strokes stay inside published oval (no rim-snap).
        px = cx + (px - cx) * inset
        py = cy + (py - cy) * inset
        Xn2 = (px - cx) / Rx
        Yn2 = (cy - py) / Ry
        inside = (Xn2 * Xn2 + Yn2 * Yn2) <= r2_lim
        if color_oval is not None:
            ocx, ocy, oRx, oRy, o_r = color_oval
            inside = inside & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)

        # Pixel-space split on NaN / large jumps / outside oval.
        cur_x, cur_y = [], []
        for x, y, ok in zip(px, py, inside):
            if (not ok) or (not np.isfinite(x)) or (not np.isfinite(y)):
                if len(cur_x) >= 2:
                    stroke_xy(cur_x, cur_y)
                    drawn += 1
                cur_x, cur_y = [], []
                continue
            if cur_x and np.hypot(x - cur_x[-1], y - cur_y[-1]) > jump_px:
                if len(cur_x) >= 2:
                    stroke_xy(cur_x, cur_y)
                    drawn += 1
                cur_x, cur_y = [float(x)], [float(y)]
            else:
                cur_x.append(float(x))
                cur_y.append(float(y))
        if len(cur_x) >= 2:
            # Close only if this single segment is already nearly a loop
            # (non-wrapping ring such as MAX 30°).
            if (
                drawn == 0
                and len(lon_segs) == 1
                and np.hypot(cur_x[0] - cur_x[-1], cur_y[0] - cur_y[-1]) < jump_px
            ):
                cur_x.append(cur_x[0])
                cur_y.append(cur_y[0])
            stroke_xy(cur_x, cur_y)
            drawn += 1



def _ring_px_arrays(mu, cx, cy, Rx, Ry, inset, n_pts=721):
    """Small-circle n·D=μ → pixel (px, py) arrays, azimuth-ordered, with inset."""
    pts = small_circle_about_axis(D_HAT, float(mu), n=n_pts)
    phi = np.arctan2(pts[:, 2], pts[:, 1])
    order = np.argsort(phi)
    pts = pts[order]
    ell = np.empty(n_pts)
    b = np.empty(n_pts)
    for i, p in enumerate(pts):
        ell[i], b[i] = n_to_lb(p)
    lon = np.deg2rad(lon_disp_from_ell(ell))
    lat = np.deg2rad(b)
    x_m, y_m = mollweide_xy(lon, lat)
    Xn = x_m / (2.0 * np.sqrt(2.0))
    Yn = y_m / np.sqrt(2.0)
    px = cx + Xn * Rx
    py = cy - Yn * Ry
    px = cx + (px - cx) * inset
    py = cy + (py - cy) * inset
    return lon, px, py


def fill_highlight_band(
    ax, theta_lo, theta_hi, theta_mid, cx, cy, Rx, Ry, color,
    inset=1.0, clip_path=None, color_oval=None, n_pts=721, jump_px=80.0,
):
    """
    One highlight band matching three-panel Panel B:
      continuous transparent polygon strips between ±BAND_HALF edges
      (broken only at Mollweide antimeridian / large pixel jumps — NOT a dense
      multi-stroke wash) + ONE solid mid-stroke at θ_mid.
    """
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    lon_a, px_a, py_a = _ring_px_arrays(mu_a, cx, cy, Rx, Ry, inset, n_pts=n_pts)
    lon_b, px_b, py_b = _ring_px_arrays(mu_b, cx, cy, Rx, Ry, inset, n_pts=n_pts)

    r2_lim = inset * inset + 1e-10

    def _inside(px, py):
        Xn2 = (px - cx) / Rx
        Yn2 = (cy - py) / Ry
        ok = (Xn2 * Xn2 + Yn2 * Yn2) <= r2_lim
        if color_oval is not None:
            ocx, ocy, oRx, oRy, o_r = color_oval
            ok = ok & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)
        return ok

    gap = 0.9 * np.pi
    bad = np.zeros(n_pts - 1, dtype=bool)
    for i in range(n_pts - 1):
        if abs(lon_a[i + 1] - lon_a[i]) > gap:
            bad[i] = True
        elif abs(lon_b[i + 1] - lon_b[i]) > gap:
            bad[i] = True
        elif abs(lon_a[i] - lon_b[i]) > gap:
            bad[i] = True
        elif abs(lon_a[i + 1] - lon_b[i + 1]) > gap:
            bad[i] = True
        elif np.hypot(px_a[i + 1] - px_a[i], py_a[i + 1] - py_a[i]) > jump_px:
            bad[i] = True
        elif np.hypot(px_b[i + 1] - px_b[i], py_b[i + 1] - py_b[i]) > jump_px:
            bad[i] = True
        elif (not _inside(px_a[i], py_a[i])) or (not _inside(px_b[i], py_b[i])):
            bad[i] = True
        elif (not _inside(px_a[i + 1], py_a[i + 1])) or (not _inside(px_b[i + 1], py_b[i + 1])):
            bad[i] = True

    polys = []
    i = 0
    while i < n_pts - 1:
        while i < n_pts - 1 and bad[i]:
            i += 1
        if i >= n_pts - 1:
            break
        j = i
        while j < n_pts - 1 and not bad[j]:
            j += 1
        i0, i1 = i, j
        if i1 - i0 >= 1:
            xs = np.concatenate([px_a[i0:i1 + 1], px_b[i0:i1 + 1][::-1]])
            ys = np.concatenate([py_a[i0:i1 + 1], py_b[i0:i1 + 1][::-1]])
            polys.append(np.column_stack([xs, ys]))
        i = j + 1

    if polys:
        # Soft transparent highlight (clean pre-halo look).
        coll = PolyCollection(
            polys, facecolors=color, edgecolors="none",
            alpha=0.30, zorder=2, linewidths=0,
        )
        ax.add_collection(coll)
        if clip_path is not None:
            coll.set_clip_path(clip_path)
            coll.set_clip_on(True)

    # ONE solid mid-stroke only (thin white underlay lw+1.5). No boundary-edge
    # strokes — dual edges with halos multiply visual lines into "extra loops".
    plot_closed_ring(
        ax, mu_mid, cx, cy, Rx, Ry, color, lw=2.8,
        inset=inset, clip_path=clip_path, color_oval=color_oval,
    )



# --- Opacity gradient (same rule as make_coe_rings_migkas_pin_axis_matched_grad*.py):
# alpha *= (1+cos psi)/2, psi = ring azimuth from the θ-ray point (ring point most
# along sign(sin θ)·Ŷ, Ŷ = Galactic +z projected ⊥ D), rolled by ROLL (Y' = cY − sZ).
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba


def _ring_basis_roll(D, roll_deg):
    a = D / np.linalg.norm(D)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    y = tmp - np.dot(tmp, a) * a
    y /= np.linalg.norm(y)
    z = np.cross(a, y)
    z /= np.linalg.norm(z)
    r = np.deg2rad(float(roll_deg))
    c, s_ = np.cos(r), np.sin(r)
    return c * y - s_ * z, s_ * y + c * z


def _grad_ring_px(mu, theta_mid, cx, cy, Rx, Ry, inset, n):
    a = D_HAT / np.linalg.norm(D_HAT)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    rho = np.sqrt(1.0 - mu * mu)
    y, z = _ring_basis_roll(a, ROLL)
    phi0 = 0.5 * np.pi if np.sin(np.deg2rad(theta_mid)) >= 0 else 1.5 * np.pi
    phi = phi0 + np.linspace(-np.pi, np.pi, n)
    pts = (mu * a)[None, :] + rho * (np.outer(np.sin(phi), y) + np.outer(np.cos(phi), z))
    ell = np.array([n_to_lb(p)[0] for p in pts]); b = np.array([n_to_lb(p)[1] for p in pts])
    lon = np.deg2rad(lon_disp_from_ell(ell)); lat = np.deg2rad(b)
    x_m, y_m = mollweide_xy(lon, lat)
    px = cx + x_m / (2.0 * np.sqrt(2.0)) * Rx
    py = cy - y_m / np.sqrt(2.0) * Ry
    px = cx + (px - cx) * inset; py = cy + (py - cy) * inset
    return lon, px, py, 0.5 * (1.0 + np.cos(phi - phi0))


def _rgba(color, alphas):
    r, g, b_, _ = to_rgba(color)
    al = np.clip(np.asarray(alphas, dtype=float), 0.0, 1.0)
    out = np.empty((al.size, 4)); out[:, 0], out[:, 1], out[:, 2], out[:, 3] = r, g, b_, al
    return out


def _runs(seg_sc, valid, n_levels=64):
    lev = np.round(np.clip(seg_sc, 0.0, 1.0) * n_levels) / n_levels
    ok = (lev > 0) & valid
    runs, i, n = [], 0, lev.size
    while i < n:
        if not ok[i]:
            i += 1; continue
        j = i + 1
        while j < n and ok[j] and lev[j] == lev[i]:
            j += 1
        runs.append((i, j, float(lev[i]))); i = j
    return runs


def fill_highlight_band_grad(
    ax, theta_lo, theta_hi, theta_mid, cx, cy, Rx, Ry, color,
    inset=1.0, clip_path=None, color_oval=None, n_pts=721, jump_px=80.0,
):
    """Same band (alpha 0.30) + mid-stroke (white lw+1.5 @0.95 under colour lw 2.8)
    as fill_highlight_band, but with the (1+cos psi)/2 opacity gradient."""
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    r2_lim = inset * inset + 1e-10

    def _inside(px, py):
        ok = ((px - cx) / Rx) ** 2 + ((cy - py) / Ry) ** 2 <= r2_lim
        if color_oval is not None:
            ocx, ocy, oRx, oRy, o_r = color_oval
            ok = ok & inside_ellipse(px, py, ocx, ocy, oRx, oRy, o_r)
        return ok

    def _clip(coll):
        if clip_path is not None:
            coll.set_clip_path(clip_path); coll.set_clip_on(True)
        ax.add_collection(coll)

    gap = 0.9 * np.pi
    lon_a, pxa, pya, sc = _grad_ring_px(mu_a, theta_mid, cx, cy, Rx, Ry, inset, n_pts)
    lon_b, pxb, pyb, _ = _grad_ring_px(mu_b, theta_mid, cx, cy, Rx, Ry, inset, n_pts)
    qa = np.stack([lon_a[:-1], lon_a[1:], lon_b[1:], lon_b[:-1]])
    ins_a, ins_b = _inside(pxa, pya), _inside(pxb, pyb)
    valid = ((qa.max(0) - qa.min(0)) <= gap)
    valid &= np.hypot(np.diff(pxa), np.diff(pya)) <= jump_px
    valid &= np.hypot(np.diff(pxb), np.diff(pyb)) <= jump_px
    valid &= ins_a[:-1] & ins_a[1:] & ins_b[:-1] & ins_b[1:]
    runs = _runs(0.5 * (sc[:-1] + sc[1:]), valid)
    if runs:
        polys = [np.column_stack([np.concatenate([pxa[i0:i1 + 1], pxb[i0:i1 + 1][::-1]]),
                                  np.concatenate([pya[i0:i1 + 1], pyb[i0:i1 + 1][::-1]])])
                 for i0, i1, _ in runs]
        _clip(PolyCollection(polys, facecolors=_rgba(color, [0.30 * r[2] for r in runs]),
                             edgecolors="none", linewidths=0, zorder=2))
    lon_m, pxm, pym, scm = _grad_ring_px(mu_mid, theta_mid, cx, cy, Rx, Ry, inset, 4001)
    ins = _inside(pxm, pym)
    valid = (np.abs(np.diff(lon_m)) <= np.pi) & (np.hypot(np.diff(pxm), np.diff(pym)) <= jump_px)
    valid &= ins[:-1] & ins[1:]
    runs = _runs(0.5 * (scm[:-1] + scm[1:]), valid)
    if runs:
        segs = [np.column_stack([pxm[i0:i1 + 1], pym[i0:i1 + 1]]) for i0, i1, _ in runs]
        _clip(LineCollection(segs, colors=_rgba("white", [0.95 * r[2] for r in runs]),
                             linewidths=2.8 + 1.5, zorder=3, capstyle="butt", joinstyle="round"))
        _clip(LineCollection(segs, colors=_rgba(color, [1.0 * r[2] for r in runs]),
                             linewidths=2.8, zorder=4, capstyle="butt", joinstyle="round"))


def draw_overlay(rgb, reg, seed, source):

    h, w, _ = rgb.shape
    cx, cy, Rx, Ry = reg["cx"], reg["cy"], reg["Rx"], reg["Ry"]
    # Published colored Mollweide (seed detect) — used to keep overlays inside
    # the Migkas map even when registration ellipse is slightly larger.
    scx = float(seed["cx"])
    scy = float(seed["cy"])
    sRx = float(seed["Rx"])
    sRy = float(seed["Ry"])

    dpi = 120
    _extra_h = 1.35 if SHOW_LEGEND else 0.35  # tighten when legend omitted
    fig, ax = plt.subplots(figsize=(w / dpi, (h / dpi) + _extra_h), dpi=dpi)
    ax.imshow(rgb, origin="upper", extent=(0, w, h, 0), interpolation="nearest")
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    # clip_path = registered oval (cx,cy,Rx,Ry). Registration itself is unchanged.
    # Because the refined ellipse can sit a bit outside the published colored
    # map, we ALSO inset strokes and cull/inset markers against the seed oval.
    clip_oval = Ellipse(
        (cx, cy), width=2.0 * Rx, height=2.0 * Ry,
        facecolor="none", edgecolor="none", transform=ax.transData,
    )
    ax.add_patch(clip_oval)

    # Stroke inset so thick white+color underlays do not AA-bleed past the rim.
    STROKE_INSET = 0.97
    # Cull points outside this fraction of the published colored seed oval.
    COLOR_R = 0.978
    color_oval = (scx, scy, sRx, sRy, COLOR_R)

    def stroke_curve(ell, b, color, lw, ls="-"):
        # White underlay then color — stands out on magma without looking oversized.
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=STROKE_INSET, clip_path=clip_oval,
            color_oval=color_oval,
            color="white", lw=lw + 1.6, ls=ls, alpha=0.95,
            solid_capstyle="round", zorder=3,
        )
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=STROKE_INSET, clip_path=clip_oval,
            color_oval=color_oval,
            color=color, lw=lw, ls=ls, alpha=1.0,
            solid_capstyle="round", zorder=4,
        )

    def mark_point(ell, b, marker, color, ms):
        """Place marker; pull limb symbols inward so halo stays in colored oval."""
        px, py = as_xy(*lb_to_px(ell, b, cx, cy, Rx, Ry))
        # Approx marker+halo radius in pixels (ms is points).
        halo_ms = ms + 1.8
        ms_px = 0.5 * halo_ms * (dpi / 72.0)
        # Normalized margin inside the published colored oval.
        margin = ms_px / min(sRx, sRy) + 0.012
        r_lim = max(0.55, COLOR_R - margin)

        Xn = (px - scx) / sRx
        Yn = (py - scy) / sRy
        r_col = float(np.hypot(Xn, Yn))

        if r_col > 1e-9:
            # Scale toward colored-oval center so the whole symbol fits.
            if r_col > r_lim:
                scale = r_lim / r_col
                px = scx + (px - scx) * scale
                py = scy + (py - scy) * scale
                Xn = (px - scx) / sRx
                Yn = (py - scy) / sRy
                r_col = float(np.hypot(Xn, Yn))
            # Mild extra inset for everyone (matches stroke inset spirit).
            mild = min(STROKE_INSET, 0.995)
            px = scx + (px - scx) * mild
            py = scy + (py - scy) * mild
            Xn = (px - scx) / sRx
            Yn = (py - scy) / sRy
            r_col = float(np.hypot(Xn, Yn))

        # Skip if still too close to the rim for this marker size.
        if r_col + margin > COLOR_R + 1e-6:
            return

        def _mark(m_size, m_color, m_edge, m_ew, z):
            (line,) = ax.plot(
                px, py, marker=marker, ms=m_size, color=m_color,
                markeredgecolor=m_edge, markeredgewidth=m_ew, linestyle="None",
                zorder=z, clip_path=clip_oval, clip_on=True,
            )
            line.set_clip_path(clip_oval)
            line.set_clip_on(True)

        _mark(ms + 1.8, "white", "white", 0, 6)
        _mark(ms, color, "k", 0.55, 7)

    # CoE axis at geometric inset=1.0. Do NOT overlay Dec=0 — Migkas already
    # printed the celestial equator on Fig. 8; a second green dashed stroke fights it.
    def stroke_curve_geom(ell, b, color, lw, ls="-"):
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=1.0, clip_path=clip_oval,
            color_oval=color_oval,
            color="white", lw=lw + 1.6, ls=ls, alpha=0.95,
            solid_capstyle="round", zorder=3,
        )
        plot_lb_curve(
            ax, ell, b, cx, cy, Rx, Ry, inset=1.0, clip_path=clip_oval,
            color_oval=color_oval,
            color=color, lw=lw, ls=ls, alpha=1.0,
            solid_capstyle="round", zorder=4,
        )

    # CoE↔anti-CoE axis = great circle through ±D (tilted +30°; not Galactic equator)
    _ax_ell, _ax_b = great_circle_through_poles(D_HAT, n=721)
    stroke_curve_geom(_ax_ell, _ax_b, C_AXIS, lw=1.15)

    # Exactly four highlight bands (match three-panel Panel B): each ±9°
    # (= 18° end-to-end = 5% of 360°) + one solid mid-stroke. No companion.
    for (lo, hi), mid, color in zip(MAX_BANDS, THETA_MAX, COLORS_MAX):
        (fill_highlight_band_grad if GRAD else fill_highlight_band)(
            ax, lo, hi, mid, cx, cy, Rx, Ry, color,
            inset=STROKE_INSET, clip_path=clip_oval, color_oval=color_oval,
        )
    for (lo, hi), mid, color in zip(MIN_BANDS, THETA_MIN, COLORS_MIN):
        (fill_highlight_band_grad if GRAD else fill_highlight_band)(
            ax, lo, hi, mid, cx, cy, Rx, Ry, color,
            inset=STROKE_INSET, clip_path=clip_oval, color_oval=color_oval,
        )

    # Markers (limb ones get stronger inset inside mark_point); per-loop colors
    mark_point(ANTI_COE_LB[0], ANTI_COE_LB[1], "o", C_ANTI, ms=5.5)
    mark_point(COE_LB[0], COE_LB[1], "*", C_COE, ms=7.5)
    for (ell, b), color in zip(MIGKAS_MAXES, COLORS_MAX):
        mark_point(ell, b, "^", color, ms=6.5)
    for (ell, b), color in zip(MIGKAS_MINS, COLORS_MIN):
        mark_point(ell, b, "s", color, ms=6.0)

    handles = [
        Patch(
            facecolor=C_MAX_30, alpha=0.55, edgecolor=C_MAX_30,
            label=(
                rf"MAX $\theta\!=\!{THETA_MAX[0]:g}^\circ$ band "
                rf"[{MAX_BANDS[0][0]:g},{MAX_BANDS[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$=5% of $360^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MAX_185, alpha=0.55, edgecolor=C_MAX_185,
            label=(
                rf"MAX $\theta\!=\!{THETA_MAX[1]:g}^\circ$ band "
                rf"[{MAX_BANDS[1][0]:g},{MAX_BANDS[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_115, alpha=0.55, edgecolor=C_MIN_115,
            label=(
                rf"MIN $\theta\!=\!{THETA_MIN[0]:g}^\circ$ band "
                rf"[{MIN_BANDS[0][0]:g},{MIN_BANDS[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_265, alpha=0.55, edgecolor=C_MIN_265,
            label=(
                rf"MIN $\theta\!=\!{THETA_MIN[1]:g}^\circ$ band "
                rf"[{MIN_BANDS[1][0]:g},{MIN_BANDS[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=$\pm9^\circ$)"
            ),
        ),
        Line2D([0], [0], color=C_AXIS, lw=2.0, label=r"CoE axis tilted $+15^\circ$"),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_30,
            markeredgecolor="k", linestyle="None", ms=8,
            label=rf"MAX $({MIGKAS_MAXES[0][0]:.0f}^\circ,{MIGKAS_MAXES[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_185,
            markeredgecolor="k", linestyle="None", ms=8,
            label=rf"MAX $({MIGKAS_MAXES[1][0]:.0f}^\circ,{MIGKAS_MAXES[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_115,
            markeredgecolor="k", linestyle="None", ms=7,
            label=rf"MIN $({MIGKAS_MINS[0][0]:.0f}^\circ,{MIGKAS_MINS[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_265,
            markeredgecolor="k", linestyle="None", ms=7,
            label=rf"MIN $({MIGKAS_MINS[1][0]:.0f}^\circ,{MIGKAS_MINS[1][1]:+.0f}^\circ)$",
        ),
        Line2D([0], [0], marker="o", color=C_ANTI, markeredgecolor="k",
               linestyle="None", ms=6,
               label=r"anti-CoE $+x$ $(0^\circ,+15^\circ)$"),
        Line2D([0], [0], marker="*", color=C_COE, markeredgecolor="k",
               linestyle="None", ms=9,
               label=r"CoE $(180^\circ,-15^\circ)$"),
    ]
    if SHOW_LEGEND:
        leg = ax.legend(
            handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.012),
            ncol=4, fontsize=7.5, frameon=True, fancybox=False,
            edgecolor="#222222", framealpha=1.0, facecolor="white",
            labelcolor="#111111", borderpad=0.6, handlelength=2.0,
        )
        leg.get_frame().set_linewidth(1.2)
    fig.suptitle(
        r"Migkas Fig. 8 $\theta=75^\circ$ $A/A_{\rm all}$"
        r" + four CoE highlight bands, each $18^\circ$ end-to-end ($\pm9^\circ$ = 5% of $360^\circ$)"
        "\n"
        r"(two MAX + two MIN; cyan/green/magenta/blue; anti-CoE at $(0^\circ,+15^\circ)$; "
        r"CoE axis tilted $+15^\circ$; $\ell$ increases LEFT"
        + (r"; gradient: darkest at $\theta$ ray, transparent opposite" if GRAD else "")
        + (rf"; ${ROLL:+g}^\circ$ CCW about CoE axis" if abs(ROLL) > 1e-9 else "")
        + ")",
        fontsize=10.5, y=0.985,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=(0.13 if SHOW_LEGEND else 0.02))
    fig.savefig(PNG_PATH, dpi=dpi, facecolor="white")
    plt.close(fig)
    print(f"wrote {PNG_PATH}")


def main():
    source = choose_source()
    print(f"source: {source}")
    rgb = np.asarray(Image.open(source).convert("RGB"))
    seed = detect_oval_bbox(rgb)
    print("oval seed:", {k: seed[k] for k in ("cx", "cy", "Rx", "Ry", "x0", "x1", "y0", "y1")})
    # Restore previous visually-good oval (Chris preferred). No medial-axis Rx stretch.
    reg = dict(GOOD_REG)
    h, w, _ = rgb.shape
    dark = rgb.mean(axis=2) < 55
    dist = chamfer_distance_to_dark(dark)
    ell_eq, b_eq = celestial_equator_lb(2500)
    lon = np.deg2rad(lon_disp_from_ell(ell_eq))
    lat = np.deg2rad(b_eq)
    mx, my = mollweide_xy(lon, lat)
    Xn = mx / (2.0 * np.sqrt(2.0))
    Yn = my / np.sqrt(2.0)
    px = reg["cx"] + Xn * reg["Rx"]
    py = reg["cy"] - Yn * reg["Ry"]
    m = (px >= 1) & (px < w - 1) & (py >= 1) & (py < h - 1)
    xi = np.clip(np.round(px[m]).astype(int), 0, w - 1)
    yi = np.clip(np.round(py[m]).astype(int), 0, h - 1)
    reg["ceq_median_px"] = float(np.median(dist[yi, xi]))
    reg["registration_source"] = "hard-set GOOD_REG (pre medial-axis retune)"
    print("restored reg:", {k: reg[k] for k in ("cx", "cy", "Rx", "Ry", "ceq_median_px")})
    checks = validate_landmarks(rgb, reg)
    print(
        "validation OK — max0 RGB", checks["migkas_max0_rgb_mean"],
        "min1 RGB", checks["migkas_min1_rgb_mean"],
    )
    global PNG_PATH, JSON_PATH, GRAD, ROLL
    for name, g, r in VERSIONS:
        PNG_PATH = SM_DIR / f"fig8_overlay_{name}.png"
        JSON_PATH = SM_DIR / f"fig8_overlay_{name}.json"
        GRAD, ROLL = g, r
        draw_overlay(rgb, reg, seed, source)
        img = Image.open(PNG_PATH).convert("RGB")
        img.resize((1700, round(img.height * 1700 / img.width)), Image.LANCZOS).save(
            SM_DIR / f"fig8_overlay_{name}_web.jpg", quality=88)

    meta = {
        "png": str(PNG_PATH),
        "script": str(Path(__file__).resolve()),
        "source_image": str(source),
        "convention": (
            "Galactic Mollweide; lon_disp=-wrap(ell); ell increases LEFT "
            "(Migkas Fig.8)"
        ),
        "axis_tilt_b_deg": AXIS_TILT_B_DEG,
        "anti_CoE_lb_deg": [float(ANTI_COE_LB[0]), float(ANTI_COE_LB[1])],
        "CoE_lb_deg": [float(COE_LB[0]), float(COE_LB[1])],
        "D_hat_galactic_xyz": D_HAT.tolist(),
        "axis_tilt_reason": (
            "Chris 2026-09-22 — markers fixed, rings rotate up 30°"
        ),
        "migkas_max_lb_deg": [list(p) for p in MIGKAS_MAXES],
        "migkas_min_lb_deg": [list(p) for p in MIGKAS_MINS],
        "theta_max_from_plus_x_deg": list(THETA_MAX),
        "theta_min_from_plus_x_deg": list(THETA_MIN),
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "band_width_note": "5% of 360° = 18° end-to-end = ±9° about each MAX/MIN center",
        "max_bands_deg": [list(b) for b in MAX_BANDS],
        "min_bands_deg": [list(b) for b in MIN_BANDS],
        "band_intervals": {
            "MAX_30": [float(MAX_BANDS[0][0]), float(MAX_BANDS[0][1])],
            "MAX_185": [float(MAX_BANDS[1][0]), float(MAX_BANDS[1][1])],
            "MIN_115": [float(MIN_BANDS[0][0]), float(MIN_BANDS[0][1])],
            "MIN_265": [float(MIN_BANDS[1][0]), float(MIN_BANDS[1][1])],
        },
        "colors": {
            "MAX_30": C_MAX_30,
            "MAX_185": C_MAX_185,
            "MIN_115": C_MIN_115,
            "MIN_265": C_MIN_265,
            "palette": "magma-popping cyan/magenta/lime/sky (no reds/oranges)",
        },
        "mu_max_bands": [list(mu_pair_for_band(lo, hi)) for lo, hi in MAX_BANDS],
        "mu_min_bands": [list(mu_pair_for_band(lo, hi)) for lo, hi in MIN_BANDS],
        "scan_stop_source": (
            "best_fit_2ell_multiplier.npz A_after extrema (two max, two min); "
            f"each drawn as one transparent {BAND_FULL:.0f}° end-to-end highlight "
            f"(±{BAND_HALF:.0f}° = 5% of 360°) + one solid mid-stroke; "
            "matches three-panel Panel B; no companion"
        ),
        "matches_three_panel": str(OUT_DIR / "make_coe_rings_migkas_pin_axis.py"),
        "oval_seed_bbox": {
            "x0": seed["x0"], "x1": seed["x1"], "y0": seed["y0"], "y1": seed["y1"],
            "cx": seed["cx"], "cy": seed["cy"], "Rx": seed["Rx"], "Ry": seed["Ry"],
        },
        "registration": {
            "cx": reg["cx"], "cy": reg["cy"], "Rx": reg["Rx"], "Ry": reg["Ry"],
            "ceq_median_px": reg["ceq_median_px"],
            "registration_source": reg.get("registration_source"),
            "ceq_residual_before_registration": dict(GOOD_REG),
            "rejected_medial_axis_reg": {
                "cx": 865.6, "cy": 766.4, "Rx": 706.19, "Ry": 393.70,
                "note": "medial-axis / achromatic DT-ridge retune; blew up Rx; reverted",
            },
            "dofs": ["cx", "cy", "Rx", "Ry"],
            "drawing": {
                "dec0_drawn": False,
                "dec0_note": "not drawn; use Migkas published celestial equator",
                "coe_axis_inset": 1.0,
                "rings_inset": 0.97,
                "note": "CoE axis geometric inset=1; rings mild inset; no Dec=0 overlay",
            },
            "method": (
                "oval seed detect + hard-set to pre-retune visually-good (cx,cy,Rx,Ry); "
                "medial-axis Rx stretch rejected; Dec=0 NOT overlaid (published CEQ kept); "
                "CoE-axis inset=1.0; rings inset~0.97"
            ),
        },
        "landmark_checks": checks,
        "drawing": (
            "Exactly 4 features. Each = ONE smooth transparent PolyCollection "
            "band (θ±9°, 18° end-to-end) + ONE solid mid-stroke at center θ. "
            "Fill alpha≈0.30. Mid-stroke: thin white underlay (lw+1.5) then color "
            "(lw=2.8). NO black outer rim; NO thick white halo stacks; NO dual "
            "boundary-edge strokes. Colors: cyan/magenta/lime/sky for magma. "
            "No 283° companion. Dec=0 NOT overlaid."
        ),
        "legend_summary": (
            "four highlight bands, each 18° end-to-end (±9° = 5% of 360%); "
            "MAX cyan+magenta, MIN lime+sky; intervals [21,39]/[176,194]/[106,124]/[256,274]"
        ),
        "notes": [
            "Exactly 4 features matching three-panel Panel B geometry/width",
            "Each = ONE smooth transparent PolyCollection ribbon (±9°) + ONE solid mid-stroke only",
            "Clean draw: fill alpha=0.30; mid-stroke thin white underlay (lw+1.5) + color lw=2.8; no black rim; no boundary-edge halos",
            "No multi-stroke scribble wash; no 283° companion; no Dec=0 overlay",
            "MAX bands [21,39] & [176,194]; MIN bands [106,124] & [256,274]",
            "Colors: MAX_30=#00e5ff, MAX_185=#ff2d95, MIN_115=#b8f200, MIN_265=#4fc3f7",
            "band_half_width_deg=9, band_full_width_deg=18 (5% of 360°)",
            "Exact n·x̂=μ small-circles; seam-split polygon fill; mid via plot_closed_ring",
            "No ellipse-fit, no interior chords, no decorative oval-rim loops",
            "No blocking text box on the map; legend below",
            "Colorbar preserved from source (overlays stay on oval)",
            "Overlays clip_path = registered Ellipse(cx,cy,Rx,Ry); rings stroke inset 0.97",
            "Dec=0 NOT drawn (Migkas Fig.8 already has celestial equator); CoE axis = great circle through ±D (tilted +30°) at inset=1.0",
            "axis_tilt_b_deg=30; anti-CoE (0°,+30°), CoE (180°,−30°); Migkas markers fixed at true sky positions",
            "Registration hard-set to pre medial-axis params (863,763.25,693.9927,392.542425); DT-ridge fit rejected",
            "Limb markers inset toward colored-oval center so halo stays inside",
            "Points culled at COLOR_R=0.978 of seed oval; limb markers inset to fit halo",
        ],
    }
    pass  # per-version JSON not written (PNG only)
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
```

### 4.2 Fig05 row generators — `make_coe_rings_*` + `make_stack_clean.py`

Exact Fig05 = copy `money_shot_four_row_only.jpg`. Below are the best available row scripts + an older poster stacker.

#### `make_coe_rings_notilt.py`

```python
# source: /workspace/cosmology/figures/make_coe_rings_notilt.py
#!/usr/bin/env python3
"""
CoE rings + Migkas A(ell) extrema loops in the Galactic frame.

Concept:
  Pure CoE axisymmetry → constant A on small-circle rings about x = CoE↔anti-CoE.
  Scan A vs angle from +x (anti-CoE) around 360°. After-2ℓ A(ℓ) has TWO maxima
  and TWO minima (from best_fit_2ell_multiplier.npz A_after):
    MAX at θ = 30° and 185° from +x,
    MIN at θ = 115° and 265° from +x.

  Exactly FOUR loops — each ONE transparent highlight band, 18° end-to-end
  (= ±9° = 5% of 360°) about the center θ:
    MAX: 30±9 → [21°, 39°] and 185±9 → [176°, 194°]
    MIN: 115±9 → [106°, 124°] and 265±9 → [256°, 274°]
  Drawing: ONE smooth transparent filled ribbon between θ−9° and θ+9°
  + ONE solid mid-stroke at the exact center θ. No multi-stroke wash,
  no dense linspace hairball, no dual strong boundaries.
  Four distinct colors (two reds + two blues). No orange, no 283° companion.

  Galactic orientation: anti-CoE (+x) at ℓ=0°, b=+15°; CoE (−x) at ℓ=180°, b=−15°
  (axis tilted +15° so features near ℓ=0 move north on Mollweide; Migkas markers fixed).
  Sky markers (ℓ, b):
    Maxes: map extrema (17°, −9°) for ~30° max; strip (185°, 0°) for second max.
    Mins:  strip (115°, 0°) for ~115° min; map extrema (274°, −22°) for ~265° min.

Outputs:
  coe_rings_migkas_pin_axis.png
  coe_rings_migkas_pin_axis.json
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle, Patch, Polygon
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

OUT_DIR = Path(__file__).resolve().parent
PNG_PATH = Path("/workspace/cosmology/papers/spreadsheet_model/loop_threepanel_notilt.png")
JSON_PATH = Path("/workspace/cosmology/papers/spreadsheet_model/loop_threepanel_notilt.json")

# 5% of 360° = 18° end-to-end → half-width ±9° about each center
BAND_HALF = 9.0
BAND_FULL = 2.0 * BAND_HALF  # 18° = 5% of 360°
BAND_TOTAL = BAND_FULL  # alias
STRIP_THETA_MAX = [29.5, 187.5]
STRIP_THETA_MIN = [122.6, 264.4]

MIGKAS_MAXES = [(17.0, -9.0), (185.0, 0.0)]
MIGKAS_MINS = [(115.0, 0.0), (274.0, -22.0)]

C_COE = "#6b3fa0"
C_ANTI = "#8e1b1b"  # dark red marker; distinct from loop reds
# Four distinct loop colors: two reds (MAX) + two blues (MIN); no orange/yellow
C_MAX_30 = "#0097b2"    # cyan — MAX at 29.5°
C_MAX_185 = "#e0107a"   # magenta — MAX at 187.5°
C_MIN_115 = "#6fa800"   # green — MIN at 122.6°
C_MIN_265 = "#2166d4"   # blue — MIN at 264.4°
C_AXIS = "#5b2c6f"
C_EQ = "#1e8449"

# Per-loop color maps (order matches STRIP_THETA_MAX / STRIP_THETA_MIN)
COLORS_MAX = [C_MAX_30, C_MAX_185]
COLORS_MIN = [C_MIN_115, C_MIN_265]


def wrap_lon_deg(lon: np.ndarray) -> np.ndarray:
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def lon_disp_from_ell(ell_deg: np.ndarray) -> np.ndarray:
    """Migkas Galactic: ℓ increases LEFT → lon_disp = −ℓ_wrapped."""
    return -wrap_lon_deg(ell_deg)


def lb_to_n(l_deg: float, b_deg: float) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    return np.array(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], dtype=float
    )


def n_to_lb(n: np.ndarray) -> tuple[float, float]:
    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    b = float(np.rad2deg(np.arcsin(np.clip(n[2], -1.0, 1.0))))
    l = float(np.rad2deg(np.arctan2(n[1], n[0])) % 360.0)
    return l, b


def celestial_equator_lb(n: int = 720) -> tuple[np.ndarray, np.ndarray]:
    """Dec=0 (ICRS) as Galactic (ℓ, b). IAU J2000 NGP convention, no astropy."""
    ra_ngp = np.deg2rad(192.85948)
    dec_ngp = np.deg2rad(27.12825)
    l_ncp = np.deg2rad(122.93192)
    ra = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    dec = 0.0
    sin_b = np.sin(dec) * np.sin(dec_ngp) + np.cos(dec) * np.cos(dec_ngp) * np.cos(ra - ra_ngp)
    sin_b = np.clip(sin_b, -1.0, 1.0)
    b = np.arcsin(sin_b)
    y = np.cos(dec) * np.sin(ra - ra_ngp)
    x = np.sin(dec) * np.cos(dec_ngp) - np.cos(dec) * np.sin(dec_ngp) * np.cos(ra - ra_ngp)
    l = (l_ncp - np.arctan2(y, x)) % (2.0 * np.pi)
    return np.rad2deg(l), np.rad2deg(b)


def great_circle_through_poles(D: np.ndarray, n: int = 361) -> tuple[np.ndarray, np.ndarray]:
    """Great circle through ±D (any meridian of the CoE axis)."""
    a = D / np.linalg.norm(D)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([1.0, 0.0, 0.0])
    e2 = np.cross(a, tmp)
    e2 /= np.linalg.norm(e2)
    t = np.linspace(0.0, 2.0 * np.pi, n)
    pts = np.outer(np.cos(t), a) + np.outer(np.sin(t), e2)
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    return ls, bs


def small_circle_about_axis(axis: np.ndarray, mu: float, n: int = 361) -> np.ndarray:
    """Unit vectors on the small circle n·axis = mu (|mu|<1)."""
    a = axis / np.linalg.norm(axis)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    tmp = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(a, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    t = np.linspace(0.0, 2.0 * np.pi, n)
    return (mu * a)[None, :] + rho * (np.outer(np.cos(t), e1) + np.outer(np.sin(t), e2))


def split_lon_segments(lon_rad: np.ndarray, lat_rad: np.ndarray, gap: float = np.pi):
    """Split polyline where Mollweide lon jumps (wrap discontinuity)."""
    d = np.diff(lon_rad)
    cuts = np.where(np.abs(d) > gap)[0] + 1
    idx = np.concatenate([[0], cuts, [len(lon_rad)]])
    segs = []
    for i0, i1 in zip(idx[:-1], idx[1:]):
        if i1 - i0 >= 2:
            segs.append((lon_rad[i0:i1], lat_rad[i0:i1]))
    return segs


def plot_lb_curve_mollweide(ax, l_deg, b_deg, **kwargs):
    lon = np.deg2rad(lon_disp_from_ell(l_deg))
    lat = np.deg2rad(b_deg)
    for xr, yr in split_lon_segments(lon, lat):
        ax.plot(xr, yr, **kwargs)


def band_edges(centers: list[float], half: float = BAND_HALF) -> list[tuple[float, float]]:
    """Return (theta_lo, theta_hi) for each band center (±half about center)."""
    return [(c - half, c + half) for c in centers]


def mu_pair_for_band(lo: float, hi: float) -> tuple[float, float]:
    """Edge μ=cos(θ) for a scan-angle band (two boundary cones only)."""
    return float(np.cos(np.deg2rad(lo))), float(np.cos(np.deg2rad(hi)))


def migkas_scan_stops() -> dict:
    th_max = [float(t) for t in STRIP_THETA_MAX]
    th_min = [float(t) for t in STRIP_THETA_MIN]
    max_bands = band_edges(th_max)
    min_bands = band_edges(th_min)
    return {
        "theta_max_from_plus_x_deg": th_max,
        "theta_min_from_plus_x_deg": th_min,
        "max_bands_deg": max_bands,
        "min_bands_deg": min_bands,
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "band_total_width_deg": BAND_TOTAL,
        "band_width_note": "5% of 360° = 18° end-to-end = ±9° about each MAX/MIN center",
        "mu_max_centers": [float(np.cos(np.deg2rad(t))) for t in th_max],
        "mu_min_centers": [float(np.cos(np.deg2rad(t))) for t in th_min],
        "mu_max_bands": [mu_pair_for_band(lo, hi) for lo, hi in max_bands],
        "mu_min_bands": [mu_pair_for_band(lo, hi) for lo, hi in min_bands],
        "source": (
            "best_fit_2ell_multiplier.npz A_after extrema (two max, two min); "
            f"each drawn as one transparent {BAND_FULL:.0f}° end-to-end highlight "
            f"(±{BAND_HALF:.0f}° = 5% of 360°) + one solid mid-stroke"
        ),
    }


def compute_geometry(st: dict) -> dict:
    anti_l, anti_b = 0.0, 0.0
    coe_l, coe_b = 180.0, 0.0

    n_coe = lb_to_n(coe_l, coe_b)
    n_anti = lb_to_n(anti_l, anti_b)
    D = n_anti  # +D = anti-CoE = +x

    th_max = [float(t) for t in st["theta_max_from_plus_x_deg"]]
    th_min = [float(t) for t in st["theta_min_from_plus_x_deg"]]
    max_bands = list(st["max_bands_deg"])
    min_bands = list(st["min_bands_deg"])

    n_maxes = [lb_to_n(*lb) for lb in MIGKAS_MAXES]
    n_mins = [lb_to_n(*lb) for lb in MIGKAS_MINS]

    sep = float(
        np.rad2deg(
            np.arccos(np.clip(np.dot(n_maxes[0], n_mins[1]), -1.0, 1.0))
        )
    )

    return {
        "D": D,
        "n_coe": n_coe,
        "n_anti": n_anti,
        "n_maxes": n_maxes,
        "n_mins": n_mins,
        "anti_lb": (anti_l, anti_b),
        "coe_lb": (coe_l, coe_b),
        "theta_max_from_plus_x_deg": th_max,
        "theta_min_from_plus_x_deg": th_min,
        "max_bands_deg": max_bands,
        "min_bands_deg": min_bands,
        "mu_max_bands": list(st["mu_max_bands"]),
        "mu_min_bands": list(st["mu_min_bands"]),
        "sep_deg": sep,
        "band_half_width_deg": st["band_half_width_deg"],
    }


def draw_wire_sphere(ax, n_mer=12, n_par=6, **kwargs):
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(-np.pi / 2, np.pi / 2, 40)
    xs = np.outer(np.cos(u), np.cos(v))
    ys = np.outer(np.sin(u), np.cos(v))
    zs = np.outer(np.ones_like(u), np.sin(v))
    ax.plot_surface(xs, ys, zs, color="#dfe6e9", alpha=0.15, linewidth=0, shade=False)
    kw = dict(color="#95a5a6", lw=0.4, alpha=0.45)
    kw.update(kwargs)
    for i in range(n_mer):
        phi = 2 * np.pi * i / n_mer
        t = np.linspace(-np.pi / 2, np.pi / 2, 80)
        ax.plot(np.cos(t) * np.cos(phi), np.cos(t) * np.sin(phi), np.sin(t), **kw)
    for j in range(1, n_par):
        b = -np.pi / 2 + np.pi * j / n_par
        t = np.linspace(0, 2 * np.pi, 80)
        r = np.cos(b)
        ax.plot(r * np.cos(t), r * np.sin(t), np.full_like(t, np.sin(b)), **kw)


def foreshortened_ellipse_xy(mu: float, n: int = 200):
    """Side-view foreshortened ellipse of a 3D ring at n·x = mu."""
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    phi = np.linspace(0, 2 * np.pi, n)
    ex = mu + 0.045 * rho * np.cos(phi)
    ey = rho * np.sin(phi)
    return ex, ey, rho


def draw_panelA_ribbon(ax, theta_lo: float, theta_hi: float, theta_mid: float, color: str):
    """
    One Panel-A highlight: clean vertical strip between x=cos(θ_hi)..cos(θ_lo)
    clipped to the unit circle (uniform alpha), plus ONE solid foreshortened
    mid-ellipse at the exact center θ. No multi-stroke wash / hairball.
    """
    mu_lo = float(np.cos(np.deg2rad(theta_lo)))
    mu_hi = float(np.cos(np.deg2rad(theta_hi)))
    x0, x1 = (min(mu_lo, mu_hi), max(mu_lo, mu_hi))
    # Vertical strip polygon (tall enough to cover the circle; clipped below)
    strip = Polygon(
        np.array([[x0, -1.15], [x1, -1.15], [x1, 1.15], [x0, 1.15]]),
        closed=True,
        facecolor=color,
        edgecolor="none",
        alpha=0.28,
        zorder=2,
        linewidth=0,
    )
    ax.add_patch(strip)
    clip_circ = Circle((0.0, 0.0), 1.0, transform=ax.transData)
    strip.set_clip_path(clip_circ)

    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    ex_m, ey_m, _ = foreshortened_ellipse_xy(mu_mid, n=256)
    ax.plot(ex_m, ey_m, color=color, lw=2.8, alpha=0.95, zorder=4, solid_capstyle="round")


def draw_panelA_theta_rays(ax, geo):
    """
    Illustrate θ: from the observer (origin), one ray per loop out to the
    ring at that loop's center θ. θ is measured from +x (anti-CoE, right)
    up to the ring. Emphasize MAX θ=30° (thick semi-transparent red to the
    top of that foreshortened ring), matching Chris's markup.
    """
    items = list(zip(geo["theta_max_from_plus_x_deg"], COLORS_MAX)) + list(
        zip(geo["theta_min_from_plus_x_deg"], COLORS_MIN)
    )
    for mid, color in items:
        t = float(np.deg2rad(mid))
        # Endpoint = top/rim point of the ring in side view = (cos θ, sin θ)
        x1, y1 = float(np.cos(t)), float(np.sin(t))
        is_max30 = abs(mid - 29.5) < 0.5
        lw = 3.2 if is_max30 else 1.8
        alpha = 0.70 if is_max30 else 0.45
        ax.plot(
            [0.0, x1], [0.0, y1],
            color=color, lw=lw, alpha=alpha, zorder=3,
            solid_capstyle="round",
        )
        # Small tick at the ring hit
        ax.plot([x1], [y1], "o", color=color, ms=4.5 if is_max30 else 3.5,
                alpha=0.9, zorder=5, markeredgecolor="none")

    # Angle arc at θ=30° from +x axis up to the ray (makes the measure obvious)
    arc_r = 0.38
    phi = np.linspace(0.0, np.deg2rad(29.5), 60)
    ax.plot(
        arc_r * np.cos(phi), arc_r * np.sin(phi),
        color=C_MAX_30, lw=1.6, alpha=0.85, zorder=3,
    )
    ax.text(
        arc_r * np.cos(np.deg2rad(15.0)) + 0.06,
        arc_r * np.sin(np.deg2rad(15.0)) + 0.04,
        r"$\theta$",
        color=C_MAX_30, fontsize=10, fontweight="bold",
        ha="left", va="bottom", zorder=6,
    )


def panel_A_model_frame(ax, geo):
    """Side view: exactly four thick foreshortened loop ribbons about x."""
    ax.set_title(
        "A. Model frame (no Galactic yet)\n"
        r"axisymmetric about $x$ = CoE$\leftrightarrow$anti-CoE"
        "\n"
        r"two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$)",
        fontsize=10,
        pad=8,
    )
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color="#7f8c8d", lw=1.2, zorder=1)
    ax.plot(0, 0, "ko", ms=6, zorder=5)
    ax.text(0.08, -0.12, "observer", fontsize=8, color="k", ha="left")

    ax.annotate(
        "",
        xy=(1.15, 0),
        xytext=(-1.15, 0),
        arrowprops=dict(arrowstyle="<->", color=C_AXIS, lw=2.5),
        zorder=4,
    )
    ax.text(1.18, 0.08, "anti-CoE\n(+x)", color=C_ANTI, fontsize=8, ha="left", va="bottom")
    ax.text(-1.18, 0.08, "CoE\n(−x)", color=C_COE, fontsize=8, ha="right", va="bottom")

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        draw_panelA_ribbon(ax, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        draw_panelA_ribbon(ax, lo, hi, mid, color)

    # Illustrate θ from +x (anti-CoE) to each loop (Chris markup)
    draw_panelA_theta_rays(ax, geo)

    label_r = 1.30
    for mid, color in list(zip(geo["theta_max_from_plus_x_deg"], COLORS_MAX)) + list(
        zip(geo["theta_min_from_plus_x_deg"], COLORS_MIN)
    ):
        t = np.deg2rad(mid)
        ax.text(
            label_r * np.cos(t),
            label_r * np.sin(t),
            rf"$\theta\!=\!{mid:g}^\circ$",
            color=color, fontsize=7.5, ha="center", va="center", fontweight="bold",
        )

    ax.text(
        0.0,
        -1.38,
        (
            r"four transparent highlight bands about $x$  ·  each $18^\circ$ end-to-end ($\pm9^\circ$ = 5% of $360^\circ$)"
            + "\n"
            + r"(after-2$\ell$ $A(\ell)$: two MAX + two MIN; one mid-stroke each)"
            + "\n"
            + r"$y,z$ free (roll) — model does not pick celestial equator"
        ),
        fontsize=7.5,
        ha="center",
        va="top",
        style="italic",
        color="#34495e",
    )
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")


def _ring_lb(D, mu, n_pts=721):
    pts = small_circle_about_axis(D, mu, n=n_pts)
    rl = np.array([n_to_lb(p)[0] for p in pts])
    rb = np.array([n_to_lb(p)[1] for p in pts])
    return rl, rb


def _ring_lonlat_rad(D, mu, n_pts=361):
    """Small-circle → Mollweide display (lon_disp, lat) in radians."""
    pts = small_circle_about_axis(D, float(mu), n=n_pts)
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    lon = np.deg2rad(lon_disp_from_ell(ls))
    lat = np.deg2rad(bs)
    return lon, lat


def fill_mollweide_ribbon(ax, D, theta_lo, theta_hi, theta_mid, color, n_pts=721):
    """
    One Mollweide highlight band:
      continuous transparent polygon strips between ±9° edges
      (broken only at Mollweide antimeridian seams — NOT a dense quad wash)
      + ONE solid mid-stroke at θ_mid.
    """
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    lon_a, lat_a = _ring_lonlat_rad(D, mu_a, n_pts=n_pts)
    lon_b, lat_b = _ring_lonlat_rad(D, mu_b, n_pts=n_pts)
    gap = 0.9 * np.pi
    # Mark steps that cross a seam on either edge (or between edges)
    bad = np.zeros(n_pts - 1, dtype=bool)
    for i in range(n_pts - 1):
        if abs(lon_a[i + 1] - lon_a[i]) > gap:
            bad[i] = True
        elif abs(lon_b[i + 1] - lon_b[i]) > gap:
            bad[i] = True
        elif abs(lon_a[i] - lon_b[i]) > gap:
            bad[i] = True
        elif abs(lon_a[i + 1] - lon_b[i + 1]) > gap:
            bad[i] = True
    # Group into continuous good runs → one polygon per run
    polys = []
    i = 0
    while i < n_pts - 1:
        while i < n_pts - 1 and bad[i]:
            i += 1
        if i >= n_pts - 1:
            break
        j = i
        while j < n_pts - 1 and not bad[j]:
            j += 1
        # indices i .. j inclusive are a continuous ribbon segment
        i0, i1 = i, j  # i1 is last good vertex index
        if i1 - i0 >= 1:
            xs = np.concatenate([lon_a[i0:i1 + 1], lon_b[i0:i1 + 1][::-1]])
            ys = np.concatenate([lat_a[i0:i1 + 1], lat_b[i0:i1 + 1][::-1]])
            polys.append(np.column_stack([xs, ys]))
        i = j + 1
    if polys:
        ax.add_collection(
            PolyCollection(
                polys, facecolors=color, edgecolors="none",
                alpha=0.26, zorder=2, linewidths=0,
            )
        )
    rl, rb = _ring_lb(D, mu_mid, n_pts=max(n_pts, 721))
    plot_lb_curve_mollweide(
        ax, rl, rb, color=color, lw=2.6, alpha=0.95, zorder=4
    )


def panel_B_mollweide(ax, geo):
    ax.set_title(
        r"B. Migkas $A(\ell)$ loops on sky (anti-CoE at $(0^\circ,0^\circ)$)"
        "\n"
        r"(Galactic Mollweide; CoE axis untilted; $\ell$ increases LEFT; "
        r"two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$))",
        fontsize=10,
        pad=8,
    )
    for b in [-60, -30, 0, 30, 60]:
        ell = np.linspace(0, 360, 361)
        plot_lb_curve_mollweide(ax, ell, np.full_like(ell, b), color="#bdc3c7", lw=0.5, zorder=0)
    for ell0 in range(0, 360, 30):
        b = np.linspace(-90, 90, 181)
        plot_lb_curve_mollweide(ax, np.full_like(b, ell0), b, color="#bdc3c7", lw=0.5, zorder=0)

    D = geo["D"]
    ls, bs = great_circle_through_poles(D, n=721)
    plot_lb_curve_mollweide(ax, ls, bs, color=C_AXIS, lw=2.0, zorder=2)

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        fill_mollweide_ribbon(ax, D, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        fill_mollweide_ribbon(ax, D, lo, hi, mid, color)

    ell_eq, b_eq = celestial_equator_lb(721)
    plot_lb_curve_mollweide(ax, ell_eq, b_eq, color=C_EQ, lw=1.8, ls="--", zorder=2)

    ax.plot(
        np.deg2rad(lon_disp_from_ell(geo["anti_lb"][0])),
        np.deg2rad(geo["anti_lb"][1]),
        marker="o", ms=8, color=C_ANTI,
        markerfacecolor=C_ANTI, markeredgecolor="k", markeredgewidth=0.6,
        zorder=6, linestyle="None",
    )
    ax.plot(
        np.deg2rad(lon_disp_from_ell(geo["coe_lb"][0])),
        np.deg2rad(geo["coe_lb"][1]),
        marker="*", ms=14, color=C_COE,
        markerfacecolor=C_COE, markeredgecolor="k", markeredgewidth=0.6,
        zorder=6, linestyle="None",
    )
    for (l, b), color in zip(MIGKAS_MAXES, COLORS_MAX):
        ax.plot(
            np.deg2rad(lon_disp_from_ell(l)),
            np.deg2rad(b),
            marker="^",
            ms=10,
            color=color,
            markeredgecolor="k",
            markeredgewidth=0.7,
            zorder=7,
            linestyle="None",
        )
    for (l, b), color in zip(MIGKAS_MINS, COLORS_MIN):
        ax.plot(
            np.deg2rad(lon_disp_from_ell(l)),
            np.deg2rad(b),
            marker="s",
            ms=10,
            color=color,
            markeredgecolor="k",
            markeredgewidth=0.7,
            zorder=7,
            linestyle="None",
        )

    xticks = np.deg2rad(lon_disp_from_ell(np.array([180, 90, 0, 270])))
    ax.set_xticks(xticks)
    ax.set_xticklabels([r"$180^\circ$", r"$90^\circ$", r"$0^\circ$", r"$270^\circ$"], fontsize=8)
    ax.set_yticks(np.deg2rad([-60, -30, 0, 30, 60]))
    ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$+30^\circ$", r"$+60^\circ$"], fontsize=8)
    ax.grid(False)
    ax.set_xlabel(r"Galactic $\ell$ (increases LEFT)", fontsize=8)
    ax.set_ylabel(r"Galactic $b$", fontsize=8)


def fill_sphere_ribbon(ax, D, theta_lo, theta_hi, theta_mid, color, n_pts=160):
    """
    One 3D highlight belt: SINGLE transparent strip between ±9° edges
    + ONE solid mid-stroke. No linspace wash / stacked alpha hairball.
    """
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    p0 = small_circle_about_axis(D, float(mu_a), n=n_pts)
    p1 = small_circle_about_axis(D, float(mu_b), n=n_pts)
    # Single belt between the two edge rings (one collection, no stacked wash)
    verts = [[p0[j], p0[j + 1], p1[j + 1], p1[j]] for j in range(n_pts - 1)]
    coll = Poly3DCollection(
        verts, facecolors=color, edgecolors="none", linewidths=0.0,
        alpha=0.25, zorder=2,
    )
    ax.add_collection3d(coll)
    pts = small_circle_about_axis(D, mu_mid, n=n_pts)
    ax.plot(
        pts[:, 0], pts[:, 1], pts[:, 2],
        color=color, lw=2.6, alpha=0.95, zorder=5,
    )


def panel_C_sphere(ax, geo):
    th_max = geo["theta_max_from_plus_x_deg"]
    th_min = geo["theta_min_from_plus_x_deg"]
    ax.set_title(
        r"C. Sphere 3D — two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$)"
        "\n"
        r"(axis line emphasized; celestial equator dashed)",
        fontsize=10,
        pad=6,
    )
    draw_wire_sphere(ax)
    D = geo["D"]
    ell_eq, b_eq = celestial_equator_lb(400)
    eq_xyz = np.array([lb_to_n(l, b) for l, b in zip(ell_eq, b_eq)])
    ax.plot(eq_xyz[:, 0], eq_xyz[:, 1], eq_xyz[:, 2], color=C_EQ, lw=1.8, ls="--", zorder=3)
    ax.plot(
        [-D[0], D[0]],
        [-D[1], D[1]],
        [-D[2], D[2]],
        color=C_AXIS,
        lw=3.0,
        zorder=5,
    )
    ax.scatter(*(-D), color=C_COE, s=120, marker="*", edgecolors="k", linewidths=0.6,
               depthshade=False, zorder=8)
    ax.scatter(*(D), color=C_ANTI, s=80, edgecolors="k", linewidths=0.6, depthshade=False, zorder=8)
    ax.text(*(-D * 1.22), "CoE", color=C_COE, fontsize=8, fontweight="bold")
    ax.text(*(D * 1.18), "anti-CoE", color=C_ANTI, fontsize=8, fontweight="bold")
    ax.scatter([0], [0], [0], color="k", s=25, depthshade=False, zorder=9)

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        fill_sphere_ribbon(ax, D, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        fill_sphere_ribbon(ax, D, lo, hi, mid, color)

    for n, color in zip(geo["n_maxes"], COLORS_MAX):
        ax.scatter(
            [n[0]], [n[1]], [n[2]],
            c=[color], marker="^", s=70,
            edgecolors="k", linewidths=0.7, depthshade=False, zorder=10,
        )
        ax.plot([0, n[0]], [0, n[1]], [0, n[2]], color=color, lw=0.8, ls=":", alpha=0.7)
    for n, color in zip(geo["n_mins"], COLORS_MIN):
        ax.scatter(
            [n[0]], [n[1]], [n[2]],
            c=[color], marker="s", s=60,
            edgecolors="k", linewidths=0.7, depthshade=False, zorder=10,
        )
        ax.plot([0, n[0]], [0, n[1]], [0, n[2]], color=color, lw=0.8, ls=":", alpha=0.7)

    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.set_zlim(-1.25, 1.25)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=18, azim=-55)
    ax.set_axis_off()
    ax.text2D(
        0.5,
        -0.02,
        rf"MAX $\theta\!=\!{th_max[0]:g}^\circ,{th_max[1]:g}^\circ$; "
        rf"MIN $\theta\!=\!{th_min[0]:g}^\circ,{th_min[1]:g}^\circ$  "
        rf"($18^\circ$ end-to-end = $\pm9^\circ$ = 5% of $360^\circ$)",
        transform=ax.transAxes,
        ha="center",
        fontsize=7.5,
        color="#34495e",
    )


def main():
    st = migkas_scan_stops()
    geo = compute_geometry(st)

    th_max = geo["theta_max_from_plus_x_deg"]
    th_min = geo["theta_min_from_plus_x_deg"]
    max_bands = geo["max_bands_deg"]
    min_bands = geo["min_bands_deg"]

    fig = plt.figure(figsize=(14.5, 5.6), dpi=160)
    fig.suptitle(
        r"Pure CoE axisymmetry: two MAX + two MIN highlight bands"
        "\n"
        rf"(each $18^\circ$ end-to-end = $\pm9^\circ$ = 5% of $360^\circ$; "
        rf"MAX $\theta\!=\!{th_max[0]:g}^\circ/{th_max[1]:g}^\circ$; "
        rf"MIN $\theta\!=\!{th_min[0]:g}^\circ/{th_min[1]:g}^\circ$; "
        r"anti-CoE $(0^\circ,0^\circ)$, CoE axis untilted, no roll)",
        fontsize=11,
        fontweight="bold",
        y=1.02,
    )

    axA = fig.add_subplot(1, 3, 1)
    panel_A_model_frame(axA, geo)

    axB = fig.add_subplot(1, 3, 2, projection="mollweide")
    panel_B_mollweide(axB, geo)

    axC = fig.add_subplot(1, 3, 3, projection="3d")
    panel_C_sphere(axC, geo)

    handles = [
        Line2D([0], [0], color=C_AXIS, lw=2.5, label=r"CoE$\leftrightarrow$anti-CoE axis"),
        Patch(
            facecolor=C_MAX_30, alpha=0.45, edgecolor=C_MAX_30,
            label=(
                rf"MAX $\theta\!=\!{th_max[0]:g}^\circ$ band "
                rf"[{max_bands[0][0]:g},{max_bands[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$=5% of $360^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MAX_185, alpha=0.45, edgecolor=C_MAX_185,
            label=(
                rf"MAX $\theta\!=\!{th_max[1]:g}^\circ$ band "
                rf"[{max_bands[1][0]:g},{max_bands[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_115, alpha=0.45, edgecolor=C_MIN_115,
            label=(
                rf"MIN $\theta\!=\!{th_min[0]:g}^\circ$ band "
                rf"[{min_bands[0][0]:g},{min_bands[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_265, alpha=0.45, edgecolor=C_MIN_265,
            label=(
                rf"MIN $\theta\!=\!{th_min[1]:g}^\circ$ band "
                rf"[{min_bands[1][0]:g},{min_bands[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Line2D([0], [0], color=C_EQ, lw=1.8, ls="--", label="celestial equator (Dec$=0$)"),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_30, markeredgecolor="k",
            ms=9,
            label=rf"MAX $({MIGKAS_MAXES[0][0]:.0f}^\circ,{MIGKAS_MAXES[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_185, markeredgecolor="k",
            ms=9,
            label=rf"MAX $({MIGKAS_MAXES[1][0]:.0f}^\circ,{MIGKAS_MAXES[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_115, markeredgecolor="k",
            ms=8,
            label=rf"MIN $({MIGKAS_MINS[0][0]:.0f}^\circ,{MIGKAS_MINS[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_265, markeredgecolor="k",
            ms=8,
            label=rf"MIN $({MIGKAS_MINS[1][0]:.0f}^\circ,{MIGKAS_MINS[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="o", color="w", markerfacecolor=C_ANTI, markeredgecolor="k",
            ms=8, label=rf"anti-CoE $({geo['anti_lb'][0]:.1f}^\circ,{geo['anti_lb'][1]:+.1f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="*", color="w", markerfacecolor=C_COE, markeredgecolor="k",
            ms=11, label=rf"CoE $({geo['coe_lb'][0]:.1f}^\circ,{geo['coe_lb'][1]:+.1f}^\circ)$",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=6,
        fontsize=7.0,
        frameon=True,
        fancybox=True,
        bbox_to_anchor=(0.5, -0.10),
    )

    fig.tight_layout()
    # --- Chris 2026-09-26: panel C sphere drawn at the same on-page diameter as
    # panel A's model-frame circle; panel C axes matched to A's vertical
    # position/height. Nothing else changes.
    from mpl_toolkits.mplot3d import proj3d
    pa, pc = axA.get_position(), axC.get_position()
    axC.set_position([pc.x0, pa.y0, pc.width, pa.height])

    def _diam_A():
        fig.canvas.draw()
        p0 = axA.transData.transform((-1.0, 0.0)); p1 = axA.transData.transform((1.0, 0.0))
        return float(p1[0] - p0[0])

    def _diam_C():
        fig.canvas.draw()
        u = np.random.default_rng(0).normal(size=(20000, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]
        xs, ys, _ = proj3d.proj_transform(u[:, 0], u[:, 1], u[:, 2], axC.get_proj())
        d = axC.transData.transform(np.column_stack([xs, ys]))
        return float(d[:, 0].max() - d[:, 0].min()), float(d[:, 1].max() - d[:, 1].min())

    zoom = 1.0
    for _ in range(8):
        dA = _diam_A(); dCx, dCy = _diam_C()
        zoom *= dA / dCx
        axC.set_box_aspect((1, 1, 1), zoom=zoom)
    dA = _diam_A(); dCx, dCy = _diam_C()
    print(f"PANEL_DIAM A={dA:.1f}px C={dCx:.1f}x{dCy:.1f}px zoom={zoom:.3f}")
    fig.savefig(PNG_PATH, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    payload = {
        "png": str(PNG_PATH),
        "script": str(OUT_DIR / "make_coe_rings_notilt.py"),
        "migkas_max_lb_deg": [list(lb) for lb in MIGKAS_MAXES],
        "migkas_min_lb_deg": [list(lb) for lb in MIGKAS_MINS],
        "strip_theta_max_deg": list(STRIP_THETA_MAX),
        "strip_theta_min_deg": list(STRIP_THETA_MIN),
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "band_total_width_deg": BAND_TOTAL,
        "band_width_note": "5% of 360° = 18° end-to-end = ±9° about each MAX/MIN center",
        "max_bands_deg": [list(b) for b in max_bands],
        "min_bands_deg": [list(b) for b in min_bands],
        "band_intervals": {
            "MAX_30": [float(max_bands[0][0]), float(max_bands[0][1])],
            "MAX_185": [float(max_bands[1][0]), float(max_bands[1][1])],
            "MIN_115": [float(min_bands[0][0]), float(min_bands[0][1])],
            "MIN_265": [float(min_bands[1][0]), float(min_bands[1][1])],
        },
        "scan_stop_source": st["source"],
        "theta_max_from_plus_x_deg": [round(t, 4) for t in th_max],
        "theta_min_from_plus_x_deg": [round(t, 4) for t in th_min],
        "mu_max_bands": [[round(a, 8), round(b, 8)] for a, b in geo["mu_max_bands"]],
        "mu_min_bands": [[round(a, 8), round(b, 8)] for a, b in geo["mu_min_bands"]],
        "CoE_lb_deg": [round(geo["coe_lb"][0], 4), round(geo["coe_lb"][1], 4)],
        "anti_CoE_lb_deg": [round(geo["anti_lb"][0], 4), round(geo["anti_lb"][1], 4)],
        "D_hat_galactic_xyz": geo["D"].tolist(),
        "angular_sep_migkas_max0_min1_deg": round(geo["sep_deg"], 6),
        "axis_tilt_b_deg": 15,
        "axis_tilt_reason": (
            "Chris 2026-09-22 — markers fixed, rings rotate up 30°"
        ),
        "orientation": (
            "Galactic orientation: anti-CoE at (ell=0°, b=+15°), CoE at "
            "(ell=180°, b=−15°); +D = anti-CoE = +x; axis tilted +15° "
            "(Chris 2026-09-22 — markers fixed, rings rotate up 30°)"
        ),
        "drawing": (
            "Exactly 4 features. Each = ONE smooth transparent highlight band "
            "between θ−9° and θ+9° (18° end-to-end = 5% of 360°) + ONE solid "
            "mid-stroke at the exact center θ. No multi-stroke wash / hairball. "
            "Four distinct colors (two reds + two blues); no orange; no 283° "
            "companion. Panel A: clean vertical strip clipped to unit circle; theta rays from observer to each loop (emphasized at MAX 30°). Mollweide "
            "seam-split honestly. Sphere: single Poly3D strip (no n_wash stack)."
        ),
        "note": (
            "Stops are after-2ℓ A(ℓ) extrema from best_fit_2ell_multiplier.npz "
            f"A_after: two MAX at θ={th_max[0]:g}°/{th_max[1]:g}° and "
            f"two MIN at θ={th_min[0]:g}°/{th_min[1]:g}°, each drawn as "
            f"one {BAND_FULL:.0f}° end-to-end highlight (±{BAND_HALF:.0f}° = 5% of 360°). "
            "Rings use μ=cos(θ); sky markers are map extrema (17°,−9°) and "
            "(274°,−22°) plus strip stand-ins (185°,0°) and (115°,0°)."
        ),
        "plot_convention": "Migkas Galactic Mollweide: lon_disp=-ell_wrapped (ell increases LEFT)",
        "legend_summary": (
            "two MAX + two MIN highlight bands, each 18° end-to-end "
            "(±9° = 5% of 360°); four distinct colors"
        ),
        "colors": {
            "MAX_30": C_MAX_30,
            "MAX_185": C_MAX_185,
            "MIN_115": C_MIN_115,
            "MIN_265": C_MIN_265,
            "no_orange": True,
        },
    }
    JSON_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    print(f"\nWrote {PNG_PATH}")
    print(f"Wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
```

#### `make_coe_rings_migkas_pin_axis_matched.py`

```python
# source: /workspace/cosmology/figures/make_coe_rings_migkas_pin_axis_matched.py
#!/usr/bin/env python3
"""
CoE rings + Migkas A(ell) extrema loops in the Galactic frame.

Concept:
  Pure CoE axisymmetry → constant A on small-circle rings about x = CoE↔anti-CoE.
  Scan A vs angle from +x (anti-CoE) around 360°. After-2ℓ A(ℓ) has TWO maxima
  and TWO minima (from best_fit_2ell_multiplier.npz A_after):
    MAX at θ = 30° and 185° from +x,
    MIN at θ = 115° and 265° from +x.

  Exactly FOUR loops — each ONE transparent highlight band, 18° end-to-end
  (= ±9° = 5% of 360°) about the center θ:
    MAX: 30±9 → [21°, 39°] and 185±9 → [176°, 194°]
    MIN: 115±9 → [106°, 124°] and 265±9 → [256°, 274°]
  Drawing: ONE smooth transparent filled ribbon between θ−9° and θ+9°
  + ONE solid mid-stroke at the exact center θ. No multi-stroke wash,
  no dense linspace hairball, no dual strong boundaries.
  Four distinct colors (two reds + two blues). No orange, no 283° companion.

  Galactic orientation: anti-CoE (+x) at ℓ=0°, b=+15°; CoE (−x) at ℓ=180°, b=−15°
  (axis tilted +15° so features near ℓ=0 move north on Mollweide; Migkas markers fixed).
  Sky markers (ℓ, b):
    Maxes: map extrema (17°, −9°) for ~30° max; strip (185°, 0°) for second max.
    Mins:  strip (115°, 0°) for ~115° min; map extrema (274°, −22°) for ~265° min.

Outputs:
  coe_rings_migkas_pin_axis.png
  coe_rings_migkas_pin_axis.json
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle, Patch, Polygon
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

OUT_DIR = Path(__file__).resolve().parent
PNG_PATH = Path("/workspace/cosmology/papers/spreadsheet_model/loop_threepanel_matched.png")
JSON_PATH = Path("/workspace/cosmology/papers/spreadsheet_model/loop_threepanel_matched.json")

# 5% of 360° = 18° end-to-end → half-width ±9° about each center
BAND_HALF = 9.0
BAND_FULL = 2.0 * BAND_HALF  # 18° = 5% of 360°
BAND_TOTAL = BAND_FULL  # alias
STRIP_THETA_MAX = [29.5, 187.5]
STRIP_THETA_MIN = [122.6, 264.4]

MIGKAS_MAXES = [(17.0, -9.0), (185.0, 0.0)]
MIGKAS_MINS = [(115.0, 0.0), (274.0, -22.0)]

C_COE = "#6b3fa0"
C_ANTI = "#8e1b1b"  # dark red marker; distinct from loop reds
# Four distinct loop colors: two reds (MAX) + two blues (MIN); no orange/yellow
C_MAX_30 = "#0097b2"    # cyan — MAX at 29.5°
C_MAX_185 = "#e0107a"   # magenta — MAX at 187.5°
C_MIN_115 = "#6fa800"   # green — MIN at 122.6°
C_MIN_265 = "#2166d4"   # blue — MIN at 264.4°
C_AXIS = "#5b2c6f"
C_EQ = "#1e8449"

# Per-loop color maps (order matches STRIP_THETA_MAX / STRIP_THETA_MIN)
COLORS_MAX = [C_MAX_30, C_MAX_185]
COLORS_MIN = [C_MIN_115, C_MIN_265]


def wrap_lon_deg(lon: np.ndarray) -> np.ndarray:
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def lon_disp_from_ell(ell_deg: np.ndarray) -> np.ndarray:
    """Migkas Galactic: ℓ increases LEFT → lon_disp = −ℓ_wrapped."""
    return -wrap_lon_deg(ell_deg)


def lb_to_n(l_deg: float, b_deg: float) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    return np.array(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], dtype=float
    )


def n_to_lb(n: np.ndarray) -> tuple[float, float]:
    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    b = float(np.rad2deg(np.arcsin(np.clip(n[2], -1.0, 1.0))))
    l = float(np.rad2deg(np.arctan2(n[1], n[0])) % 360.0)
    return l, b


def celestial_equator_lb(n: int = 720) -> tuple[np.ndarray, np.ndarray]:
    """Dec=0 (ICRS) as Galactic (ℓ, b). IAU J2000 NGP convention, no astropy."""
    ra_ngp = np.deg2rad(192.85948)
    dec_ngp = np.deg2rad(27.12825)
    l_ncp = np.deg2rad(122.93192)
    ra = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    dec = 0.0
    sin_b = np.sin(dec) * np.sin(dec_ngp) + np.cos(dec) * np.cos(dec_ngp) * np.cos(ra - ra_ngp)
    sin_b = np.clip(sin_b, -1.0, 1.0)
    b = np.arcsin(sin_b)
    y = np.cos(dec) * np.sin(ra - ra_ngp)
    x = np.sin(dec) * np.cos(dec_ngp) - np.cos(dec) * np.sin(dec_ngp) * np.cos(ra - ra_ngp)
    l = (l_ncp - np.arctan2(y, x)) % (2.0 * np.pi)
    return np.rad2deg(l), np.rad2deg(b)


def great_circle_through_poles(D: np.ndarray, n: int = 361) -> tuple[np.ndarray, np.ndarray]:
    """Great circle through ±D (any meridian of the CoE axis)."""
    a = D / np.linalg.norm(D)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([1.0, 0.0, 0.0])
    e2 = np.cross(a, tmp)
    e2 /= np.linalg.norm(e2)
    t = np.linspace(0.0, 2.0 * np.pi, n)
    pts = np.outer(np.cos(t), a) + np.outer(np.sin(t), e2)
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    return ls, bs


def small_circle_about_axis(axis: np.ndarray, mu: float, n: int = 361) -> np.ndarray:
    """Unit vectors on the small circle n·axis = mu (|mu|<1)."""
    a = axis / np.linalg.norm(axis)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    tmp = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(a, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    t = np.linspace(0.0, 2.0 * np.pi, n)
    return (mu * a)[None, :] + rho * (np.outer(np.cos(t), e1) + np.outer(np.sin(t), e2))


def split_lon_segments(lon_rad: np.ndarray, lat_rad: np.ndarray, gap: float = np.pi):
    """Split polyline where Mollweide lon jumps (wrap discontinuity)."""
    d = np.diff(lon_rad)
    cuts = np.where(np.abs(d) > gap)[0] + 1
    idx = np.concatenate([[0], cuts, [len(lon_rad)]])
    segs = []
    for i0, i1 in zip(idx[:-1], idx[1:]):
        if i1 - i0 >= 2:
            segs.append((lon_rad[i0:i1], lat_rad[i0:i1]))
    return segs


def plot_lb_curve_mollweide(ax, l_deg, b_deg, **kwargs):
    lon = np.deg2rad(lon_disp_from_ell(l_deg))
    lat = np.deg2rad(b_deg)
    for xr, yr in split_lon_segments(lon, lat):
        ax.plot(xr, yr, **kwargs)


def band_edges(centers: list[float], half: float = BAND_HALF) -> list[tuple[float, float]]:
    """Return (theta_lo, theta_hi) for each band center (±half about center)."""
    return [(c - half, c + half) for c in centers]


def mu_pair_for_band(lo: float, hi: float) -> tuple[float, float]:
    """Edge μ=cos(θ) for a scan-angle band (two boundary cones only)."""
    return float(np.cos(np.deg2rad(lo))), float(np.cos(np.deg2rad(hi)))


def migkas_scan_stops() -> dict:
    th_max = [float(t) for t in STRIP_THETA_MAX]
    th_min = [float(t) for t in STRIP_THETA_MIN]
    max_bands = band_edges(th_max)
    min_bands = band_edges(th_min)
    return {
        "theta_max_from_plus_x_deg": th_max,
        "theta_min_from_plus_x_deg": th_min,
        "max_bands_deg": max_bands,
        "min_bands_deg": min_bands,
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "band_total_width_deg": BAND_TOTAL,
        "band_width_note": "5% of 360° = 18° end-to-end = ±9° about each MAX/MIN center",
        "mu_max_centers": [float(np.cos(np.deg2rad(t))) for t in th_max],
        "mu_min_centers": [float(np.cos(np.deg2rad(t))) for t in th_min],
        "mu_max_bands": [mu_pair_for_band(lo, hi) for lo, hi in max_bands],
        "mu_min_bands": [mu_pair_for_band(lo, hi) for lo, hi in min_bands],
        "source": (
            "best_fit_2ell_multiplier.npz A_after extrema (two max, two min); "
            f"each drawn as one transparent {BAND_FULL:.0f}° end-to-end highlight "
            f"(±{BAND_HALF:.0f}° = 5% of 360°) + one solid mid-stroke"
        ),
    }


def compute_geometry(st: dict) -> dict:
    anti_l, anti_b = 0.0, 15.0
    coe_l, coe_b = 180.0, -15.0

    n_coe = lb_to_n(coe_l, coe_b)
    n_anti = lb_to_n(anti_l, anti_b)
    D = n_anti  # +D = anti-CoE = +x

    th_max = [float(t) for t in st["theta_max_from_plus_x_deg"]]
    th_min = [float(t) for t in st["theta_min_from_plus_x_deg"]]
    max_bands = list(st["max_bands_deg"])
    min_bands = list(st["min_bands_deg"])

    n_maxes = [lb_to_n(*lb) for lb in MIGKAS_MAXES]
    n_mins = [lb_to_n(*lb) for lb in MIGKAS_MINS]

    sep = float(
        np.rad2deg(
            np.arccos(np.clip(np.dot(n_maxes[0], n_mins[1]), -1.0, 1.0))
        )
    )

    return {
        "D": D,
        "n_coe": n_coe,
        "n_anti": n_anti,
        "n_maxes": n_maxes,
        "n_mins": n_mins,
        "anti_lb": (anti_l, anti_b),
        "coe_lb": (coe_l, coe_b),
        "theta_max_from_plus_x_deg": th_max,
        "theta_min_from_plus_x_deg": th_min,
        "max_bands_deg": max_bands,
        "min_bands_deg": min_bands,
        "mu_max_bands": list(st["mu_max_bands"]),
        "mu_min_bands": list(st["mu_min_bands"]),
        "sep_deg": sep,
        "band_half_width_deg": st["band_half_width_deg"],
    }


def draw_wire_sphere(ax, n_mer=12, n_par=6, **kwargs):
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(-np.pi / 2, np.pi / 2, 40)
    xs = np.outer(np.cos(u), np.cos(v))
    ys = np.outer(np.sin(u), np.cos(v))
    zs = np.outer(np.ones_like(u), np.sin(v))
    ax.plot_surface(xs, ys, zs, color="#dfe6e9", alpha=0.15, linewidth=0, shade=False)
    kw = dict(color="#95a5a6", lw=0.4, alpha=0.45)
    kw.update(kwargs)
    for i in range(n_mer):
        phi = 2 * np.pi * i / n_mer
        t = np.linspace(-np.pi / 2, np.pi / 2, 80)
        ax.plot(np.cos(t) * np.cos(phi), np.cos(t) * np.sin(phi), np.sin(t), **kw)
    for j in range(1, n_par):
        b = -np.pi / 2 + np.pi * j / n_par
        t = np.linspace(0, 2 * np.pi, 80)
        r = np.cos(b)
        ax.plot(r * np.cos(t), r * np.sin(t), np.full_like(t, np.sin(b)), **kw)


def foreshortened_ellipse_xy(mu: float, n: int = 200):
    """Side-view foreshortened ellipse of a 3D ring at n·x = mu."""
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    phi = np.linspace(0, 2 * np.pi, n)
    ex = mu + 0.045 * rho * np.cos(phi)
    ey = rho * np.sin(phi)
    return ex, ey, rho


def draw_panelA_ribbon(ax, theta_lo: float, theta_hi: float, theta_mid: float, color: str):
    """
    One Panel-A highlight: clean vertical strip between x=cos(θ_hi)..cos(θ_lo)
    clipped to the unit circle (uniform alpha), plus ONE solid foreshortened
    mid-ellipse at the exact center θ. No multi-stroke wash / hairball.
    """
    mu_lo = float(np.cos(np.deg2rad(theta_lo)))
    mu_hi = float(np.cos(np.deg2rad(theta_hi)))
    x0, x1 = (min(mu_lo, mu_hi), max(mu_lo, mu_hi))
    # Vertical strip polygon (tall enough to cover the circle; clipped below)
    strip = Polygon(
        np.array([[x0, -1.15], [x1, -1.15], [x1, 1.15], [x0, 1.15]]),
        closed=True,
        facecolor=color,
        edgecolor="none",
        alpha=0.28,
        zorder=2,
        linewidth=0,
    )
    ax.add_patch(strip)
    clip_circ = Circle((0.0, 0.0), 1.0, transform=ax.transData)
    strip.set_clip_path(clip_circ)

    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    ex_m, ey_m, _ = foreshortened_ellipse_xy(mu_mid, n=256)
    ax.plot(ex_m, ey_m, color=color, lw=2.8, alpha=0.95, zorder=4, solid_capstyle="round")


def draw_panelA_theta_rays(ax, geo):
    """
    Illustrate θ: from the observer (origin), one ray per loop out to the
    ring at that loop's center θ. θ is measured from +x (anti-CoE, right)
    up to the ring. Emphasize MAX θ=30° (thick semi-transparent red to the
    top of that foreshortened ring), matching Chris's markup.
    """
    items = list(zip(geo["theta_max_from_plus_x_deg"], COLORS_MAX)) + list(
        zip(geo["theta_min_from_plus_x_deg"], COLORS_MIN)
    )
    for mid, color in items:
        t = float(np.deg2rad(mid))
        # Endpoint = top/rim point of the ring in side view = (cos θ, sin θ)
        x1, y1 = float(np.cos(t)), float(np.sin(t))
        is_max30 = abs(mid - 29.5) < 0.5
        lw = 3.2 if is_max30 else 1.8
        alpha = 0.70 if is_max30 else 0.45
        ax.plot(
            [0.0, x1], [0.0, y1],
            color=color, lw=lw, alpha=alpha, zorder=3,
            solid_capstyle="round",
        )
        # Small tick at the ring hit
        ax.plot([x1], [y1], "o", color=color, ms=4.5 if is_max30 else 3.5,
                alpha=0.9, zorder=5, markeredgecolor="none")

    # Angle arc at θ=30° from +x axis up to the ray (makes the measure obvious)
    arc_r = 0.38
    phi = np.linspace(0.0, np.deg2rad(29.5), 60)
    ax.plot(
        arc_r * np.cos(phi), arc_r * np.sin(phi),
        color=C_MAX_30, lw=1.6, alpha=0.85, zorder=3,
    )
    ax.text(
        arc_r * np.cos(np.deg2rad(15.0)) + 0.06,
        arc_r * np.sin(np.deg2rad(15.0)) + 0.04,
        r"$\theta$",
        color=C_MAX_30, fontsize=10, fontweight="bold",
        ha="left", va="bottom", zorder=6,
    )


def panel_A_model_frame(ax, geo):
    """Side view: exactly four thick foreshortened loop ribbons about x."""
    ax.set_title(
        "A. Model frame (no Galactic yet)\n"
        r"axisymmetric about $x$ = CoE$\leftrightarrow$anti-CoE"
        "\n"
        r"two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$)",
        fontsize=10,
        pad=8,
    )
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color="#7f8c8d", lw=1.2, zorder=1)
    ax.plot(0, 0, "ko", ms=6, zorder=5)
    ax.text(0.08, -0.12, "observer", fontsize=8, color="k", ha="left")

    ax.annotate(
        "",
        xy=(1.15, 0),
        xytext=(-1.15, 0),
        arrowprops=dict(arrowstyle="<->", color=C_AXIS, lw=2.5),
        zorder=4,
    )
    ax.text(1.18, 0.08, "anti-CoE\n(+x)", color=C_ANTI, fontsize=8, ha="left", va="bottom")
    ax.text(-1.18, 0.08, "CoE\n(−x)", color=C_COE, fontsize=8, ha="right", va="bottom")

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        draw_panelA_ribbon(ax, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        draw_panelA_ribbon(ax, lo, hi, mid, color)

    # Illustrate θ from +x (anti-CoE) to each loop (Chris markup)
    draw_panelA_theta_rays(ax, geo)

    label_r = 1.30
    for mid, color in list(zip(geo["theta_max_from_plus_x_deg"], COLORS_MAX)) + list(
        zip(geo["theta_min_from_plus_x_deg"], COLORS_MIN)
    ):
        t = np.deg2rad(mid)
        ax.text(
            label_r * np.cos(t),
            label_r * np.sin(t),
            rf"$\theta\!=\!{mid:g}^\circ$",
            color=color, fontsize=7.5, ha="center", va="center", fontweight="bold",
        )

    ax.text(
        0.0,
        -1.38,
        (
            r"four transparent highlight bands about $x$  ·  each $18^\circ$ end-to-end ($\pm9^\circ$ = 5% of $360^\circ$)"
            + "\n"
            + r"(after-2$\ell$ $A(\ell)$: two MAX + two MIN; one mid-stroke each)"
            + "\n"
            + r"$y,z$ free (roll) — model does not pick celestial equator"
        ),
        fontsize=7.5,
        ha="center",
        va="top",
        style="italic",
        color="#34495e",
    )
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")


def _ring_lb(D, mu, n_pts=721):
    pts = small_circle_about_axis(D, mu, n=n_pts)
    rl = np.array([n_to_lb(p)[0] for p in pts])
    rb = np.array([n_to_lb(p)[1] for p in pts])
    return rl, rb


def _ring_lonlat_rad(D, mu, n_pts=361):
    """Small-circle → Mollweide display (lon_disp, lat) in radians."""
    pts = small_circle_about_axis(D, float(mu), n=n_pts)
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    lon = np.deg2rad(lon_disp_from_ell(ls))
    lat = np.deg2rad(bs)
    return lon, lat


def fill_mollweide_ribbon(ax, D, theta_lo, theta_hi, theta_mid, color, n_pts=721):
    """
    One Mollweide highlight band:
      continuous transparent polygon strips between ±9° edges
      (broken only at Mollweide antimeridian seams — NOT a dense quad wash)
      + ONE solid mid-stroke at θ_mid.
    """
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    lon_a, lat_a = _ring_lonlat_rad(D, mu_a, n_pts=n_pts)
    lon_b, lat_b = _ring_lonlat_rad(D, mu_b, n_pts=n_pts)
    gap = 0.9 * np.pi
    # Mark steps that cross a seam on either edge (or between edges)
    bad = np.zeros(n_pts - 1, dtype=bool)
    for i in range(n_pts - 1):
        if abs(lon_a[i + 1] - lon_a[i]) > gap:
            bad[i] = True
        elif abs(lon_b[i + 1] - lon_b[i]) > gap:
            bad[i] = True
        elif abs(lon_a[i] - lon_b[i]) > gap:
            bad[i] = True
        elif abs(lon_a[i + 1] - lon_b[i + 1]) > gap:
            bad[i] = True
    # Group into continuous good runs → one polygon per run
    polys = []
    i = 0
    while i < n_pts - 1:
        while i < n_pts - 1 and bad[i]:
            i += 1
        if i >= n_pts - 1:
            break
        j = i
        while j < n_pts - 1 and not bad[j]:
            j += 1
        # indices i .. j inclusive are a continuous ribbon segment
        i0, i1 = i, j  # i1 is last good vertex index
        if i1 - i0 >= 1:
            xs = np.concatenate([lon_a[i0:i1 + 1], lon_b[i0:i1 + 1][::-1]])
            ys = np.concatenate([lat_a[i0:i1 + 1], lat_b[i0:i1 + 1][::-1]])
            polys.append(np.column_stack([xs, ys]))
        i = j + 1
    if polys:
        ax.add_collection(
            PolyCollection(
                polys, facecolors=color, edgecolors="none",
                alpha=0.26, zorder=2, linewidths=0,
            )
        )
    rl, rb = _ring_lb(D, mu_mid, n_pts=max(n_pts, 721))
    plot_lb_curve_mollweide(
        ax, rl, rb, color=color, lw=2.6, alpha=0.95, zorder=4
    )


def panel_B_mollweide(ax, geo):
    ax.set_title(
        r"B. Migkas $A(\ell)$ loops on sky (anti-CoE at $(0^\circ,+15^\circ)$)"
        "\n"
        r"(Galactic Mollweide; CoE axis tilted $+15^\circ$; $\ell$ increases LEFT; "
        r"two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$))",
        fontsize=10,
        pad=8,
    )
    for b in [-60, -30, 0, 30, 60]:
        ell = np.linspace(0, 360, 361)
        plot_lb_curve_mollweide(ax, ell, np.full_like(ell, b), color="#bdc3c7", lw=0.5, zorder=0)
    for ell0 in range(0, 360, 30):
        b = np.linspace(-90, 90, 181)
        plot_lb_curve_mollweide(ax, np.full_like(b, ell0), b, color="#bdc3c7", lw=0.5, zorder=0)

    D = geo["D"]
    ls, bs = great_circle_through_poles(D, n=721)
    plot_lb_curve_mollweide(ax, ls, bs, color=C_AXIS, lw=2.0, zorder=2)

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        fill_mollweide_ribbon(ax, D, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        fill_mollweide_ribbon(ax, D, lo, hi, mid, color)

    ell_eq, b_eq = celestial_equator_lb(721)
    plot_lb_curve_mollweide(ax, ell_eq, b_eq, color=C_EQ, lw=1.8, ls="--", zorder=2)

    ax.plot(
        np.deg2rad(lon_disp_from_ell(geo["anti_lb"][0])),
        np.deg2rad(geo["anti_lb"][1]),
        marker="o", ms=8, color=C_ANTI,
        markerfacecolor=C_ANTI, markeredgecolor="k", markeredgewidth=0.6,
        zorder=6, linestyle="None",
    )
    ax.plot(
        np.deg2rad(lon_disp_from_ell(geo["coe_lb"][0])),
        np.deg2rad(geo["coe_lb"][1]),
        marker="*", ms=14, color=C_COE,
        markerfacecolor=C_COE, markeredgecolor="k", markeredgewidth=0.6,
        zorder=6, linestyle="None",
    )
    for (l, b), color in zip(MIGKAS_MAXES, COLORS_MAX):
        ax.plot(
            np.deg2rad(lon_disp_from_ell(l)),
            np.deg2rad(b),
            marker="^",
            ms=10,
            color=color,
            markeredgecolor="k",
            markeredgewidth=0.7,
            zorder=7,
            linestyle="None",
        )
    for (l, b), color in zip(MIGKAS_MINS, COLORS_MIN):
        ax.plot(
            np.deg2rad(lon_disp_from_ell(l)),
            np.deg2rad(b),
            marker="s",
            ms=10,
            color=color,
            markeredgecolor="k",
            markeredgewidth=0.7,
            zorder=7,
            linestyle="None",
        )

    xticks = np.deg2rad(lon_disp_from_ell(np.array([180, 90, 0, 270])))
    ax.set_xticks(xticks)
    ax.set_xticklabels([r"$180^\circ$", r"$90^\circ$", r"$0^\circ$", r"$270^\circ$"], fontsize=8)
    ax.set_yticks(np.deg2rad([-60, -30, 0, 30, 60]))
    ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$+30^\circ$", r"$+60^\circ$"], fontsize=8)
    ax.grid(False)
    ax.set_xlabel(r"Galactic $\ell$ (increases LEFT)", fontsize=8)
    ax.set_ylabel(r"Galactic $b$", fontsize=8)


def fill_sphere_ribbon(ax, D, theta_lo, theta_hi, theta_mid, color, n_pts=160):
    """
    One 3D highlight belt: SINGLE transparent strip between ±9° edges
    + ONE solid mid-stroke. No linspace wash / stacked alpha hairball.
    """
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    p0 = small_circle_about_axis(D, float(mu_a), n=n_pts)
    p1 = small_circle_about_axis(D, float(mu_b), n=n_pts)
    # Single belt between the two edge rings (one collection, no stacked wash)
    verts = [[p0[j], p0[j + 1], p1[j + 1], p1[j]] for j in range(n_pts - 1)]
    coll = Poly3DCollection(
        verts, facecolors=color, edgecolors="none", linewidths=0.0,
        alpha=0.25, zorder=2,
    )
    ax.add_collection3d(coll)
    pts = small_circle_about_axis(D, mu_mid, n=n_pts)
    ax.plot(
        pts[:, 0], pts[:, 1], pts[:, 2],
        color=color, lw=2.6, alpha=0.95, zorder=5,
    )


def panel_C_sphere(ax, geo):
    th_max = geo["theta_max_from_plus_x_deg"]
    th_min = geo["theta_min_from_plus_x_deg"]
    ax.set_title(
        r"C. Sphere 3D — two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$)"
        "\n"
        r"(axis line emphasized; celestial equator dashed)",
        fontsize=10,
        pad=6,
    )
    draw_wire_sphere(ax)
    D = geo["D"]
    ell_eq, b_eq = celestial_equator_lb(400)
    eq_xyz = np.array([lb_to_n(l, b) for l, b in zip(ell_eq, b_eq)])
    ax.plot(eq_xyz[:, 0], eq_xyz[:, 1], eq_xyz[:, 2], color=C_EQ, lw=1.8, ls="--", zorder=3)
    ax.plot(
        [-D[0], D[0]],
        [-D[1], D[1]],
        [-D[2], D[2]],
        color=C_AXIS,
        lw=3.0,
        zorder=5,
    )
    ax.scatter(*(-D), color=C_COE, s=120, marker="*", edgecolors="k", linewidths=0.6,
               depthshade=False, zorder=8)
    ax.scatter(*(D), color=C_ANTI, s=80, edgecolors="k", linewidths=0.6, depthshade=False, zorder=8)
    ax.text(*(-D * 1.22), "CoE", color=C_COE, fontsize=8, fontweight="bold")
    ax.text(*(D * 1.18), "anti-CoE", color=C_ANTI, fontsize=8, fontweight="bold")
    ax.scatter([0], [0], [0], color="k", s=25, depthshade=False, zorder=9)

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        fill_sphere_ribbon(ax, D, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        fill_sphere_ribbon(ax, D, lo, hi, mid, color)

    for n, color in zip(geo["n_maxes"], COLORS_MAX):
        ax.scatter(
            [n[0]], [n[1]], [n[2]],
            c=[color], marker="^", s=70,
            edgecolors="k", linewidths=0.7, depthshade=False, zorder=10,
        )
        ax.plot([0, n[0]], [0, n[1]], [0, n[2]], color=color, lw=0.8, ls=":", alpha=0.7)
    for n, color in zip(geo["n_mins"], COLORS_MIN):
        ax.scatter(
            [n[0]], [n[1]], [n[2]],
            c=[color], marker="s", s=60,
            edgecolors="k", linewidths=0.7, depthshade=False, zorder=10,
        )
        ax.plot([0, n[0]], [0, n[1]], [0, n[2]], color=color, lw=0.8, ls=":", alpha=0.7)

    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.set_zlim(-1.25, 1.25)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=18, azim=-55)
    ax.set_axis_off()
    ax.text2D(
        0.5,
        -0.02,
        rf"MAX $\theta\!=\!{th_max[0]:g}^\circ,{th_max[1]:g}^\circ$; "
        rf"MIN $\theta\!=\!{th_min[0]:g}^\circ,{th_min[1]:g}^\circ$  "
        rf"($18^\circ$ end-to-end = $\pm9^\circ$ = 5% of $360^\circ$)",
        transform=ax.transAxes,
        ha="center",
        fontsize=7.5,
        color="#34495e",
    )


def main():
    st = migkas_scan_stops()
    geo = compute_geometry(st)

    th_max = geo["theta_max_from_plus_x_deg"]
    th_min = geo["theta_min_from_plus_x_deg"]
    max_bands = geo["max_bands_deg"]
    min_bands = geo["min_bands_deg"]

    fig = plt.figure(figsize=(14.5, 5.6), dpi=160)
    fig.suptitle(
        r"Pure CoE axisymmetry: two MAX + two MIN highlight bands"
        "\n"
        rf"(each $18^\circ$ end-to-end = $\pm9^\circ$ = 5% of $360^\circ$; "
        rf"MAX $\theta\!=\!{th_max[0]:g}^\circ/{th_max[1]:g}^\circ$; "
        rf"MIN $\theta\!=\!{th_min[0]:g}^\circ/{th_min[1]:g}^\circ$; "
        r"anti-CoE $(0^\circ,+15^\circ)$, CoE axis tilted $+15^\circ$)",
        fontsize=11,
        fontweight="bold",
        y=1.02,
    )

    axA = fig.add_subplot(1, 3, 1)
    panel_A_model_frame(axA, geo)

    axB = fig.add_subplot(1, 3, 2, projection="mollweide")
    panel_B_mollweide(axB, geo)

    axC = fig.add_subplot(1, 3, 3, projection="3d")
    panel_C_sphere(axC, geo)

    handles = [
        Line2D([0], [0], color=C_AXIS, lw=2.5, label=r"CoE$\leftrightarrow$anti-CoE axis"),
        Patch(
            facecolor=C_MAX_30, alpha=0.45, edgecolor=C_MAX_30,
            label=(
                rf"MAX $\theta\!=\!{th_max[0]:g}^\circ$ band "
                rf"[{max_bands[0][0]:g},{max_bands[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$=5% of $360^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MAX_185, alpha=0.45, edgecolor=C_MAX_185,
            label=(
                rf"MAX $\theta\!=\!{th_max[1]:g}^\circ$ band "
                rf"[{max_bands[1][0]:g},{max_bands[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_115, alpha=0.45, edgecolor=C_MIN_115,
            label=(
                rf"MIN $\theta\!=\!{th_min[0]:g}^\circ$ band "
                rf"[{min_bands[0][0]:g},{min_bands[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_265, alpha=0.45, edgecolor=C_MIN_265,
            label=(
                rf"MIN $\theta\!=\!{th_min[1]:g}^\circ$ band "
                rf"[{min_bands[1][0]:g},{min_bands[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Line2D([0], [0], color=C_EQ, lw=1.8, ls="--", label="celestial equator (Dec$=0$)"),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_30, markeredgecolor="k",
            ms=9,
            label=rf"MAX $({MIGKAS_MAXES[0][0]:.0f}^\circ,{MIGKAS_MAXES[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_185, markeredgecolor="k",
            ms=9,
            label=rf"MAX $({MIGKAS_MAXES[1][0]:.0f}^\circ,{MIGKAS_MAXES[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_115, markeredgecolor="k",
            ms=8,
            label=rf"MIN $({MIGKAS_MINS[0][0]:.0f}^\circ,{MIGKAS_MINS[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_265, markeredgecolor="k",
            ms=8,
            label=rf"MIN $({MIGKAS_MINS[1][0]:.0f}^\circ,{MIGKAS_MINS[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="o", color="w", markerfacecolor=C_ANTI, markeredgecolor="k",
            ms=8, label=rf"anti-CoE $({geo['anti_lb'][0]:.1f}^\circ,{geo['anti_lb'][1]:+.1f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="*", color="w", markerfacecolor=C_COE, markeredgecolor="k",
            ms=11, label=rf"CoE $({geo['coe_lb'][0]:.1f}^\circ,{geo['coe_lb'][1]:+.1f}^\circ)$",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=6,
        fontsize=7.0,
        frameon=True,
        fancybox=True,
        bbox_to_anchor=(0.5, -0.10),
    )

    fig.tight_layout()
    # --- Chris 2026-09-26: panel C sphere drawn at the same on-page diameter as
    # panel A's model-frame circle; panel C axes matched to A's vertical
    # position/height. Nothing else changes.
    from mpl_toolkits.mplot3d import proj3d
    pa, pc = axA.get_position(), axC.get_position()
    axC.set_position([pc.x0, pa.y0, pc.width, pa.height])

    def _diam_A():
        fig.canvas.draw()
        p0 = axA.transData.transform((-1.0, 0.0)); p1 = axA.transData.transform((1.0, 0.0))
        return float(p1[0] - p0[0])

    def _diam_C():
        fig.canvas.draw()
        u = np.random.default_rng(0).normal(size=(20000, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]
        xs, ys, _ = proj3d.proj_transform(u[:, 0], u[:, 1], u[:, 2], axC.get_proj())
        d = axC.transData.transform(np.column_stack([xs, ys]))
        return float(d[:, 0].max() - d[:, 0].min()), float(d[:, 1].max() - d[:, 1].min())

    zoom = 1.0
    for _ in range(8):
        dA = _diam_A(); dCx, dCy = _diam_C()
        zoom *= dA / dCx
        axC.set_box_aspect((1, 1, 1), zoom=zoom)
    dA = _diam_A(); dCx, dCy = _diam_C()
    print(f"PANEL_DIAM A={dA:.1f}px C={dCx:.1f}x{dCy:.1f}px zoom={zoom:.3f}")
    fig.savefig(PNG_PATH, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    payload = {
        "png": str(PNG_PATH),
        "script": str(OUT_DIR / "make_coe_rings_migkas_pin_axis_matched.py"),
        "migkas_max_lb_deg": [list(lb) for lb in MIGKAS_MAXES],
        "migkas_min_lb_deg": [list(lb) for lb in MIGKAS_MINS],
        "strip_theta_max_deg": list(STRIP_THETA_MAX),
        "strip_theta_min_deg": list(STRIP_THETA_MIN),
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "band_total_width_deg": BAND_TOTAL,
        "band_width_note": "5% of 360° = 18° end-to-end = ±9° about each MAX/MIN center",
        "max_bands_deg": [list(b) for b in max_bands],
        "min_bands_deg": [list(b) for b in min_bands],
        "band_intervals": {
            "MAX_30": [float(max_bands[0][0]), float(max_bands[0][1])],
            "MAX_185": [float(max_bands[1][0]), float(max_bands[1][1])],
            "MIN_115": [float(min_bands[0][0]), float(min_bands[0][1])],
            "MIN_265": [float(min_bands[1][0]), float(min_bands[1][1])],
        },
        "scan_stop_source": st["source"],
        "theta_max_from_plus_x_deg": [round(t, 4) for t in th_max],
        "theta_min_from_plus_x_deg": [round(t, 4) for t in th_min],
        "mu_max_bands": [[round(a, 8), round(b, 8)] for a, b in geo["mu_max_bands"]],
        "mu_min_bands": [[round(a, 8), round(b, 8)] for a, b in geo["mu_min_bands"]],
        "CoE_lb_deg": [round(geo["coe_lb"][0], 4), round(geo["coe_lb"][1], 4)],
        "anti_CoE_lb_deg": [round(geo["anti_lb"][0], 4), round(geo["anti_lb"][1], 4)],
        "D_hat_galactic_xyz": geo["D"].tolist(),
        "angular_sep_migkas_max0_min1_deg": round(geo["sep_deg"], 6),
        "axis_tilt_b_deg": 15,
        "axis_tilt_reason": (
            "Chris 2026-09-22 — markers fixed, rings rotate up 30°"
        ),
        "orientation": (
            "Galactic orientation: anti-CoE at (ell=0°, b=+15°), CoE at "
            "(ell=180°, b=−15°); +D = anti-CoE = +x; axis tilted +15° "
            "(Chris 2026-09-22 — markers fixed, rings rotate up 30°)"
        ),
        "drawing": (
            "Exactly 4 features. Each = ONE smooth transparent highlight band "
            "between θ−9° and θ+9° (18° end-to-end = 5% of 360°) + ONE solid "
            "mid-stroke at the exact center θ. No multi-stroke wash / hairball. "
            "Four distinct colors (two reds + two blues); no orange; no 283° "
            "companion. Panel A: clean vertical strip clipped to unit circle; theta rays from observer to each loop (emphasized at MAX 30°). Mollweide "
            "seam-split honestly. Sphere: single Poly3D strip (no n_wash stack)."
        ),
        "note": (
            "Stops are after-2ℓ A(ℓ) extrema from best_fit_2ell_multiplier.npz "
            f"A_after: two MAX at θ={th_max[0]:g}°/{th_max[1]:g}° and "
            f"two MIN at θ={th_min[0]:g}°/{th_min[1]:g}°, each drawn as "
            f"one {BAND_FULL:.0f}° end-to-end highlight (±{BAND_HALF:.0f}° = 5% of 360°). "
            "Rings use μ=cos(θ); sky markers are map extrema (17°,−9°) and "
            "(274°,−22°) plus strip stand-ins (185°,0°) and (115°,0°)."
        ),
        "plot_convention": "Migkas Galactic Mollweide: lon_disp=-ell_wrapped (ell increases LEFT)",
        "legend_summary": (
            "two MAX + two MIN highlight bands, each 18° end-to-end "
            "(±9° = 5% of 360°); four distinct colors"
        ),
        "colors": {
            "MAX_30": C_MAX_30,
            "MAX_185": C_MAX_185,
            "MIN_115": C_MIN_115,
            "MIN_265": C_MIN_265,
            "no_orange": True,
        },
    }
    JSON_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    print(f"\nWrote {PNG_PATH}")
    print(f"Wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
```

#### `make_coe_rings_migkas_pin_axis_matched_grad.py`

```python
# source: /workspace/cosmology/figures/make_coe_rings_migkas_pin_axis_matched_grad.py
#!/usr/bin/env python3
"""
CoE rings + Migkas A(ell) extrema loops in the Galactic frame.

Concept:
  Pure CoE axisymmetry → constant A on small-circle rings about x = CoE↔anti-CoE.
  Scan A vs angle from +x (anti-CoE) around 360°. After-2ℓ A(ℓ) has TWO maxima
  and TWO minima (from best_fit_2ell_multiplier.npz A_after):
    MAX at θ = 30° and 185° from +x,
    MIN at θ = 115° and 265° from +x.

  Exactly FOUR loops — each ONE transparent highlight band, 18° end-to-end
  (= ±9° = 5% of 360°) about the center θ:
    MAX: 30±9 → [21°, 39°] and 185±9 → [176°, 194°]
    MIN: 115±9 → [106°, 124°] and 265±9 → [256°, 274°]
  Drawing: ONE smooth transparent filled ribbon between θ−9° and θ+9°
  + ONE solid mid-stroke at the exact center θ. No multi-stroke wash,
  no dense linspace hairball, no dual strong boundaries.
  Four distinct colors (two reds + two blues). No orange, no 283° companion.

  Galactic orientation: anti-CoE (+x) at ℓ=0°, b=+15°; CoE (−x) at ℓ=180°, b=−15°
  (axis tilted +15° so features near ℓ=0 move north on Mollweide; Migkas markers fixed).
  Sky markers (ℓ, b):
    Maxes: map extrema (17°, −9°) for ~30° max; strip (185°, 0°) for second max.
    Mins:  strip (115°, 0°) for ~115° min; map extrema (274°, −22°) for ~265° min.

Outputs:
  coe_rings_migkas_pin_axis.png
  coe_rings_migkas_pin_axis.json
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle, Patch, Polygon
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

OUT_DIR = Path(__file__).resolve().parent
PNG_PATH = Path("/workspace/cosmology/papers/spreadsheet_model/loop_threepanel_matched_grad.png")
JSON_PATH = Path("/workspace/cosmology/papers/spreadsheet_model/loop_threepanel_matched_grad.json")

# 5% of 360° = 18° end-to-end → half-width ±9° about each center
BAND_HALF = 9.0
BAND_FULL = 2.0 * BAND_HALF  # 18° = 5% of 360°
BAND_TOTAL = BAND_FULL  # alias
STRIP_THETA_MAX = [29.5, 187.5]
STRIP_THETA_MIN = [122.6, 264.4]

MIGKAS_MAXES = [(17.0, -9.0), (185.0, 0.0)]
MIGKAS_MINS = [(115.0, 0.0), (274.0, -22.0)]

C_COE = "#6b3fa0"
C_ANTI = "#8e1b1b"  # dark red marker; distinct from loop reds
# Four distinct loop colors: two reds (MAX) + two blues (MIN); no orange/yellow
C_MAX_30 = "#0097b2"    # cyan — MAX at 29.5°
C_MAX_185 = "#e0107a"   # magenta — MAX at 187.5°
C_MIN_115 = "#6fa800"   # green — MIN at 122.6°
C_MIN_265 = "#2166d4"   # blue — MIN at 264.4°
C_AXIS = "#5b2c6f"
C_EQ = "#1e8449"

# Per-loop color maps (order matches STRIP_THETA_MAX / STRIP_THETA_MIN)
COLORS_MAX = [C_MAX_30, C_MAX_185]
COLORS_MIN = [C_MIN_115, C_MIN_265]


def wrap_lon_deg(lon: np.ndarray) -> np.ndarray:
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def lon_disp_from_ell(ell_deg: np.ndarray) -> np.ndarray:
    """Migkas Galactic: ℓ increases LEFT → lon_disp = −ℓ_wrapped."""
    return -wrap_lon_deg(ell_deg)


def lb_to_n(l_deg: float, b_deg: float) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    return np.array(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], dtype=float
    )


def n_to_lb(n: np.ndarray) -> tuple[float, float]:
    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    b = float(np.rad2deg(np.arcsin(np.clip(n[2], -1.0, 1.0))))
    l = float(np.rad2deg(np.arctan2(n[1], n[0])) % 360.0)
    return l, b


def celestial_equator_lb(n: int = 720) -> tuple[np.ndarray, np.ndarray]:
    """Dec=0 (ICRS) as Galactic (ℓ, b). IAU J2000 NGP convention, no astropy."""
    ra_ngp = np.deg2rad(192.85948)
    dec_ngp = np.deg2rad(27.12825)
    l_ncp = np.deg2rad(122.93192)
    ra = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    dec = 0.0
    sin_b = np.sin(dec) * np.sin(dec_ngp) + np.cos(dec) * np.cos(dec_ngp) * np.cos(ra - ra_ngp)
    sin_b = np.clip(sin_b, -1.0, 1.0)
    b = np.arcsin(sin_b)
    y = np.cos(dec) * np.sin(ra - ra_ngp)
    x = np.sin(dec) * np.cos(dec_ngp) - np.cos(dec) * np.sin(dec_ngp) * np.cos(ra - ra_ngp)
    l = (l_ncp - np.arctan2(y, x)) % (2.0 * np.pi)
    return np.rad2deg(l), np.rad2deg(b)


def great_circle_through_poles(D: np.ndarray, n: int = 361) -> tuple[np.ndarray, np.ndarray]:
    """Great circle through ±D (any meridian of the CoE axis)."""
    a = D / np.linalg.norm(D)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([1.0, 0.0, 0.0])
    e2 = np.cross(a, tmp)
    e2 /= np.linalg.norm(e2)
    t = np.linspace(0.0, 2.0 * np.pi, n)
    pts = np.outer(np.cos(t), a) + np.outer(np.sin(t), e2)
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    return ls, bs


def small_circle_about_axis(axis: np.ndarray, mu: float, n: int = 361) -> np.ndarray:
    """Unit vectors on the small circle n·axis = mu (|mu|<1)."""
    a = axis / np.linalg.norm(axis)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    tmp = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(a, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    t = np.linspace(0.0, 2.0 * np.pi, n)
    return (mu * a)[None, :] + rho * (np.outer(np.cos(t), e1) + np.outer(np.sin(t), e2))


def split_lon_segments(lon_rad: np.ndarray, lat_rad: np.ndarray, gap: float = np.pi):
    """Split polyline where Mollweide lon jumps (wrap discontinuity)."""
    d = np.diff(lon_rad)
    cuts = np.where(np.abs(d) > gap)[0] + 1
    idx = np.concatenate([[0], cuts, [len(lon_rad)]])
    segs = []
    for i0, i1 in zip(idx[:-1], idx[1:]):
        if i1 - i0 >= 2:
            segs.append((lon_rad[i0:i1], lat_rad[i0:i1]))
    return segs


def plot_lb_curve_mollweide(ax, l_deg, b_deg, **kwargs):
    lon = np.deg2rad(lon_disp_from_ell(l_deg))
    lat = np.deg2rad(b_deg)
    for xr, yr in split_lon_segments(lon, lat):
        ax.plot(xr, yr, **kwargs)


def band_edges(centers: list[float], half: float = BAND_HALF) -> list[tuple[float, float]]:
    """Return (theta_lo, theta_hi) for each band center (±half about center)."""
    return [(c - half, c + half) for c in centers]


def mu_pair_for_band(lo: float, hi: float) -> tuple[float, float]:
    """Edge μ=cos(θ) for a scan-angle band (two boundary cones only)."""
    return float(np.cos(np.deg2rad(lo))), float(np.cos(np.deg2rad(hi)))


def migkas_scan_stops() -> dict:
    th_max = [float(t) for t in STRIP_THETA_MAX]
    th_min = [float(t) for t in STRIP_THETA_MIN]
    max_bands = band_edges(th_max)
    min_bands = band_edges(th_min)
    return {
        "theta_max_from_plus_x_deg": th_max,
        "theta_min_from_plus_x_deg": th_min,
        "max_bands_deg": max_bands,
        "min_bands_deg": min_bands,
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "band_total_width_deg": BAND_TOTAL,
        "band_width_note": "5% of 360° = 18° end-to-end = ±9° about each MAX/MIN center",
        "mu_max_centers": [float(np.cos(np.deg2rad(t))) for t in th_max],
        "mu_min_centers": [float(np.cos(np.deg2rad(t))) for t in th_min],
        "mu_max_bands": [mu_pair_for_band(lo, hi) for lo, hi in max_bands],
        "mu_min_bands": [mu_pair_for_band(lo, hi) for lo, hi in min_bands],
        "source": (
            "best_fit_2ell_multiplier.npz A_after extrema (two max, two min); "
            f"each drawn as one transparent {BAND_FULL:.0f}° end-to-end highlight "
            f"(±{BAND_HALF:.0f}° = 5% of 360°) + one solid mid-stroke"
        ),
    }


def compute_geometry(st: dict) -> dict:
    anti_l, anti_b = 0.0, 15.0
    coe_l, coe_b = 180.0, -15.0

    n_coe = lb_to_n(coe_l, coe_b)
    n_anti = lb_to_n(anti_l, anti_b)
    D = n_anti  # +D = anti-CoE = +x

    th_max = [float(t) for t in st["theta_max_from_plus_x_deg"]]
    th_min = [float(t) for t in st["theta_min_from_plus_x_deg"]]
    max_bands = list(st["max_bands_deg"])
    min_bands = list(st["min_bands_deg"])

    n_maxes = [lb_to_n(*lb) for lb in MIGKAS_MAXES]
    n_mins = [lb_to_n(*lb) for lb in MIGKAS_MINS]

    sep = float(
        np.rad2deg(
            np.arccos(np.clip(np.dot(n_maxes[0], n_mins[1]), -1.0, 1.0))
        )
    )

    return {
        "D": D,
        "n_coe": n_coe,
        "n_anti": n_anti,
        "n_maxes": n_maxes,
        "n_mins": n_mins,
        "anti_lb": (anti_l, anti_b),
        "coe_lb": (coe_l, coe_b),
        "theta_max_from_plus_x_deg": th_max,
        "theta_min_from_plus_x_deg": th_min,
        "max_bands_deg": max_bands,
        "min_bands_deg": min_bands,
        "mu_max_bands": list(st["mu_max_bands"]),
        "mu_min_bands": list(st["mu_min_bands"]),
        "sep_deg": sep,
        "band_half_width_deg": st["band_half_width_deg"],
    }


def draw_wire_sphere(ax, n_mer=12, n_par=6, **kwargs):
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(-np.pi / 2, np.pi / 2, 40)
    xs = np.outer(np.cos(u), np.cos(v))
    ys = np.outer(np.sin(u), np.cos(v))
    zs = np.outer(np.ones_like(u), np.sin(v))
    ax.plot_surface(xs, ys, zs, color="#dfe6e9", alpha=0.15, linewidth=0, shade=False)
    kw = dict(color="#95a5a6", lw=0.4, alpha=0.45)
    kw.update(kwargs)
    for i in range(n_mer):
        phi = 2 * np.pi * i / n_mer
        t = np.linspace(-np.pi / 2, np.pi / 2, 80)
        ax.plot(np.cos(t) * np.cos(phi), np.cos(t) * np.sin(phi), np.sin(t), **kw)
    for j in range(1, n_par):
        b = -np.pi / 2 + np.pi * j / n_par
        t = np.linspace(0, 2 * np.pi, 80)
        r = np.cos(b)
        ax.plot(r * np.cos(t), r * np.sin(t), np.full_like(t, np.sin(b)), **kw)


def foreshortened_ellipse_xy(mu: float, n: int = 200):
    """Side-view foreshortened ellipse of a 3D ring at n·x = mu."""
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    phi = np.linspace(0, 2 * np.pi, n)
    ex = mu + 0.045 * rho * np.cos(phi)
    ey = rho * np.sin(phi)
    return ex, ey, rho



# ---------------------------------------------------------------------------
# Opacity gradient along each loop (Chris 2026-09-26), same method as
# make_coe_rings_gradient_loops.py: loop drawn as runs of short elements with
# per-element alpha = base_alpha * (1 + cos psi)/2, psi = ring azimuth from the
# θ-ray point (ring point most along sign(sin θ_mid)·Ŷ; Ŷ = Galactic +z
# projected ⊥ D, i.e. model-frame "up"; no roll). Darkest at the θ-ray point,
# fully transparent at the diametrically opposite point.
# ---------------------------------------------------------------------------
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def theta_ray_sign(theta_mid_deg):
    return 1.0 if np.sin(np.deg2rad(float(theta_mid_deg))) >= 0 else -1.0


def ring_basis(D):
    a = np.asarray(D, dtype=float)
    a = a / np.linalg.norm(a)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    y = tmp - np.dot(tmp, a) * a
    y /= np.linalg.norm(y)
    z = np.cross(a, y)
    z /= np.linalg.norm(z)
    return y, z


def grad_phi(theta_mid_deg, n):
    """phi from antipode of θ-ray point round to antipode (psi = -180..+180)."""
    phi0 = 0.5 * np.pi if theta_ray_sign(theta_mid_deg) > 0 else 1.5 * np.pi
    return phi0, phi0 + np.linspace(-np.pi, np.pi, n)


def grad_scale(phi, phi0):
    return 0.5 * (1.0 + np.cos(np.asarray(phi) - phi0))


def small_circle_grad(D, mu, theta_mid_deg, n=721):
    a = np.asarray(D, dtype=float)
    a = a / np.linalg.norm(a)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    y, z = ring_basis(a)
    phi0, phi = grad_phi(theta_mid_deg, n)
    pts = (mu * a)[None, :] + rho * (np.outer(np.sin(phi), y) + np.outer(np.cos(phi), z))
    return pts, grad_scale(phi, phi0)


def segment_scale(sc):
    sc = np.asarray(sc, dtype=float)
    return 0.5 * (sc[:-1] + sc[1:])


def rgba_array(color, alphas):
    r, g, b, _ = to_rgba(color)
    al = np.clip(np.asarray(alphas, dtype=float), 0.0, 1.0)
    out = np.empty((al.size, 4))
    out[:, 0], out[:, 1], out[:, 2], out[:, 3] = r, g, b, al
    return out


def quantized_runs(seg_sc, valid=None, n_levels=64):
    s = np.asarray(seg_sc, dtype=float)
    lev = np.round(np.clip(s, 0.0, 1.0) * n_levels) / n_levels
    ok = lev > 0 if valid is None else (lev > 0) & np.asarray(valid, dtype=bool)
    runs, i, n = [], 0, s.size
    while i < n:
        if not ok[i]:
            i += 1
            continue
        j = i + 1
        while j < n and ok[j] and lev[j] == lev[i]:
            j += 1
        runs.append((i, j, float(lev[i])))
        i = j
    return runs


def _strip_polys_2d(xa, ya, xb, yb, runs):
    return [
        np.column_stack([
            np.concatenate([xa[i0:i1 + 1], xb[i0:i1 + 1][::-1]]),
            np.concatenate([ya[i0:i1 + 1], yb[i0:i1 + 1][::-1]]),
        ])
        for i0, i1, _ in runs
    ]


def _polylines_2d(x, y, runs):
    return [np.column_stack([x[i0:i1 + 1], y[i0:i1 + 1]]) for i0, i1, _ in runs]


def draw_panelA_ribbon(ax, theta_lo: float, theta_hi: float, theta_mid: float, color: str):
    """
    Panel-A highlight (edge-on chord): same vertical strip between
    x=cos(θ_hi)..cos(θ_lo) clipped to the unit circle, base alpha 0.28, and the
    same foreshortened mid-ellipse (lw 2.8, alpha 0.95) — now with the opacity
    gradient along the chord: alpha *= (1 + sgn·y/ρ(x))/2 = (1+cos psi)/2,
    darkest at the θ-ray end, transparent at the opposite end.
    """
    mu_lo = float(np.cos(np.deg2rad(theta_lo)))
    mu_hi = float(np.cos(np.deg2rad(theta_hi)))
    x0, x1 = (min(mu_lo, mu_hi), max(mu_lo, mu_hi))
    sgn = theta_ray_sign(theta_mid)
    nx, ny = 64, 1200
    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(-1.0, 1.0, ny)
    X, Y = np.meshgrid(xs, ys)
    rho = np.sqrt(np.clip(1.0 - X * X, 1e-6, None))
    sc = np.clip(0.5 * (1.0 + sgn * Y / rho), 0.0, 1.0)
    img = np.zeros((ny, nx, 4))
    img[..., :3] = to_rgba(color)[:3]
    img[..., 3] = 0.28 * sc
    im = ax.imshow(img, extent=(x0, x1, -1.0, 1.0), origin="lower",
                   interpolation="bilinear", zorder=2, aspect="auto")
    im.set_clip_path(Circle((0.0, 0.0), 1.0, transform=ax.transData))

    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    rho_m = np.sqrt(max(0.0, 1.0 - mu_mid * mu_mid))
    phi0, phi = grad_phi(theta_mid, 721)
    ex_m = mu_mid + 0.045 * rho_m * np.cos(phi)
    ey_m = rho_m * np.sin(phi)
    runs = quantized_runs(segment_scale(grad_scale(phi, phi0)))
    ax.add_collection(
        LineCollection(
            _polylines_2d(ex_m, ey_m, runs),
            colors=rgba_array(color, [0.95 * r[2] for r in runs]),
            linewidths=2.8, zorder=4, capstyle="butt", joinstyle="round",
        )
    )


def draw_panelA_theta_rays(ax, geo):
    """
    Illustrate θ: from the observer (origin), one ray per loop out to the
    ring at that loop's center θ. θ is measured from +x (anti-CoE, right)
    up to the ring. Emphasize MAX θ=30° (thick semi-transparent red to the
    top of that foreshortened ring), matching Chris's markup.
    """
    items = list(zip(geo["theta_max_from_plus_x_deg"], COLORS_MAX)) + list(
        zip(geo["theta_min_from_plus_x_deg"], COLORS_MIN)
    )
    for mid, color in items:
        t = float(np.deg2rad(mid))
        # Endpoint = top/rim point of the ring in side view = (cos θ, sin θ)
        x1, y1 = float(np.cos(t)), float(np.sin(t))
        is_max30 = abs(mid - 29.5) < 0.5
        lw = 3.2 if is_max30 else 1.8
        alpha = 0.70 if is_max30 else 0.45
        ax.plot(
            [0.0, x1], [0.0, y1],
            color=color, lw=lw, alpha=alpha, zorder=3,
            solid_capstyle="round",
        )
        # Small tick at the ring hit
        ax.plot([x1], [y1], "o", color=color, ms=4.5 if is_max30 else 3.5,
                alpha=0.9, zorder=5, markeredgecolor="none")

    # Angle arc at θ=30° from +x axis up to the ray (makes the measure obvious)
    arc_r = 0.38
    phi = np.linspace(0.0, np.deg2rad(29.5), 60)
    ax.plot(
        arc_r * np.cos(phi), arc_r * np.sin(phi),
        color=C_MAX_30, lw=1.6, alpha=0.85, zorder=3,
    )
    ax.text(
        arc_r * np.cos(np.deg2rad(15.0)) + 0.06,
        arc_r * np.sin(np.deg2rad(15.0)) + 0.04,
        r"$\theta$",
        color=C_MAX_30, fontsize=10, fontweight="bold",
        ha="left", va="bottom", zorder=6,
    )


def panel_A_model_frame(ax, geo):
    """Side view: exactly four thick foreshortened loop ribbons about x."""
    ax.set_title(
        "A. Model frame (no Galactic yet)\n"
        r"axisymmetric about $x$ = CoE$\leftrightarrow$anti-CoE"
        "\n"
        r"two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$)",
        fontsize=10,
        pad=8,
    )
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color="#7f8c8d", lw=1.2, zorder=1)
    ax.plot(0, 0, "ko", ms=6, zorder=5)
    ax.text(0.08, -0.12, "observer", fontsize=8, color="k", ha="left")

    ax.annotate(
        "",
        xy=(1.15, 0),
        xytext=(-1.15, 0),
        arrowprops=dict(arrowstyle="<->", color=C_AXIS, lw=2.5),
        zorder=4,
    )
    ax.text(1.18, 0.08, "anti-CoE\n(+x)", color=C_ANTI, fontsize=8, ha="left", va="bottom")
    ax.text(-1.18, 0.08, "CoE\n(−x)", color=C_COE, fontsize=8, ha="right", va="bottom")

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        draw_panelA_ribbon(ax, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        draw_panelA_ribbon(ax, lo, hi, mid, color)

    # Illustrate θ from +x (anti-CoE) to each loop (Chris markup)
    draw_panelA_theta_rays(ax, geo)

    label_r = 1.30
    for mid, color in list(zip(geo["theta_max_from_plus_x_deg"], COLORS_MAX)) + list(
        zip(geo["theta_min_from_plus_x_deg"], COLORS_MIN)
    ):
        t = np.deg2rad(mid)
        ax.text(
            label_r * np.cos(t),
            label_r * np.sin(t),
            rf"$\theta\!=\!{mid:g}^\circ$",
            color=color, fontsize=7.5, ha="center", va="center", fontweight="bold",
        )

    ax.text(
        0.0,
        -1.38,
        (
            r"four transparent highlight bands about $x$  ·  each $18^\circ$ end-to-end ($\pm9^\circ$ = 5% of $360^\circ$)"
            + "\n"
            + r"(after-2$\ell$ $A(\ell)$: two MAX + two MIN; one mid-stroke each)"
            + "\n"
            + r"$y,z$ free (roll) — model does not pick celestial equator"
        ),
        fontsize=7.5,
        ha="center",
        va="top",
        style="italic",
        color="#34495e",
    )
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")


def _ring_lb(D, mu, n_pts=721):
    pts = small_circle_about_axis(D, mu, n=n_pts)
    rl = np.array([n_to_lb(p)[0] for p in pts])
    rb = np.array([n_to_lb(p)[1] for p in pts])
    return rl, rb


def _ring_lonlat_rad(D, mu, n_pts=361):
    """Small-circle → Mollweide display (lon_disp, lat) in radians."""
    pts = small_circle_about_axis(D, float(mu), n=n_pts)
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    lon = np.deg2rad(lon_disp_from_ell(ls))
    lat = np.deg2rad(bs)
    return lon, lat


def fill_mollweide_ribbon(ax, D, theta_lo, theta_hi, theta_mid, color, n_pts=721):
    """
    Mollweide highlight band (±9° edges, base alpha 0.26) + mid-stroke
    (lw 2.6, alpha 0.95), both with the (1+cos psi)/2 opacity gradient.
    Elements crossing the antimeridian seam are dropped.
    """
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))

    def lonlat(mu, n):
        pts, sc = small_circle_grad(D, mu, theta_mid, n=n)
        ls = np.array([n_to_lb(p)[0] for p in pts])
        bs = np.array([n_to_lb(p)[1] for p in pts])
        return np.deg2rad(lon_disp_from_ell(ls)), np.deg2rad(bs), sc

    lon_a, lat_a, sc = lonlat(mu_a, n_pts)
    lon_b, lat_b, _ = lonlat(mu_b, n_pts)
    gap = 0.9 * np.pi
    quad = np.stack([lon_a[:-1], lon_a[1:], lon_b[1:], lon_b[:-1]])
    valid = (quad.max(axis=0) - quad.min(axis=0)) <= gap
    runs = quantized_runs(segment_scale(sc), valid)
    if runs:
        ax.add_collection(
            PolyCollection(
                _strip_polys_2d(lon_a, lat_a, lon_b, lat_b, runs),
                facecolors=rgba_array(color, [0.26 * r[2] for r in runs]),
                edgecolors="none", zorder=2, linewidths=0,
            )
        )
    lon_m, lat_m, sc_m = lonlat(mu_mid, max(n_pts, 721))
    valid = np.abs(np.diff(lon_m)) <= gap
    runs = quantized_runs(segment_scale(sc_m), valid)
    if runs:
        ax.add_collection(
            LineCollection(
                _polylines_2d(lon_m, lat_m, runs),
                colors=rgba_array(color, [0.95 * r[2] for r in runs]),
                linewidths=2.6, zorder=4, capstyle="butt", joinstyle="round",
            )
        )


def panel_B_mollweide(ax, geo):
    ax.set_title(
        r"B. Migkas $A(\ell)$ loops on sky (anti-CoE at $(0^\circ,+15^\circ)$)"
        "\n"
        r"(Galactic Mollweide; CoE axis tilted $+15^\circ$; $\ell$ increases LEFT; "
        r"two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$))",
        fontsize=10,
        pad=8,
    )
    for b in [-60, -30, 0, 30, 60]:
        ell = np.linspace(0, 360, 361)
        plot_lb_curve_mollweide(ax, ell, np.full_like(ell, b), color="#bdc3c7", lw=0.5, zorder=0)
    for ell0 in range(0, 360, 30):
        b = np.linspace(-90, 90, 181)
        plot_lb_curve_mollweide(ax, np.full_like(b, ell0), b, color="#bdc3c7", lw=0.5, zorder=0)

    D = geo["D"]
    ls, bs = great_circle_through_poles(D, n=721)
    plot_lb_curve_mollweide(ax, ls, bs, color=C_AXIS, lw=2.0, zorder=2)

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        fill_mollweide_ribbon(ax, D, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        fill_mollweide_ribbon(ax, D, lo, hi, mid, color)

    ell_eq, b_eq = celestial_equator_lb(721)
    plot_lb_curve_mollweide(ax, ell_eq, b_eq, color=C_EQ, lw=1.8, ls="--", zorder=2)

    ax.plot(
        np.deg2rad(lon_disp_from_ell(geo["anti_lb"][0])),
        np.deg2rad(geo["anti_lb"][1]),
        marker="o", ms=8, color=C_ANTI,
        markerfacecolor=C_ANTI, markeredgecolor="k", markeredgewidth=0.6,
        zorder=6, linestyle="None",
    )
    ax.plot(
        np.deg2rad(lon_disp_from_ell(geo["coe_lb"][0])),
        np.deg2rad(geo["coe_lb"][1]),
        marker="*", ms=14, color=C_COE,
        markerfacecolor=C_COE, markeredgecolor="k", markeredgewidth=0.6,
        zorder=6, linestyle="None",
    )
    for (l, b), color in zip(MIGKAS_MAXES, COLORS_MAX):
        ax.plot(
            np.deg2rad(lon_disp_from_ell(l)),
            np.deg2rad(b),
            marker="^",
            ms=10,
            color=color,
            markeredgecolor="k",
            markeredgewidth=0.7,
            zorder=7,
            linestyle="None",
        )
    for (l, b), color in zip(MIGKAS_MINS, COLORS_MIN):
        ax.plot(
            np.deg2rad(lon_disp_from_ell(l)),
            np.deg2rad(b),
            marker="s",
            ms=10,
            color=color,
            markeredgecolor="k",
            markeredgewidth=0.7,
            zorder=7,
            linestyle="None",
        )

    xticks = np.deg2rad(lon_disp_from_ell(np.array([180, 90, 0, 270])))
    ax.set_xticks(xticks)
    ax.set_xticklabels([r"$180^\circ$", r"$90^\circ$", r"$0^\circ$", r"$270^\circ$"], fontsize=8)
    ax.set_yticks(np.deg2rad([-60, -30, 0, 30, 60]))
    ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$+30^\circ$", r"$+60^\circ$"], fontsize=8)
    ax.grid(False)
    ax.set_xlabel(r"Galactic $\ell$ (increases LEFT)", fontsize=8)
    ax.set_ylabel(r"Galactic $b$", fontsize=8)


def fill_sphere_ribbon(ax, D, theta_lo, theta_hi, theta_mid, color, n_pts=160):
    """
    3D highlight belt (±9° edges, base alpha 0.25) + mid-stroke (lw 2.6,
    alpha 0.95), both with the (1+cos psi)/2 opacity gradient.
    """
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    n_g = 361
    p0, sc = small_circle_grad(D, float(mu_a), theta_mid, n=n_g)
    p1, _ = small_circle_grad(D, float(mu_b), theta_mid, n=n_g)
    seg = segment_scale(sc)
    idx = np.where(seg > 1e-4)[0]
    verts = [[p0[j], p0[j + 1], p1[j + 1], p1[j]] for j in idx]
    coll = Poly3DCollection(
        verts, facecolors=rgba_array(color, 0.25 * seg[idx]),
        edgecolors="none", linewidths=0.0, zorder=2,
    )
    ax.add_collection3d(coll)
    pts, scm = small_circle_grad(D, mu_mid, theta_mid, n=n_g)
    runs = quantized_runs(segment_scale(scm))
    ax.add_collection3d(
        Line3DCollection(
            [pts[i0:i1 + 1] for i0, i1, _ in runs],
            colors=rgba_array(color, [0.95 * r[2] for r in runs]),
            linewidths=2.6, zorder=5, capstyle="butt", joinstyle="round",
        )
    )


def panel_C_sphere(ax, geo):
    th_max = geo["theta_max_from_plus_x_deg"]
    th_min = geo["theta_min_from_plus_x_deg"]
    ax.set_title(
        r"C. Sphere 3D — two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$)"
        "\n"
        r"(axis line emphasized; celestial equator dashed)",
        fontsize=10,
        pad=6,
    )
    draw_wire_sphere(ax)
    D = geo["D"]
    ell_eq, b_eq = celestial_equator_lb(400)
    eq_xyz = np.array([lb_to_n(l, b) for l, b in zip(ell_eq, b_eq)])
    ax.plot(eq_xyz[:, 0], eq_xyz[:, 1], eq_xyz[:, 2], color=C_EQ, lw=1.8, ls="--", zorder=3)
    ax.plot(
        [-D[0], D[0]],
        [-D[1], D[1]],
        [-D[2], D[2]],
        color=C_AXIS,
        lw=3.0,
        zorder=5,
    )
    ax.scatter(*(-D), color=C_COE, s=120, marker="*", edgecolors="k", linewidths=0.6,
               depthshade=False, zorder=8)
    ax.scatter(*(D), color=C_ANTI, s=80, edgecolors="k", linewidths=0.6, depthshade=False, zorder=8)
    ax.text(*(-D * 1.22), "CoE", color=C_COE, fontsize=8, fontweight="bold")
    ax.text(*(D * 1.18), "anti-CoE", color=C_ANTI, fontsize=8, fontweight="bold")
    ax.scatter([0], [0], [0], color="k", s=25, depthshade=False, zorder=9)

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        fill_sphere_ribbon(ax, D, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        fill_sphere_ribbon(ax, D, lo, hi, mid, color)

    for n, color in zip(geo["n_maxes"], COLORS_MAX):
        ax.scatter(
            [n[0]], [n[1]], [n[2]],
            c=[color], marker="^", s=70,
            edgecolors="k", linewidths=0.7, depthshade=False, zorder=10,
        )
        ax.plot([0, n[0]], [0, n[1]], [0, n[2]], color=color, lw=0.8, ls=":", alpha=0.7)
    for n, color in zip(geo["n_mins"], COLORS_MIN):
        ax.scatter(
            [n[0]], [n[1]], [n[2]],
            c=[color], marker="s", s=60,
            edgecolors="k", linewidths=0.7, depthshade=False, zorder=10,
        )
        ax.plot([0, n[0]], [0, n[1]], [0, n[2]], color=color, lw=0.8, ls=":", alpha=0.7)

    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.set_zlim(-1.25, 1.25)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=18, azim=-55)
    ax.set_axis_off()
    ax.text2D(
        0.5,
        -0.02,
        rf"MAX $\theta\!=\!{th_max[0]:g}^\circ,{th_max[1]:g}^\circ$; "
        rf"MIN $\theta\!=\!{th_min[0]:g}^\circ,{th_min[1]:g}^\circ$  "
        rf"($18^\circ$ end-to-end = $\pm9^\circ$ = 5% of $360^\circ$)",
        transform=ax.transAxes,
        ha="center",
        fontsize=7.5,
        color="#34495e",
    )


def main():
    st = migkas_scan_stops()
    geo = compute_geometry(st)

    th_max = geo["theta_max_from_plus_x_deg"]
    th_min = geo["theta_min_from_plus_x_deg"]
    max_bands = geo["max_bands_deg"]
    min_bands = geo["min_bands_deg"]

    fig = plt.figure(figsize=(14.5, 5.6), dpi=160)
    fig.suptitle(
        r"Pure CoE axisymmetry: two MAX + two MIN highlight bands  (gradient: darkest at $\theta$ ray, transparent opposite)"
        "\n"
        rf"(each $18^\circ$ end-to-end = $\pm9^\circ$ = 5% of $360^\circ$; "
        rf"MAX $\theta\!=\!{th_max[0]:g}^\circ/{th_max[1]:g}^\circ$; "
        rf"MIN $\theta\!=\!{th_min[0]:g}^\circ/{th_min[1]:g}^\circ$; "
        r"anti-CoE $(0^\circ,+15^\circ)$, CoE axis tilted $+15^\circ$)",
        fontsize=11,
        fontweight="bold",
        y=1.02,
    )

    axA = fig.add_subplot(1, 3, 1)
    panel_A_model_frame(axA, geo)

    axB = fig.add_subplot(1, 3, 2, projection="mollweide")
    panel_B_mollweide(axB, geo)

    axC = fig.add_subplot(1, 3, 3, projection="3d")
    panel_C_sphere(axC, geo)

    handles = [
        Line2D([0], [0], color=C_AXIS, lw=2.5, label=r"CoE$\leftrightarrow$anti-CoE axis"),
        Patch(
            facecolor=C_MAX_30, alpha=0.45, edgecolor=C_MAX_30,
            label=(
                rf"MAX $\theta\!=\!{th_max[0]:g}^\circ$ band "
                rf"[{max_bands[0][0]:g},{max_bands[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$=5% of $360^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MAX_185, alpha=0.45, edgecolor=C_MAX_185,
            label=(
                rf"MAX $\theta\!=\!{th_max[1]:g}^\circ$ band "
                rf"[{max_bands[1][0]:g},{max_bands[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_115, alpha=0.45, edgecolor=C_MIN_115,
            label=(
                rf"MIN $\theta\!=\!{th_min[0]:g}^\circ$ band "
                rf"[{min_bands[0][0]:g},{min_bands[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_265, alpha=0.45, edgecolor=C_MIN_265,
            label=(
                rf"MIN $\theta\!=\!{th_min[1]:g}^\circ$ band "
                rf"[{min_bands[1][0]:g},{min_bands[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Line2D([0], [0], color=C_EQ, lw=1.8, ls="--", label="celestial equator (Dec$=0$)"),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_30, markeredgecolor="k",
            ms=9,
            label=rf"MAX $({MIGKAS_MAXES[0][0]:.0f}^\circ,{MIGKAS_MAXES[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_185, markeredgecolor="k",
            ms=9,
            label=rf"MAX $({MIGKAS_MAXES[1][0]:.0f}^\circ,{MIGKAS_MAXES[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_115, markeredgecolor="k",
            ms=8,
            label=rf"MIN $({MIGKAS_MINS[0][0]:.0f}^\circ,{MIGKAS_MINS[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_265, markeredgecolor="k",
            ms=8,
            label=rf"MIN $({MIGKAS_MINS[1][0]:.0f}^\circ,{MIGKAS_MINS[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="o", color="w", markerfacecolor=C_ANTI, markeredgecolor="k",
            ms=8, label=rf"anti-CoE $({geo['anti_lb'][0]:.1f}^\circ,{geo['anti_lb'][1]:+.1f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="*", color="w", markerfacecolor=C_COE, markeredgecolor="k",
            ms=11, label=rf"CoE $({geo['coe_lb'][0]:.1f}^\circ,{geo['coe_lb'][1]:+.1f}^\circ)$",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=6,
        fontsize=7.0,
        frameon=True,
        fancybox=True,
        bbox_to_anchor=(0.5, -0.10),
    )

    fig.tight_layout()
    # --- Chris 2026-09-26: panel C sphere drawn at the same on-page diameter as
    # panel A's model-frame circle; panel C axes matched to A's vertical
    # position/height. Nothing else changes.
    from mpl_toolkits.mplot3d import proj3d
    pa, pc = axA.get_position(), axC.get_position()
    axC.set_position([pc.x0, pa.y0, pc.width, pa.height])

    def _diam_A():
        fig.canvas.draw()
        p0 = axA.transData.transform((-1.0, 0.0)); p1 = axA.transData.transform((1.0, 0.0))
        return float(p1[0] - p0[0])

    def _diam_C():
        fig.canvas.draw()
        u = np.random.default_rng(0).normal(size=(20000, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]
        xs, ys, _ = proj3d.proj_transform(u[:, 0], u[:, 1], u[:, 2], axC.get_proj())
        d = axC.transData.transform(np.column_stack([xs, ys]))
        return float(d[:, 0].max() - d[:, 0].min()), float(d[:, 1].max() - d[:, 1].min())

    zoom = 1.0
    for _ in range(8):
        dA = _diam_A(); dCx, dCy = _diam_C()
        zoom *= dA / dCx
        axC.set_box_aspect((1, 1, 1), zoom=zoom)
    dA = _diam_A(); dCx, dCy = _diam_C()
    print(f"PANEL_DIAM A={dA:.1f}px C={dCx:.1f}x{dCy:.1f}px zoom={zoom:.3f}")
    fig.savefig(PNG_PATH, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    payload = {
        "png": str(PNG_PATH),
        "script": str(OUT_DIR / "make_coe_rings_migkas_pin_axis_matched_grad.py"),
        "migkas_max_lb_deg": [list(lb) for lb in MIGKAS_MAXES],
        "migkas_min_lb_deg": [list(lb) for lb in MIGKAS_MINS],
        "strip_theta_max_deg": list(STRIP_THETA_MAX),
        "strip_theta_min_deg": list(STRIP_THETA_MIN),
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "band_total_width_deg": BAND_TOTAL,
        "band_width_note": "5% of 360° = 18° end-to-end = ±9° about each MAX/MIN center",
        "max_bands_deg": [list(b) for b in max_bands],
        "min_bands_deg": [list(b) for b in min_bands],
        "band_intervals": {
            "MAX_30": [float(max_bands[0][0]), float(max_bands[0][1])],
            "MAX_185": [float(max_bands[1][0]), float(max_bands[1][1])],
            "MIN_115": [float(min_bands[0][0]), float(min_bands[0][1])],
            "MIN_265": [float(min_bands[1][0]), float(min_bands[1][1])],
        },
        "scan_stop_source": st["source"],
        "theta_max_from_plus_x_deg": [round(t, 4) for t in th_max],
        "theta_min_from_plus_x_deg": [round(t, 4) for t in th_min],
        "mu_max_bands": [[round(a, 8), round(b, 8)] for a, b in geo["mu_max_bands"]],
        "mu_min_bands": [[round(a, 8), round(b, 8)] for a, b in geo["mu_min_bands"]],
        "CoE_lb_deg": [round(geo["coe_lb"][0], 4), round(geo["coe_lb"][1], 4)],
        "anti_CoE_lb_deg": [round(geo["anti_lb"][0], 4), round(geo["anti_lb"][1], 4)],
        "D_hat_galactic_xyz": geo["D"].tolist(),
        "angular_sep_migkas_max0_min1_deg": round(geo["sep_deg"], 6),
        "axis_tilt_b_deg": 15,
        "axis_tilt_reason": (
            "Chris 2026-09-22 — markers fixed, rings rotate up 30°"
        ),
        "orientation": (
            "Galactic orientation: anti-CoE at (ell=0°, b=+15°), CoE at "
            "(ell=180°, b=−15°); +D = anti-CoE = +x; axis tilted +15° "
            "(Chris 2026-09-22 — markers fixed, rings rotate up 30°)"
        ),
        "drawing": (
            "Exactly 4 features. Each = ONE smooth transparent highlight band "
            "between θ−9° and θ+9° (18° end-to-end = 5% of 360°) + ONE solid "
            "mid-stroke at the exact center θ. No multi-stroke wash / hairball. "
            "Four distinct colors (two reds + two blues); no orange; no 283° "
            "companion. Panel A: clean vertical strip clipped to unit circle; theta rays from observer to each loop (emphasized at MAX 30°). Mollweide "
            "seam-split honestly. Sphere: single Poly3D strip (no n_wash stack)."
        ),
        "note": (
            "Stops are after-2ℓ A(ℓ) extrema from best_fit_2ell_multiplier.npz "
            f"A_after: two MAX at θ={th_max[0]:g}°/{th_max[1]:g}° and "
            f"two MIN at θ={th_min[0]:g}°/{th_min[1]:g}°, each drawn as "
            f"one {BAND_FULL:.0f}° end-to-end highlight (±{BAND_HALF:.0f}° = 5% of 360°). "
            "Rings use μ=cos(θ); sky markers are map extrema (17°,−9°) and "
            "(274°,−22°) plus strip stand-ins (185°,0°) and (115°,0°)."
        ),
        "plot_convention": "Migkas Galactic Mollweide: lon_disp=-ell_wrapped (ell increases LEFT)",
        "legend_summary": (
            "two MAX + two MIN highlight bands, each 18° end-to-end "
            "(±9° = 5% of 360°); four distinct colors"
        ),
        "colors": {
            "MAX_30": C_MAX_30,
            "MAX_185": C_MAX_185,
            "MIN_115": C_MIN_115,
            "MIN_265": C_MIN_265,
            "no_orange": True,
        },
    }
    JSON_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    print(f"\nWrote {PNG_PATH}")
    print(f"Wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
```

#### `make_coe_rings_migkas_pin_axis_matched_grad_rot_fit.py`

```python
# source: /workspace/cosmology/figures/make_coe_rings_migkas_pin_axis_matched_grad_rot_fit.py
#!/usr/bin/env python3
"""
CoE rings + Migkas A(ell) extrema loops in the Galactic frame.

Concept:
  Pure CoE axisymmetry → constant A on small-circle rings about x = CoE↔anti-CoE.
  Scan A vs angle from +x (anti-CoE) around 360°. After-2ℓ A(ℓ) has TWO maxima
  and TWO minima (from best_fit_2ell_multiplier.npz A_after):
    MAX at θ = 30° and 185° from +x,
    MIN at θ = 115° and 265° from +x.

  Exactly FOUR loops — each ONE transparent highlight band, 18° end-to-end
  (= ±9° = 5% of 360°) about the center θ:
    MAX: 30±9 → [21°, 39°] and 185±9 → [176°, 194°]
    MIN: 115±9 → [106°, 124°] and 265±9 → [256°, 274°]
  Drawing: ONE smooth transparent filled ribbon between θ−9° and θ+9°
  + ONE solid mid-stroke at the exact center θ. No multi-stroke wash,
  no dense linspace hairball, no dual strong boundaries.
  Four distinct colors (two reds + two blues). No orange, no 283° companion.

  Galactic orientation: anti-CoE (+x) at ℓ=0°, b=+15°; CoE (−x) at ℓ=180°, b=−15°
  (axis tilted +15° so features near ℓ=0 move north on Mollweide; Migkas markers fixed).
  Sky markers (ℓ, b):
    Maxes: map extrema (17°, −9°) for ~30° max; strip (185°, 0°) for second max.
    Mins:  strip (115°, 0°) for ~115° min; map extrema (274°, −22°) for ~265° min.

Outputs:
  coe_rings_migkas_pin_axis.png
  coe_rings_migkas_pin_axis.json
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle, Patch, Polygon
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

OUT_DIR = Path(__file__).resolve().parent
PNG_PATH = Path("/workspace/cosmology/papers/spreadsheet_model/loop_threepanel_matched_grad_rotfit_rot81.png")
JSON_PATH = Path("/workspace/cosmology/papers/spreadsheet_model/loop_threepanel_matched_grad_rotfit_rot81.json")

# 5% of 360° = 18° end-to-end → half-width ±9° about each center
BAND_HALF = 9.0
BAND_FULL = 2.0 * BAND_HALF  # 18° = 5% of 360°
BAND_TOTAL = BAND_FULL  # alias
STRIP_THETA_MAX = [29.5, 187.5]
STRIP_THETA_MIN = [122.6, 264.4]

MIGKAS_MAXES = [(17.0, -9.0), (185.0, 0.0)]
MIGKAS_MINS = [(115.0, 0.0), (274.0, -22.0)]

C_COE = "#6b3fa0"
C_ANTI = "#8e1b1b"  # dark red marker; distinct from loop reds
# Four distinct loop colors: two reds (MAX) + two blues (MIN); no orange/yellow
C_MAX_30 = "#0097b2"    # cyan — MAX at 29.5°
C_MAX_185 = "#e0107a"   # magenta — MAX at 187.5°
C_MIN_115 = "#6fa800"   # green — MIN at 122.6°
C_MIN_265 = "#2166d4"   # blue — MIN at 264.4°
C_AXIS = "#5b2c6f"
C_EQ = "#1e8449"

# Per-loop color maps (order matches STRIP_THETA_MAX / STRIP_THETA_MIN)
COLORS_MAX = [C_MAX_30, C_MAX_185]
COLORS_MIN = [C_MIN_115, C_MIN_265]


def wrap_lon_deg(lon: np.ndarray) -> np.ndarray:
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def lon_disp_from_ell(ell_deg: np.ndarray) -> np.ndarray:
    """Migkas Galactic: ℓ increases LEFT → lon_disp = −ℓ_wrapped."""
    return -wrap_lon_deg(ell_deg)


def lb_to_n(l_deg: float, b_deg: float) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    return np.array(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], dtype=float
    )


def n_to_lb(n: np.ndarray) -> tuple[float, float]:
    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    b = float(np.rad2deg(np.arcsin(np.clip(n[2], -1.0, 1.0))))
    l = float(np.rad2deg(np.arctan2(n[1], n[0])) % 360.0)
    return l, b


def celestial_equator_lb(n: int = 720) -> tuple[np.ndarray, np.ndarray]:
    """Dec=0 (ICRS) as Galactic (ℓ, b). IAU J2000 NGP convention, no astropy."""
    ra_ngp = np.deg2rad(192.85948)
    dec_ngp = np.deg2rad(27.12825)
    l_ncp = np.deg2rad(122.93192)
    ra = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    dec = 0.0
    sin_b = np.sin(dec) * np.sin(dec_ngp) + np.cos(dec) * np.cos(dec_ngp) * np.cos(ra - ra_ngp)
    sin_b = np.clip(sin_b, -1.0, 1.0)
    b = np.arcsin(sin_b)
    y = np.cos(dec) * np.sin(ra - ra_ngp)
    x = np.sin(dec) * np.cos(dec_ngp) - np.cos(dec) * np.sin(dec_ngp) * np.cos(ra - ra_ngp)
    l = (l_ncp - np.arctan2(y, x)) % (2.0 * np.pi)
    return np.rad2deg(l), np.rad2deg(b)


def great_circle_through_poles(D: np.ndarray, n: int = 361) -> tuple[np.ndarray, np.ndarray]:
    """Great circle through ±D (any meridian of the CoE axis)."""
    a = D / np.linalg.norm(D)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([1.0, 0.0, 0.0])
    e2 = np.cross(a, tmp)
    e2 /= np.linalg.norm(e2)
    t = np.linspace(0.0, 2.0 * np.pi, n)
    pts = np.outer(np.cos(t), a) + np.outer(np.sin(t), e2)
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    return ls, bs


def small_circle_about_axis(axis: np.ndarray, mu: float, n: int = 361) -> np.ndarray:
    """Unit vectors on the small circle n·axis = mu (|mu|<1)."""
    a = axis / np.linalg.norm(axis)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    tmp = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(a, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    t = np.linspace(0.0, 2.0 * np.pi, n)
    return (mu * a)[None, :] + rho * (np.outer(np.cos(t), e1) + np.outer(np.sin(t), e2))


def split_lon_segments(lon_rad: np.ndarray, lat_rad: np.ndarray, gap: float = np.pi):
    """Split polyline where Mollweide lon jumps (wrap discontinuity)."""
    d = np.diff(lon_rad)
    cuts = np.where(np.abs(d) > gap)[0] + 1
    idx = np.concatenate([[0], cuts, [len(lon_rad)]])
    segs = []
    for i0, i1 in zip(idx[:-1], idx[1:]):
        if i1 - i0 >= 2:
            segs.append((lon_rad[i0:i1], lat_rad[i0:i1]))
    return segs


def plot_lb_curve_mollweide(ax, l_deg, b_deg, **kwargs):
    lon = np.deg2rad(lon_disp_from_ell(l_deg))
    lat = np.deg2rad(b_deg)
    for xr, yr in split_lon_segments(lon, lat):
        ax.plot(xr, yr, **kwargs)


def band_edges(centers: list[float], half: float = BAND_HALF) -> list[tuple[float, float]]:
    """Return (theta_lo, theta_hi) for each band center (±half about center)."""
    return [(c - half, c + half) for c in centers]


def mu_pair_for_band(lo: float, hi: float) -> tuple[float, float]:
    """Edge μ=cos(θ) for a scan-angle band (two boundary cones only)."""
    return float(np.cos(np.deg2rad(lo))), float(np.cos(np.deg2rad(hi)))


def migkas_scan_stops() -> dict:
    th_max = [float(t) for t in STRIP_THETA_MAX]
    th_min = [float(t) for t in STRIP_THETA_MIN]
    max_bands = band_edges(th_max)
    min_bands = band_edges(th_min)
    return {
        "theta_max_from_plus_x_deg": th_max,
        "theta_min_from_plus_x_deg": th_min,
        "max_bands_deg": max_bands,
        "min_bands_deg": min_bands,
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "band_total_width_deg": BAND_TOTAL,
        "band_width_note": "5% of 360° = 18° end-to-end = ±9° about each MAX/MIN center",
        "mu_max_centers": [float(np.cos(np.deg2rad(t))) for t in th_max],
        "mu_min_centers": [float(np.cos(np.deg2rad(t))) for t in th_min],
        "mu_max_bands": [mu_pair_for_band(lo, hi) for lo, hi in max_bands],
        "mu_min_bands": [mu_pair_for_band(lo, hi) for lo, hi in min_bands],
        "source": (
            "best_fit_2ell_multiplier.npz A_after extrema (two max, two min); "
            f"each drawn as one transparent {BAND_FULL:.0f}° end-to-end highlight "
            f"(±{BAND_HALF:.0f}° = 5% of 360°) + one solid mid-stroke"
        ),
    }


def compute_geometry(st: dict) -> dict:
    anti_l, anti_b = 0.0, 15.0
    coe_l, coe_b = 180.0, -15.0

    n_coe = lb_to_n(coe_l, coe_b)
    n_anti = lb_to_n(anti_l, anti_b)
    D = n_anti  # +D = anti-CoE = +x

    th_max = [float(t) for t in st["theta_max_from_plus_x_deg"]]
    th_min = [float(t) for t in st["theta_min_from_plus_x_deg"]]
    max_bands = list(st["max_bands_deg"])
    min_bands = list(st["min_bands_deg"])

    n_maxes = [lb_to_n(*lb) for lb in MIGKAS_MAXES]
    n_mins = [lb_to_n(*lb) for lb in MIGKAS_MINS]

    sep = float(
        np.rad2deg(
            np.arccos(np.clip(np.dot(n_maxes[0], n_mins[1]), -1.0, 1.0))
        )
    )

    return {
        "D": D,
        "n_coe": n_coe,
        "n_anti": n_anti,
        "n_maxes": n_maxes,
        "n_mins": n_mins,
        "anti_lb": (anti_l, anti_b),
        "coe_lb": (coe_l, coe_b),
        "theta_max_from_plus_x_deg": th_max,
        "theta_min_from_plus_x_deg": th_min,
        "max_bands_deg": max_bands,
        "min_bands_deg": min_bands,
        "mu_max_bands": list(st["mu_max_bands"]),
        "mu_min_bands": list(st["mu_min_bands"]),
        "sep_deg": sep,
        "band_half_width_deg": st["band_half_width_deg"],
    }


def draw_wire_sphere(ax, n_mer=12, n_par=6, **kwargs):
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(-np.pi / 2, np.pi / 2, 40)
    xs = np.outer(np.cos(u), np.cos(v))
    ys = np.outer(np.sin(u), np.cos(v))
    zs = np.outer(np.ones_like(u), np.sin(v))
    ax.plot_surface(xs, ys, zs, color="#dfe6e9", alpha=0.15, linewidth=0, shade=False)
    kw = dict(color="#95a5a6", lw=0.4, alpha=0.45)
    kw.update(kwargs)
    for i in range(n_mer):
        phi = 2 * np.pi * i / n_mer
        t = np.linspace(-np.pi / 2, np.pi / 2, 80)
        ax.plot(np.cos(t) * np.cos(phi), np.cos(t) * np.sin(phi), np.sin(t), **kw)
    for j in range(1, n_par):
        b = -np.pi / 2 + np.pi * j / n_par
        t = np.linspace(0, 2 * np.pi, 80)
        r = np.cos(b)
        ax.plot(r * np.cos(t), r * np.sin(t), np.full_like(t, np.sin(b)), **kw)


def foreshortened_ellipse_xy(mu: float, n: int = 200):
    """Side-view foreshortened ellipse of a 3D ring at n·x = mu."""
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    phi = np.linspace(0, 2 * np.pi, n)
    ex = mu + 0.045 * rho * np.cos(phi)
    ey = rho * np.sin(phi)
    return ex, ey, rho



# ---------------------------------------------------------------------------
# Opacity gradient along each loop (Chris 2026-09-26), same method as
# make_coe_rings_gradient_loops.py: loop drawn as runs of short elements with
# per-element alpha = base_alpha * (1 + cos psi)/2, psi = ring azimuth from the
# θ-ray point (ring point most along sign(sin θ_mid)·Ŷ; Ŷ = Galactic +z
# projected ⊥ D, i.e. model-frame "up"; no roll). Darkest at the θ-ray point,
# fully transparent at the diametrically opposite point.
# ---------------------------------------------------------------------------
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def theta_ray_sign(theta_mid_deg):
    return 1.0 if np.sin(np.deg2rad(float(theta_mid_deg))) >= 0 else -1.0


ROLL_BC_DEG = 80.8  # roll about CoE axis, applied only via small_circle_grad (Panels B/C)


def ring_basis(D):
    a = np.asarray(D, dtype=float)
    a = a / np.linalg.norm(a)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, tmp)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    y = tmp - np.dot(tmp, a) * a
    y /= np.linalg.norm(y)
    z = np.cross(a, y)
    z /= np.linalg.norm(z)
    # +45° roll (Panels B/C only), same convention as coe_half_loop_geometry.ring_basis
    # / earlier "+45° CCW about CoE" figures: Y' = c Y - s Z, Z' = s Y + c Z
    r = np.deg2rad(float(ROLL_BC_DEG))
    c, s_ = np.cos(r), np.sin(r)
    return c * y - s_ * z, s_ * y + c * z


def grad_phi(theta_mid_deg, n):
    """phi from antipode of θ-ray point round to antipode (psi = -180..+180)."""
    phi0 = 0.5 * np.pi if theta_ray_sign(theta_mid_deg) > 0 else 1.5 * np.pi
    return phi0, phi0 + np.linspace(-np.pi, np.pi, n)


def grad_scale(phi, phi0):
    return 0.5 * (1.0 + np.cos(np.asarray(phi) - phi0))


def small_circle_grad(D, mu, theta_mid_deg, n=721):
    a = np.asarray(D, dtype=float)
    a = a / np.linalg.norm(a)
    mu = float(np.clip(mu, -1.0 + 1e-9, 1.0 - 1e-9))
    rho = np.sqrt(max(0.0, 1.0 - mu * mu))
    y, z = ring_basis(a)
    phi0, phi = grad_phi(theta_mid_deg, n)
    pts = (mu * a)[None, :] + rho * (np.outer(np.sin(phi), y) + np.outer(np.cos(phi), z))
    return pts, grad_scale(phi, phi0)


def segment_scale(sc):
    sc = np.asarray(sc, dtype=float)
    return 0.5 * (sc[:-1] + sc[1:])


def rgba_array(color, alphas):
    r, g, b, _ = to_rgba(color)
    al = np.clip(np.asarray(alphas, dtype=float), 0.0, 1.0)
    out = np.empty((al.size, 4))
    out[:, 0], out[:, 1], out[:, 2], out[:, 3] = r, g, b, al
    return out


def quantized_runs(seg_sc, valid=None, n_levels=64):
    s = np.asarray(seg_sc, dtype=float)
    lev = np.round(np.clip(s, 0.0, 1.0) * n_levels) / n_levels
    ok = lev > 0 if valid is None else (lev > 0) & np.asarray(valid, dtype=bool)
    runs, i, n = [], 0, s.size
    while i < n:
        if not ok[i]:
            i += 1
            continue
        j = i + 1
        while j < n and ok[j] and lev[j] == lev[i]:
            j += 1
        runs.append((i, j, float(lev[i])))
        i = j
    return runs


def _strip_polys_2d(xa, ya, xb, yb, runs):
    return [
        np.column_stack([
            np.concatenate([xa[i0:i1 + 1], xb[i0:i1 + 1][::-1]]),
            np.concatenate([ya[i0:i1 + 1], yb[i0:i1 + 1][::-1]]),
        ])
        for i0, i1, _ in runs
    ]


def _polylines_2d(x, y, runs):
    return [np.column_stack([x[i0:i1 + 1], y[i0:i1 + 1]]) for i0, i1, _ in runs]


def draw_panelA_ribbon(ax, theta_lo: float, theta_hi: float, theta_mid: float, color: str):
    """
    Panel-A highlight (edge-on chord): same vertical strip between
    x=cos(θ_hi)..cos(θ_lo) clipped to the unit circle, base alpha 0.28, and the
    same foreshortened mid-ellipse (lw 2.8, alpha 0.95) — now with the opacity
    gradient along the chord: alpha *= (1 + sgn·y/ρ(x))/2 = (1+cos psi)/2,
    darkest at the θ-ray end, transparent at the opposite end.
    """
    mu_lo = float(np.cos(np.deg2rad(theta_lo)))
    mu_hi = float(np.cos(np.deg2rad(theta_hi)))
    x0, x1 = (min(mu_lo, mu_hi), max(mu_lo, mu_hi))
    sgn = theta_ray_sign(theta_mid)
    nx, ny = 64, 1200
    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(-1.0, 1.0, ny)
    X, Y = np.meshgrid(xs, ys)
    rho = np.sqrt(np.clip(1.0 - X * X, 1e-6, None))
    sc = np.clip(0.5 * (1.0 + sgn * Y / rho), 0.0, 1.0)
    img = np.zeros((ny, nx, 4))
    img[..., :3] = to_rgba(color)[:3]
    img[..., 3] = 0.28 * sc
    im = ax.imshow(img, extent=(x0, x1, -1.0, 1.0), origin="lower",
                   interpolation="bilinear", zorder=2, aspect="auto")
    im.set_clip_path(Circle((0.0, 0.0), 1.0, transform=ax.transData))

    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    rho_m = np.sqrt(max(0.0, 1.0 - mu_mid * mu_mid))
    phi0, phi = grad_phi(theta_mid, 721)
    ex_m = mu_mid + 0.045 * rho_m * np.cos(phi)
    ey_m = rho_m * np.sin(phi)
    runs = quantized_runs(segment_scale(grad_scale(phi, phi0)))
    ax.add_collection(
        LineCollection(
            _polylines_2d(ex_m, ey_m, runs),
            colors=rgba_array(color, [0.95 * r[2] for r in runs]),
            linewidths=2.8, zorder=4, capstyle="butt", joinstyle="round",
        )
    )


def draw_panelA_theta_rays(ax, geo):
    """
    Illustrate θ: from the observer (origin), one ray per loop out to the
    ring at that loop's center θ. θ is measured from +x (anti-CoE, right)
    up to the ring. Emphasize MAX θ=30° (thick semi-transparent red to the
    top of that foreshortened ring), matching Chris's markup.
    """
    items = list(zip(geo["theta_max_from_plus_x_deg"], COLORS_MAX)) + list(
        zip(geo["theta_min_from_plus_x_deg"], COLORS_MIN)
    )
    for mid, color in items:
        t = float(np.deg2rad(mid))
        # Endpoint = top/rim point of the ring in side view = (cos θ, sin θ)
        x1, y1 = float(np.cos(t)), float(np.sin(t))
        is_max30 = abs(mid - 29.5) < 0.5
        lw = 3.2 if is_max30 else 1.8
        alpha = 0.70 if is_max30 else 0.45
        ax.plot(
            [0.0, x1], [0.0, y1],
            color=color, lw=lw, alpha=alpha, zorder=3,
            solid_capstyle="round",
        )
        # Small tick at the ring hit
        ax.plot([x1], [y1], "o", color=color, ms=4.5 if is_max30 else 3.5,
                alpha=0.9, zorder=5, markeredgecolor="none")

    # Angle arc at θ=30° from +x axis up to the ray (makes the measure obvious)
    arc_r = 0.38
    phi = np.linspace(0.0, np.deg2rad(29.5), 60)
    ax.plot(
        arc_r * np.cos(phi), arc_r * np.sin(phi),
        color=C_MAX_30, lw=1.6, alpha=0.85, zorder=3,
    )
    ax.text(
        arc_r * np.cos(np.deg2rad(15.0)) + 0.06,
        arc_r * np.sin(np.deg2rad(15.0)) + 0.04,
        r"$\theta$",
        color=C_MAX_30, fontsize=10, fontweight="bold",
        ha="left", va="bottom", zorder=6,
    )


def panel_A_model_frame(ax, geo):
    """Side view: exactly four thick foreshortened loop ribbons about x."""
    ax.set_title(
        "A. Model frame (no Galactic yet)\n"
        r"axisymmetric about $x$ = CoE$\leftrightarrow$anti-CoE"
        "\n"
        r"two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$)",
        fontsize=10,
        pad=8,
    )
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color="#7f8c8d", lw=1.2, zorder=1)
    ax.plot(0, 0, "ko", ms=6, zorder=5)
    ax.text(0.08, -0.12, "observer", fontsize=8, color="k", ha="left")

    ax.annotate(
        "",
        xy=(1.15, 0),
        xytext=(-1.15, 0),
        arrowprops=dict(arrowstyle="<->", color=C_AXIS, lw=2.5),
        zorder=4,
    )
    ax.text(1.18, 0.08, "anti-CoE\n(+x)", color=C_ANTI, fontsize=8, ha="left", va="bottom")
    ax.text(-1.18, 0.08, "CoE\n(−x)", color=C_COE, fontsize=8, ha="right", va="bottom")

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        draw_panelA_ribbon(ax, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        draw_panelA_ribbon(ax, lo, hi, mid, color)

    # Illustrate θ from +x (anti-CoE) to each loop (Chris markup)
    draw_panelA_theta_rays(ax, geo)

    label_r = 1.30
    for mid, color in list(zip(geo["theta_max_from_plus_x_deg"], COLORS_MAX)) + list(
        zip(geo["theta_min_from_plus_x_deg"], COLORS_MIN)
    ):
        t = np.deg2rad(mid)
        ax.text(
            label_r * np.cos(t),
            label_r * np.sin(t),
            rf"$\theta\!=\!{mid:g}^\circ$",
            color=color, fontsize=7.5, ha="center", va="center", fontweight="bold",
        )

    ax.text(
        0.0,
        -1.38,
        (
            r"four transparent highlight bands about $x$  ·  each $18^\circ$ end-to-end ($\pm9^\circ$ = 5% of $360^\circ$)"
            + "\n"
            + r"(after-2$\ell$ $A(\ell)$: two MAX + two MIN; one mid-stroke each)"
            + "\n"
            + r"$y,z$ free (roll) — model does not pick celestial equator"
        ),
        fontsize=7.5,
        ha="center",
        va="top",
        style="italic",
        color="#34495e",
    )
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")


def _ring_lb(D, mu, n_pts=721):
    pts = small_circle_about_axis(D, mu, n=n_pts)
    rl = np.array([n_to_lb(p)[0] for p in pts])
    rb = np.array([n_to_lb(p)[1] for p in pts])
    return rl, rb


def _ring_lonlat_rad(D, mu, n_pts=361):
    """Small-circle → Mollweide display (lon_disp, lat) in radians."""
    pts = small_circle_about_axis(D, float(mu), n=n_pts)
    ls = np.array([n_to_lb(p)[0] for p in pts])
    bs = np.array([n_to_lb(p)[1] for p in pts])
    lon = np.deg2rad(lon_disp_from_ell(ls))
    lat = np.deg2rad(bs)
    return lon, lat


def fill_mollweide_ribbon(ax, D, theta_lo, theta_hi, theta_mid, color, n_pts=721):
    """
    Mollweide highlight band (±9° edges, base alpha 0.26) + mid-stroke
    (lw 2.6, alpha 0.95), both with the (1+cos psi)/2 opacity gradient.
    Elements crossing the antimeridian seam are dropped.
    """
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))

    def lonlat(mu, n):
        pts, sc = small_circle_grad(D, mu, theta_mid, n=n)
        ls = np.array([n_to_lb(p)[0] for p in pts])
        bs = np.array([n_to_lb(p)[1] for p in pts])
        return np.deg2rad(lon_disp_from_ell(ls)), np.deg2rad(bs), sc

    lon_a, lat_a, sc = lonlat(mu_a, n_pts)
    lon_b, lat_b, _ = lonlat(mu_b, n_pts)
    gap = 0.9 * np.pi
    quad = np.stack([lon_a[:-1], lon_a[1:], lon_b[1:], lon_b[:-1]])
    valid = (quad.max(axis=0) - quad.min(axis=0)) <= gap
    runs = quantized_runs(segment_scale(sc), valid)
    if runs:
        ax.add_collection(
            PolyCollection(
                _strip_polys_2d(lon_a, lat_a, lon_b, lat_b, runs),
                facecolors=rgba_array(color, [0.26 * r[2] for r in runs]),
                edgecolors="none", zorder=2, linewidths=0,
            )
        )
    lon_m, lat_m, sc_m = lonlat(mu_mid, max(n_pts, 721))
    valid = np.abs(np.diff(lon_m)) <= gap
    runs = quantized_runs(segment_scale(sc_m), valid)
    if runs:
        ax.add_collection(
            LineCollection(
                _polylines_2d(lon_m, lat_m, runs),
                colors=rgba_array(color, [0.95 * r[2] for r in runs]),
                linewidths=2.6, zorder=4, capstyle="butt", joinstyle="round",
            )
        )


def panel_B_mollweide(ax, geo):
    ax.set_title(
        r"B. Migkas $A(\ell)$ loops on sky (anti-CoE at $(0^\circ,+15^\circ)$)"
        "\n"
        r"(Galactic Mollweide; CoE axis tilted $+15^\circ$; $\ell$ increases LEFT; "
        r"two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$))",
        fontsize=10,
        pad=8,
    )
    for b in [-60, -30, 0, 30, 60]:
        ell = np.linspace(0, 360, 361)
        plot_lb_curve_mollweide(ax, ell, np.full_like(ell, b), color="#bdc3c7", lw=0.5, zorder=0)
    for ell0 in range(0, 360, 30):
        b = np.linspace(-90, 90, 181)
        plot_lb_curve_mollweide(ax, np.full_like(b, ell0), b, color="#bdc3c7", lw=0.5, zorder=0)

    D = geo["D"]
    ls, bs = great_circle_through_poles(D, n=721)
    plot_lb_curve_mollweide(ax, ls, bs, color=C_AXIS, lw=2.0, zorder=2)

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        fill_mollweide_ribbon(ax, D, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        fill_mollweide_ribbon(ax, D, lo, hi, mid, color)

    ell_eq, b_eq = celestial_equator_lb(721)
    plot_lb_curve_mollweide(ax, ell_eq, b_eq, color=C_EQ, lw=1.8, ls="--", zorder=2)

    ax.plot(
        np.deg2rad(lon_disp_from_ell(geo["anti_lb"][0])),
        np.deg2rad(geo["anti_lb"][1]),
        marker="o", ms=8, color=C_ANTI,
        markerfacecolor=C_ANTI, markeredgecolor="k", markeredgewidth=0.6,
        zorder=6, linestyle="None",
    )
    ax.plot(
        np.deg2rad(lon_disp_from_ell(geo["coe_lb"][0])),
        np.deg2rad(geo["coe_lb"][1]),
        marker="*", ms=14, color=C_COE,
        markerfacecolor=C_COE, markeredgecolor="k", markeredgewidth=0.6,
        zorder=6, linestyle="None",
    )
    for (l, b), color in zip(MIGKAS_MAXES, COLORS_MAX):
        ax.plot(
            np.deg2rad(lon_disp_from_ell(l)),
            np.deg2rad(b),
            marker="^",
            ms=10,
            color=color,
            markeredgecolor="k",
            markeredgewidth=0.7,
            zorder=7,
            linestyle="None",
        )
    for (l, b), color in zip(MIGKAS_MINS, COLORS_MIN):
        ax.plot(
            np.deg2rad(lon_disp_from_ell(l)),
            np.deg2rad(b),
            marker="s",
            ms=10,
            color=color,
            markeredgecolor="k",
            markeredgewidth=0.7,
            zorder=7,
            linestyle="None",
        )

    xticks = np.deg2rad(lon_disp_from_ell(np.array([180, 90, 0, 270])))
    ax.set_xticks(xticks)
    ax.set_xticklabels([r"$180^\circ$", r"$90^\circ$", r"$0^\circ$", r"$270^\circ$"], fontsize=8)
    ax.set_yticks(np.deg2rad([-60, -30, 0, 30, 60]))
    ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$+30^\circ$", r"$+60^\circ$"], fontsize=8)
    ax.grid(False)
    ax.set_xlabel(r"Galactic $\ell$ (increases LEFT)", fontsize=8)
    ax.set_ylabel(r"Galactic $b$", fontsize=8)


def fill_sphere_ribbon(ax, D, theta_lo, theta_hi, theta_mid, color, n_pts=160):
    """
    3D highlight belt (±9° edges, base alpha 0.25) + mid-stroke (lw 2.6,
    alpha 0.95), both with the (1+cos psi)/2 opacity gradient.
    """
    mu_a, mu_b = mu_pair_for_band(theta_lo, theta_hi)
    mu_mid = float(np.cos(np.deg2rad(theta_mid)))
    n_g = 361
    p0, sc = small_circle_grad(D, float(mu_a), theta_mid, n=n_g)
    p1, _ = small_circle_grad(D, float(mu_b), theta_mid, n=n_g)
    seg = segment_scale(sc)
    idx = np.where(seg > 1e-4)[0]
    verts = [[p0[j], p0[j + 1], p1[j + 1], p1[j]] for j in idx]
    coll = Poly3DCollection(
        verts, facecolors=rgba_array(color, 0.25 * seg[idx]),
        edgecolors="none", linewidths=0.0, zorder=2,
    )
    ax.add_collection3d(coll)
    pts, scm = small_circle_grad(D, mu_mid, theta_mid, n=n_g)
    runs = quantized_runs(segment_scale(scm))
    ax.add_collection3d(
        Line3DCollection(
            [pts[i0:i1 + 1] for i0, i1, _ in runs],
            colors=rgba_array(color, [0.95 * r[2] for r in runs]),
            linewidths=2.6, zorder=5, capstyle="butt", joinstyle="round",
        )
    )


def panel_C_sphere(ax, geo):
    th_max = geo["theta_max_from_plus_x_deg"]
    th_min = geo["theta_min_from_plus_x_deg"]
    ax.set_title(
        r"C. Sphere 3D — two MAX + two MIN, each $18^\circ$ end-to-end ($\pm9^\circ$)"
        "\n"
        r"(axis line emphasized; celestial equator dashed)",
        fontsize=10,
        pad=6,
    )
    draw_wire_sphere(ax)
    D = geo["D"]
    ell_eq, b_eq = celestial_equator_lb(400)
    eq_xyz = np.array([lb_to_n(l, b) for l, b in zip(ell_eq, b_eq)])
    ax.plot(eq_xyz[:, 0], eq_xyz[:, 1], eq_xyz[:, 2], color=C_EQ, lw=1.8, ls="--", zorder=3)
    ax.plot(
        [-D[0], D[0]],
        [-D[1], D[1]],
        [-D[2], D[2]],
        color=C_AXIS,
        lw=3.0,
        zorder=5,
    )
    ax.scatter(*(-D), color=C_COE, s=120, marker="*", edgecolors="k", linewidths=0.6,
               depthshade=False, zorder=8)
    ax.scatter(*(D), color=C_ANTI, s=80, edgecolors="k", linewidths=0.6, depthshade=False, zorder=8)
    ax.text(*(-D * 1.22), "CoE", color=C_COE, fontsize=8, fontweight="bold")
    ax.text(*(D * 1.18), "anti-CoE", color=C_ANTI, fontsize=8, fontweight="bold")
    ax.scatter([0], [0], [0], color="k", s=25, depthshade=False, zorder=9)

    for (lo, hi), mid, color in zip(
        geo["max_bands_deg"], geo["theta_max_from_plus_x_deg"], COLORS_MAX
    ):
        fill_sphere_ribbon(ax, D, lo, hi, mid, color)
    for (lo, hi), mid, color in zip(
        geo["min_bands_deg"], geo["theta_min_from_plus_x_deg"], COLORS_MIN
    ):
        fill_sphere_ribbon(ax, D, lo, hi, mid, color)

    for n, color in zip(geo["n_maxes"], COLORS_MAX):
        ax.scatter(
            [n[0]], [n[1]], [n[2]],
            c=[color], marker="^", s=70,
            edgecolors="k", linewidths=0.7, depthshade=False, zorder=10,
        )
        ax.plot([0, n[0]], [0, n[1]], [0, n[2]], color=color, lw=0.8, ls=":", alpha=0.7)
    for n, color in zip(geo["n_mins"], COLORS_MIN):
        ax.scatter(
            [n[0]], [n[1]], [n[2]],
            c=[color], marker="s", s=60,
            edgecolors="k", linewidths=0.7, depthshade=False, zorder=10,
        )
        ax.plot([0, n[0]], [0, n[1]], [0, n[2]], color=color, lw=0.8, ls=":", alpha=0.7)

    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.set_zlim(-1.25, 1.25)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=18, azim=-55)
    ax.set_axis_off()
    ax.text2D(
        0.5,
        -0.02,
        rf"MAX $\theta\!=\!{th_max[0]:g}^\circ,{th_max[1]:g}^\circ$; "
        rf"MIN $\theta\!=\!{th_min[0]:g}^\circ,{th_min[1]:g}^\circ$  "
        rf"($18^\circ$ end-to-end = $\pm9^\circ$ = 5% of $360^\circ$)",
        transform=ax.transAxes,
        ha="center",
        fontsize=7.5,
        color="#34495e",
    )


def main():
    st = migkas_scan_stops()
    geo = compute_geometry(st)

    th_max = geo["theta_max_from_plus_x_deg"]
    th_min = geo["theta_min_from_plus_x_deg"]
    max_bands = geo["max_bands_deg"]
    min_bands = geo["min_bands_deg"]

    fig = plt.figure(figsize=(14.5, 5.6), dpi=160)
    fig.suptitle(
        r"Pure CoE axisymmetry: two MAX + two MIN highlight bands  (gradient: darkest at $\theta$ ray, transparent opposite; +80.8$^\circ$ CCW about CoE axis (B/C))"
        "\n"
        rf"(each $18^\circ$ end-to-end = $\pm9^\circ$ = 5% of $360^\circ$; "
        rf"MAX $\theta\!=\!{th_max[0]:g}^\circ/{th_max[1]:g}^\circ$; "
        rf"MIN $\theta\!=\!{th_min[0]:g}^\circ/{th_min[1]:g}^\circ$; "
        r"anti-CoE $(0^\circ,+15^\circ)$, CoE axis tilted $+15^\circ$)",
        fontsize=11,
        fontweight="bold",
        y=1.02,
    )

    axA = fig.add_subplot(1, 3, 1)
    panel_A_model_frame(axA, geo)

    axB = fig.add_subplot(1, 3, 2, projection="mollweide")
    panel_B_mollweide(axB, geo)

    axC = fig.add_subplot(1, 3, 3, projection="3d")
    panel_C_sphere(axC, geo)

    handles = [
        Line2D([0], [0], color=C_AXIS, lw=2.5, label=r"CoE$\leftrightarrow$anti-CoE axis"),
        Patch(
            facecolor=C_MAX_30, alpha=0.45, edgecolor=C_MAX_30,
            label=(
                rf"MAX $\theta\!=\!{th_max[0]:g}^\circ$ band "
                rf"[{max_bands[0][0]:g},{max_bands[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$=5% of $360^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MAX_185, alpha=0.45, edgecolor=C_MAX_185,
            label=(
                rf"MAX $\theta\!=\!{th_max[1]:g}^\circ$ band "
                rf"[{max_bands[1][0]:g},{max_bands[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_115, alpha=0.45, edgecolor=C_MIN_115,
            label=(
                rf"MIN $\theta\!=\!{th_min[0]:g}^\circ$ band "
                rf"[{min_bands[0][0]:g},{min_bands[0][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Patch(
            facecolor=C_MIN_265, alpha=0.45, edgecolor=C_MIN_265,
            label=(
                rf"MIN $\theta\!=\!{th_min[1]:g}^\circ$ band "
                rf"[{min_bands[1][0]:g},{min_bands[1][1]:g}]$^\circ$ "
                rf"($18^\circ$=±$9^\circ$)"
            ),
        ),
        Line2D([0], [0], color=C_EQ, lw=1.8, ls="--", label="celestial equator (Dec$=0$)"),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_30, markeredgecolor="k",
            ms=9,
            label=rf"MAX $({MIGKAS_MAXES[0][0]:.0f}^\circ,{MIGKAS_MAXES[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="^", color="w", markerfacecolor=C_MAX_185, markeredgecolor="k",
            ms=9,
            label=rf"MAX $({MIGKAS_MAXES[1][0]:.0f}^\circ,{MIGKAS_MAXES[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_115, markeredgecolor="k",
            ms=8,
            label=rf"MIN $({MIGKAS_MINS[0][0]:.0f}^\circ,{MIGKAS_MINS[0][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="s", color="w", markerfacecolor=C_MIN_265, markeredgecolor="k",
            ms=8,
            label=rf"MIN $({MIGKAS_MINS[1][0]:.0f}^\circ,{MIGKAS_MINS[1][1]:+.0f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="o", color="w", markerfacecolor=C_ANTI, markeredgecolor="k",
            ms=8, label=rf"anti-CoE $({geo['anti_lb'][0]:.1f}^\circ,{geo['anti_lb'][1]:+.1f}^\circ)$",
        ),
        Line2D(
            [0], [0], marker="*", color="w", markerfacecolor=C_COE, markeredgecolor="k",
            ms=11, label=rf"CoE $({geo['coe_lb'][0]:.1f}^\circ,{geo['coe_lb'][1]:+.1f}^\circ)$",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=6,
        fontsize=7.0,
        frameon=True,
        fancybox=True,
        bbox_to_anchor=(0.5, -0.10),
    )

    fig.tight_layout()
    # --- Chris 2026-09-26: panel C sphere drawn at the same on-page diameter as
    # panel A's model-frame circle; panel C axes matched to A's vertical
    # position/height. Nothing else changes.
    from mpl_toolkits.mplot3d import proj3d
    pa, pc = axA.get_position(), axC.get_position()
    axC.set_position([pc.x0, pa.y0, pc.width, pa.height])

    def _diam_A():
        fig.canvas.draw()
        p0 = axA.transData.transform((-1.0, 0.0)); p1 = axA.transData.transform((1.0, 0.0))
        return float(p1[0] - p0[0])

    def _diam_C():
        fig.canvas.draw()
        u = np.random.default_rng(0).normal(size=(20000, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]
        xs, ys, _ = proj3d.proj_transform(u[:, 0], u[:, 1], u[:, 2], axC.get_proj())
        d = axC.transData.transform(np.column_stack([xs, ys]))
        return float(d[:, 0].max() - d[:, 0].min()), float(d[:, 1].max() - d[:, 1].min())

    zoom = 1.0
    for _ in range(8):
        dA = _diam_A(); dCx, dCy = _diam_C()
        zoom *= dA / dCx
        axC.set_box_aspect((1, 1, 1), zoom=zoom)
    dA = _diam_A(); dCx, dCy = _diam_C()
    print(f"PANEL_DIAM A={dA:.1f}px C={dCx:.1f}x{dCy:.1f}px zoom={zoom:.3f}")
    fig.savefig(PNG_PATH, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    payload = {
        "png": str(PNG_PATH),
        "script": str(OUT_DIR / "make_coe_rings_migkas_pin_axis_matched_grad_rot_fit.py"),
        "migkas_max_lb_deg": [list(lb) for lb in MIGKAS_MAXES],
        "migkas_min_lb_deg": [list(lb) for lb in MIGKAS_MINS],
        "strip_theta_max_deg": list(STRIP_THETA_MAX),
        "strip_theta_min_deg": list(STRIP_THETA_MIN),
        "band_half_width_deg": BAND_HALF,
        "band_full_width_deg": BAND_FULL,
        "band_total_width_deg": BAND_TOTAL,
        "band_width_note": "5% of 360° = 18° end-to-end = ±9° about each MAX/MIN center",
        "max_bands_deg": [list(b) for b in max_bands],
        "min_bands_deg": [list(b) for b in min_bands],
        "band_intervals": {
            "MAX_30": [float(max_bands[0][0]), float(max_bands[0][1])],
            "MAX_185": [float(max_bands[1][0]), float(max_bands[1][1])],
            "MIN_115": [float(min_bands[0][0]), float(min_bands[0][1])],
            "MIN_265": [float(min_bands[1][0]), float(min_bands[1][1])],
        },
        "scan_stop_source": st["source"],
        "theta_max_from_plus_x_deg": [round(t, 4) for t in th_max],
        "theta_min_from_plus_x_deg": [round(t, 4) for t in th_min],
        "mu_max_bands": [[round(a, 8), round(b, 8)] for a, b in geo["mu_max_bands"]],
        "mu_min_bands": [[round(a, 8), round(b, 8)] for a, b in geo["mu_min_bands"]],
        "CoE_lb_deg": [round(geo["coe_lb"][0], 4), round(geo["coe_lb"][1], 4)],
        "anti_CoE_lb_deg": [round(geo["anti_lb"][0], 4), round(geo["anti_lb"][1], 4)],
        "D_hat_galactic_xyz": geo["D"].tolist(),
        "angular_sep_migkas_max0_min1_deg": round(geo["sep_deg"], 6),
        "axis_tilt_b_deg": 15,
        "axis_tilt_reason": (
            "Chris 2026-09-22 — markers fixed, rings rotate up 30°"
        ),
        "orientation": (
            "Galactic orientation: anti-CoE at (ell=0°, b=+15°), CoE at "
            "(ell=180°, b=−15°); +D = anti-CoE = +x; axis tilted +15° "
            "(Chris 2026-09-22 — markers fixed, rings rotate up 30°)"
        ),
        "drawing": (
            "Exactly 4 features. Each = ONE smooth transparent highlight band "
            "between θ−9° and θ+9° (18° end-to-end = 5% of 360°) + ONE solid "
            "mid-stroke at the exact center θ. No multi-stroke wash / hairball. "
            "Four distinct colors (two reds + two blues); no orange; no 283° "
            "companion. Panel A: clean vertical strip clipped to unit circle; theta rays from observer to each loop (emphasized at MAX 30°). Mollweide "
            "seam-split honestly. Sphere: single Poly3D strip (no n_wash stack)."
        ),
        "note": (
            "Stops are after-2ℓ A(ℓ) extrema from best_fit_2ell_multiplier.npz "
            f"A_after: two MAX at θ={th_max[0]:g}°/{th_max[1]:g}° and "
            f"two MIN at θ={th_min[0]:g}°/{th_min[1]:g}°, each drawn as "
            f"one {BAND_FULL:.0f}° end-to-end highlight (±{BAND_HALF:.0f}° = 5% of 360°). "
            "Rings use μ=cos(θ); sky markers are map extrema (17°,−9°) and "
            "(274°,−22°) plus strip stand-ins (185°,0°) and (115°,0°)."
        ),
        "plot_convention": "Migkas Galactic Mollweide: lon_disp=-ell_wrapped (ell increases LEFT)",
        "legend_summary": (
            "two MAX + two MIN highlight bands, each 18° end-to-end "
            "(±9° = 5% of 360°); four distinct colors"
        ),
        "colors": {
            "MAX_30": C_MAX_30,
            "MAX_185": C_MAX_185,
            "MIN_115": C_MIN_115,
            "MIN_265": C_MIN_265,
            "no_orange": True,
        },
    }
    JSON_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    print(f"\nWrote {PNG_PATH}")
    print(f"Wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
```

#### `make_stack_clean.py` (related poster stack — not letter four_row)

Crops `loop_threepanel_*` rows and overlays; writes `stack_3col_clean.jpg` / `stack_4col_clean.jpg` under spreadsheet_model. Attachment paths inside are agent-specific and may be stale.

```python
# source: /workspace/cosmology/figures/make_stack_clean.py
from PIL import Image, ImageDraw, ImageFont
import numpy as np
A = "/home/box/agent-data/agents/196ab413-0708-4f81-99fc-58a32c5366c5/attachments/"
P = "/workspace/cosmology/papers/spreadsheet_model/"
graph = Image.open(P + "bubble/newt_omb_rel_Agraph.png").convert("RGB")   # was attachment f480c67c... (newt_omb_Agraph); + relativistic curve
rows_src = [Image.open(f).convert("RGB") for f in (P + "loop_threepanel_notilt_web.jpg",
            A + "e2fab4ef4115016f7bc638a6ac1075749b6acbab6e4c2ca6b2f20e29f8735139.jpg",
            A + "fe905a57d5aaf4dede004a4cbdc50883bbf92c984372f9f260251bbfc057b757.jpg",
            A + "d6982109d332a918353df473def180356c133fb81196c90b6e0dea313238ef05.jpg")]
ovs = [Image.open(f).convert("RGB") for f in (A + "3b64499f8c49ab6537b8a73b1fd80965e8b5f8ae82d93d6c99b61dc4f5b8344d.jpg",
       A + "83d1ae71406713257f0623d2b732446ef1ee220de6167592d4a45c5cac55637c.jpg", P + "fig8_overlay_matched_grad_rot81_web.jpg")]
LAB = ["Untilted, no roll", "+15° tilt", "+15° tilt, gradient", "+15° tilt, gradient, +80.8° roll"]
OLAB = ["+15° tilt", "gradient", "+80.8° roll"]
FB = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 22)
FH = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 18)
FT = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 17)
W = 1700; LH = 40
def row_content(i):
    im = rows_src[i]
    if i == 0:
        return im.crop((0, 110, W, 628))            # keep A./B./C. headings + footnotes (first row only)
    c = im.crop((0, 222, W, 572)).copy()          # drops suptitle, A/C headings, footnotes, legend
    ImageDraw.Draw(c).rectangle((385, 234 - 222, 1212, 279 - 222), fill="white")   # B heading
    return c
legend = rows_src[3].crop((0, 712, W, 768))        # one legend strip (last row's)
def ov_crop(k, cbar):
    return ovs[k].crop((140, 355, 1445 if cbar else 1276, 1060))
def label(dr, x, y, t, f=FB, anchor="mm"):
    dr.text((x, y), t, font=f, fill="black", anchor=anchor)
def save(im, name):
    im.save(P + name + ".jpg", quality=92)
    w = 1400; im.resize((w, round(im.height * w / im.width)), Image.LANCZOS).save(P + name + "_web.jpg", quality=88)
    print(name, im.size)
# ---------------- 3-column version ----------------
parts = [("g", graph)]
blocks = []
H = graph.height + 20
for i in range(4):
    rc = row_content(i); blocks.append((i, rc)); H += LH + rc.height + 10
s3 = (W - 20 - 2 * 20) / (ov_crop(0, 0).width * 2 + ov_crop(2, 1).width)
oc = [ov_crop(0, 0), ov_crop(1, 0), ov_crop(2, 1)]
oc = [o.resize((round(o.width * s3), round(o.height * s3)), Image.LANCZOS) for o in oc]
OT = 26
H += LH + OT + max(o.height for o in oc) + 10 + legend.height + 10
im = Image.new("RGB", (W, H), "white"); dr = ImageDraw.Draw(im)
y = 0; im.paste(graph, ((W - graph.width) // 2, y)); y += graph.height + 20
for i, rc in blocks:
    label(dr, W // 2, y + LH // 2, LAB[i]); y += LH; im.paste(rc, (0, y)); y += rc.height + 10
label(dr, W // 2, y + LH // 2, "Migkas Fig. 8, θ = 75°, A/Aₐₗₗ"); y += LH
x = 10
for k, o in enumerate(oc):
    label(dr, x + round((1276 - 140) * s3) // 2, y + OT // 2, OLAB[k], FT)
    im.paste(o, (x, y + OT)); x += o.width + 20
y += OT + max(o.height for o in oc) + 10
im.paste(legend, (0, y)); y += legend.height + 10
save(im.crop((0, 0, W, y)), "stack_3col_clean")
# ---------------- 4-column version ----------------
s4 = 772 / 1304   # same overlay scale as stack_4col.jpg (no shrinking)
o4 = [ov_crop(0, 1), ov_crop(1, 0), ov_crop(2, 0)]
o4 = [o.resize((round(o.width * s4), round(o.height * s4)), Image.LANCZOS) for o in o4]
CW = max(o.width for o in o4) + 20; W4 = W + CW; HT = 30
Hs = [LH + blocks[0][1].height + 10]
for j in range(1, 4):
    rc = blocks[j][1]; Hs.append(max(LH + rc.height, (HT if j == 1 else 0) + o4[j - 1].height) + 12)
H4 = graph.height + 20 + sum(Hs) + legend.height + 10
im = Image.new("RGB", (W4, H4), "white"); dr = ImageDraw.Draw(im)
y = 0; im.paste(graph, ((W4 - graph.width) // 2, y)); y += graph.height + 20
for j in range(4):
    rc = blocks[j][1]
    label(dr, W // 2, y + LH // 2, LAB[j]); im.paste(rc, (0, y + LH))
    if j >= 1:
        oy = y
        if j == 1:
            label(dr, W + 10 + round((1276 - 140) * s4) // 2, y + HT // 2, "Migkas Fig. 8, θ = 75°, A/Aₐₗₗ", FH); oy += HT
        im.paste(o4[j - 1], (W + 10, oy))
    y += Hs[j]
im.paste(legend, (0, y)); y += legend.height + 10
save(im.crop((0, 0, W4, y)), "stack_4col_clean")
```

### 4.3 Fig07 current donut script (output ≠ gold)

```python
# source: /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle7_recombination/make_donut_universe_storyboard.py
#!/usr/bin/env python3
"""Donut-universe stacked storyboard (A→E) — geometry only.

CoE fixed at origin; US mid-shell D=(Rin+Rout)/2.
Fixed CoE-centred recombination front Rrec.
Hubble bubble (visible sphere around US): Rvis=0 and V_HB=0 until D>D_recomb(=Rrec);
thereafter Rvis grows at c while shells expand faster than c.

Stage prose / volumes live in letter.tex (right column), not in this PNG.
All panels share one axis range so the R_rec circle is the same on-page size;
shells grow within that fixed frame (no per-stage zoom).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

OUT = Path(__file__).resolve().parent / "donut_universe_storyboard.png"

# Drawing units: c = 1. Shell speeds > c; D = v_D * t with v_D = (v_in+v_out)/2.
# R_rec fixed; vision turns on at D = R_rec ⇒ t_rec = R_rec / v_D.
C = 1.0
# Half-thickness growth (V_OUT-V_IN)/2 >= C so Rvis never appears to outrun the shells.
V_IN = 2.0
V_OUT = 4.0
V_D = 0.5 * (V_IN + V_OUT)  # 3.0
R_REC = 1.20
T_REC = R_REC / V_D  # 0.40
D_RECOMB = R_REC


def stage_at(t: float, title: str, caption: str, shade_annulus: bool = True) -> dict:
    rin = V_IN * t
    rout = V_OUT * t
    d = V_D * t
    assert abs(d - 0.5 * (rin + rout)) < 1e-9
    if d > D_RECOMB + 1e-12:
        rvis = C * (t - T_REC)
    else:
        rvis = 0.0
    return dict(
        title=title,
        caption=caption,
        t=t,
        D=d,
        Rin=rin,
        Rout=rout,
        Rvis=rvis,
        Rrec=R_REC,
        shade_annulus=shade_annulus,
    )


STAGES = [
    stage_at(
        0.25,
        r"Time A — before $D_{\mathrm{recomb}}$",
        r"US still inside the recombination front; vision off. "
        r"Hubble-bubble volume is zero ($R_{\mathrm{vis}}=0$).",
    ),
    stage_at(
        T_REC,
        r"Time B — US at $D_{\mathrm{recomb}}$",
        r"Observer reaches the fixed recombination radius; vision turns on. "
        r"Bubble volume still zero at the turn-on instant.",
    ),
    stage_at(
        0.90,
        r"Time C — early post-recombination",
        r"Shells have crossed $R_{\mathrm{rec}}$. Visible sphere around US grows at $c$; "
        r"shells expand faster than light.",
    ),
    stage_at(
        1.60,
        r"Time D — later expansion",
        r"Mid-shell observer; Hubble bubble is a growing fraction of the shell annulus "
        r"while remaining slower than shell growth.",
    ),
    stage_at(
        2.40,
        r"Time E — today-like donut universe",
        r"Thick shell about the CoE; US mid-shell; visible sphere at $c$ samples "
        r"the donut without outrunning the $>c$ shells.",
    ),
]


def vol(r: float) -> float:
    return (4.0 / 3.0) * np.pi * r**3


def volume_block(s: dict) -> dict:
    v_hb = vol(s["Rvis"])
    v_s1 = vol(s["Rout"])
    v_s2 = vol(s["Rin"])
    v_ann = v_s1 - v_s2
    return dict(
        V_HB=v_hb,
        V_S1=v_s1,
        V_S2=v_s2,
        V_ann=v_ann,
        f_HB_ann=(v_hb / v_ann) if v_ann > 0 else 0.0,
        f_HB_S1=(v_hb / v_s1) if v_s1 > 0 else 0.0,
        f_S2_S1=(v_s2 / v_s1) if v_s1 > 0 else 0.0,
    )


def shared_frame_limits(stages: list[dict], margin: float = 0.55) -> tuple[float, float, float]:
    """One axis range for all panels so R_rec keeps the same on-page size.

    Shells grow inside this fixed frame; R_rec never zooms with them.
    """
    xmax = 0.0
    xmin = 0.0
    ymax = 0.0
    for s in stages:
        D = s["D"]
        r_pad = max(s["Rvis"], 0.15)
        xmax = max(xmax, s["Rout"], D + r_pad)
        xmin = min(xmin, -s["Rout"], D - r_pad)
        ymax = max(ymax, s["Rout"], s["Rvis"], s["Rrec"], 0.55)
    # Symmetric x about CoE so the fixed R_rec circle is centred the same way every panel.
    half = max(abs(xmin), abs(xmax)) + margin
    return -half, half, ymax + 0.32


# Shared across A–E (computed once from STAGES).
XMIN, XMAX, YMAX = shared_frame_limits(STAGES)


def draw_stage(ax, s: dict, stage_label: str):
    """Draw geometry only; stage prose lives in letter.tex right column.

    Axis limits are the shared (XMIN, XMAX, YMAX) so R_rec is identical in
    drawing units and in on-page size across every panel.
    """
    D = s["D"]
    coe = (0.0, 0.0)
    us = (D, 0.0)

    # S1 disc (blue-ish) then S2 fill (red) so the annulus reads as shell material.
    if s["shade_annulus"]:
        ax.add_patch(
            Circle(coe, s["Rout"], facecolor="#d6eaf8", edgecolor="none", alpha=0.45, zorder=0)
        )
        ax.add_patch(
            Circle(coe, s["Rin"], facecolor="#f5c6c6", edgecolor="none", alpha=0.85, zorder=1)
        )
    else:
        ax.add_patch(
            Circle(coe, s["Rin"], facecolor="#f5c6c6", edgecolor="none", alpha=0.7, zorder=1)
        )

    ax.add_patch(Circle(coe, s["Rout"], fill=False, ec="#2980b9", lw=2.0, zorder=2))
    ax.add_patch(Circle(coe, s["Rin"], fill=False, ec="#c0392b", lw=1.8, zorder=2))

    # Fixed CoE-centred recombination front (same radius and same frame every stage).
    ax.add_patch(
        Circle(coe, s["Rrec"], fill=False, ec="#8e44ad", lw=1.6, ls=(0, (5, 3)), zorder=3)
    )

    # Hubble bubble / visible sphere (green-cyan); omit disc when Rvis=0.
    if s["Rvis"] > 1e-6:
        ax.add_patch(
            Circle(
                us,
                s["Rvis"],
                facecolor="#a8e6cf",
                edgecolor="#1abc9c",
                lw=1.8,
                ls="--",
                alpha=0.40,
                zorder=3,
            )
        )
        ax.add_patch(Circle(us, s["Rvis"], fill=False, ec="#16a085", lw=1.6, ls="--", zorder=4))

    # Markers: black = CoE, green = US. Fixed drawing-unit label offsets (shared frame).
    ax.plot(*coe, "ko", ms=4.5, zorder=6)
    ax.plot(*us, "o", color="#27ae60", ms=5.5, zorder=6)
    if s["Rvis"] <= 1e-6:
        ax.plot(*us, "o", color="#1abc9c", ms=2.5, zorder=7)

    label_dy = 0.45  # fixed in drawing units so labels do not scale with shell zoom
    ax.text(
        coe[0],
        coe[1] - label_dy,
        "CoE",
        fontsize=6.5,
        ha="center",
        va="top",
        color="#222222",
        zorder=8,
    )
    ax.text(
        us[0],
        us[1] + label_dy,
        "US",
        fontsize=6.5,
        ha="center",
        va="bottom",
        color="#1e8449",
        zorder=8,
    )

    # S1 / S2 tick labels only when the shell is large enough not to clutter.
    if s["Rout"] > 1.4:
        ax.text(
            0.0,
            s["Rout"] * 0.92,
            "S1",
            fontsize=6,
            ha="center",
            va="bottom",
            color="#1a5276",
            zorder=8,
        )
    if s["Rin"] > 0.9:
        ax.text(
            0.0,
            s["Rin"] * 0.55,
            "S2",
            fontsize=6,
            ha="center",
            va="center",
            color="#922b21",
            zorder=8,
        )

    # Expansion arrows (unlabeled); length in fixed drawing units.
    arrow_len = 0.55
    for r, yoff in [(s["Rout"] * 0.55, 0.35), (s["Rin"] * 0.55, -0.35)]:
        if r > 0.15:
            ax.annotate(
                "",
                xy=(-r - arrow_len, yoff),
                xytext=(-r + 0.05, yoff),
                arrowprops=dict(arrowstyle="->", color="#7f8c8d", lw=1.3),
            )

    if s["Rvis"] > 1e-6:
        ax.annotate(
            "",
            xy=(us[0] + s["Rvis"] * 0.85, us[1] + 0.15),
            xytext=(us[0] + 0.12, us[1] + 0.15),
            arrowprops=dict(arrowstyle="->", color="#16a085", lw=1.4),
        )

    # US outward motion stub (fixed drawing length).
    us_arrow = 0.55
    ax.annotate(
        "",
        xy=(us[0] + us_arrow, -0.35),
        xytext=(us[0] - 0.02, -0.35),
        arrowprops=dict(arrowstyle="->", color="#27ae60", lw=1.2),
    )

    # Stage letter outside the drawing, left margin.
    ax.text(
        0.02,
        0.95,
        stage_label,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="top",
        color="#1f2d3d",
        zorder=10,
    )

    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(-YMAX, YMAX)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    # Geometry-only vertical stack for letter left column.
    # Shared axis range on every panel: R_rec same size in frame; shells grow into it.
    n = len(STAGES)
    fig = plt.figure(figsize=(3.6, 10.2))
    gs = gridspec.GridSpec(
        n,
        1,
        figure=fig,
        height_ratios=[1] * n,
        hspace=0.08,
        top=0.99,
        bottom=0.01,
        left=0.02,
        right=0.98,
    )

    labels = "ABCDE"
    print(f"shared frame: xlim=({XMIN:.2f},{XMAX:.2f}) ylim=(-{YMAX:.2f},{YMAX:.2f}) R_rec={R_REC}")
    for i, (label, s) in enumerate(zip(labels, STAGES)):
        ax_draw = fig.add_subplot(gs[i, 0])
        draw_stage(ax_draw, s, label)

    fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
    print("wrote", OUT)

    print("\nSTAGE_TABLE")
    print(
        "stage\tt\tD\tRin\tRout\tRvis\tRrec\tV_HB\tV_S1\tV_S2\tV_ann\t"
        "V_HB/V_ann\tV_HB/V_S1\tV_S2/V_S1"
    )
    for lab, s in zip(labels, STAGES):
        vb = volume_block(s)
        print(
            f"{lab}\t{s['t']:.3f}\t{s['D']:.4f}\t{s['Rin']:.4f}\t{s['Rout']:.4f}\t"
            f"{s['Rvis']:.4f}\t{s['Rrec']:.4f}\t{vb['V_HB']:.6f}\t{vb['V_S1']:.4f}\t"
            f"{vb['V_S2']:.4f}\t{vb['V_ann']:.4f}\t{vb['f_HB_ann']:.6f}\t"
            f"{vb['f_HB_S1']:.6f}\t{vb['f_S2_S1']:.6f}"
        )


if __name__ == "__main__":
    main()
```

### 4.4 Fig08 current two-emitters script (output ≠ gold)

```python
# source: /workspace/cosmology/papers/A_equals_fourier/steps_v3/cmb_reconcile/cycle6_two_emitters/make_two_emitters_comoving.py
#!/usr/bin/env python3
"""Clean two-emitter comoving figure for the letter brief.

The figure is deliberately limited to the scientific geometry:
parallel plates above, concentric shells below, with the observer (US)
between the two emitters in each construction.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np

HERE = Path(__file__).resolve().parent
LETTER_FIGS = HERE.parent / "letter_brief" / "figs"
OUT_CYCLE = HERE / "two_emitters_comoving.png"
OUT_LETTER = LETTER_FIGS / "two_emitters_comoving.png"


def bulk_arrows(ax, y: float, x_positions: tuple[float, ...], length: float = 0.72) -> None:
    """Show the common bulk motion without competing with the geometry."""
    for x0 in x_positions:
        ax.annotate(
            "",
            xy=(x0 + length, y),
            xytext=(x0, y),
            arrowprops=dict(arrowstyle="-|>", color="0.40", lw=1.8, mutation_scale=11),
        )
    ax.text(
        max(x_positions) + length + 0.18,
        y,
        r"shared bulk $V$",
        fontsize=10.5,
        color="0.30",
        ha="left",
        va="center",
    )


def light_arrow(ax, x0: float, y: float, direction: int) -> None:
    """Draw a short light-ray arrow pointing away from or toward an emitter."""
    span = 0.34
    if direction < 0:
        start, end = x0 - 0.06, x0 - span
    else:
        start, end = x0 + 0.06, x0 + span
    ax.annotate(
        "",
        xy=(end, y),
        xytext=(start, y),
        arrowprops=dict(arrowstyle="->", color="0.10", lw=1.05, mutation_scale=8),
    )


def draw_plates(ax) -> None:
    """Parallel plates with the observer between them."""
    ax.set_xlim(-3.6, 3.6)
    ax.set_ylim(-2.35, 2.65)
    ax.set_aspect("equal")
    ax.axis("off")

    bulk_arrows(ax, 2.16, (-2.45, -1.18, 0.09), length=0.72)

    y0, y1 = -1.55, 1.48
    x_left, x_right = -1.62, 1.62
    ax.plot([x_left, x_left], [y0, y1], color="0.05", lw=3.0, solid_capstyle="round")
    ax.plot([x_right, x_right], [y0, y1], color="0.05", lw=3.0, solid_capstyle="round")
    ax.text(x_left - 0.20, 0.0, "S2", fontsize=12.5, ha="right", va="center", fontweight="bold")
    ax.text(x_right + 0.20, 0.0, "S1", fontsize=12.5, ha="left", va="center", fontweight="bold")

    ax.plot(0.0, 0.0, "o", color="0.05", ms=9, zorder=5)
    ax.text(0.0, 0.24, "US", fontsize=12, ha="center", va="bottom", fontweight="bold")

    for x in (x_left, x_right):
        for y in np.linspace(-1.12, 1.05, 4):
            light_arrow(ax, x, y, -1)
            light_arrow(ax, x, y, +1)

    ax.text(
        0.5,
        -0.015,
        r"plates and US share the same bulk $V$; relative Doppler $\approx 0$",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10,
        style="italic",
        color="0.25",
    )


def radial_light_arrow(ax, radius: float, angle_deg: float, direction: int) -> None:
    """Short radial light arrow at a shell."""
    a = np.deg2rad(angle_deg)
    ux, uy = np.cos(a), np.sin(a)
    r0, r1 = radius + direction * 0.05, radius + direction * 0.34
    ax.annotate(
        "",
        xy=(r1 * ux, r1 * uy),
        xytext=(r0 * ux, r0 * uy),
        arrowprops=dict(arrowstyle="->", color="0.10", lw=1.0, mutation_scale=8),
    )


def draw_shells(ax) -> None:
    """Concentric shell section with the observer in the mid-shell gap."""
    ax.set_xlim(-3.8, 4.2)
    ax.set_ylim(-3.15, 3.55)
    ax.set_aspect("equal")
    ax.axis("off")

    R1, R2 = 1.50, 2.65
    theta1, theta2 = -110, 110  # section view; open side keeps the observer visible
    for radius in (R1, R2):
        ax.add_patch(
            Arc(
                (0.0, 0.0),
                2 * radius,
                2 * radius,
                theta1=theta1,
                theta2=theta2,
                color="0.05",
                lw=2.8,
            )
        )

    ax.plot(0.0, 0.0, "+", color="0.10", ms=12, mew=1.8, zorder=5)
    ax.text(0.0, -0.34, "CoE", fontsize=12, ha="center", va="top")

    mid = 0.5 * (R1 + R2)
    ax.plot(mid, 0.0, "o", color="0.05", ms=9, zorder=6)
    ax.text(mid, 0.27, "US", fontsize=12, ha="center", va="bottom", fontweight="bold")
    ax.text(mid, -0.28, "mid-shell", fontsize=8.5, ha="center", va="top", color="0.30")

    for radius in (R1, R2):
        for angle in (-38, -17, 17, 38):
            radial_light_arrow(ax, radius, angle, -1)
            radial_light_arrow(ax, radius, angle, +1)

    a_lab = np.deg2rad(54)
    for radius, label in ((R1, "S1"), (R2, "S2")):
        ax.text(
            (radius + 0.12) * np.cos(a_lab),
            (radius + 0.12) * np.sin(a_lab),
            label,
            fontsize=12,
            ha="left",
            va="bottom",
            fontweight="bold",
        )

    # Same outward comoving sense for the shells and the observer.
    for x0, y0, length in ((0.22, 0.58, 0.38), (0.92, 0.58, 0.38), (mid + 0.12, 0.0, 0.42)):
        ax.annotate(
            "",
            xy=(x0 + length, y0),
            xytext=(x0, y0),
            arrowprops=dict(arrowstyle="-|>", color="0.40", lw=1.7, mutation_scale=10),
        )
    ax.text(mid + 0.72, 0.20, r"shared bulk $V$", fontsize=10.5, color="0.30", ha="left", va="bottom")

    ax.text(
        0.5,
        -0.015,
        r"shells and US coast outward together; relative Doppler $\approx 0$",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10,
        style="italic",
        color="0.25",
    )


def main() -> None:
    fig = plt.figure(figsize=(7.4, 8.8), facecolor="white")
    fig.suptitle(
        "US between two emitters — comoving",
        fontsize=16,
        fontweight="bold",
        y=0.975,
    )
    gs = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=(1.0, 1.22),
        left=0.08,
        right=0.94,
        top=0.92,
        bottom=0.07,
        hspace=0.34,
    )

    ax_plates = fig.add_subplot(gs[0, 0])
    ax_shells = fig.add_subplot(gs[1, 0])
    draw_plates(ax_plates)
    draw_shells(ax_shells)

    for out in (OUT_CYCLE, OUT_LETTER):
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=180, bbox_inches="tight", facecolor="white")
        print("wrote", out)
    plt.close(fig)


if __name__ == "__main__":
    main()
```

## 5. Full sources — best-available approximates (do not pixel-match)

Letter-era one-offs missing for Fig03/04/09/10/A1/A2/A2b. These `coe-figures` scripts regenerate *something* visually related; delivery INVENTORY explicitly did **not** ship them as exact.

### 5.1 Fig03 approx — `fig03_migkas_slices.py`

```python
# source: /workspace/coe-figures/src/fig03_migkas_slices.py
#!/usr/bin/env python3
"""Fig 3 — Migkas A(ℓ) in redshift slices + D–core-mass degeneracy."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from physics import A_F, TH1, load_fit_curves, load_migkas_fig6, save_fig, R_VISION, H_MEAN
from paths import CMB, MIGKAS_FIG7_CSV, BUBBLE


def _load_slice_digitised():
    """Fall back: use fig6/fig7 digitised as all-z; synthesize slice-like variants."""
    # Prefer cycle3 slices.json if present
    p = CMB / "cycle3/slices.json"
    if p.exists():
        return json.loads(p.read_text())
    return None


def _synthetic_slices_from_model():
    """Build three panels using digitised all-z + model at different Rs depths.

    Slice medians approximate Migkas redshift bins via Hubble law.
    """
    ell, A, e = load_migkas_fig6()
    th, A_newt, A_free, A_rel, meta = load_fit_curves()
    # Approximate median distances (Gpc) for z~0.05, z~0.1, all-z~0.06
    # Rs = c z / H / 1000
    c = 299792.458
    slices = {
        "0.03–0.07": dict(zmed=0.05, Rs=c * 0.05 / H_MEAN / 1e3, data=(ell, A * 0.98 + 0.02, e)),
        r"$z>0.07$": dict(zmed=0.10, Rs=c * 0.10 / H_MEAN / 1e3, data=(ell, A * 1.01 - 0.01, e)),
        "all-$z$": dict(zmed=0.06, Rs=R_VISION, data=(ell, A, e)),
    }
    # Try to recompute model at each Rs with fixed geometry
    models = {}
    try:
        import os, sys
        cwd = os.getcwd()
        sys.path.insert(0, str(BUBBLE))
        os.chdir(BUBBLE)
        import newt_omb as N
        from bubble_core import A_from, shifted
        mn = meta["newt"]
        # recover kin,kout from R_in,R_out
        kin = (mn["D"] - mn["R_in"]) / R_VISION
        kout = (mn["R_out"] - mn["D"]) / R_VISION
        for name, s in slices.items():
            Rs = min(max(s["Rs"], 0.2), R_VISION)
            a = N.geo(mn["D"], mn["mcc"], kin, kout, Rs)
            out = N.af_of(a)
            if out is None:
                models[name] = A_newt
            else:
                af = out[0]
                models[name] = A_from(shifted(af, mn["l0"]), TH1 - mn["l0"], mn["amp"], mn["phi"])
        os.chdir(cwd)
    except Exception as exc:
        print("slice model recompute failed, using headline newt:", exc)
        for name in slices:
            models[name] = A_newt
    return slices, models, meta


def _degeneracy_panel(ax, meta):
    """D vs core mass (mcc) trade-off schematic for fixed cluster a2."""
    D = np.linspace(15, 100, 200)
    # Schematic degeneracy: mcc ~ const / D^α to hold a2; use letter anchor mcc~0.0017 at D=49
    mcc0, D0 = meta["newt"]["mcc"], meta["newt"]["D"]
    # Approximate: Mc fixed ⇒ mcc ∝ 1/D^3; also shell cavity scales
    mcc_fixedMc = mcc0 * (D0 / D) ** 3
    mcc_a2 = mcc0 * (D0 / D) ** 1.2  # softer track for a2-held
    ax.plot(D, mcc_fixedMc, "k--", lw=1.2, label=r"fixed $M_c$ ($m_{cc}\propto D^{-3}$)")
    ax.plot(D, mcc_a2, color="C3", lw=1.8, label=r"$a_2=0.084$ held (schematic)")
    ax.plot(meta["newt"]["D"], meta["newt"]["mcc"], "o", color="C3", ms=8, label=r"Newt $\Omega_s=0.049$")
    ax.plot(meta["rel"]["D"], meta["rel"]["mcc"], "s", color="C1", ms=7, label="Rel. LTB")
    ax.set_xlabel(r"$D$ [Gpc]")
    ax.set_ylabel(r"$m_{cc}=M_c/M_{\mathrm{crit}}(D)$")
    ax.set_yscale("log")
    ax.set_title(r"$D$–core-mass degeneracy", loc="left", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False)
    ax.grid(True, alpha=0.3)


def main():
    slices_json = _load_slice_digitised()
    fig = plt.figure(figsize=(10.5, 7.2), dpi=160)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0], hspace=0.35, wspace=0.3)

    if slices_json and "slices" in slices_json:
        # Real per-slice profiles if available
        order = [("0.03-0.07", "0.03–0.07"), ("z>0.07", r"$z>0.07$"), ("all", "all-$z$")]
        th, A_newt, _, _, meta = load_fit_curves()
        for i, (key, title) in enumerate(order):
            ax = fig.add_subplot(gs[0, i])
            e = slices_json["slices"][key]
            A = np.asarray(e["A"], float)
            # LC from migkas_A style: assume 5-deg steps 0..360
            lc = np.arange(0, 360, 5.0)[: len(A)]
            err = np.asarray(e.get("A_err", np.full_like(A, 0.02)), float)
            ok = np.isfinite(A)
            ax.errorbar(lc[ok], A[ok], yerr=err[ok], fmt="o", ms=2.5, color="0.4", ecolor="0.7", lw=0.7)
            ax.plot(th, A_newt, color="C3", lw=1.4, label="CoE model @ slice depth")
            ax.set_ylim(0.7, 1.25)
            ax.set_xlim(0, 360)
            ax.set_title(f"{title}\n$z_{{\\mathrm{{med}}}}={e.get('z_median', float('nan')):.3f}$", fontsize=10)
            ax.set_xlabel(r"$\ell$ (deg)")
            if i == 0:
                ax.set_ylabel(r"$A/A_{\mathrm{all}}$")
            ax.legend(fontsize=7, frameon=False)
        axd = fig.add_subplot(gs[1, :])
        _degeneracy_panel(axd, meta)
    else:
        slices, models, meta = _synthetic_slices_from_model()
        for i, (name, s) in enumerate(slices.items()):
            ax = fig.add_subplot(gs[0, i])
            ell, A, e = s["data"]
            ax.errorbar(ell, A, yerr=e, fmt="o", ms=2.8, color="0.4", ecolor="0.7", lw=0.7, label="Migkas (digitised)")
            ax.plot(TH1, models[name], color="C3", lw=1.5, label=rf"CoE @ $R_s$={s['Rs']:.2f} Gpc")
            ax.plot(TH1, A_F(TH1), color="C0", lw=1.0, alpha=0.7, label=r"$A_F$")
            ax.set_ylim(0.65, 1.25)
            ax.set_xlim(0, 360)
            ax.set_title(f"{name}\n$z_{{\\mathrm{{med}}}}\\approx{s['zmed']:.2f}$", fontsize=10)
            ax.set_xlabel(r"$\ell$ (deg)")
            if i == 0:
                ax.set_ylabel(r"$A/A_{\mathrm{all}}$")
            ax.legend(fontsize=7, frameon=False)
        axd = fig.add_subplot(gs[1, :])
        _degeneracy_panel(axd, meta)

    fig.suptitle(
        "Migkas $A(\\ell)$ in redshift slices + $D$–core-mass degeneracy",
        fontsize=12,
    )
    p = save_fig(fig, "fig03_migkas_slices", also_pdf=True)
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()
```
### 5.2 Fig04 approx — `fig04_amp_full.py`

```python
# source: /workspace/coe-figures/src/fig04_amp_full.py
#!/usr/bin/env python3
"""Fig 4 — Full Migkas A/Aall profile + residual panel for Newt and rel fits."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from physics import (
    A_F,
    BAND_CENTRES,
    BAND_COLORS,
    BAND_HALF,
    load_fit_curves,
    load_migkas_fig6,
    r2_score,
    save_fig,
)


def main():
    ell_d, A_d, e_d = load_migkas_fig6()
    th, A_newt, A_free, A_rel, meta = load_fit_curves()
    th2 = np.append(th, 360.0)
    cl = lambda y: np.append(y, y[0])
    Af = A_F(th)

    fig, (a, b) = plt.subplots(
        2, 1, figsize=(8.0, 6.2), sharex=True, gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.07}
    )
    for L, col in zip(BAND_CENTRES, BAND_COLORS):
        for ax in (a, b):
            ax.axvspan(L - BAND_HALF, L + BAND_HALF, color=col, alpha=0.20, lw=0, zorder=0)

    a.errorbar(ell_d, A_d, yerr=e_d, fmt="o", ms=4, color="0.35", ecolor="0.65", lw=0.9, zorder=3, label="Migkas (digitised)")
    a.plot(th2, A_F(th2), color="C0", lw=1.7, label=r"$A_F$", zorder=2)
    mn, mr = meta["newt"], meta["rel"]
    a.plot(
        th2,
        cl(A_newt),
        color="#1a7f37",
        lw=2.0,
        label=rf"Newtonian $\Omega_s$={mn['Os']}: $D$={mn['D']:.1f}, $R_{{\rm in}}$={mn['R_in']:.1f}, $R_{{\rm out}}$={mn['R_out']:.1f} Gpc, $\ell_0$={mn['l0']:.0f}° ($R^2$={mn['R2']:.3f})",
    )
    a.plot(
        th2,
        cl(A_rel),
        color="C1",
        lw=1.8,
        ls="--",
        label=rf"Relativistic LTB $\Omega_s$={mr['Os']}: $D$={mr['D']:.1f}, $R_{{\rm in}}$={mr['R_in']:.1f}, $R_{{\rm out}}$={mr['R_out']:.1f} Gpc, $\ell_0$={mr['l0']:.1f}° ($R^2$={mr['R2']:.3f})",
    )
    a.set_ylabel(r"$A/A_{\mathrm{all}}$", fontsize=12)
    a.set_ylim(0.60, 1.24)
    a.legend(fontsize=7.2, loc="lower left", frameon=False)
    a.set_title("Full Migkas $A/A_{\\mathrm{all}}$ profile", loc="left", fontsize=12)

    # Residuals vs data (interp)
    def res(model):
        return 100 * (A_d - np.interp(ell_d, th, model, period=360)) / A_d

    b.axhline(0, color="k", lw=0.8)
    b.plot(ell_d, res(A_newt), "o-", ms=3.5, color="#1a7f37", lw=1.0, label="data − Newt")
    b.plot(ell_d, res(A_rel), "s--", ms=3.5, color="C1", lw=1.0, label="data − Rel")
    b.set_ylabel("residual [%]", fontsize=11)
    b.set_xlabel(r"Galactic longitude $\ell$ (deg)", fontsize=12)
    b.set_xlim(0, 360)
    b.set_xticks(range(0, 361, 45))
    b.legend(fontsize=8, frameon=False, loc="upper right")
    b.grid(True, alpha=0.25)

    p = save_fig(fig, "fig04_amp_full", also_pdf=True)
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()
```
### 5.3 Fig05 analytic approx (no Migkas underlay) — `fig05_sky_orientation.py`

```python
# source: /workspace/coe-figures/src/fig05_sky_orientation.py
#!/usr/bin/env python3
"""Fig 5 — Sky orientation sequence: untilted → tilt → gradient → roll."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from physics import (
    BAND_CENTRES,
    BAND_COLORS,
    BAND_HALF,
    BAND_KINDS,
    ROLL_DEG,
    TILT_DEG,
    band_circle_xyz,
    mollweide_xy,
    rotate_xyz,
    save_fig,
    xyz_to_galactic,
)


def _draw_model_frame(ax, tilt, roll, gradient=False):
    ax.set_aspect("equal")
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.axis("off")
    circ = Circle((0, 0), 1.0, fill=False, ec="0.4", lw=1.4)
    ax.add_patch(circ)
    ax.plot(0, 0, "ko", ms=5)
    # CoE–anti-CoE axis in the drawing plane (x horizontal)
    ax.annotate("", xy=(1.15, 0), xytext=(-1.15, 0), arrowprops=dict(arrowstyle="<->", color="#6c3483", lw=1.5))
    ax.text(-1.28, 0.12, "CoE (−x)", color="#6c3483", fontsize=7, ha="left")
    ax.text(1.28, 0.12, "anti-CoE (+x)", color="#6c3483", fontsize=7, ha="right")

    # After tilt/roll, project band midplanes onto the x–z plane of the *drawing*
    # (model frame before Galactic embedding): show θ wedges as chords.
    for th, col, kind in zip(BAND_CENTRES, BAND_COLORS, BAND_KINDS):
        # chord at polar angle th from +x
        t = np.radians(th)
        # In 2D side view (x horizontal, z vertical): band appears as two points at ±sin
        # Draw as thick arc segments on the circle at that polar angle from +x.
        # Polar angle from +x: point (cos th, sin th) in x–z.
        # For a belt, draw a short radial tick band across the circle.
        alpha = 0.35
        if gradient:
            # fade on CoE side (θ near 180)
            dist = min(abs(th - 180), abs(th - 180 + 360), abs(th - 180 - 360))
            alpha = 0.15 + 0.55 * (1 - dist / 180.0)
        # draw two symmetric marks at ±y
        for sgn in (+1, -1):
            x = np.cos(t)
            z = sgn * np.sin(t)
            # rotate mark for tilt visualization in side view: tilt lifts anti-CoE
            tt = np.radians(tilt)
            xr = x * np.cos(tt) - z * np.sin(tt)
            zr = x * np.sin(tt) + z * np.cos(tt)
            # roll mixes y into z — approximate by rotating the ± marks
            rr = np.radians(roll)
            # treat current (0, z) as yz and roll
            yr = 0.0 * np.cos(rr) - zr * np.sin(rr)  # unused in xz
            zr2 = 0.0 * np.sin(rr) + zr * np.cos(rr)
            ax.plot([0.82 * xr, 1.05 * xr], [0.82 * zr2, 1.05 * zr2], color=col, lw=4.5, alpha=alpha, solid_capstyle="round")
        # label once
        xl = 1.18 * np.cos(t)
        zl = 1.18 * np.sin(t)
        ax.text(xl, zl, rf"$\theta$={th:g}°", color=col, fontsize=6.5, ha="center", va="center", alpha=min(1, alpha + 0.3))


def _draw_mollweide(ax, tilt, roll, gradient=False):
    # outline
    phi = np.linspace(-np.pi, np.pi, 400)
    # Mollweide boundary
    bx = 2 * np.sqrt(2) / np.pi * phi  # will flip
    by = np.zeros_like(phi)  # equator placeholder
    # proper ellipse outline
    t = np.linspace(0, 2 * np.pi, 400)
    ox = 2 * np.sqrt(2) * np.cos(t)
    oy = np.sqrt(2) * np.sin(t)
    ax.plot(-ox, oy, "k-", lw=0.8)  # lon increases left → already in -x of our mollweide_xy

    # celestial equator (b=0)
    lon = np.linspace(0, 360, 400)
    x, y = mollweide_xy(lon, np.zeros_like(lon))
    ax.plot(x, y, color="#1a7f37", lw=0.9, ls="--", alpha=0.8)

    for th, col, kind in zip(BAND_CENTRES, BAND_COLORS, BAND_KINDS):
        band = band_circle_xyz(th, BAND_HALF, n=720)
        for edge_name, lw in (("mid", 2.2), ("lo", 0.9), ("hi", 0.9)):
            xyz = rotate_xyz(band[edge_name], tilt_deg=tilt, roll_deg=roll)
            ell, b = xyz_to_galactic(xyz)
            # split at longitude wrap
            x, y = mollweide_xy(ell, b)
            # break segments where jumps
            d = np.hypot(np.diff(x), np.diff(y))
            breaks = np.where(d > 0.35)[0]
            idx = np.r_[0, breaks + 1, len(x)]
            for i0, i1 in zip(idx[:-1], idx[1:]):
                if i1 - i0 < 2:
                    continue
                xs, ys = x[i0:i1], y[i0:i1]
                if gradient:
                    # fade by Galactic longitude near CoE (ℓ~180)
                    # use mean ell of segment
                    alpha = 0.25 + 0.6 * (1 - abs(((ell[i0:i1].mean() - 0 + 180) % 360) - 180) / 180)
                    # stronger near anti-CoE ℓ~0
                    mid = ell[i0:i1].mean()
                    # distance to ℓ=0 wrap
                    d0 = min(mid % 360, 360 - (mid % 360))
                    alpha = 0.2 + 0.75 * (1 - d0 / 180.0)
                else:
                    alpha = 0.9 if edge_name == "mid" else 0.45
                ax.plot(xs, ys, color=col, lw=lw if edge_name == "mid" else 0.8, alpha=alpha)

    # anti-CoE marker
    ac = rotate_xyz(np.array([[1.0, 0.0, 0.0]]), tilt_deg=tilt, roll_deg=roll)
    ell0, b0 = xyz_to_galactic(ac)
    x0, y0 = mollweide_xy(ell0, b0)
    ax.plot(x0, y0, "o", color="#7b241c", ms=6, zorder=5)
    ax.text(x0[0], y0[0] + 0.12, "anti-CoE", color="#7b241c", fontsize=7, ha="center")

    ax.set_aspect("equal")
    ax.set_xlim(-3.0, 3.0)
    ax.set_ylim(-1.7, 1.7)
    ax.axis("off")
    # lon labels
    for L, lab in ((180, "180"), (270, "270"), (0, "0"), (90, "90"), (180, "180")):
        pass
    for L, xpos in ((180, -2.0), (270, -1.0), (0, 0.0), (90, 1.0), (180, 2.0)):
        ax.text(xpos, -1.55, lab if L != 180 or xpos < 0 else "180", fontsize=7, ha="center", color="0.35")
    ax.text(0, -1.72, r"Galactic $\ell$ (increases ←)", fontsize=7, ha="center", color="0.35")


def _draw_sphere3d(ax, tilt, roll, gradient=False):
    # wire sphere
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 20)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, color="0.75", lw=0.3, alpha=0.5)
    # axis
    ax.plot([-1.2, 1.2], [0, 0], [0, 0], color="#6c3483", lw=1.5)
    for th, col in zip(BAND_CENTRES, BAND_COLORS):
        band = band_circle_xyz(th, 0.0, n=200)["mid"]
        xyz = rotate_xyz(band, tilt_deg=tilt, roll_deg=roll)
        alpha = 0.9
        if gradient:
            d0 = min(abs(th - 0), abs(th - 360))
            # th from anti-CoE: small th = anti-CoE side
            alpha = 0.25 + 0.7 * (1 - min(th, 360 - th) / 180.0) if th > 90 else 0.95
            if th > 90:
                alpha = 0.35
        ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=col, lw=2.0, alpha=alpha)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(-1.2, 1.2)
    ax.set_axis_off()
    ax.view_init(elev=18, azim=35)


def main():
    rows = [
        (0.0, 0.0, False, "Untilted, no roll"),
        (TILT_DEG, 0.0, False, rf"+{TILT_DEG:g}° tilt"),
        (TILT_DEG, 0.0, True, rf"+{TILT_DEG:g}° tilt, gradient"),
        (TILT_DEG, ROLL_DEG, True, rf"+{TILT_DEG:g}° tilt, gradient, +{ROLL_DEG:g}° roll"),
    ]
    fig = plt.figure(figsize=(12.5, 14.0), dpi=150)
    for i, (tilt, roll, grad, title) in enumerate(rows):
        ax0 = fig.add_subplot(4, 3, 3 * i + 1)
        _draw_model_frame(ax0, tilt, roll, gradient=grad)
        ax0.set_title("Model frame" if i == 0 else "", fontsize=10)
        ax0.text(0.0, 1.22, title, ha="center", fontsize=9, transform=ax0.transData, clip_on=False)

        ax1 = fig.add_subplot(4, 3, 3 * i + 2)
        _draw_mollweide(ax1, tilt, roll, gradient=grad)
        ax1.set_title("Mollweide (Galactic)" if i == 0 else "", fontsize=10)

        ax2 = fig.add_subplot(4, 3, 3 * i + 3, projection="3d")
        _draw_sphere3d(ax2, tilt, roll, gradient=grad)
        ax2.set_title("Sphere 3D" if i == 0 else "", fontsize=10)

    fig.suptitle(
        "Sky orientation sequence: untilted → +15° tilt → gradient → +80.8° roll",
        fontsize=13,
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    p = save_fig(fig, "fig05_sky_orientation", also_pdf=False, dpi=160)
    # also jpg-friendly alias naming: fig05 already png
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()
```
### 5.4 Fig09 approx — `fig09_kill.py`

```python
# source: /workspace/coe-figures/src/fig09_kill.py
#!/usr/bin/env python3
"""Fig 9 — Unit-normalized kill function vs plate/shell/shear residuals at three D."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from physics import (
    extend_0_360,
    kill_from_doppler,
    load_dtt,
    load_vs_kill,
    save_fig,
    unit_norm,
)


def main():
    dtt = load_dtt()
    rows = dtt["rows"]
    vs = load_vs_kill()

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2), dpi=160, sharey=True)
    for ax, row in zip(axes, rows):
        th, kill = kill_from_doppler(row)
        th360, kill360 = extend_0_360(th, kill)
        kn = unit_norm(kill360)
        ax.plot(th360, kn, color="k", lw=2.2, label="kill (unit-norm)")

        # exploratory shapes from vs_kill if available
        comps = vs.get("comparisons", [])
        D = row["D"]
        # find nearest comparison or synthesize plate/shell/shear toys
        phi = np.radians(th360)
        # plate-like step residual (odd about 90)
        plate = unit_norm(np.cos(phi))
        shell = unit_norm(np.cos(2 * phi))
        shear = unit_norm(np.sin(2 * phi))
        ax.plot(th360, plate, color="C0", lw=1.3, ls="--", alpha=0.85, label="plate-like")
        ax.plot(th360, shell, color="C1", lw=1.3, ls=":", alpha=0.85, label="shell/P2-like")
        ax.plot(th360, shear, color="C3", lw=1.2, ls="-.", alpha=0.85, label="shear-like")

        ax.set_title(rf"$D={D:.2f}$ Gpc", fontsize=11)
        ax.set_xlabel(r"$\phi$ (deg) on great circle")
        ax.set_xlim(0, 360)
        ax.set_xticks([0, 90, 180, 270, 360])
        ax.axhline(0, color="0.5", lw=0.6)
        ax.grid(True, alpha=0.25)
        if ax is axes[0]:
            ax.set_ylabel("unit-normalized amplitude")
            ax.legend(fontsize=7.5, frameon=False, loc="lower left")

    fig.suptitle(
        "Kill function vs exploratory plate / shell / shear residuals",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    p = save_fig(fig, "fig09_kill", also_pdf=True)
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()
```
### 5.5 Fig10 approx — `fig10_cancel_partner.py`

```python
# source: /workspace/coe-figures/src/fig10_cancel_partner.py
#!/usr/bin/env python3
"""Fig 10 — CoE Doppler + multipole partner residuals (dipole+P2 and ℓ≤4)."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from physics import kill_from_doppler, load_dtt, partner_cancel, save_fig, unit_norm


def main():
    dtt = load_dtt()
    # use D≈49 Gpc row (headline-like); also show D=22.28
    prefer = sorted(dtt["rows"], key=lambda r: abs(r["D"] - 49.0))[:2]
    if len(prefer) < 2:
        prefer = dtt["rows"][:2]

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.5), dpi=160, sharex="col")
    for col, row in enumerate(prefer):
        th = np.asarray(row["theta_toward_deg"], float)
        dop = np.asarray(row["dtt_doppler"], float)
        dop0 = dop - dop.mean()
        _, res_p2, a = partner_cancel(dop0, th, mode="dip_p2")
        _, res_l4, _ = partner_cancel(dop0, th, mode="l_le_4")

        ax0 = axes[0, col]
        ax0.plot(th, dop0, color="C0", lw=1.8, label="Doppler (demeaned)")
        ax0.plot(th, -dop0, color="0.35", lw=1.2, ls="--", label="−Doppler (= kill shape)")
        ax0.set_title(rf"$D={row['D']:.2f}$ Gpc", fontsize=11)
        ax0.set_ylabel(r"$\Delta T/T$")
        ax0.legend(fontsize=8, frameon=False)
        ax0.grid(True, alpha=0.25)
        ax0.axhline(0, color="k", lw=0.5)

        ax1 = axes[1, col]
        ax1.plot(th, res_p2, color="C3", lw=1.8, label=r"residual after dipole+$P_2$")
        ax1.plot(th, res_l4, color="C2", lw=1.6, ls="--", label=r"residual after $\ell\leq 4$")
        ax1.axhline(0, color="k", lw=0.5)
        ax1.set_xlabel(r"$\theta$ from CoE (deg)")
        ax1.set_ylabel(r"$\Delta T/T$ residual")
        ax1.legend(fontsize=8, frameon=False)
        ax1.grid(True, alpha=0.25)
        # annotate a2
        ax0.text(
            0.02,
            0.05,
            rf"$a_2^{{\mathrm{{Dop}}}}={row.get('doppler_a2', a[2]):.3f}$",
            transform=ax0.transAxes,
            fontsize=8,
            color="0.3",
        )

    fig.suptitle("CoE Doppler + multipole partner cancellation", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    p = save_fig(fig, "fig10_cancel_partner", also_pdf=True)
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()
```
### 5.6 FigA1 approx — `figA1_ltb_shell.py`

```python
# source: /workspace/coe-figures/src/figA1_ltb_shell.py
#!/usr/bin/env python3
"""Fig A1 — LTB shell evolution: flat m+r; E=0; E>0; E<0."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from physics import save_fig


def evolve_dust(E: float, t_max: float = 8.0, a0: float = 1e-4):
    """Integrate ȧ² = 2E + 2GM/a with GM=1 units; a(0)=a0, ȧ>0."""
    GM = 1.0

    def f(t, y):
        a, = y
        rhs = 2 * E + 2 * GM / max(a, 1e-12)
        if rhs < 0:
            return [0.0]
        return [np.sqrt(rhs)]

    sol = solve_ivp(f, (0, t_max), [a0], dense_output=True, max_step=0.01, rtol=1e-6)
    t = np.linspace(0, min(t_max, sol.t[-1]), 600)
    a = sol.sol(t)[0]
    # stop if recollapse starts
    if E < 0:
        # analytic cycloid max a = GM/|E|
        amax = GM / abs(E)
        # continue through recollapse with signed root
        def f2(t, y):
            a, = y
            rhs = 2 * E + 2 * GM / max(a, 1e-12)
            # use phase from energy: after turnaround ȧ < 0
            return [-np.sqrt(max(rhs, 0.0))]

        # find turnaround
        i_turn = int(np.argmax(a))
        if a[i_turn] > 0.5 * amax:
            t_turn = t[i_turn]
            sol2 = solve_ivp(f2, (t_turn, t_max), [a[i_turn]], dense_output=True, max_step=0.01)
            t2 = np.linspace(t_turn, min(t_max, sol2.t[-1]), 400)
            a2 = sol2.sol(t2)[0]
            t = np.concatenate([t[: i_turn + 1], t2[1:]])
            a = np.concatenate([a[: i_turn + 1], a2[1:]])
    return t, a


def evolve_flat_mr(t_max=8.0):
    """Schematic flat matter+radiation: early a~t^{1/2}, late a~t^{2/3}."""
    t = np.linspace(1e-3, t_max, 800)
    # match at teq=1: continuous a and prefer radiation then matter
    teq = 1.0
    a_eq = teq ** 0.5
    a = np.where(t < teq, t ** 0.5, a_eq * (t / teq) ** (2 / 3))
    return t, a


def main():
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.5), dpi=160)
    panels = [
        (axes[0, 0], "Flat matter + radiation", evolve_flat_mr()),
        (axes[0, 1], r"Dust $E=0$ (marginally bound)", evolve_dust(0.0)),
        (axes[1, 0], r"Dust $E>0$ (expand forever)", evolve_dust(+0.15)),
        (axes[1, 1], r"Dust $E<0$ (recollapse)", evolve_dust(-0.25, t_max=10.0)),
    ]
    for ax, title, (t, a) in panels:
        ax.plot(t, a, color="C0", lw=2.0)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r"$t$ (code units)")
        ax.set_ylabel(r"$a(t)$ / areal factor")
        ax.grid(True, alpha=0.3)
        if "radiation" in title.lower():
            ax.axvline(1.0, color="0.5", ls=":", lw=1)
            ax.text(1.05, ax.get_ylim()[1] * 0.1 if False else 0.3, r"$t_{\mathrm{eq}}$", fontsize=8, color="0.4")

    # annotate asymptotics on first panel
    axes[0, 0].text(0.35, 0.55, r"$a\sim t^{1/2}$", fontsize=9, color="0.35", transform=axes[0, 0].transAxes)
    axes[0, 0].text(0.62, 0.75, r"$a\sim t^{2/3}$", fontsize=9, color="0.35", transform=axes[0, 0].transAxes)

    fig.suptitle("LTB shell evolution (Einstein/Friedmann shell constraint)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    p = save_fig(fig, "figA1_ltb_shell", also_pdf=True)
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()
```
### 5.7 FigA2 approx — `figA2_coe_sun.py`

```python
# source: /workspace/coe-figures/src/figA2_coe_sun.py
#!/usr/bin/env python3
"""Fig A2 — CoE and Sun timing from Plano through 2026."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time
import astropy.units as u
from astropy.utils import iers

iers.conf.auto_download = False

from physics import COE_GAL_L, COE_GAL_B, PLANO_LAT, PLANO_LON, PLANO_H, save_fig


def find_events(coord, site, year, n_days=365):
    """Approximate rise / culmination / set local hours for each day."""
    # sample every 4 minutes over each day in UTC
    doys = np.arange(1, n_days + 1)
    rise_h, culm_h, set_h = [], [], []
    culm_alt = []
    lon_hours = site.lon.to(u.hourangle).value
    t0 = Time(f"{year}-01-01 00:00:00", scale="utc")
    for d in doys:
        hours = np.arange(0, 24, 4 / 60)
        # interpret hours as Local Mean Time, convert to UTC
        utc_days = (d - 1) + hours / 24.0 - lon_hours / 24.0
        times = t0 + utc_days * u.day
        aa = coord.transform_to(AltAz(obstime=times, location=site))
        alt = aa.alt.deg
        az = aa.az.deg
        # culmination = max alt
        i_c = int(np.argmax(alt))
        culm_h.append(hours[i_c])
        culm_alt.append(alt[i_c])
        # rise: first upward zero crossing
        rise = set_ = np.nan
        for i in range(1, len(alt)):
            if alt[i - 1] < 0 <= alt[i]:
                # linear interp
                f = -alt[i - 1] / (alt[i] - alt[i - 1] + 1e-12)
                rise = hours[i - 1] + f * (hours[i] - hours[i - 1])
                break
        for i in range(1, len(alt)):
            if alt[i - 1] > 0 >= alt[i]:
                f = alt[i - 1] / (alt[i - 1] - alt[i] + 1e-12)
                set_ = hours[i - 1] + f * (hours[i] - hours[i - 1])
        rise_h.append(rise)
        set_h.append(set_)
    return (
        doys,
        np.asarray(rise_h),
        np.asarray(culm_h),
        np.asarray(set_h),
        np.asarray(culm_alt),
    )


def main():
    year = 2026
    site = EarthLocation(lat=PLANO_LAT * u.deg, lon=PLANO_LON * u.deg, height=PLANO_H * u.m)
    coe = SkyCoord(l=COE_GAL_L * u.deg, b=COE_GAL_B * u.deg, frame="galactic").icrs

    # CoE is essentially fixed RA/Dec → rise/set nearly constant in LMT
    doy, cr, cc, cs, ca = find_events(coe, site, year, n_days=365)

    # Sun: sample weekly for speed then interpolate display
    sun_doys = np.arange(1, 366, 3)
    sr, sc, ss = [], [], []
    t0 = Time(f"{year}-01-01 00:00:00", scale="utc")
    lon_hours = site.lon.to(u.hourangle).value
    for d in sun_doys:
        hours = np.arange(0, 24, 4 / 60)
        utc_days = (d - 1) + hours / 24.0 - lon_hours / 24.0
        times = t0 + utc_days * u.day
        sun = get_sun(times)
        aa = sun.transform_to(AltAz(obstime=times, location=site))
        alt = aa.alt.deg
        i_c = int(np.argmax(alt))
        sc.append(hours[i_c])
        rise = set_ = np.nan
        for i in range(1, len(alt)):
            if alt[i - 1] < 0 <= alt[i]:
                f = -alt[i - 1] / (alt[i] - alt[i - 1] + 1e-12)
                rise = hours[i - 1] + f * (hours[i] - hours[i - 1])
                break
        for i in range(1, len(alt)):
            if alt[i - 1] > 0 >= alt[i]:
                f = alt[i - 1] / (alt[i - 1] - alt[i] + 1e-12)
                set_ = hours[i - 1] + f * (hours[i] - hours[i - 1])
        sr.append(rise)
        ss.append(set_)
    sr, sc, ss = map(np.asarray, (sr, sc, ss))

    fig, ax = plt.subplots(figsize=(11.0, 5.2), dpi=160)
    ax.plot(doy, cr, color="#6c3483", lw=1.8, label="CoE rise")
    ax.plot(doy, cc, color="#8e44ad", lw=1.8, ls="--", label="CoE culmination")
    ax.plot(doy, cs, color="#5b2c6f", lw=1.8, ls=":", label="CoE set")
    ax.plot(sun_doys, sr, color="#e67e22", lw=1.6, label="Sun rise")
    ax.plot(sun_doys, sc, color="#d35400", lw=1.6, ls="--", label="Sun culmination")
    ax.plot(sun_doys, ss, color="#a04000", lw=1.6, ls=":", label="Sun set")
    ax.set_xlabel("Day of year 2026")
    ax.set_ylabel("Local mean time (hours)")
    ax.set_xlim(1, 365)
    ax.set_ylim(0, 24)
    ax.set_yticks(range(0, 25, 3))
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2, frameon=False, loc="upper right")
    ax.set_title(
        "CoE and Sun timing from Plano, TX (2026)\n"
        rf"CoE at Galactic $(\ell,b)=({COE_GAL_L:.0f}^\circ,{COE_GAL_B:.0f}^\circ)$",
        loc="left",
        fontsize=11,
    )
    p = save_fig(fig, "figA2_coe_sun", also_pdf=True)
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()
```
### 5.8 FigA2b / crossings approx — `figA_table_crossings.py`

```python
# source: /workspace/coe-figures/src/figA_table_crossings.py
#!/usr/bin/env python3
"""Table A1 visual — Sun/CoE event crossings Central Texas 2027."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time
import astropy.units as u
from astropy.utils import iers

iers.conf.auto_download = False

from physics import COE_GAL_L, COE_GAL_B, PLANO_LAT, PLANO_LON, PLANO_H, save_fig


def daily_events(coord_fn, site, year):
    """Return rise, culm, set LMT hours for each DOY (coord_fn(times)->SkyCoord)."""
    lon_h = site.lon.to(u.hourangle).value
    t0 = Time(f"{year}-01-01 00:00:00", scale="utc")
    doys = np.arange(1, 366)
    out = {k: np.full(len(doys), np.nan) for k in ("rise", "culm", "set")}
    hours = np.arange(0, 24, 5 / 60)
    for i, d in enumerate(doys):
        utc_days = (d - 1) + hours / 24.0 - lon_h / 24.0
        times = t0 + utc_days * u.day
        aa = coord_fn(times).transform_to(AltAz(obstime=times, location=site))
        alt = aa.alt.deg
        out["culm"][i] = hours[int(np.argmax(alt))]
        for j in range(1, len(alt)):
            if alt[j - 1] < 0 <= alt[j] and np.isnan(out["rise"][i]):
                f = -alt[j - 1] / (alt[j] - alt[j - 1] + 1e-12)
                out["rise"][i] = hours[j - 1] + f * (hours[j] - hours[j - 1])
            if alt[j - 1] > 0 >= alt[j] and np.isnan(out["set"][i]):
                f = alt[j - 1] / (alt[j - 1] - alt[j] + 1e-12)
                out["set"][i] = hours[j - 1] + f * (hours[j] - hours[j - 1])
    return doys, out


def crossings(doys, sun, coe, tol_h=0.15):
    """Find DOYs where Sun and CoE event times match within tol."""
    labels = ["rise", "culm", "set"]
    pairs = []
    for i, s_ev in enumerate(labels):
        for j, c_ev in enumerate(labels):
            ds = sun[s_ev] - coe[c_ev]
            # unwrap near day boundary not needed for LMT hours usually
            hit = np.where(np.isfinite(ds) & (np.abs(ds) < tol_h))[0]
            # also pick closest approach even if > tol (min |ds|)
            if len(hit) == 0:
                k = int(np.nanargmin(np.abs(ds)))
                if np.isfinite(ds[k]) and abs(ds[k]) < 0.75:
                    hit = np.array([k])
            for k in hit:
                pairs.append(
                    dict(
                        doy=int(doys[k]),
                        sun=s_ev,
                        coe=c_ev,
                        sun_h=float(sun[s_ev][k]),
                        coe_h=float(coe[c_ev][k]),
                        diag=(i == j),
                    )
                )
    # unique by (sun,coe) keeping closest
    best = {}
    for p in pairs:
        key = (p["sun"], p["coe"])
        err = abs(p["sun_h"] - p["coe_h"])
        if key not in best or err < abs(best[key]["sun_h"] - best[key]["coe_h"]):
            best[key] = p
    return sorted(best.values(), key=lambda p: (not p["diag"], p["doy"]))


def main():
    year = 2027
    site = EarthLocation(lat=PLANO_LAT * u.deg, lon=PLANO_LON * u.deg, height=PLANO_H * u.m)
    coe_icrs = SkyCoord(l=COE_GAL_L * u.deg, b=COE_GAL_B * u.deg, frame="galactic").icrs

    def coe_fn(times):
        return coe_icrs

    def sun_fn(times):
        return get_sun(times)

    # subsample for speed: every 2 days then refine
    doys, sun = daily_events(sun_fn, site, year)
    _, coe = daily_events(coe_fn, site, year)
    rows = crossings(doys, sun, coe, tol_h=0.25)
    # keep up to 9
    rows = rows[:9]

    fig, ax = plt.subplots(figsize=(10.5, 3.6), dpi=160)
    ax.axis("off")
    ax.set_title(
        f"Sun / CoE event crossings — Central Texas {year}\n"
        r"Highlighted diagonal: rise↔rise, culm↔culm, set↔set",
        fontsize=11,
        loc="left",
        pad=12,
    )

    col_labels = ["#", "DOY", "Date", "Sun event", "CoE event", "Sun LMT", "CoE LMT", "Δ (h)"]
    cell = []
    colors = []
    for n, r in enumerate(rows, 1):
        t = Time(f"{year}-01-01") + (r["doy"] - 1) * u.day
        date = t.strftime("%Y-%m-%d")
        cell.append(
            [
                str(n),
                str(r["doy"]),
                date,
                r["sun"],
                r["coe"],
                f"{r['sun_h']:.2f}",
                f"{r['coe_h']:.2f}",
                f"{r['sun_h']-r['coe_h']:+.2f}",
            ]
        )
        colors.append(["#fdebd0" if r["diag"] else "white"] * len(col_labels))

    if not cell:
        cell = [["—"] * len(col_labels)]
        colors = [["white"] * len(col_labels)]

    table = ax.table(cellText=cell, colLabels=col_labels, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.15, 1.45)
    for (i, j), cell_obj in table.get_celld().items():
        if i == 0:
            cell_obj.set_facecolor("#2c3e50")
            cell_obj.set_text_props(color="white", weight="bold")
        elif i >= 1 and i - 1 < len(colors):
            cell_obj.set_facecolor(colors[i - 1][j])

    p = save_fig(fig, "figA_table_crossings", also_pdf=True)
    # also alias requested name
    from paths import OUT
    import shutil

    shutil.copy(p, OUT / "figA1_table_crossings.png")
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    main()
```

### Related kill/cancel/amp/slices scripts found (not letter-exact)

| Path | Role |
|------|------|
| `.../cycle4/dtt/plot_kill_term.py` | kill term vs θ |
| `.../cycle4/dtt/plot_kill_term_0_360.py` | 0–360 kill |
| `.../cycle4/dtt/plot_kill_energy_pct_0_360.py` | energy % |
| `.../cycle4/dtt_rec_track/plot_kill_0_360.py` | rec-track kill |
| `.../cycle6_two_emitters/simple_cmb_two_quad_cancel.json` | data behind Fig10 asset |
| `.../cycle6_two_emitters/vs_kill_term_0_360.json` | shapes used in kill comparisons |

No letter-era `make_*amp*` / `make_*slices*` / `make_*kill*` / `make_*cancel*` script that emits the gold dims was found on disk.

---

## 6. Dependencies

| Package | Used by |
|---------|---------|
| `numpy` | nearly all |
| `matplotlib` | nearly all |
| `pillow` (PIL) | universe_thin web jpg; stacking; verify |
| `scipy` | universe_thin (`brentq`); some coe-figures |
| `astropy` | FigA3 sky path; FigA2/crossings approximates |
| `healpy` | **not required** for the exact delivery set (Mollweide in ring/overlay scripts is matplotlib-based) |

Suggested envs already on box:
- `/workspace/coe-figures/.venv` (matplotlib, numpy, scipy, astropy, pillow)
- system / `/workspace/.venv` as available

Fig01/02 additionally need the local bubble module tree and JSON fit stores under  
`/workspace/cosmology/papers/spreadsheet_model/bubble/`.

Fig06 overlay script needs the Migkas Fig.8 raster (see script / run log: `map_compare/FULL_migkas_theta75.png`).

---

## 7. Gold embed map

| Page | Extract | Dims | Delivery file |
|------|---------|------|---------------|
| 1 | gold-000.png | 1460×488 | Fig01_geometry.png |
| 2 | gold-002.png | 1182×1287 | Fig02_Achart.png |
| 3 | gold-004.png | 1123×1014 | Fig03_migkas_slices.png |
| 4 | gold-006.png | 1800×1385 | Fig04_amp_full.png |
| 5 | gold-007.jpg | 1370×1420 | Fig05_sky_orientation_four_row.jpg |
| 6 | gold-008.png | 2164×3222 | Fig06_overlay_stack.png |
| 7 | gold-009.png | 1045×1795 | Fig07_storyboard.png |
| 8 | gold-011.png | 1181×1908 | Fig08_two_emitters.png |
| 9 | gold-013.jpg | 1206×1280 | Fig09_kill.jpg |
| 10 | gold-014.png | 2084×2212 | Fig10_cancel_partner.png |
| 12 | gold-016.png | 1444×969 | FigA1_ltb_shell.png |
| 13 | gold-018.jpg | 1280×636 | FigA2_coe_sun.jpg |
| 13 | gold-019.jpg | 1239×570 | FigA2b_sun_coe_crossings.jpg |
| 14 | gold-020.png | 2415×1432 | FigA3_coe_sky.png |
| 14 | gold-022.jpg | 2080×944 | FigA_comic_HANDDRAWN_from_PDF.jpg |
| — | (n/a) | 883×945 | eclipse_inner_shell.png |

Odd-numbered `gold-00{1,3,5,…}.png` extracts are soft masks / tiny embeds — skipped in delivery.

---

*Generated from `/workspace/figs_exact_from_code/delivery/INVENTORY.md` + on-disk scripts. No fabricated matching generators.*
