# figs_exact_from_code — inventory

Authoritative letter PDF: `attachments/02b6cafe…pdf` (14 pages, Sep 30 2026).
Gold targets: `pdfimages -all` extracts in `gold_extracts/` (smasks skipped as soft masks).
Delivery folder: `/workspace/figs_exact_from_code/delivery/`

## Hard rules followed
- Exact visual match to PDF embeds (pixel-identical; bitwise where JPEG/PNG encoding matches).
- Prefer re-run of original script, else copy on-disk asset that already matches embed.
- Did **not** ship prior `figs_regenerated_from_code/` approximates (esp. wrong Fig3 / coe-figures recomputes).
- No whiteouts. No Amp×sin squiggly figure (`money_shot_amp_full_SINE_BAD.png` avoided).
- Hand-drawn comic: PDF original only; labeled HANDDRAWN.
- Eclipse: ran attachment `.py`.
- Did **not** modify/delete CoE_calculations.md or reference attachments.

## Results

| File | Label | Method | Pixel=gold | Bitwise=gold | Status |
|------|-------|--------|------------|--------------|--------|
| `Fig01_geometry.png` | Fig1 Geometry CoE/shell/vision | rerun:universe_thin.py | yes | no (RGBA/PNG enc) | **PASS** |
| `Fig02_Achart.png` | Fig2 Migkas A/A_all + model curves | rerun:universe_thin.py | yes | no (RGBA/PNG enc) | **PASS** |
| `Fig03_migkas_slices.png` | Fig3 Migkas redshift slices + D–m_cc degeneracy | copy:exact_on_disk_asset | yes | no (RGBA/PNG enc) | **PASS** |
| `Fig04_amp_full.png` | Fig4 Full amp + residual | copy:exact_on_disk_asset | yes | no (RGBA/PNG enc) | **PASS** |
| `Fig05_sky_orientation_four_row.jpg` | Fig5 Four-row sky orientation | copy:bitwise_match_asset | yes | yes | **PASS** |
| `Fig06_overlay_stack.png` | Fig6 Migkas Fig.8 overlay stack | copy:bitwise_match_asset | yes | yes | **PASS** |
| `Fig07_storyboard.png` | Fig7 Donut-universe storyboard A–D | copy:exact_on_disk_asset | yes | no (RGBA/PNG enc) | **PASS** |
| `Fig08_two_emitters.png` | Fig8 Two-emitter baseline | copy:exact_on_disk_asset | yes | no (RGBA/PNG enc) | **PASS** |
| `Fig09_kill.jpg` | Fig9 Kill function vs plate/shell/shear | copy:bitwise_match_asset | yes | yes | **PASS** |
| `Fig10_cancel_partner.png` | Fig10 Cancel partner residuals | copy:exact_on_disk_asset | yes | no (RGBA/PNG enc) | **PASS** |
| `FigA1_ltb_shell.png` | FigA1 LTB shell evolution | copy:exact_on_disk_asset | yes | no (RGBA/PNG enc) | **PASS** |
| `FigA2_coe_sun.jpg` | FigA2 CoE/Sun timing 2026 | copy:bitwise_match_asset | yes | yes | **PASS** |
| `FigA2b_sun_coe_crossings.jpg` | TableA1 Sun/CoE crossings 2027 | copy:bitwise_match_asset | yes | yes | **PASS** |
| `FigA3_coe_sky.png` | FigA3 Annual CoE sky path from Plano | rerun:_make_coe_sky_path_dates.py | yes | no (RGBA/PNG enc) | **PASS** |
| `FigA_comic_HANDDRAWN_from_PDF.jpg` | p.14 hand-drawn comic | extract/copy:handdrawn_not_from_code | yes | yes | **PASS** |
| `eclipse_inner_shell.png` | Inner-shell eclipse geometry (supplemental; not a letter-page embed) | rerun:attachment 7ef6f478…py | yes | n/a | **PASS** |

## Per-figure notes

### Fig01_geometry.png
- **Status:** PASS
- **Gold:** gold-000.png md5=`ef40ce3d854574591d09d1edb96a38cc`
- **Delivered md5:** `d9fd6acd60b35a72f61515b87e18ee95`
- **Dims/size:** [1460, 488] / 138940 bytes
- **Method:** rerun:universe_thin.py
- **Notes:** Pixel-identical to PDF embed; bitwise-identical to spreadsheet_model/universe_thin_drawing.png

### Fig02_Achart.png
- **Status:** PASS
- **Gold:** gold-002.png md5=`5178d14e2707186ad6926d3a92ef20f7`
- **Delivered md5:** `6dac4950b3b907e5b421cd1a256df2f2`
- **Dims/size:** [1182, 1287] / 206551 bytes
- **Method:** rerun:universe_thin.py
- **Notes:** Pixel-identical to PDF embed; bitwise-identical to spreadsheet_model/universe_thin_Achart.png

### Fig03_migkas_slices.png
- **Status:** PASS
- **Gold:** gold-004.png md5=`aee1382160adba01ab10ec5cc9c37323`
- **Delivered md5:** `6e97e4c31d507d3921aa8b702adefae6`
- **Dims/size:** [1123, 1014] / 245357 bytes
- **Method:** copy:exact_on_disk_asset
- **Notes:** Original one-off plot script not found; migkas_slices_degeneracy_clean.png pixel-identical to PDF embed. Prior figs_regenerated_from_code Fig3 was WRONG.

### Fig04_amp_full.png
- **Status:** PASS
- **Gold:** gold-006.png md5=`771d0afc75dd3aa8a128de519df58216`
- **Delivered md5:** `8ecabb7ddbd69a9c605cd45c459b9c6c`
- **Dims/size:** [1800, 1385] / 554781 bytes
- **Method:** copy:exact_on_disk_asset
- **Notes:** money_shot_amp_full.png pixel-identical; NOT SINE_BAD. One-off script not found.

### Fig05_sky_orientation_four_row.jpg
- **Status:** PASS
- **Gold:** gold-007.jpg md5=`ed41b905e2563efb4df2c92b1e13e3b4`
- **Delivered md5:** `ed41b905e2563efb4df2c92b1e13e3b4`
- **Dims/size:** [1370, 1420] / 404512 bytes
- **Method:** copy:bitwise_match_asset
- **Notes:** money_shot_four_row_only.jpg bitwise-identical to PDF JPEG embed.

### Fig06_overlay_stack.png
- **Status:** PASS
- **Gold:** gold-008.png md5=`fd9f1c1c95194a3fb7864de08698704f`
- **Delivered md5:** `fd9f1c1c95194a3fb7864de08698704f`
- **Dims/size:** [2164, 3222] / 2664167 bytes
- **Method:** copy:bitwise_match_asset
- **Notes:** att14_overlay_stack.png bitwise-identical to PDF embed.

### Fig07_storyboard.png
- **Status:** PASS
- **Gold:** gold-009.png md5=`05a980ff995d5a43f75de00c14fc750e`
- **Delivered md5:** `f2b0d04b62f80b9cb32c41a1e0d397ac`
- **Dims/size:** [1045, 1795] / 211403 bytes
- **Method:** copy:exact_on_disk_asset
- **Notes:** att14_storyboard.png pixel-identical. Current donut script output differs — not used.

### Fig08_two_emitters.png
- **Status:** PASS
- **Gold:** gold-011.png md5=`6953c5ae7902824f2fa24c854884cef2`
- **Delivered md5:** `41ddd33dbefce0695c7e225c72f26b58`
- **Dims/size:** [1181, 1908] / 129240 bytes
- **Method:** copy:exact_on_disk_asset
- **Notes:** att14_two_emitters.png pixel-identical. Current two_emitters script output differs — not used.

### Fig09_kill.jpg
- **Status:** PASS
- **Gold:** gold-013.jpg md5=`15818532df3ff7b9bde79a200e779402`
- **Delivered md5:** `15818532df3ff7b9bde79a200e779402`
- **Dims/size:** [1206, 1280] / 199281 bytes
- **Method:** copy:bitwise_match_asset
- **Notes:** kill_vs_cycle6_shapes.jpg bitwise-identical to PDF JPEG embed.

### Fig10_cancel_partner.png
- **Status:** PASS
- **Gold:** gold-014.png md5=`7f6701fcf63d81511470aea0c3626dc0`
- **Delivered md5:** `e8de308591d54d8c6de1b5481729b658`
- **Dims/size:** [2084, 2212] / 652876 bytes
- **Method:** copy:exact_on_disk_asset
- **Notes:** cancel_partner_residual.png pixel-identical (= simple_cmb_two_quad_cancel.png).

### FigA1_ltb_shell.png
- **Status:** PASS
- **Gold:** gold-016.png md5=`717250da18b391623eb4618f233969fe`
- **Delivered md5:** `09462730cede0bb8dbb0052aef3eda1a`
- **Dims/size:** [1444, 969] / 161023 bytes
- **Method:** copy:exact_on_disk_asset
- **Notes:** ltb_a_vs_t_all.png pixel-identical.

### FigA2_coe_sun.jpg
- **Status:** PASS
- **Gold:** gold-018.jpg md5=`9cb249b71dbb3ea478495deb47837921`
- **Delivered md5:** `9cb249b71dbb3ea478495deb47837921`
- **Dims/size:** [1280, 636] / 89336 bytes
- **Method:** copy:bitwise_match_asset
- **Notes:** appendix_coe_sun_timing_2026.jpg bitwise-identical to PDF JPEG embed.

### FigA2b_sun_coe_crossings.jpg
- **Status:** PASS
- **Gold:** gold-019.jpg md5=`6cf8f10855878c559c54e97ea06cb029`
- **Delivered md5:** `6cf8f10855878c559c54e97ea06cb029`
- **Dims/size:** [1239, 570] / 62616 bytes
- **Method:** copy:bitwise_match_asset
- **Notes:** appendix_sun_coe_crossings_2027.jpg bitwise-identical to PDF JPEG embed.

### FigA3_coe_sky.png
- **Status:** PASS
- **Gold:** gold-020.png md5=`1c0042622ee084daa22f2dfaf7f4a086`
- **Delivered md5:** `6d68e94d82d64d87d409f67620259ce4`
- **Dims/size:** [2415, 1432] / 222910 bytes
- **Method:** rerun:_make_coe_sky_path_dates.py
- **Notes:** Re-ran _make_coe_sky_path_dates.py; pixel-identical to gold.

### FigA_comic_HANDDRAWN_from_PDF.jpg
- **Status:** PASS
- **Gold:** gold-022.jpg md5=`262874d6412bd4b1d38b7905e237e72b`
- **Delivered md5:** `262874d6412bd4b1d38b7905e237e72b`
- **Dims/size:** [2080, 944] / 237379 bytes
- **Method:** extract/copy:handdrawn_not_from_code
- **Notes:** HAND-DRAWN — cannot invent from code. Bitwise-identical to PDF JPEG embed.

### eclipse_inner_shell.png
- **Status:** PASS
- **Gold:** (not a letter embed)
- **Delivered md5:** `f76d8526ed668a5df5855f9cb66314cc`
- **Dims/size:** [883, 945] / 110236 bytes
- **Method:** rerun:attachment 7ef6f478…py
- **Notes:** Ran provided attachment eclipse script; bitwise-identical to letter_brief/figs/inner_shell_eclipse_sketch.png.

## Gold embed map (page → extract)

| Page | Extract | Dims | Letter asset |
|------|---------|------|--------------|
| 1 | gold-000.png | 1460×488 | Fig01 / universe_thin_drawing |
| 2 | gold-002.png | 1182×1287 | Fig02 / universe_thin_Achart |
| 3 | gold-004.png | 1123×1014 | Fig03 / migkas_slices_degeneracy_clean |
| 4 | gold-006.png | 1800×1385 | Fig04 / money_shot_amp_full |
| 5 | gold-007.jpg | 1370×1420 | Fig05 / money_shot_four_row_only |
| 6 | gold-008.png | 2164×3222 | Fig06 / att14_overlay_stack |
| 7 | gold-009.png | 1045×1795 | Fig07 / att14_storyboard |
| 8 | gold-011.png | 1181×1908 | Fig08 / att14_two_emitters |
| 9 | gold-013.jpg | 1206×1280 | Fig09 / kill_vs_cycle6_shapes |
| 10 | gold-014.png | 2084×2212 | Fig10 / cancel_partner_residual |
| 12 | gold-016.png | 1444×969 | FigA1 / ltb_a_vs_t_all |
| 13 | gold-018.jpg | 1280×636 | FigA2 / appendix_coe_sun_timing_2026 |
| 13 | gold-019.jpg | 1239×570 | crossings / appendix_sun_coe_crossings_2027 |
| 14 | gold-020.png | 2415×1432 | FigA3 / coe_sky_path_dates_texas |
| 14 | gold-022.jpg | 2080×944 | comic HANDDRAWN / appendix_coe_comic |

## FAIL list

None — all letter embeds pixel-matched; eclipse matched reference PNG.

## Zip
`/workspace/figs_exact_from_code/figs_exact_from_code.zip`
