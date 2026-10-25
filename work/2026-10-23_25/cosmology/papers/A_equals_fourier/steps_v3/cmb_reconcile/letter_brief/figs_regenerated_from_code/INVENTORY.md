# CoE letter figures — regenerated from code

Authoritative PDF: attachment `02b6cafe…pdf` (14 pages, LaTeX, Wed Sep 30 2026 CDT) ≡ `letter_layout_source.pdf`.

Output folder: `letter_brief/figs_regenerated_from_code/`

## Inventory (PDF figure → script → output → status)

| PDF fig | Caption (short) | Script used | Output | Status |
|---------|-----------------|-------------|--------|--------|
| Fig 1 | Geometry CoE / shell / vision | `spreadsheet_model/bubble/universe_thin.py` | `Fig01_geometry.png` | **OK — pixel match to PDF embed** |
| Fig 2 | Migkas A/A_all + model curves | same `universe_thin.py` | `Fig02_Achart.png` | **OK — pixel match** |
| Fig 3 | Redshift slices + D–m_cc degeneracy | `coe-figures/src/fig03_migkas_slices.py` (recompute; original one-off not found) | `Fig03_migkas_slices.png` | OK from code; layout differs from PDF (~13% px) |
| Fig 4 | Full amp + residual (no Amp×sin squiggles) | `coe-figures/src/fig04_amp_full.py` (recompute; `money_shot_amp_full` one-off not found) | `Fig04_amp_full.png` | OK from code; **not** the forbidden SINE_BAD figure; visual differs from PDF money_shot |
| Fig 5 | Four-row sky orientation | `make_coe_rings_notilt.py` + `make_coe_rings_migkas_pin_axis_matched.py` + `_matched_grad.py` + `_matched_grad_rot_fit.py`, composed with `make_stack_clean`-style crops | `Fig05_sky_orientation_four_row.jpg` | OK from code; ~same size as PDF (1370×1441 vs 1370×1420) |
| Fig 6 | Migkas Fig.8 overlays (tilt / grad / roll) | **`figures/make_migkas_fig8_overlay_matched.py`** with `SHOW_LEGEND=True` (restored after run); stacked | `Fig06_overlay_matched.png`, `_grad.png`, `_grad_rot81.png`, `Fig06_overlay_stack.png` | **OK — full panels 1975×1515 with legend**; stack taller than PDF crop |
| Fig 7 | Donut storyboard A–D/E | `cycle7_recombination/make_donut_universe_storyboard.py` | `Fig07_storyboard.png` | OK from code; current script is narrower than PDF-era storyboard |
| Fig 8 | Two-emitter baseline | `cycle6_two_emitters/make_two_emitters_comoving.py` | `Fig08_two_emitters.png` | OK from code |
| Fig 9 | Kill vs plate/shell/shear | `coe-figures/src/fig09_kill.py` (from `cycle4/dtt/dtt.json`; letter jpg one-off not found) | `Fig09_kill.png` | OK from code; layout differs from `kill_vs_cycle6_shapes.jpg` |
| Fig 10 | Cancel partner residuals | `coe-figures/src/fig10_cancel_partner.py` (from `simple_cmb_two_quad_cancel.json`) | `Fig10_cancel_partner.png` | OK from code; PDF embed matches older `simple_cmb_two_quad_cancel.png` |
| Fig A1 | LTB shell evolution | `coe-figures/src/figA1_ltb_shell.py` | `FigA1_ltb_shell.png` | OK from code |
| Fig A2 | CoE / Sun timing 2026 | `coe-figures/src/figA2_coe_sun.py` | `FigA2_coe_sun.png` | OK from code (PDF jpg had no surviving plot script) |
| Table A1 | Sun/CoE crossings 2027 | `coe-figures/src/figA_table_crossings.py` | `FigA_table_crossings.png` | OK from code |
| Fig A3 | Annual CoE sky path | `letter_brief/figs/_make_coe_sky_path_dates.py` | `FigA3_coe_sky.png` | **OK — pixel match** |
| p.14 comic | “Mommy look…” hand-drawn | *none* (not a plot) | `FigA_comic.png` = schematic redraw; `FigA_comic_HANDDRAWN_not_from_code.jpg` = PDF original | **FAIL exact** — hand-drawn; schematic only from code |
| *(extra)* | Inner-shell eclipse geometry | `cycle6_two_emitters/make_inner_shell_eclipse_sketch.py` (**recreated**; original inline session not on disk) | `eclipse_inner_shell.png` | OK from code; matches reference geometry (Rin=1,Rout=3,D=1.8,α≈33.7°); not pixel-identical to attachment |

## Hard-rule compliance

1. Regenerated with Python plotting code — no whiteouts/screenshot pastes of damaged assets.
2. Delivered figure files (not a rebuilt letter PDF).
3. Fig 6 used `make_migkas_fig8_overlay_matched.py`; full undamaged panels with legend at 1975×1515.
4. Did **not** regenerate the Amp×sin figure with squiggly connectors (`money_shot_amp_full_SINE_BAD.png` avoided).
5. Matched letter `\includegraphics` provenance where original scripts exist (`universe_thin.py`, overlay matched, donut, two-emitters, CoE sky).

## Zip

See sibling `figs_regenerated_from_code.zip`.
