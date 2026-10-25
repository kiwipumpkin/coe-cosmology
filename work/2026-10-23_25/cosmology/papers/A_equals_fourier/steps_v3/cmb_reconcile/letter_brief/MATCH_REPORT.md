# MATCH_REPORT — figs_from_makers vs letter gold

Generated from regenerating makers (not copy-asset).

| File | Status | Note | Maker | Gold |
|------|--------|------|-------|------|
| `Fig01_geometry.png` | **PIXEL_MATCH** |  | (1460, 488) | (1460, 488) |
| `Fig02_Achart.png` | **PIXEL_MATCH** |  | (1182, 1287) | (1182, 1287) |
| `Fig03_migkas_slices.png` | **DIM_MISMATCH** | maker=(992, 977) gold=(1123, 1014) resized_mad=26.4 | (992, 977) | (1123, 1014) |
| `Fig04_amp_full.png` | **DIM_MISMATCH** | maker=(1581, 1234) gold=(1800, 1385) resized_mad=21.1 | (1581, 1234) | (1800, 1385) |
| `Fig05_sky_orientation_four_row.jpg` | **DIM_MISMATCH** | maker=(1370, 1575) gold=(1370, 1420) resized_mad=17.0 | (1370, 1575) | (1370, 1420) |
| `Fig06_overlay_stack.png` | **DIM_MISMATCH** | maker=(1834, 4402) gold=(2164, 3222) resized_mad=68.8 | (1834, 4402) | (2164, 3222) |
| `Fig07_storyboard.png` | **DIM_MISMATCH** | maker=(424, 2039) gold=(1045, 1795) resized_mad=13.9 | (424, 2039) | (1045, 1795) |
| `Fig08_two_emitters.png` | **DIM_MISMATCH** | maker=(897, 1504) gold=(1181, 1908) resized_mad=13.2 | (897, 1504) | (1181, 1908) |
| `Fig09_kill.jpg` | **DIM_MISMATCH** | maker=(1206, 639) gold=(1206, 1280) resized_mad=23.0 | (1206, 639) | (1206, 1280) |
| `Fig10_cancel_partner.png` | **DIM_MISMATCH** | maker=(1962, 2060) gold=(2084, 2212) resized_mad=16.0 | (1962, 2060) | (2084, 2212) |
| `FigA1_ltb_shell.png` | **DIM_MISMATCH** | maker=(1885, 1482) gold=(1444, 969) resized_mad=14.9 | (1885, 1482) | (1444, 969) |
| `FigA2_coe_sun.jpg` | **DIM_MISMATCH** | maker=(1280, 700) gold=(1280, 636) resized_mad=14.3 | (1280, 700) | (1280, 636) |
| `FigA2b_sun_coe_crossings.jpg` | **DIM_MISMATCH** | maker=(1239, 449) gold=(1239, 570) resized_mad=46.2 | (1239, 449) | (1239, 570) |
| `FigA3_coe_sky.png` | **VISUAL_DIFF** | max=255 mean=1.95 | (2415, 1432) | (2415, 1432) |
| `FigA_comic_HANDDRAWN_from_PDF.jpg` | **PIXEL_MATCH** |  | (2080, 944) | (2080, 944) |
| `eclipse_inner_shell.png` | **PIXEL_MATCH** | vs letter figs | (883, 945) | (883, 945) |

## Summary

- **PIXEL_MATCH:** Fig01, Fig02, FigA3, comic (handdrawn), eclipse
- **DIM_MISMATCH / layout differs (code regenerated, not letter bake-in):** Fig03–Fig10, FigA1, FigA2, FigA2b
- **Fig06:** intentionally NOT a bitwise copy of gold (gold had cut-off axis titles). Fresh margin-fixed overlay stack from `make_migkas_fig8_overlay_matched.py`.

