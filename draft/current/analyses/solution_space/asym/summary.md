# Rerun with no mirror-symmetric layouts (masses + deep voids)

Chris: "Rerun with the rule that it can't be exactly symmetric."
Same model, target and floor as `../` (ss_core.py unchanged). Leftover is given in units of the 1e-5 CMB-flat floor, with a target of 1×.

## The rule (applied in every fit, then checked again on the result)
"Mirror" means reflection across the plane that contains us, the CoE and the galactic pole (ℓ → 2×15.16° − ℓ). The ring can't tell b from −b, so a mirror image with its latitude flipped also counts.
- **No mirror pairs:** two objects of the same kind are forbidden if all three hold at once:
  - one sits within **20°** of the other's mirror image,
  - their distances from us differ by less than **15%**,
  - their sizes differ by less than **25%** (mass for masses, radius for voids).
- **Nothing on the mirror plane:** no object may be within **10°** of the plane, because such an object is its own mirror image.
- **Why these tolerances:** they rule out every earlier good fit. The old 0.58× floor layout has 5 exact mirror pairs plus 1 object on the axis. The old 4+4 best (4.46×) still had pairs only 16° off-mirror.
- **How it's applied:** a smooth penalty during the fit, then a hard check with the rule above. Only fits that pass are counted, and all of them did.
- **Starts:** 16 per mix (24 for 6+5 and 6+6). About 60% were warm starts from the old fits, perturbed independently object by object so the starts are not symmetric. The rest were cold random starts.

## Result (`a_best_vs_N_asym_logy.png`)
|N|before (mirror allowed)|now (no mirror)|best mix now|
|-|-|-|-|
|8|4.46|6.19|4 masses + 4 voids|
|9|2.96|3.90|6 + 3|
|10|1.26|2.30|6 + 4|
|11|**0.58**|1.25|6 + 5|
|12|0.72|**0.73**|7 + 5 (6 + 6: 0.78)|

- **Cost of the rule:** fits get about 1.3–2× worse at a given N.
- **Still reaches 1×:** it now takes 12 objects instead of 11. 4 of the 24 starts at 6 + 6 and 1 of 16 at 7 + 5 go below 1×.

## Best layout without mirror images (`d_asym_solution_map.png`): 7 masses + 5 voids, 0.73×
- **Mass bands:** masses sit at ℓ ≈ 27°, 150°, 164°, 177°, 215°, 230° and 305°.
  - The 162° and 228° bands and the 305–315° band are still there.
  - The 75° band moved to about 27° in this layout. The 6 + 5 best (1.25×) still has it, at ℓ 82°.
- **Voids:** 4 sit on the far side at 24–32 Gpc (ℓ 26°, 51°, 120°, 357°), plus one near the Migkas side (ℓ 267°).
- **Mass balance:** 5.82e52 kg in masses against −5.81e52 kg in voids, which agree to 0.17%. Nothing forced this.
- **Objects at limits:** 1 void at the 4 Gpc cap, no mass at the 2e52 cap, and nothing touching our sphere or the CoE interior. Several objects sit exactly 10.5° from the mirror plane, which is the edge of the allowed zone.

## What it tells us
- **The overall pattern is still symmetric.** The residual itself is symmetric about the ℓ 15°/195° axis, so the fit still builds a pattern that is mirror-symmetric to 0.1–0.2%. It now does this with objects that are not mirror copies of each other.
- **The fit uses the gaps in the rule.** Several pairs sit at mirror positions (2–11° off) but differ in mass or radius by just over 25%. Others are a single object facing two smaller ones.
- **A stricter rule costs much more.** As a quick check, any same-kind object within 20° of another's mirror position was forbidden, whatever its mass or distance:
  - 6 + 5 best is **4.47×** (16 starts).
  - 6 + 6 best is 4.64×, but only 3 starts finished (cut for time).
  - So without mirror-placed pairs, about 4.5× is the best found at N = 11–12.

Files:
- Code: `asym_core.py`, `run_asym.py`, `run_all.py`, `run_pos.py`, `report_asym.py`, `check_old.py`
- Results: `raw/`, `raw_pos/`, `asym_results.json`
