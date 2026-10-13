# All-sky acceleration-sign test, CoE shell Paper I (Newtonian, all at one instant)

Script: `accel_sign_test.py` · numbers: `results.json` · maps: `maps/` (Mollweide, CoE ★ at map centre, red +, blue −, white 0, black line = sign flip)

## Inputs taken from Paper I (`/workspace/cosmo_preprint/coe_shell_preprint.tex`)
- Geometry: D = 58.2 Gpc, shell 43.9–72.5 Gpc, visible radius 14.3 Gpc. Every sightline with s ≤ 14.3 Gpc stays inside the shell (r_g runs from 43.9 to 72.5), so no outside-shell case comes up.
- Velocity law: the paper's only explicit law, from the SR section (eq. v2, run 2): v² = β1² r0² + 2μ(1/r − 1/r0), with β1 = 0.5866, μ = 0.09332, s = 3, and d1 = D/3 = 19.4 Gpc. That gives v_obs = 0.469c, and v runs from 0.176c at 43.9 Gpc to 0.659c at 72.5 Gpc. Expanding uses this outward. Collapsing uses the same speed inward.
- The paper has no single M(r), so I ran three mass profiles:
  - **P:** central mass only, GM = μc²d1, so M = 7.5×10⁵² kg. This is the mass that sets the paper's velocity law.
  - **PS:** P plus the exp-lock shell, ρ = 2.9×10⁻²⁷ exp(−(r−D)/922.65 Gpc). The shell holds 1.06×10⁵⁶ kg.
  - **U:** a uniform shell at the same density with an empty interior, for contrast.
- Check: the vector projection matches the reviewer's formula a_los = −f(r_g)(s + D cosθ) + f(D) D cosθ to within 3×10⁻¹⁵.

## Key structural result
On a single snapshot, both a_rel and the transverse term |v_rel|² − v_los² are unchanged when every velocity is reversed. So **d²s/dt² is identical for the expanding and collapsing cases**. Reversing the motion only flips the sign of v_los. Whether the recession or approach speed is growing therefore flips too: if a fraction f of the sky shows growing speed while expanding, then 1 − f does while collapsing.

## Results (fractions of the sky; angles measured from the CoE direction)
| Profile | Case | v_los | d²s/dt² > 0 | Speed growing | Sign-flip rings |
|---|---|---|---|---|---|
| P | expanding | recedes over the whole sky (2.4–4.5k km/s at s=1; 34–85k km/s at s=14) | 0.63 | **0.63** (two caps: within about 66–68° of the CoE and about 69–73° of the anti-CoE) | 68° and 111° at s=1; 63° and 107° at s=14 |
| P | collapsing | approaches over the whole sky | 0.63 | 0.37 (the belt between the rings) | same rings |
| PS / U | expanding | recedes over the whole sky | 0.00 | 0.00 (everywhere decelerating) | none |
| PS / U | collapsing | approaches over the whole sky | 0.00 | 1.00 (approach speeds up everywhere) | none |

The results barely depend on s; the rings move by only about 4–5° between s = 1 and 14 Gpc. Magnitudes:
- **P:** |s̈| ≲ 6×10⁻¹⁴ m/s² at s = 1 and ≲ 1.8×10⁻¹² m/s² at s = 14.
- **PS and U:** s̈ ≈ −(1.4–4.7)×10⁻¹¹ m/s² at s = 1 and −(2–8)×10⁻¹⁰ m/s² at s = 14. The shell's own mass outweighs the central mass about 1400×.
- The transverse term reaches at most 1.2×10⁻¹² m/s². It matters for P but is negligible for PS and U.

## Direct answers
1. **Does expansion give apparent acceleration in every direction?** No, for all three profiles. With the paper's own central mass, about 63% of the sky shows it, in two polar caps around the CoE and anti-CoE directions, with a decelerating belt about 68–111° from the CoE. Once the shell's own matter at 2.9×10⁻²⁷ kg/m³ is included, recession decelerates in every direction.
2. **Does collapse keep apparent acceleration only in a restricted region?** Under this velocity law, nothing recedes during collapse: v(r) rises outward, so every sightline approaches. "Acceleration" can then only mean the approach speeding up. That happens in the complementary region: the 37% belt for P, or the whole sky for PS and U. The claimed pattern, receding with apparent acceleration in a limited cone, does not appear with this velocity law and mass. Getting receding objects during collapse would need a velocity law where inward speed falls with r.

## Caveats
- Newtonian dynamics, with speeds up to 0.66c.
- All positions and velocities are taken at one instant. There is no light-cone or lookback correction, and no photon redshift. The observable is z(t_obs) along the past light cone, which this does not compute.
- The paper's velocity law comes from energy conservation in the central mass only. PS and U therefore combine a velocity field and a gravity field that are not mutually consistent. A consistent LTB model with E(r) from the full M(r) is the next step.
- Rotation, the +15° tilt, the four clumps and dark matter are all left out.

## Files
- `maps/sdd_{P,PS,U}_{expanding,collapsing}_s{1,3,7,14}.png`: d²s/dt² (24 maps).
- `maps/speedup_*`: whether the recession or approach speed is growing (sign only, 24 maps).
- `maps/alos_{P,PS,U}_s*.png`: gravitational a_los, which is the same in both cases (12 maps).
- `maps/vlos_P_*`: v_los (8 maps).

---
# Redo (self-consistent, Chris's spreadsheet definition) — `selfconsistent.py`, `selfconsistent.json`
This redo uses the definition from `cosmology/papers/spreadsheet_model/columns.tex`:
- V = k·r about the CoE (Sedov–Taylor, homologous).
- ρ ∝ e^(−r/L), with the mass inside D held fixed.
- Relative acceleration = a(P) − a(observer).
- The "recession v with gravity" is |v_rel + a_rel dt|.
- No "megus" definition was found anywhere in /workspace.

Because V ∝ r, v_rel = k·s·n points straight away from us. That makes the transverse term exactly zero. Two consequences:
- Change of recession speed = sign(k) × a_los.
- d²s/dt² = a_los, the same whether the shells expand or collapse.

## Fraction of sky with relative acceleration away from us (a_los > 0)
| ρ profile | Fraction (s = 1 / 7 / 14 Gpc) | Inward-pulling band, ° from CoE (s = 7) |
|---|---|---|
| Point mass, or L ≤ 6.07 Gpc (paper's free-L fit) | 42% | 52° to 123° |
| L = 10 Gpc | 39% | 49° to 125° |
| L = 20 Gpc | 25% | 38° to 135° |
| L ≥ 50 Gpc, incl. the CatWISE lock L = 923 (≈ uniform) | 0% | whole sky |

## Plain answers
1. **Claim (1), acceleration away from us across the whole sky while expanding, is impossible for any positive spherical mass.** At 90° from the CoE axis the observer feels no sideways pull, but the galaxy is pulled toward the CoE. That gives a_los = −f(r_g)·s < 0 exactly. A non-homologous v(r) adds a transverse term, but for nearby objects that term grows like s³ while the pull grows like s, so it cannot rescue small s.
2. **Claim (2), most of the sky accelerating away with a band near 90° pulled toward us, does appear, but only for centrally concentrated mass (L ≲ 10 Gpc or a point mass).** About 40% of the sky (two caps) accelerates away; the band from about 50° to 123° accelerates toward us. With the paper's locked, nearly uniform shell, everything is pulled toward us.
3. **The map is the same for expanding and collapsing at a given moment**, because gravity does not depend on velocity. Expansion versus collapse changes only whether that acceleration speeds up or slows down the recession. Expanding: recession speeds up in the caps (~42%). Collapsing: approach speeds up in the band (~58%) and slows in the caps.
4. **A way to get a real expanding/collapsing difference:** model the light cone, because the objects are seen at earlier times when the shell speeds were different. That is not done here.

Maps (L = 6.07, s = 7; identical by construction): `sc_expanding_L6.07_s7.png`, `sc_collapsing_L6.07_s7.png`.

## Point mass at the CoE, with velocities consistent with its gravity
- **Velocities:** v² = β1²r0² + 2μ(1/r − 1/r0), the paper's run-2 energy law. It is exact energy conservation for a point mass of 7.5×10⁵² kg, so velocities and gravity are consistent here.
- **Transverse term:** the flow is not ∝ r, so this case includes the velocity-dependent term.
- **Acceleration away from us (d²s/dt² > 0):** 63% of the sky, in two caps. One cap is within about 66° of the CoE; the other is the matching region around the anti-CoE.
- **Band pulled toward us:** about 66° to 109° from the CoE at s = 7 Gpc.
- **Band at other distances:** 68°–111° at s = 1 and 63°–107° at s = 14.
- **Expanding vs collapsing:** the same d²s/dt² map in both. Expanding, recession speeds up in the caps (63%). Collapsing, everything approaches and the approach slows in the caps (63%) but speeds up in the band (37%).
- **Maps:** `pm_expanding_s7.png` and `pm_collapsing_s7.png`.
