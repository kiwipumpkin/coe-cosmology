# Sequential placement rerun with shallow, observation-based voids (Clampitt & Jain 2015, arXiv:1404.1834)

## Void model
- **Core (r < Rv):** δ(r) = −0.5 + 0.25 (r/Rv)². The centre is −0.5 and the mean inside Rv is −0.35 (the old model was −0.70).
- **Compensation:** no compact rim clusters. The missing mass sits in a smooth shell from Rv to 3Rv. The enclosed deficit falls smoothly, is still half-uncompensated at 2Rv, and is fully compensated at 3Rv. Net mass is zero per void.
- **Tidal pull:** each void is spherically symmetric, so its pull at every ring point and at us comes exactly from the enclosed mass M_enc(r). This is the same tidal formula as before, with M replaced by M_enc(r).
  - A void more than 3Rv from both the CMB ring and us has **exactly zero** effect (shell theorem).
- **Overlap rule:** void cores (radius Rv) never overlap (gap ≥ 0.02 Gpc, checked on every layout). Compensation shells may overlap neighbours and other cores; their masses simply add (linear superposition).
- **Shells near us:** a shell may also reach inside the 14.3 Gpc sphere, which is where its pull on the ring comes from.
- **Same as before:** region, scoring, floor and void sizes (1.5–3.75 Gpc near, 1.5 Gpc elsewhere). There are no clusters now; voids only.

## 1. Near-region fix (14.3–30 Gpc)
|case|leftover|
|-|-|
|start|237.3×|
|four clumps|88.8×|
|old deep fit (24 voids + 24 compact clusters)|27.4×|
|old deep layout, voids swapped to shallow, clusters removed (no refit)|237.3× (no change at all)|
|old deep layout, voids swapped to shallow, old clusters kept|245.5×|
|**shallow voids refitted** (4 seeds, 28–36 voids each)|**135.3×** (all seeds 135–137×)|

- **Why the swap does nothing:** every old void sits 3.2–10 Rv from the ring, so a self-compensated version is invisible.
- **Where the old fix came from:** the compact clusters and the unbalanced deep voids.
- **Shallow refit:** voids alone get about 40% of the way down (237× → 135×). They stall well above the four clumps' 89×.

## 2. Outer fill (30–42.9 Gpc) around the frozen shallow near layout (135.3×)
|fill|layouts|voids added|median leftover|change|helped / hurt|
|-|-|-|-|-|-|
|foam-like|5|~4,540|135.3×|0.0×|0 / 0|
|random|5|~1,700|135.3×|0.0×|0 / 0|

- With shallow voids the outer shell has **no effect at all**: every outer void is far more than 3Rv (4.5 Gpc) from the ring and from us.
- For comparison, the old deep model's outer fills added +1.3× (foam) and +1.8× (random).
- The random fill stopped earlier (it ends after 3,000 total failed tries), so it has fewer voids. That doesn't matter because each void's effect is exactly zero.

## 3. Growth from empty, one void at a time (Rv = 1.5 Gpc, 2 seeds, capped at 6,000 voids ≈ 30% void volume)
|rule|lowest while growing (deep → shallow)|when full (deep → shallow)|
|-|-|-|
|foam (grown next to neighbours)|223× → 237×|610× → 334×|
|random|215× → 229×|1,370× → 422×|
|CMB-aware|37× → 166×|186× → 166×|

- Only about 45–55 of the 6,000 voids sit close enough to the ring (within 3Rv) to matter.
- Shallow voids hurt much less when the shell fills, because there are no compact clusters to throw off the pattern. They also can't help much.
- Even the CMB-aware growth only reaches 166×, and it never gets near 89×.

## Charts
- `near_fix_deep_vs_shallow.png`
- `outer_deep_vs_shallow.png`
- `growth_deep_vs_shallow.png`

## Cuts
- The outer fills ran 5 layouts per fill type, not 30.
- Growth ran 2 seeds per rule and was capped at 6,000 voids; the jamming tail was too slow.
- The old deep values are from ../results.json.
