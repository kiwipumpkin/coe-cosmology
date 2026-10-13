# Radial blast-wave shells, only inside our visible sphere

**Setup**
- Concentric over- and under-dense shells centred on the CoE (58.2 Gpc away toward l = 195°), present only inside our 14.3 Gpc sphere. Outside it there is no wave.
- Contrast: δ = A√2·cos(2π r_CoE / λ), with A = 0.10 rms, the same as the checkerboard runs.
- Wavelengths: λ = 120 h⁻¹ Mpc = **171 Mpc** (Einasto 1997, with h = 0.7), plus a BAO-like **150 Mpc** for comparison.
- Same ring kernel, sampling and body-fitting method as `../checkerboard/`.

**Geometry (Fig 1)**
- Along the CoE axis our sphere spans 43.9 to 72.5 Gpc from the CoE. That is 28.6 Gpc / 171 Mpc ≈ **167 crests across our sphere** (190 at 150 Mpc).
- Perpendicular to the axis the shells curve away, so there are only about 10 crests from our centre to the edge.
- Inside our sphere the shells therefore look like nearly parallel, gently curved bands.
- The ring runs *along* the shells near l = 195° and 15°, and cuts *across* them near l = 105° and 285°.

**Results (A = 0.10, phase 0)**

| | 171 Mpc (Einasto) | 150 Mpc (BAO) |
|---|---|---|
| Waves alone: total | 7.5× | 5.8× |
| … m = 1–5 / 6–12 | 1.6× / 1.9× | 0.6× / 1.1× |
| … m ≥ 13 | **7.0×** | **5.7×** |
| … symmetric share, 1° sampling (all m / m 1–12) | 60% / 86% | 57% / 97% |
| Other shell phases (7): m 1–12 | 2.4–17.8× | 1.2–13.9× |
| Other shell phases (7): m ≥ 13 | 7.05–7.07× | 5.6–5.7× |
| With the 237× residual, before bodies | 237× | 238× |
| **After 270 bodies:** m 1–12 | **0.004×** | **0.001×** |
| **After 270 bodies:** m ≥ 13 left | **7.0×** | **5.7×** |
| After bodies, symmetric share | 57% | 54% |

- **The bodies.** 144 masses and 126 voids (148 and 122 in the BAO case), each at most 2×10⁵² kg and all at least 3.7 Gpc from the ring (4.6 Gpc in the BAO case), against a required 1.5 Gpc. 150 of them are the same set used for the checkerboard to cancel the 237× residual; 120 cancel the waves' own low modes. Together they add about 1× of fine texture.
- **Main point.** As with the checkerboard, the bodies remove m = 1–12 completely, but the shells leave a noise-like floor at m ≥ 13. That floor is 7× at 171 Mpc and 5.7× at 150 Mpc, and no distant body touches it (fitting the waves alone only takes 7.05 → 6.94×). It is the same for every shell phase.
- **Low modes depend on phase.** m 1–12 varies strongly with where the crests fall: it is set by how the ring runs along a single shell near l = 15° and 195°. These modes are always mirror-symmetric and always removable by the bodies. I fitted bodies only for phase 0, not for the worst phase (17.8×).

**Caveat on symmetry (important)**
Perfect spherical shells centred on the CoE axis are *exactly* mirror-symmetric about that axis on the ring. With 3600 ring points the symmetric share comes out at **100%**. The roughly 55–60% seen at the standard 1° sampling is an artefact: the 360 sample points are not placed symmetrically about l = 15.16°, and the shells are finer than one sample. So clean shells give mirror-symmetric noise. It takes irregularities to make the noise asymmetric, such as the jitter used in the earlier blast-wave run or the web on top. The checkerboard does not have this property, because its grid is not aligned with the CoE axis.

**Scaling and convergence**
- **Contrast and density.** The pattern is linear in A (Fig 2c): at A = 0.03, multiply by 0.3. It is also linear in matter density, so with ordinary matter only (no dark matter) it is about **6× lower**: m ≥ 13 ≈ 1.2× at 171 Mpc and about 1× at 150 Mpc, near the CMB level.
- **Distance from the ring.** Summing out to 2 Gpc instead of 1 Gpc changes results by about 1%.
- **Cell size.** Finer sampling of the waves (19 instead of 28.6 Mpc) raises m ≥ 13 by about 1.3×.
- **Ring sampling.** With 3600 ring points, m ≥ 13 in total is unchanged (7.3× and 5.8×), but only 4.7× and 3.7× of it falls at m 13–180; the rest is above m = 180. m 1–12 at 3600 points is 2.0× and 1.1×.
- Treat all numbers as good to about a factor of 1.5.

**Files**
- `fig1_side_view_radial_shells.png`: side view with zooms of 3 Gpc, 1 Gpc and 300 Mpc on the ring point toward the CoE.
- `fig2_radial_ring_alone.png`: ring residual with the green ±1× band, spectrum by mode, and scaling.
- `fig3_bodies_cancel_noise_left.png`: before and after the bodies, the spectrum after, and a side view of the bodies.
- Code: `ri_core.py`, `run_fit.py`, `conv.py`, `fig_*.py`. Data: `final.json`, `final.npz`, `convergence.json`.
