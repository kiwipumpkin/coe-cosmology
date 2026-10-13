# Einasto-style checkerboard through our visible universe

## The period, checked against the paper
Einasto et al. 1997, Nature 385, 139, [arXiv:astro-ph/9701018](https://arxiv.org/abs/astro-ph/9701018):
- Power-spectrum peak at k0 = 0.052 h/Mpc, so λ0 = 2π/k0 = **120 ± 15 h⁻¹ Mpc**. The grid step is 120 ± 20 h⁻¹ Mpc, which they describe as "distances between high-density regions across voids". The paper says the whole distribution "resembles a three-dimensional chessboard".
- The figure is in **h⁻¹ Mpc**. With h = 0.7 that is **171 Mpc**, not 120 Mpc. (The earlier lattice run in sequential_open used 120 Mpc without converting.)
- 120 h⁻¹ Mpc is the **full repeat**: one supercluster region to the next one across a void. So:
  - in a chessboard of alternating cubes, each cube is half a period, **85.7 Mpc**;
  - in Einasto's own picture (one void per cell, superclusters on the walls and corners), each cell is one full period, **171 Mpc**.
- Counts from us to the ring (14.3 Gpc): **167 chess cells (334 across)**, or **83 full periods (167 across)**. The "about 115 cells to the ring" figure matches 120 Mpc cells, i.e. without the h conversion. With h = 0.7 the correct numbers are 167 cells or 83 periods.
- At the ring, one chess cell spans 0.34° and one period spans 0.69°. Going round the ring, one period corresponds to mode **m ≈ 524**, far above m = 1–12.

## Two versions of the background
- **chess**: cubes of 85.7 Mpc alternating +A and −A.
- **Einasto cell**: 171 Mpc cells with a void at each centre and overdense walls, edges and corners: δ = −A(cos kx + cos ky + cos kz)/√1.5.
- Both use the same rms A. The grid axes follow the supergalactic axes, as in Einasto's Fig 1. A further 8 runs used random orientations and phases.

**Contrast.** I used A = 0.10. Linear CDM with σ8 = 0.81 gives an rms contrast of 0.19 for a cell-sized region (37 h⁻¹ Mpc). Einasto's peak is only about 1.4× above CDM in power, so about 30% of that power is the regular part, which gives √0.3 × 0.19 ≈ 0.10. The pattern scales exactly linearly with A (Fig 4c):
- A = 0.03 (a BAO-level ripple) gives 0.3 times these numbers.
- A = 0.19 (every fluctuation at this scale arranged as a checkerboard) gives 1.9 times.

The density used is the same total-matter density as all the earlier runs. With baryons only, as we're ignoring dark matter, everything is about 6× smaller.

## Results (main run: supergalactic orientation, A = 0.10)
| | chessboard | Einasto cell |
|---|---|---|
| Checkerboard alone, total | 4.9× | 10.2× |
| … m = 1–12 | 4.0× (almost all m = 1, from the cells right around us) | 3.6× |
| … m ≥ 13 | **2.7×** | **9.6×** |
| … symmetric share | 39% | 53% |
| 9 orientations/phases: m ≥ 13 | 2.6–3.6× | 9.6–10.1× |
| Combined with the 237× residual, before bodies | 237× (m ≥ 13: 2.7×) | 237× (m ≥ 13: 9.6×) |
| **After 270 bodies**, m = 1–5 / 6–12 | 0.000× / 0.000× | 0.000× / 0.006× |
| After bodies, m ≥ 13 (what's left) | **2.7×** | **9.3×** |
| After bodies, symmetric share | 52% | 55% |

- **Bodies.** 270 bodies: 144 point masses and 126 voids (negative masses). Each is at most 2×10⁵² kg, kept at least 1.5 Gpc off the ring (in practice all are at least 4.7 Gpc away) and outside the CoE core. 150 of them cancel the 237× residual and 120 cancel the checkerboard's own low modes. The bodies add about 0.6–2.4× of their own fine-scale texture.
- **The point.** Placed bodies can drive m = 1–12 to essentially zero. The checkerboard still leaves a flat floor across all m ≥ 13 that no distant body can remove: fitting the checkerboard alone only takes m ≥ 13 from 2.7 to 2.6× (chess) and from 9.6 to 9.3× (Einasto). That floor shows no mirror symmetry (a share of about 50% is what pure noise gives) and looks like noise, but it is fixed entirely by the grid. A grid is about as orderly as structure gets, and it still shows up on the ring as noise.

## Caveats (order of magnitude)
- **Ring sampling.** The cells are much smaller than the 1° ring sampling (m ≤ 180), so their power folds back into all the sampled modes. With 3600 ring points instead of 360, m = 1–12 is unchanged for chess (4.1×), and m = 13–180 drops to 1.4× (chess) and 6.4× (Einasto), with the rest landing above m = 180.
- **Resolution.** Sampling each cell more finely (19 instead of 28.6 Mpc) raises m ≥ 13 by 1.3–1.5×. Treat the floor values as good to about a factor of 1.5.
- **Our own cell.** The chess m = 1 term is the pull of the few cells around us. In the real CMB a dipole is removed as our own motion.
- **Distance range.** Only cells within 1 Gpc of the ring are summed, plus 1 Gpc around us. Going to 0.5 or 1.5 Gpc changes the result by less than 2%, because a periodic field cancels at larger distances.
- **Figure 1** is a recreation from the numbers stated in the paper (k0, λ0, the 1.4× excess, slopes n = 1 and −1.8). It contains no catalogue data, and the choice of which grid lines carry supercluster chains is illustrative.
- **runs.json** holds the 9-orientation ensemble for the checkerboard alone. Its "after_*" fields come from an earlier one-step body fit and are superseded by final.json.

## Files
- fig1_einasto_recreation.png: power spectrum, correlation function and grid picture in the style of Einasto 1997.
- fig2_cell_and_scaling.png: one Einasto cell (void at the centre, superclusters on the corners and edges), a chessboard block, a length ladder, and how a cell maps onto the ring.
- fig3_side_view_checkerboard.png: the full 14.3 Gpc sphere with zooms of 3 Gpc, 1 Gpc and 300 Mpc on the ring point toward the CoE.
- fig4_checkerboard_ring_alone.png: ring residual with the green ±1× band, spectrum by mode, and scaling with contrast.
- fig5_bodies_cancel_noise_left.png: before and after the bodies, the spectrum after, and a side view of the bodies.
- Code: cb_core.py, conv.py, run_cb.py, main_fit.py, combine.py, fig_*.py. Data: final.json/npz, convergence.json, runs.json.
