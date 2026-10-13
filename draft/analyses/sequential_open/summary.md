# Sequential placement with the CoE core open and our sphere open — status 9 Oct 2026, ~4:55 PM CT

All leftovers are RMS around the CMB ring as multiples of the 1× CMB level (the "floor"). SR + central mass alone leaves 237×.


## 0. CORRECTED fill (Chris, 9 Oct): our sphere + shell out to 30 Gpc (tidal depth), CoE core NOT filled — supersedes the core-filled versions in §1

Same B growth rules and packing; region = inside 14.3 Gpc (no straddling the ring, nothing within 0.5 Gpc of us) plus the shell out to 30 Gpc from us, core excluded (100,443 Gpc³). Runs fill until no room (B) or to the scaled cap (shallow: 2,241 voids).

| run (figure seed) | old final / best | tidal-depth final / best | bodies | inside our sphere | 5-seed range of final | symmetric share (5 seeds) |
|---|---|---|---|---|---|---|
| B random (s2) | 1045 / 213 | 1744 / 305 | 14,987 | 1,542 | 1744–2830 | 16–86% |
| B foam rule (s3) | 610 / 221 | 1136 / 219 | 15,253 | 1,561 | 979–1406 | 26–79% |
| B CMB-aware (s1; 2 seeds) | 169 / 36 | 344 / 37 | 16,744 | 1,645 | 344–639 | 30–52% |
| B CMB-aware shallow (s2) | 163 / 163 | 133 / 133 | 2,241 | 77 | 133–195 | 79–87% |

Reading: filling our sphere roughly doubles the leftover of the blind rules (bodies inside our sphere alone pull 900–1230× for random/foam), the CMB-aware run ends 2× higher with the same best point (~37×), and only the shallow-void run is a bit lower (133 vs 163×). Symmetry: random/foam scatter widely around 50% (none ≥90%); shallow sits at ~80% because it fits the symmetric 237× target. Figures: fig7_growth_tidal.png, fig8_growth_tidal.png (residual panels use a y-axis that is linear inside ±1 and logarithmic beyond, so the green ±1× band is always visible; fig8_growth_variants.png re-rendered the same way), symmetry_tidal_vs_leftover.png, symmetry_tidal_by_mode.png.

## 0b. 30 Gpc fill INCLUDING the part of the CoE core inside 30 Gpc (rest of core empty)

The core cap within 30 Gpc of us is 12,689 Gpc³, which is 11% of the 113,132 Gpc³ fill region. It is packed with the same rules, and there is no wall at the core boundary.

| run (figure seed) | core-empty tidal final / best | incl. core part final / best | bodies (in core part / inside our sphere) | 5-seed final range (core-empty → incl. core) | symmetric share over seeds |
|---|---|---|---|---|---|
| B random (s2) | 1744 / 305 | 2106 / 285 | 16,905 (1,932 / 1,441) | 1744–2830 → 1958–3167 | 16–86% → 19–76% |
| B foam rule (s3) | 1136 / 219 | 2164 / 219 | 17,367 (2,024 / 1,547) | 979–1406 → 794–2164 | 26–79% → 3–66% |
| B CMB-aware (s1; 2 seeds) | 344 / 37 | 601 / 23 | 19,257 (2,292 / 1,652) | 344–639 → 359–601 | 30–52% → 65–83% |
| B CMB-aware shallow (s2) | 133 / 133 | 134 / 134 | 2,524 (270 / 103) | 133–195 → 103–193 | 79–87% → 67–84% |

Reading: including the nearby core cap does not help. Final leftovers are similar or higher (foam s3 doubles). The CMB-aware best point dips (37 → 23×) but the run ends higher. The ~2,000 core-cap bodies have a ring pattern of their own of 160–490×. Symmetry stays broadly scattered for the blind rules, and foam even turns strongly antisymmetric in 3 of 5 seeds. Figures: fig7_growth_tidalcore.png (old dashed grey, core-empty dashed colour, incl. core solid), fig8_growth_tidalcore.png (core cap inside 30 Gpc shaded pink; symlog residual axes with the green ±1× band), symmetry_tidalcore_vs_leftover.png, symmetry_tidalcore_by_mode.png. Bodies-on-top on these webs (≥1.5 Gpc off the ring) are queued.

## 3b. Bodies on top with added bodies kept ≥1.5 Gpc off the ring (running; supersedes §3/§4 ring-hugging fits, moved to web_add/old_ring_hugging/)
- no background: 237× → 0.98× with 68 bodies (34 masses + 34 voids) — reaches 1×.
- CMB-aware shallow tidal web (s1): 157× → 0.94× with 77 bodies.
- quasi-lattice s1: → 1.83× with 130 bodies, then no allowed body helps (the lattice's 2.2× fine texture cannot be removed from 1.5 Gpc away).
- radial wave, random/foam/CMB-aware tidal webs and remaining seeds: running (web_add/*_omp0_ex1.5.json; render with `python3 report_web.py` -> web_*_ex1.5.png).

Standalone wave/lattice figure (no bodies): fig_wave_lattice_alone.png.

## 1. (superseded) The four B growth runs, three filled versions each (done)

"Core filled" = bodies may go anywhere in the CoE interior, our sphere stays empty. "Our sphere filled" = bodies may go inside 14.3 Gpc, core blocked as before. "Both" = both open. Same packing rule as the old runs (void volume fraction ~0.28–0.34), so filled runs place about twice as many bodies.

| run | old (core blocked, 14.3–42.9) final / best | core filled final / best | our sphere filled final / best | both filled final / best |
|---|---|---|---|---|
| B random | 1045 / 213 (40,474 bodies) | 1081 / 214 (87,654) | 1188 / 230 (40,600) | 1980 / 237 (88,228) |
| B foam rule | 610 / 221 (41,818) | 630 / 234 (91,539) | 742 / 218 (42,861) | 1063 / 234 (92,106) |
| B CMB-aware | 169 / 36 (45,283) | 438 / 26 (100,044) | 370 / 36 (47,110) | 378 / 22 (102,214) |
| B CMB-aware, shallow voids | 163 / 163 (6,000) | 168 / 168 (13,912) | 182 / 182 (6,267) | 158 / 158 (14,183) |

About half the bodies (~48–55 k) end up in the core when it is open. When our sphere is open, ~1,350–1,650 go inside it (the shallow run puts almost none there).

Does filling help? No. For the blind rules (random, foam) every filled version ends higher than the old run. Filling both is worst (random 1980× vs 1045×), because bodies inside our sphere sit right next to the ring and their pull alone is ~360–790×. Core bodies are farther away and pull less (core group on its own 200–280×). But they are not a smooth background either: they add their own random pattern. In the CMB-aware run the core lets the low point drop a bit (36 → 22–26×), but the run still ends 2–2.5× higher than the old one. The shallow-void run barely changes. So core matter does not tidy up the ring; it adds an unorganized contribution, as Chris suspected.

Figures: fig7_growth_variants.png (curves vs the old dashed runs, log axes, green 1× line, dotted 237×/89×), fig8_growth_variants.png (side views, core shaded and marked open or blocked, green ±1× band). The original six-run figures with the A near fits are fig7_open_curves.png and fig8_open_universes.png (A near fit, deep voids: 27.4 → 14.9×; A near fit, shallow voids: 135 → 16.6×; both with 48 bodies, some placed inside our sphere and in the core).

## 2. Mirror symmetry of the combined ring pattern about the CoE axis (first pass, symmetry.json, symmetry_*.png)

Seeds so far: random 5 per version; foam 3–4; CMB-aware 1–2; shallow 5. The symmetric share is the fraction of m=1–12 power that is mirror-symmetric about ℓ 15°/195° (50% = no preference).
- Random growth: 12–69%, typically ~50%. With our sphere filled it is often strongly *anti*-symmetric (12–26%). The pattern is hundreds to thousands × the 1× level.
- Foam growth: 41–73%, mostly ~60–70%. 560–1600×.
- CMB-aware: 46–71%, 370–440×.
- Shallow CMB-aware: 67–97%, mostly ~80%, 117–182×. This rule leans symmetric because it is fitting the symmetric 237× target.
- Near-symmetric (≥90%) happened in 1 of 44 filled runs (a shallow, core-filled run at 97%). It is not typical for a self-assembled web. More foam seeds (s4, s5) are still running.

## 3. Fixed web background + any bodies on top (first pass)

Method (omp_add.py): the background is fixed. Bodies are drawn from a bank of 80,000 candidate positions: half within 2 Gpc of the ring, a quarter inside our sphere, a quarter anywhere including the core. Each body's mass is fitted as a cluster (≤2e52 kg) or a deep void (0.5–4 Gpc). No body may straddle the ring. The method adds the best next body and refits all masses together. A body that would need less than 1e49 kg is removed. It stops at <1× or when no allowed body helps. Bodies are not yet checked for overlap with each other or with the background.

Result so far (CMB-aware web, our sphere filled): 370× → 5.35× with 170 added bodies (135 masses + 35 voids, 167 of them within 2 Gpc of the ring). m=1–5 is 0.11× and m=6–12 is 0.26×, so **the large scales reach CMB level easily**. What is left is all fine-scale (m ≥ 13: 5.3×) and has no symmetry preference (56%). See the ring plot in web_side_views.png: spiky, noise-like.

Reference with no background (SR + central mass only): 237× → 2.5× with 91 bodies; m ≥ 13 = 2.4×. So with this solver, even an empty background stalls at ~2.5× of sample-to-sample texture. That is the solver/body-rule limit. A body ≥1e49 kg near the ring is too coarse to trim single-degree spikes. The 1× target was not reached by any background yet.

The random/foam webs (hundreds to thousands × before) are queued now (log_omp.txt, web_add/*_omp0.json). An earlier joint-fit version took 16 bodies to get the foam web from 1017× to 333×.

## 4. Radial blast-wave and quasi-lattice backgrounds (done, 5 seeds each)

Construction (wave_bg.py):
- **Radial wave:** δ = A cos(2π r_CoE / 150 Mpc + jitter), concentric shells about the CoE through the whole filled region (core + shell + our sphere).
- **Quasi-lattice:** δ = A/3 (cos kx′ + cos ky′ + cos kz′) with 120 Mpc spacing, random orientation and jittered coordinates.
- **Amplitude and jitter:** peak contrast A = ±0.15 (typical contrast of the real web smoothed to ~60–75 Mpc). Jitter is a smooth random field (0.5–2 Gpc wavelengths, rms ¼ wavelength).
- **How the pull is computed:** the coherent spherical part of the wave is exact (shell theorem). The jittered field within 1 Gpc of the ring is summed cell by cell (5.4 M cells of 37.5 Mpc). The lattice's far field cancels by periodicity.
- **Scaling:** the pattern is linear in A, so A = 0.015 gives one tenth of every number below.
- **Caveat:** 150 Mpc shells are finer than the ring sampling (1° ≈ 250 Mpc at the ring). The m ≥ 13 part partly aliases fine structure onto the 360 samples, so treat the 13+ amplitude as order-of-magnitude.

Background alone (the background's own ring pattern, without SR):

| background | own rms | m=1–5 | m=6–12 | m≥13 | symmetric share (m=1–12) |
|---|---|---|---|---|---|
| radial wave, 150 Mpc (5 seeds) | 15.6–17.7× | 1.6–2.9× | 3.1–4.2× | 15–17× | **78–95%** |
| quasi-lattice, 120 Mpc (5 seeds) | 2.2–2.3× | 0.3× | 0.3–0.6× | 2.1–2.3× | 41–70% |
| random web (filled, s2) | 1000–2000× | 700–1800× | 400–590× | 410–650× | 37–63% |
| foam web (filled, s3) | 600–980× | 285–760× | 340–440× | 300–440× | 56–61% |

The radial wave is the only background with a clearly mirror-symmetric large-scale pattern (78–95%). Its coherent shells are centred on the CoE, which lies on the axis, so they are exactly symmetric, and only the jitter breaks that. Its texture is almost all fine-scale. The lattice is tiny (2×) and unbiased. Both are 30–1000× quieter than the coarse 1.5 Gpc webs. Combined with SR + central mass, all of them still read ~237× before adding bodies, because SR dominates.

Bodies added on top (same method as §3):

| background | bodies added | masses / voids | left | m=1–5 | m=6–12 | m≥13 | symmetric share after |
|---|---|---|---|---|---|---|---|
| none | 91 | 46 / 45 | 2.5× | 0.15 | 0.65 | 2.4 | 52% |
| quasi-lattice (5 seeds) | 78–93 | 29–48 / 39–53 | 2.5–4.1× | 0.04–0.4 | 0.06–0.8 | 2.6–4.0 | 33–54% |
| radial wave (5 seeds) | 96–119 | 60–76 / 36–43 | 10.1–12.0× | 0.09–0.23 | 0.11–0.37 | 10–12 | 23–61% |

Reading:
- In every case m=1–12 drops well under 1× with ~100 bodies.
- The lattice is indistinguishable from no background: its fine texture is below what the bodies leave anyway.
- The radial wave leaves ~11× of fine-scale texture that the added bodies cannot remove, about 4× more than the others. The bodies trim the large scales but not the 150 Mpc shell texture.
- After adding bodies, the leftover's symmetric share scatters around 50%. Even the symmetric wave ends with a non-symmetric, noise-like leftover.
- Charts: web_bodies_needed.png, web_spectrum.png, web_symmetry.png, web_side_views.png. The side views show the added bodies on top of the background, with a 3 Gpc zoom of the shells or lattice in the ring plane and a green ±1× band on every ring plot.

## Still running / not done
- Web + bodies for the random/foam/CMB-aware/shallow filled webs (queued, ~5–10 min each, 3 at a time); report_web.py re-renders all four web_* charts when they finish.
- Extra foam seeds s4, s5 for the symmetry statistics (running, checkpointed); then rerun symmetry.py.
- No background has reached under 1× yet. The stall at ~2.5× even with no background means the next step is a better fine-scale solver (continuous position refinement of the near-ring bodies) or a decision that 13+ texture below a few × is acceptable.

## 5. Bodies on top with a uniform-in-volume placement prior (prior_add.py; standard region = our sphere + shell to 30 Gpc incl. core cap)

- **Prior and penalty:** candidate positions are uniform in volume over the modeled region (80,000 per run, no extra candidates near the ring). For a body d from the ring, F(d) is the fraction of the region's volume within d of the ring. Measured from 4 M uniform points, F ≈ 0.025 d² near the ring, so F reaches 1 at d0 ≈ 6.3 Gpc. The penalty is −log F(d) = 2 log(d0/d). That is the (d0/d)² law in log form, with d0 set by the geometry rather than chosen.
- **What the penalty measures:** under the prior each body's −log F is Exp(1), so the layout's total is Gamma(N,1). Σ log F is the layout's log-probability of being that close to the ring, and the prior expects −N.
- **Objective:** χ²/2 + λ·Σ(−log F), with χ² over the 360 ring samples and 1× per sample as the noise. λ=1 is the plain Bayesian posterior. λ = 0, 20, 300, 3000 are comparisons.
- **Body rules:** clusters ≤2e52 kg, deep voids R 0.5–4 Gpc, nothing under 1e49 kg, no straddling the ring.
- **Fit procedure:** add the candidate with the best penalized gain, then refit all masses together. Stop at <1× or when no candidate improves the penalized objective.

Results: leftover (× 1×), bodies, bodies within 3 Gpc of the ring vs the uniform expectation, and log-probability vs the prior's expectation.

| background | earlier ring-weighted fit (≥1.5 Gpc) | λ=0 uniform bank | λ=1 | λ=20 | λ=300 | λ=3000 |
|---|---|---|---|---|---|---|
| none | 0.97×, 68, 61 in 3 Gpc (exp 1.5), logp −308 (exp −68) | 4.8×, 63, 41 (1.4), −225 | 4.2×, 61, 41 (1.4), −227 | 5.5×, 53, 37, −200 | 13.0×, 33, 20, −114 | 22.2×, 26, 8 (0.6), −59 (exp −26), p=4e-7 |
| CMB-aware shallow web | 0.90×, 92, 89 (2.1), −435 | 13.9×, 68, 51, −291 | 13.9×, 65, 53, −296 | 15.8×, 52, 45, −249 | 19.4×, 37, 30, −168 | 35.8×, 22, 14, −79 |
| random web | running | 138.6×, 182, 117, −700 | 139.0×, 143, 117, −693 | 140.2×, 124, 110, −649 | 153.8×, 159, 88, −532 | 187.4×, 108, 60, −375 |
| foam web | 1.74×, 311, 310 (7.0), −1383 | 30.6×, 112, 88, −525 | 30.8×, 100, 87, −518 | 33.5×, 96, 84, −490 | 41.9×, 84, 68, −406 | 73.9×, 51, 40, −242 |
| radial wave | 12.9×, 133, 132 (3.0), −631 | 15.3×, 66, 41, −232 | 15.5×, 51, 41, −227 | 16.3×, 47, 34, −190 | 19.7×, 32, 19, −113 | 31.2×, 27, 6, −49, p=3e-4 |
| quasi-lattice | 1.83×, 130, 123 (2.9), −602 | 4.3×, 68, 42, −237 | 4.9×, 57, 40, −227 | 5.5×, 49, 37, −210 | 13.7×, 35, 19, −113 | 22.4×, 26, 8, −59 |

Reading:
- **The prior-weighted penalty barely moves the fit (λ=1).** One extra body near the ring buys tens to thousands of χ² units, while the penalty is at most a few units. So the "honest" posterior still piles bodies near the ring: 60–80% of them within 3 Gpc, where the prior expects about 2%. The closeness log-probability is 3–5× more negative than the prior's −N, with p-values of 1e-36 to 1e-143.
- **With a uniform-volume bank, no background reaches 1×.** Even with no background the fit stalls at about 4–5×, and the random/foam webs stall at 139× and 31×. Earlier, 1× was reached (no background, shallow web) only because the candidate bank put half its candidates within 1.5–4.5 Gpc of the ring.
- **Forcing prior-consistency costs leftover.** At λ=3000 the no-background layout is still prior-unlikely (p=4e-7) and the leftover rises to 22×. For the wave it reaches p=3e-4 at 31×. Every background shows the same trade-off: the closer the layout is to "placed at random in volume", the further the leftover is from 1×.
- **What's left is always fine-scale.** In all fits m=1–12 is ≤0.6× and the remainder is m≥13. Its symmetric share is 34–80%, scattered with no consistent value.

Figures: prior_distance_cdf.png (body distance from the ring vs the uniform-volume expectation, per background and λ), prior_left_vs_logp.png, prior_side_views.png (λ=1; modeled region only, ±32 Gpc, arrow to the CoE, symlog residuals with the green ±1× band). Table: web_add/prior_table.json.
