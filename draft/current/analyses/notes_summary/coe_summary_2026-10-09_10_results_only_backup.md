# CoE shell cosmology: what we covered, Oct 9–10, 2026

## 1. The brief (LaTeX)
- Cleared up the 0.58× (brief, 11 objects) vs 0.34× (Writer's in-plane sketch, 12 objects) mismatch. The brief's runs are the source of record; Writer's tiled shell stays as an illustration.
- Fig 2 reworked: panel (a) is the leftover after each step; panel (b) is now a plain step-by-step demo of placing masses in the side view.
- Pages 8–10: all nine examples sorted globally from largest to smallest leftover (numbered 1–9).
- Added: the patch-ratio result, the asymmetry test, a light-through-structure section (Sachs–Wolfe, ISW, Rees–Sciama, lensing), and the faster-than-light recession note (caveats plus one clause in the conclusions).
- Layout: page 13 charts stacked vertically, legends below graphs, taller charts, references not forced onto a new page, arXiv links.

## 2. Figure standards adopted
- Every residual-around-the-ring plot shades the ±1× CMB band.
- Side views show only the modeled region (our 14.3 Gpc sphere plus the shell to ~30 Gpc, including the core cap).
- The CoE direction is a vector from us toward the CoE.
- Before/after views sit side by side with no arrows: clusters red, voids blue, anything that moved green in both views.
- Every before/after view has the before and after residual underneath.
- One run per image with moderate fonts; don't cram figures.
- Area-change maps use Doppler colours: blue = compressed, red = expanded, white = no change.
- Exception: on force maps each body is drawn as a body-sized arrow along its net force.

## 3. Placement fits and solution space
- Joint inside-plus-outside fits; runs with the core open and bodies allowed inside our sphere.
- Growth runs: random and foam rule, with the core filled, the visible sphere filled, or both.
- Statistical framing: draw random universes from real statistics, estimate the rare tail, measure the "good" region; entropy as the surface to optimize.
- Symmetry: real CMB anomalies are mostly asymmetric, so we measure the symmetric share of every leftover.
- Large scales (m=1–12) can always be cancelled; the web leaves fine-scale (m13+) texture.

## 4. Void–cluster network rules
- Node-link census: cluster–cluster, cluster–void and void–void links; connections per node, link lengths, angles between links; no-crossing rule.
- Voids effectively "attract": they push matter out, drain walls and merge.
- Literature: Clampitt & Jain ([arXiv:1404.1834](https://arxiv.org/abs/1404.1834)), whose real voids are about half as deep as ours.

## 5. Periodic structure and blast-wave energy
- Einasto et al. 1997 periodicity, ~120 h⁻¹ Mpc ≈ 171 Mpc ([Paper I, astro-ph/9701018](https://arxiv.org/abs/astro-ph/9701018); [Paper II Fig 2, astro-ph/9704127](https://arxiv.org/abs/astro-ph/9704127)): checkerboard and radial-shell backgrounds across the visible sphere.
- Wave energy (order of magnitude): ~10⁶³–10⁶⁵ J in our sphere, 10⁶⁵–10⁶⁷ J in the shell, ~10⁷⁴ J total from Sedov–Taylor. This implies a faster-than-light edge, flagged as next-stage relativistic LTB work.
- Model ignores dark matter and compares observables only.

## 6. Minimal-edit tests (start from a random filled web, edit to reach CMB level)
- First the fitter piled masses on the CMB ring. We reverted to the joint fitter with a 3 Gpc ring backstop and filled the core cap.
- With deletion allowed: reached 0.87–1.00×, but only by deleting bodies near the ring.
- No deletion (move and resize only): two foam webs reached 1×. Foam s3 needed 489 edits (2.8% of the web, 0.995×) and foam s1 needed 614 (3.6%, 0.95×). Near-edge bodies moved a median ~2 Gpc, others under 1 Gpc; most masses changed less than 2×.
- The 3 Gpc wall was a judgment call, not physics. It was replaced by an entropy rule: no edit may make the web more ordered.
- Under the entropy rule nothing reaches 1×. Best is foam s2: 1223× to 2.4× with 638 edits (3.8%). m1–12 drops to ~0.01×; the rest is m13+ ripple. No pile-up at the ring (the edits partly refill the gap left by the growth rule).

## 7. Mesh, area, entropy and energy (foam s2)
- Non-crossing neighbour mesh: cells near the ring grow (median 0.60 to 1.20 Gpc² within 3 Gpc).
- Full Delaunay triangulation (no cutoff, every spot covered): near the ring most spots are compressed (big gap triangles split), but the population evens out (median 0.79 to 1.35 Gpc²; area entropy 4.70 to 5.03). Whole slab: median 0.59 to 0.66 Gpc², about 6 connections per node, median link 1.77 to 1.87 Gpc.
- Energy: near the ring it rises (~+23%); farther out it falls (~−14% at 3–6 Gpc); overall the web is ~3% more bound. The entropy term is negligible next to gravity.
- Density field view: overall clumpiness rises only 2–5%; void-trough area near the ring falls from 19% to 12%.

## 8. Crystal analogy
- FCC lattice of clusters with voids in tetrahedral holes (4-cluster pyramids), 171 Mpc cell; octahedral holes for comparison.
- Ladder from perfect to random: defects, jitter, polycrystal. Our foam web sits *past* random in disorder and is the least bound; the edits nudge it slightly toward the lattice.
- 2D triangular crystal, equal numbers, voids as bare negative mass: the lowest-energy state is full segregation. That's unphysical, because it ignores where the void's missing matter goes.
- Balance (compensation) rule: each void's missing mass M_def = |δ|ρ(4/3)πR³ sits in clusters around it, in a wall at R to 1.3R; N ≈ 0.9 + M_def / 2.1×10⁵⁰ kg. At lowest energy, clusters bunch into knots where 2–3 voids meet, like the cosmic web.

## 9. Forces
- Pair energies (foam s2): cluster–void repulsion +4.7×10⁶⁹ J dominates cluster–cluster −2.1×10⁶⁹ J and void–void −0.7×10⁶⁹ J; edits change each by 1–3%.
- Moved bodies end ~20% closer to force balance and are not pushed toward the ring.

## 10. In progress
- New random web built with the balance rule plus the project's size/abundance rules, then move/resize under the entropy rule toward CMB level (3 seeds). Deliverables: side views, coloured mesh before/after with residuals, and stacked histograms.

## Next steps raised
- Share walls between neighbouring voids in the balance rule.
- Test how results depend on the ring rule (2 or 5 Gpc, smooth penalty).
- Full 3D (tetrahedral) version of the mesh statistics.
- Relativistic LTB treatment for the faster-than-light edge.
