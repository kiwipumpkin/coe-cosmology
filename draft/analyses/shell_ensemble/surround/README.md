# Surround variant of the shell ensemble (Chris: "points surrounding the visible sphere ... add some beyond D+d" + "we have a density gradient")

Scripts: `ensemble_surround.py` (run: `python3 ensemble_surround.py 1000`), `counts_surround.py`, `figs_surround.py`.
Physics, residual and scoring are imported verbatim from `../ensemble.py` (field3, amps, summarize, LEVELS, KG, residual R from scan.py). Originals untouched.

## What the original placement was
Uniform in the volume of the CoE-centred shell 43.9–72.5 Gpc (D−d..D+d, D=58.2, d=R_vis=14.3), over all directions around the CoE (full 3D, not planar), rejecting centres inside the 14.3 Gpc vision sphere.
Seen from us that is lopsided: 97% of objects lie in the hemisphere facing the CoE, ~0% in the cap facing away, only 5% within 28.6 Gpc of us, median distance from us 84 Gpc, and nothing beyond D+d.
Reason: the visible sphere exactly spans the shell (it touches D−d on the CoE side and D+d on the far side), so the shell cannot wrap the far side of the visible sphere.

## What changed (only placement)
- Region: all space between 14.3 and 42.9 Gpc from us (3 R_vis, i.e. two visible-radii of space around the visible sphere), every direction, with CoE distance ≥ D−d.
  This covers the shell D−d..D+d **and beyond D+d**. The farthest point from the CoE is 58.2+42.9 = 101.1 Gpc = D+3d.
  The CoE interior (< D−d) stays empty, as in the original, so the sky patch pointing straight at the CoE has no nearby objects. That is geometry: the visible sphere touches D−d there.
- Density gradient: exp-lock law ρ = ρ0 e^(−r/L), L = 922.65 Gpc. L comes from the CatWISE lock δ = [ρ(D−d)−ρ(D+d)]/ρ(D) = 0.031, L = d/asinh(0.0155).
  Under the exp law, L does not depend on D, so it is the same for D = 48.3 (exp-lock Migkas fit) and 58.2 (ensemble geometry). All δ-locked laws have the same local slope there (ρ'/ρ ≈ −0.00108 per Gpc).
  Cluster number density per unit volume ∝ ρ(r_CoE); void number density ∝ 1/ρ(r_CoE). In plain terms, clusters prefer the denser, CoE-ward side and voids the thinner, outer side.
  The beyond-D+d share therefore follows volume × density instead of being forced: ≈45% of objects (46% with the buffer).
  Forcing 25–30% would have created an artificial ×0.5 jump in object density at D+d, which contradicts a smooth 3% gradient.
- Unchanged: N = 200/500; mass levels a/b/c (log-uniform |M|); half clusters / half voids with zero net contrast; `bare` (centre > 14.3) and `buf` (whole object outside the vision sphere); 1,000 realizations; residual = R − Σ k_i field3_i; best-scaling (one free multiplier); all count definitions.
- Also run: the same surround region with no gradient (`unif_*`) for comparison.
- The tidal model `field3` was already 3D, so no extension was needed.
- Superseded first attempt (equal counts per sky direction; this piles objects into thin slivers next to us and so violates "density follows ρ"): `isotropic_per_direction_check/`. It gives even more "hurt".

## Result (bare, as drawn; out of 1,000)
Helped / hurt: tiny 464/536 (N200), 516/484 (N500); medium 363/637, 243/757; large 0/1000, 0/1000.
Original: 521/479, 509/491, 474/526, 438/562, 36/964, 6/994.
The gradient changes nothing measurable (a 3% density variation), and the surround+gradient vs surround-uniform differences are within noise.
Full table: `summary_table.md`. Figures: `helped_vs_hurt.png`, `example_realization.png`.
