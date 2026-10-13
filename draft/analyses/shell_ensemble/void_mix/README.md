# Void share of the random web varied (Chris: "I don't think it should be half voids and half clusters. I think you should range that.")

Scripts: `ensemble_void_mix.py` (run: `python3 ensemble_void_mix.py 1000 <k list> 200,500 abc <tag>`; it was run as 5 parallel jobs, k=0..4), and `report.py` (produces `counts.json`, `summary_table.md`, `helped_vs_hurt_void_mix.png`).
Everything else matches `../with_clumps`: surround region, exp-lock gradient placement, clump compensation, overlap redraw, scoring, and seeds. Originals are untouched.

- **Void share:** f_v = 10, 25, 50, 75 or 90% of the N objects. In the "random" case each layout draws f_v uniformly from 0–1, with at least 1 void and 1 cluster.
- **Zero net mass (same rule as the original code):** each object's |M| is drawn log-uniform in its level. The void deficits are then rescaled together so that total void deficit = total cluster mass.
  The smooth exp-lock profile is therefore unchanged.
- **δ cap:** voids nominally have δ = −0.5. Rescaled by s, they deepen to δ = −0.5·s, capped at δ = −1 (completely empty). Past the cap they get larger instead (radius = deficit volume at δ = −1).
  At 10% and 25% voids every layout hits the cap. At 50% and above, none does. In the random mix, about 33% of layouts hit it.
  The tidal pattern treats objects as points, so depth and radius only matter for the clump-overlap redraw.
- **Edge cases:**
  - Random mix: a few layouts with only 1–2 voids get very large voids (largest: 17.8 Gpc at N=500 with large masses), which is bigger than our visible sphere. This is unphysical, but it is rare.
  - f_v = 50%: reproduces `../with_clumps` exactly for k=0. For k≥1 there are tiny differences (for example 173 vs 180 at k=4, N=200, medium) because void radii now follow the rescaled deficit.
- **Why more voids helps more often:** the web gets weaker. The total object strength sqrt(Σ M²) at N=200, medium is 6.2e51 kg at 10% voids, 2.1e51 at 50% and 0.7e51 at 90%. That is because a few deep voids are heavier objects than many shallow ones.
  A weaker web disturbs the leftover less (closer to a coin flip). It does not reach lower: the best single layouts come from the strong, few-void webs.
