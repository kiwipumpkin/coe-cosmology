# void_mix follow-up (Chris: results for 1–3 masses; "one to four masses and voids"; how much the helpful layouts helped)

- `followup_report.py` → `helped_vs_hurt_by_k.png`, `baseline_masses_vs_masses_voids.png`, `how_much_helped_k0.png`, `counts.json`. `summary_table.md` was generated inline from the same JSONs.
- `fit_masses_voids.py` → `fit_masses_voids.json`: the k=1–4 refits. Two variants: masses only (M>0) and masses or voids (sign free). Both use the same procedure:
  - objects in the ring plane, with free angle and free distance 15–42.9 Gpc from us;
  - 150 random starts plus subsets of the fitted four;
  - compensation as in `with_clumps` (net mass taken out of, or put back into, the surround region in proportion to the exp-lock density);
  - no two placed objects may overlap.
  - Size cap: |M| ≤ ρ̄·(4/3)π(s − 14.3 Gpc)³. A void can be at most completely empty (δ = −1), and its radius grows with its deficit until it just reaches our visible sphere. Clusters get the same cap.
  - All fitted objects are 1–4% of the cap, so the cap never binds.
  - Without the no-overlap rule, k=3 found a 94× "solution" that was a 1.7e54 kg empty void with an equal cluster sitting inside it. It was rejected.
- `ensemble_fitted.py` → `ensemble_fitted.json`, `raw_*.npz`: void_mix random web (seed as void_mix) on top of the refit baselines, at N=200, level b, 50% and 90% voids.
- Masses-only refits beat the earlier fitted-four subsets (k=1–4: 211.5/168.6/131.8/68.9 vs 221.8/185.8/182.1/88.8). That is because each object now has its own distance, where the earlier fit used one shared distance.
