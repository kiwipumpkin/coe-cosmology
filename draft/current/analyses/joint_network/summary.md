# Joint inside + outside fit with boundary continuity and cosmic-web rules

All files are in `/workspace/clump_scan/joint_network/`. I didn't touch the brief, `notes/` or `core_open/`.

## What was fitted
- **Bodies:** masses (1e49–2e52 kg) and voids (radius 0.5–4 Gpc, depth up to −0.7). They can go anywhere: inside our 14.3 Gpc sphere, outside it, and inside the explosion-centre (CoE) core. There's no wall and no density cap.
- **Body counts:** N = 8, 12, 16, 20, 24 and 32, half masses and half voids.
- **Starts:** 8 per case, with mirror layouts allowed and with mirror layouts forbidden (the same no-mirror rule as before).
- **Kernel:** the same Newtonian ring pull as every earlier run. Net mass is removed smoothly over our sphere, as in the all-open runs. The score is the leftover after solving for the residual with modes m=1–12, measured against the 1× target (1e-5).
- **Hard rules:** no body straddles our sphere's edge, and mass cores never overlap.
- **Soft rules:** every soft rule is a penalty with a tunable strength. The levels are:

| level | what's added | rule strength |
|-|-|-|
| no rules | nothing | – |
| continuity | inside vs outside must match in body count per volume, void share, mass sizes, void sizes and net density contrast (to within 0.5%). It's checked twice: in the 10–14.3 vs 14.3–19 Gpc shells, and for our whole sphere vs the equal-volume 14.3–18 Gpc shell outside it.* | 10 |
| web, loose | continuity (wider tolerances) + web rules (wider tolerances) + density trend allowed up to 100× | 3 |
| web, strict | continuity + web rules at the values you set: masses on void walls and between 2 voids; void overlap ≤50%; a mass at every void–void contact; 5–7 void neighbours per void; ≥3 links per mass; neighbouring masses 1–2 void diameters apart; each void's deficit matched by nearby masses within 30%; neighbourhoods net to ~0 | 30 |
| web + trend cap | as above, plus each void's depth and each mass's excess (over its own neighbourhood) held to the observed maximum-contrast trend (fits_v2: over_noDef / under_noDisp), +0.29 dex | 30 |

\*I added the whole-sphere check because, with the 10–14.3 / 14.3–19 shells alone, the fit just put every body at 5–8 Gpc and passed the rule trivially.

## Results (best_vs_N_logy.png)
| | best leftover, mirror allowed | best leftover, no mirror |
|-|-|-|
| earlier: inside only / outside only / all open (smooth limit) | 0.025 / 0.094 / 0.073 (N=12) | outside only, no mirror: 0.34 (N=12) |
| joint, no rules | **0.034** (N=24), 0.07 (N=12) | 1.6 (N=32)† |
| joint + continuity (rules met) | **0.77** (N=12) | **0.97** (N=20) |
| joint + loose web rules | 0.18–3.3, never with every rule met | 1.9–7, never with every rule met |
| joint + strict web rules | 10–80, never with every rule met | 26–170, never with every rule met |
| + density-trend cap | 7–80 while still sitting ~1000× (3 dex) above the trend. Actually holding to the cap: **260–310×**, i.e. no reduction from the 237× start (trend_cap_check.json) | 23–96 (same caveat) |

†The joint no-mirror search is under-explored: it does worse than the earlier dedicated outside no-mirror run (0.34×). Read the no-mirror joint row as an upper bound.

### What the rules cost
- **Continuity costs about 10–20×** (0.03–0.07× → 0.77×) and changes the kind of solution. The best joint solution with no rules splits 12 bodies inside our sphere and 12 outside. As soon as inside and outside must look alike, every rule-meeting solution moves **all** its bodies outside, at a typical 23–28 Gpc from us, past both check shells. The rule is satisfied because both sides of the boundary are empty, not because structure carries on across it. So continuity, in the sense you meant, is never actually met by structure on both sides.
- **The web rules cost much more.** No fit met every strict web rule. The rule that keeps breaking is "5–7 touching void neighbours per void" (11–12 of 12 voids fall short). The next two are "masses on void walls between 2 voids" and "void deficit balanced by neighbouring masses" (rule_breakdown.json). The best compromises leave 17–45×.
- **The density-trend cap leaves nothing.** At the observed maximum contrast for Gpc scales, a mass can hold only ~1e48–1e49 kg of excess over its neighbourhood. The fits need 1e51–2e52 kg, about 3 dex more, and voids have to be ~5–10× deeper than the trend allows. Holding to the cap gives back the starting leftover.

## Does the fitted web look typical? (network_vs_random.png, network_compare.json)
- **Null comparison:** for each key case, I made 100 random networks from the same bodies, placed isotropically and evenly by volume, then relaxed under the same rules but with no CMB term.
- **Geometry looks ordinary.** The fitted link lengths (~1.0–1.3 void diameters), links per void (3–5) and links per mass sit inside the random distributions (5th–95th percentile).
- **The random networks fail the strict web rules too.** Their median rule penalty is 0.17, against 0.12 for the fit. At N ≤ 32, these rules are hard to meet even with no CMB goal.
- **Two mild outliers:** the web-rule fits put more bodies inside our sphere than random networks do (top 1–2%). The trend-cap fit is aligned with the CoE axis (|cos| 0.97, top 2%). Orientation was left free.
- **The leftover is the extreme outlier.** All 600 random rule-obeying networks left 250–3,000×. None got within a factor of ~250 of the fitted layouts, so the chance by placement alone is below 1 in 100 for each case.
- **The amplitudes are far rarer than the geometry.** Each body's contrast over its own region is a median 350–2,800 ΛCDM standard deviations (σ) from the mean, up to ~7,000σ. That's a combined log10 probability of about −3,600 to −7,200. What's unusual isn't the shape of the web but how strong its bodies are.

## Bottom line
- When inside and outside are fitted together with no rules, the CMB flattens easily (0.03×).
- Requiring the structure to carry across our boundary pushes the solution entirely outside (0.8–1×).
- Adding realistic web wiring stops the fit reaching 1× (best compromise 17×+).
- Holding the bodies to the observed density-vs-size trend removes the effect completely.
- The arrangement of the fitted network looks like an ordinary web. Its body strengths don't: they're thousands of σ beyond ΛCDM and ~1000× beyond the observed trend.

## Files
- **Fits and drivers:** `jn_core.py` (model and rules, numpy), `jn_fast.py` (compiled rules and Jacobian, checked against jn_core), `run_jn.py` (phases A/B/C), `seeds.py`.
- **Null networks:** `null_nets.py`, output in `null/`.
- **Charts:** `jn_report.py` (best_vs_N_logy.png, results.json), `make_views.py` + `views.py` (side_views_mirror.png, side_views_nomirror.png), `network_compare.py`.
- **Checks:** `rule_breakdown.py` (rule_breakdown.json), `trend_cap_check.py` (trend_cap_check.json).
- **Raw data:** every fit is in `raw/<level>_<m|n>_N<N>.json`.
- **Caveats:**
  - Overlapping voids (allowed up to 50% at the rule levels) double-count the deficit where they overlap.
  - With 8 starts per case, some curves don't fall steadily with N (multi-start noise).
