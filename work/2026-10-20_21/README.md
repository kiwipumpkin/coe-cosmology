# 2026-10-20 / 2026-10-21 — CoE red universe + realistic web: CMB ℓ=2–30 and zero-centre Migkas fits

Scripts, notes (HEURISTICS.md), small JSON summaries (<2 MB) and figures (PNG) from `/workspace/redweb2/` on the analysis box. Paths inside scripts point there. Large arrays (.npy/.npz) and logs are not included.

Reference: ΛCDM χ² = 20.1 vs Planck at ℓ 2–30; flat-sky Migkas A/A_all χ² = 45.1 (71 bins).

## Key results (from the JSON files)
| Run | Result |
|---|---|
| Fine R scan, universe #113, flow build (`realweb/los_flow_R/flowR_bestR_scan.json`) | best R = 4.0 Gpc, E = 0.01305, χ² = 38.9, Pearson r² = 0.0015, log-slope 0.018 (R=3.5: χ² 50.9; R=4.5: χ² 64.1) |
| Flow build, 61 guesses, R = 4 Gpc growth version (`flowR_summary.json`) | median E* = 0.0137, 61/61 guesses, median χ² = 49.3, median r² = 0.006 |
| Monte Carlo, 1000 universes (`realweb/mc/samples/bracket_guesses.json`) | bracket-rule passes: R=4 Gpc 1000, R=8 Gpc 999, R=14 Gpc 849; all three 848; R=0.5/0.9/2 Gpc: 0 |
| Monte Carlo ΛCDM realisations (`realweb/mc/mc_results.json`) | median χ² 42.1, Pearson r² median 0.037 |
| Zero-centre Migkas fit (`zerocore/migkas_zero_fits.json`) | best D = 2 Gpc: χ² 13.19, Pearson r² 0.646; best D ≤ R_cut (1.92): D = 1 Gpc, χ² 13.37, r² 0.668; red reference χ² 34.3, r² 0.362 |
| Zero-centre scans (`zerocore/scan.json`, `scanB.json`), best χ² per D | D=0.5: 15.2; 1: 13.4; 2: 13.2; 5: 14.5; 20: 21.4; 42: 20.8 (r² 0.534); 60: 20.9 |

## Scripts
### Top level (thin-shell / jump runs, network plots, figure stacking)
- `run.py`, `run3.py`–`run11.py`: full-sky ring-point generalisation of field3, parameter variants.
- `thin_run.py`, `thin_density_run.py`, `thin4_run.py`, `thin4_add.py`: same, thin-shell and density variants; `thin_*plot*.py`, `thin4_lin*.py`, `*_bodies.py`: their plots.
- `jump_run.py`, `jump_step(s).py`, `jump_an.py`, `jump_plot.py`: per-body χ² change (add/remove body) analysis.
- `bodies.py`, `bodies2.py`, `bodies3.py`: side view of one generated network (slab, voids, clusters).
- `stack*.py`, `lin1.py`, `lin3.py`, `mk4.py`, `panel3.py`, `plot4.py`: level-set ladder stacking and linear Planck panels.
### densweb/
- `gen.py`: density-following random web (seed-2000 style); `webv.py`: web peculiar-velocity field; `cmb.py`, `mfit.py`, `check.py`: sky spectrum, Migkas plane fit, checks.
### quadrant/
- `quad.py`, `wedge.py`: one quadrant of the red CoE universe, 0..D+R from the CoE.
### realweb/ (realistic web: Press–Schechter halos + compensated voids)
- `build.py`: builds the realistic web groups; `checks.py`: homogeneity checks; `test.py`: FFT vs direct-sum validation; `cosmo.py`: cosmology helpers; `analyze.py`, `interim/analyze_interim.py`, `lin_planck.py`: plots.
### realweb/mc/ (Monte Carlo, 1000 universes × 6 R)
- `mc.py`: MC of full universes at strength 1; `mc_an*.py`: summaries; `mc_bracket.py`, `mc_samples.py`: bracket-rule guesses; `mc_ladder*.py`: level-set ladders (log/linear/fine/4-panel); `mc_mean*.py`, `mc_Ravg.py`: mean spectra of passing universes (per R and averaged over R); `mc_domain_plot.py`, `mc_overlay.py`, `mc_br_plot.py`: shaded min–max domain and overlap overlay; `mc_spec.py`, `mc_all.py`.
### realweb/los/ (static light-path build)
- `los.py`: line-of-sight CMB, components stored separately; `an.py`: summary. Notes in `HEURISTICS.md`.
### realweb/los_flow/ (web riding the red LTB flow)
- `ltb.py`: red-universe LTB expansion and photon paths; `flow.py`: light-path CMB with flow; `deform.py`, `lvl.py`, `an.py`: deformation, levels, component rms. Notes in `HEURISTICS.md`.
### realweb/los_flow_R/ (flow build with observation distance R as input)
- `flowR.py`: main build; `rscan.py`, `rbest.py`: fine R scan / best R; `fine.py`, `hi.py`: fine and higher-strength ladders; `dec.py`, `rot.py`: ℓ=5 decomposition and smooth-red axis tests; `an.py`, `tab.py`, `fix.py`: strength charts, tables, linear-panel fix (`lin4_fix.txt`). Notes in `HEURISTICS.md`.
### zerocore/ (zero-centre density law ρ ∝ (x/xp)^a exp(−(x/xp)^b))
- `zc.py`: build + Migkas fit; `scan.py`, `scanB.py`: parameter scans; `famZ0.py`: CMB trace for the CatWISE-locked family; `refA.py`, `overlay.py`: refit and overlay plots.

## Not included (size)
.npy/.npz arrays (204 .npz + 30 .npy, e.g. `realweb/mc/mc_sky.npy` 295 MB, `mc_sky_snap.npy` 217 MB, `realweb/bodies_g0_coarse.npz` 156 MB, `realweb/groups.npz` 118 MB, `jump_data.npz` 69 MB, `los_flow_R/u113_Rscan.npz` 64 MB, per-guess `*/comp/u*.npz`), 22 `.log` and run `log*.txt` files, `__pycache__`.
