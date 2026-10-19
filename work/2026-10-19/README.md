# 2026-10-19 — CoE CMB / anisotropy analysis scripts

Code only (66 Python files). Data outputs (.npz/.npy/.png/PDF) are not included; they live on the analysis box under `cosmology-course/notes/accel/gr/rings_cmb/v11/`. Paths inside scripts point there.

## Ray geometry and emission crossing (top level)
`surface*.py`, `radial*.py`, `emitcross*.py`, `grid*.py`, `dtt*.py`, `sweep_*.py`, `scan.py`, `simple.py`, `dscaled*.py`, `contained.py`, `col2x2.py`, `sim_trace.py`, `zsave.py`
CoE ray diagrams. Light is emitted when expanding matter crosses a recombination radius centred on the CoE. Computes dT/T on the sky with and without mass, and compares simultaneous vs sequential recombination cases. With no mass the sky is flat; mass makes the expansion uneven.

## static/
Static-frame feasibility checks against the Migkas cluster anisotropy (including a Galilean check), plus ray and dT/T plots.

## distributed/
Runs with mass distributed by a density law (inside and outside the visible universe) instead of a point core.

## sky/
Sky-direction calculations and sky map / angle plots.

## pecmotion/
Effect of our 370 km/s motion relative to the CMB frame (the green/orange curves).

## rotfit/ (free/, spread/)
Sky-orientation fits to Migkas: fixed axis (`rotfit/`) and free axis (`free/`), with spread of good fits (`spread/`).
**Correction:** the earlier "roll is null" result came from a run with no radial gradient. Once the radial density gradient breaks axial symmetry, roll about the CoE axis IS a real fit parameter. The "roll is null" claim is withdrawn.

## dipole_amp/
Simultaneity lever: the ~89x dipole (vs Planck) mostly comes from the spread in emission times, and it falls toward the Planck level as recombination is made more simultaneous.

## shrink_rec/ (motion/, floor/)
Recombination radius that shrinks with a critical-density threshold (shrink rate set by the CatWISE-locked outward density slope), giving a dipole of ~1.72x Planck. `motion/` folds in the 370 km/s motion. `floor/` analyses the geometric floor (~1.73x), which comes entirely from mass making the expansion uneven.

## clusters_voids/ (paper_style/, leftover_vs_n/, m345/)
Compensated random web of voids and clusters (Jennings-Li-Hu void sizes, Tinker cluster masses, ~29% void volume), fed through our density-at-recombination rule dT/T = ln(1+delta)/3. `paper_style/` draws the web like the paper's figures (CoE on the left). `leftover_vs_n/` plots leftover anisotropy vs number of bodies placed. `m345/` targets the m = 3,4,5 multipoles after the dipole is absorbed by our motion relative to the CMB frame.
Result so far: the dipole is cheap to cancel (1 void or 2 masses), but cancelling it raises the non-dipole large-scale RMS (5.44x -> ~12.8x). Under our rule a body only matters if it cuts the recombination surface, so bodies wholly inside or outside our sphere are both inert (unlike the paper's gravity kernel, where inside bodies are most effective). The greedy search was time-capped; the "~33+ bodies to flatten l=1-12" figure is an extrapolation, not a result. The m345 run may not have finished when this was committed.
