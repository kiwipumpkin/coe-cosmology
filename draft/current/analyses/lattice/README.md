# Crystal-lattice clusters + voids around the visible sphere (Chris: "I'm envisioning crystal lattices")

The question: if the outside structure is a repeating cluster/void lattice instead of random filler, how far does it bring down the leftover CMB pattern?
- Start: 237.3× the 1e-5 floor (after special relativity plus the central mass).
- Four fitted clumps: 88.8× (compensated, net mass zero).

Results: `summary_table.md` (all tables) and `results.json` (everything, including every spacing).
Figures: `best_by_lattice.png`, `leftover_vs_spacing.png`, `jitter.png`, `example_lattice.png`.
Earlier run folders and the course notes are untouched. The code imports from them read-only.

## Setup (same as the random surround runs unless stated)
- **Region:** 14.3–42.9 Gpc from us in every direction, with CoE distance ≥ D−d = 43.9 Gpc. That covers the shell and everything beyond D+d, out to D+3d. The CoE interior stays empty. D = 58.2 Gpc, d = 14.3 Gpc, CoE at the ptshell axis.
- **Tidal model, residual, floor:** `field3`, `R`, `FLOOR` from `../shell_ensemble`, unchanged. Objects are point masses. Lattice + clumps uses the compensated four-clump baseline from `../shell_ensemble/with_clumps` (88.8×).
- **Lattices** (sign = cluster + / void −). The lattice is anchored at the CoE, and its orientation is measured from the CoE axis.
  - (a) **Simple cubic, rock-salt signs:** every neighbour has the opposite sign. a = neighbour spacing.
  - (b) **Body-centred cubic (CsCl):** corners +, centres −. The amplitude sign is free, so "the reverse" is included.
  - (c) **Face-centred cubic:** alternating +/− layers (CuAu L1₀), 2+ and 2− per cell.
  - (d) **Hexagonal close-packed:** A layers +, B layers −, c = 1.633a.
  - (e) **Rings around the CoE axis:** concentric spheres around the CoE with radial step a. On each sphere there are rings around the CoE–us axis, with objects about a apart along each ring. Signs alternate along the ring, ring to ring, and sphere to sphere. The axis is fixed, so only 3 offsets are free.
- **Density gradient and net zero:** object mass = amp × ρ_exp-lock(r_CoE) × (volume per site), with L = 922.65 Gpc. amp is the fraction of its cell's mass gathered into a cluster or removed to make a void.
  Void deficits are then rescaled so total void deficit = total cluster excess (net mass 0, same rule as the random runs). |amp| is capped so no void goes below δ = −1.
- **Vision-sphere rule (main):** every object lies entirely outside our 14.3 Gpc sphere. Clusters are treated as δ=+1 spheres and voids as δ=−1 spheres holding their |M|. Offending sites are dropped and the fit is redone.
  Without this rule the fit parks objects right on the 14.3 Gpc edge (first test: closest object 14.30 Gpc). The "bare" rule (centre outside only, as in the random bare runs) is reported as a check.
- **Fit:** spacing scanned over 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 25, 28.6 (= 2d), 32, 36 and 40 Gpc.
  At each spacing there is a global search over uniform random rotations × phases (Sobol: 2048 samples; 1024 for a < 10 Gpc; 512 for rings), then a Nelder-Mead polish of the best 3–6.
  The amplitude is solved exactly (least squares, bounded). The sampling was not reduced. The first launch crashed on a bookkeeping bug and was rerun in full.

## Parameter counts (fitted)
- **3D lattices:** 8 = spacing + 3 orientation angles + 3 phase offsets + 1 amplitude. Lattice type is also a discrete choice, from 5.
- **Rings:** 5 = spacing + 3 offsets + amplitude.
- **Four clumps:** 5 = shared distance + 2 masses + 2 angles, with mirror symmetry about the axis.
- **Lattice + clumps:** 8 + 5 = 13.
- **Best random layout:** 1 (strength), but it is the best of thousands of draws.

## Results
Lattice alone (start 237×):

| lattice | best | spacing | params | objects in region | mass per object |
|-|-|-|-|-|-|
| Hexagonal close-packed | **102.6×** | 28.6 Gpc | 8 | 17 | 1.5e52 kg |
| Simple cubic (rock-salt) | 116.0× | 28.6 | 8 | 11 | 3.7e52 |
| Body-centred cubic | 120.8× | 25 | 8 | 30 | 3.6e51 |
| Face-centred cubic | 121.0× | 32 | 8 | 33 | 9.6e51 |
| Rings around CoE axis | 140.8× | 25 | 5 | 22 | 4.7e51 |

- **No lattice alone beats the four clumps (88.8×).** Every type beats the best random layout (153× tuned, 175× as drawn).
- The best lattices need heavy objects: 4e51–4e52 kg each, which is level c or above and up to ~5× a fitted clump.
- At the best spacings only 11–33 objects fall in the region.
- Orientation matters a lot. A typical (median) orientation of the same perfect lattice gives 231–237×, which is almost nothing. The good results come from the fitted orientation and phase.
- The effect is not just the nearest object or two. For hexagonal, the 4 nearest alone give 124× (8 nearest: 129×). For simple cubic, the 4 nearest alone make it worse (437×), so the cancellation comes from the whole set.
- **Bare rule** (centres outside 14.3 Gpc only): face-centred 91.6× and hexagonal 96.9×. These are still not below 88.8×.

Lattice on top of the four compensated clumps (start 88.8×):

| lattice | best | spacing | params (lattice + clumps) | mass per object |
|-|-|-|-|-|
| Body-centred cubic | **39.3×** | 14 Gpc | 8 + 5 | 7.5e50 kg |
| Hexagonal close-packed | 42.3× | 20 | 8 + 5 | 1.5e51 |
| Simple cubic | 45.2× | 14 | 8 + 5 | 8.2e50 |
| Face-centred cubic | 49.1× | 20 | 8 + 5 | 1.4e51 |
| Rings around CoE axis | 64.9× | 14 | 5 + 5 | 5.2e52 |

- **A lattice on top of the clumps improves on 88.8× for every type.** The best is 39.3×, below the best random + clumps layout (61×).
- The objects are lighter here (7.5e50–1.5e51 kg, between levels b and c), and the best spacing is 14–20 Gpc, close to d and the ~16 Gpc clump spacing.
- The rings case is the exception: it needs 5e52 kg objects.

Spacing: there is no single preferred spacing. The curves are jagged.
- **Lattice alone:** the best spacings are 25–32 Gpc, around 2d = 28.6. Hexagonal is also good at 14–16 Gpc (104–105×).
- **With clumps:** the best spacings are 14–20 Gpc.
- **Search noise:** the jaggedness is partly real (which few objects land near us) and partly search noise. Treat individual numbers as about ±10–15×. The search is global but not exhaustive.

Jitter: the best lattice, with each object shifted by a random 3D offset of RMS f × spacing. 1,000 layouts per level, strength kept as fitted.

| | 10% | 25% | 50% |
|-|-|-|-|
| Hexagonal alone (perfect 102.6×): helped / hurt vs 237× | 666 / 334 | 112 / 888 | 51 / 949 |
| median leftover | 215× | 342× | 416× |
| BCC + clumps (perfect 39.3×): helped / hurt vs 88.8× | 881 / 119 | 178 / 822 | 16 / 984 |
| median leftover | 70× | 117× | 151× |

None of the 6,000 jittered layouts beat its perfect lattice. A small shake keeps most of the benefit, and by half a spacing it behaves like the random runs (mostly hurt).
This is the lattice → random bridge: the useful range is ordered structure with modest disorder.

## Caveats
- Point masses and a Newtonian tidal pattern on the in-plane CMB ring, the same as all earlier runs.
- The leftover is the residual RMS. Mode amplitudes per best fit are in `results.json` (`amps_over_floor`, m = 1–6).
- The CoE-anchored lattice origin is arbitrary. The phase offsets cover it.

## Files
- `lattice_core.py`: lattice generators, masses and pattern.
- `run_search.py`: search. `python3 run_search.py alone,with4` runs the buffer rule; `python3 run_search.py alone bare` runs the bare rule.
- `jitter.py`, `report.py`.
- `search_alone_with4_buf.json`, `search_alone_bare.json`, `partial_*.jsonl`, `jitter.json`, `log_*.txt`.
