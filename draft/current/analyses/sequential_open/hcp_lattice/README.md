# HCP lattice analogue (11 Oct 2026)

Scripts (run from sequential_open/, in this order): `hcp_lattice.py` (lattice, holes, neighbour shells, ladder, compensated web s3) ->
`hcp_seeds.py` (16 random void fills of the perfect lattices) -> `hcp_ewald.py` (periodic all-pairs energies) -> `hcp_figs.py` (figures).

Setup: HCP, ABAB stacking, nearest-neighbour spacing d = 171/sqrt2 = 120.9 Mpc (same as the FCC version, cube 171 Mpc), c = sqrt(8/3) d = 197.5 Mpc.
Per cluster: 2 tetrahedral holes (4 clusters at 0.612 d = 74 Mpc) and 1 octahedral hole (6 clusters at 0.707 d = 85 Mpc), the same as FCC.
Masses: cluster 1.40e50 kg, void -8.39e50 kg. 'Web mix' = voids in 1/12 of tetra holes (or 1/6 of octa holes) = 1 void per 6 clusters, which is almost exactly compensated (6 x 1.40 - 8.39 = +0.01e50 kg).
Ladder recipe identical to fcc_ladder.py (14-cell box, 2-cell margin; polycrystal 18-cell box, 10 grains).

## Neighbour shells (count at distance / d)
| | 1.000 | 1.414 | 1.633 | 1.732 | 1.915 | 2.000 | 2.236 |
|---|---|---|---|---|---|---|---|
| FCC | 12 | 6 | - | 24 | - | 12 | 24 |
| HCP | 12 | 6 | 2 | 18 | 12 | 6 | 12 |
First AND second shells identical; they differ from the third shell on (HCP's 2 clusters straight above/below at 1.633 d).
Holes: HCP tetra holes come in face-sharing pairs only 0.408 d (49 Mpc) apart (FCC nearest tetra-tetra 0.707 d); HCP octa holes stack in straight columns 0.816 d apart (FCC nearest octa-octa 1.0 d). A tetra void in HCP also has a 5th cluster at 1.02 d (FCC 1.17 d).

## Delaunay statistics (ladder convention)
| structure | nodes | S_hist | cv vol | mean degree | link/spacing | u = U/N/(G M²/spacing) | links c-c/c-v/v-v |
|---|---|---|---|---|---|---|---|
| HCP perfect | 4462 | 1.367 | 0.411 | 15.21 | 1.188 | +2.57 | 0.74/0.25/0.01 |
| HCP jitter 5% | 4472 | 1.902 | 0.423 | 15.26 | 1.189 | +2.31 | 0.74/0.25/0.01 |
| HCP jitter 15% | 4526 | 2.082 | 0.412 | 14.93 | 1.196 | +1.52 | 0.76/0.23/0.01 |
| HCP jitter 30% | 4584 | 2.290 | 0.488 | 14.97 | 1.221 | +1.35 | 0.77/0.22/0.01 |
| HCP defects 3% | 4468 | 1.636 | 0.438 | 15.16 | 1.194 | +2.14 | 0.75/0.24/0.01 |
| HCP defects 10% | 4444 | 2.058 | 0.529 | 15.21 | 1.208 | +1.54 | 0.75/0.24/0.01 |
| HCP polycrystal | 12735 | 2.268 | 0.536 | 15.29 | 1.227 | +1.60 | 0.75/0.24/0.01 |
| HCP perfect, octa-hole voids | 4484 | 0.670 | 0.279 | 14.20 | 1.184 | +0.56 | 0.81/0.17/0.01 |
| HCP clusters only | 3808 | 0 | 0 | 13.99 | 1.166 | -5.46 | 1/0/0 |
| FCC clusters only | 4000 | 0 | 0 | 14.00 | 1.189 | -5.39 | 1/0/0 |
| FCC perfect (web mix) | 4648 | 1.245 | 0.417 | 15.34 | 1.210 | +1.93 | 0.72/0.26/0.02 |
| random (Poisson) | 4686 | 2.806 | 0.865 | 15.55 | 1.269 | -0.11 | 0.72/0.26/0.02 |
| foam s2 web BEFORE / AFTER | 13516 | 3.064 / 3.061 | 1.31 / 1.27 | 15.06 / 15.08 | 1.259 | +4.93 / +4.72 | 0.67/0.32/0.01 |
| balanced web s3 BEFORE / AFTER | 7941 | 2.713 / 2.726 | 0.81 / 0.82 | 15.36 / 15.37 | 1.277 | +4.05 / +3.95 | 0.64/0.34/0.02 |

16 random void fills of the perfect lattice: HCP u = +2.19 ± 0.22 (s.d.), FCC +2.01 ± 0.24; octa variant HCP +0.78 ± 0.16, FCC +0.69 ± 0.13. Not separable in this convention.
Clusters only: every tetrahedron has the same volume and links are d (86%) or sqrt2 d (14%, the octahedron diagonals) in BOTH lattices; Delaunay cannot tell them apart.

## Periodic all-pairs energy (Ewald, uniform background = compensation picture), u per node in G M_cl²/d
Madelung check (unit charges): FCC 0.895874, HCP 0.895838 (literature 0.895874 / 0.895838).
| composition | FCC | HCP | lower |
|---|---|---|---|
| clusters only | +1.62099 | +1.62093 | HCP by 6e-5 |
| web mix, tetra holes (40 fills) | +2.98 ± 0.10 | +2.49 ± 0.09 | HCP by 0.50 |
| web mix, octa holes (40 fills) | +3.14 ± 0.09 | +3.12 ± 0.08 | tie |
| all tetra, compensated void = -0.5 | +0.686 | +0.648 | HCP |
| all octa, compensated void = -1 (rock salt vs NiAs) | +1.236 | +1.222 | HCP |
| all tetra, real masses | +46.3 | +40.6 | HCP |
| all octa, real masses | +27.6 | +27.2 | HCP |
Gravity is Coulomb with the sign flipped, so u_grav = -E_Coulomb; every one of these is positive (the ordered lattice is an energy maximum under gravity).
