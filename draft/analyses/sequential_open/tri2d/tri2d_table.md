# 2D triangular crystal, equal numbers of voids and clusters, U = −G Σ MᵢMⱼ/rᵢⱼ (voids negative mass)

Lattice: triangular, periodic box 44x46 sites = 2024 (box 44.0 a x 39.8 a). All-pairs: 2D-periodic Ewald sum of 1/r with neutralizing uniform background (k=0 term); includes self-image terms (constant at fixed composition). MC: nonlocal Kawasaki swaps (random cluster <-> random void), geometric anneal; NN: T 3->0.01 (x scale), 60 steps x 40N attempts; all pairs: T 5->0.01, 40 x 15N attempts; B temperatures scaled by ((1-Mv/Mc)/2)^2=12.23; NN result replaced by the best of 6 longer anneals (T 3->0.05, 300 steps x 400N attempts).
Energy unit G·M_c²/a with M_c = 1.40e50 kg: a = 1.89 Gpc: 2.243e+64 J, a = 2.03 Gpc (2D-matched): 2.088e+64 J, a = 171 Mpc: 2.479e+65 J.
Perfect lattice geometry (every configuration, positions fixed): triangle area √3/4 a² = 1.547 Gpc² (a = 1.89 Gpc), 1.784 Gpc² (a = 2.03 Gpc (2D-matched)), 0.01266 Gpc² (a = 171 Mpc); link length a; 6 connections per node.

| case | configuration | like / unlike links per node | U per node, neighbours (G M_c²/a) | U per node, all pairs (G M_c²/a) | neighbours, J per node at a = 1.89 / 2.03 / 0.171 Gpc |
|---|---|---|---|---|---|
| A | random 50/50 | 2.98 / 3.02 | +0.018 | +0.050 | +3.99e+62 / +3.71e+62 / +4.41e+63 |
| A | alternating rows (= best alternation) | 2.00 / 4.00 | +1.000 | +0.702 | +2.24e+64 / +2.09e+64 / +2.48e+65 |
| A | double-row stripes | 4.09 / 1.91 | -1.087 | -0.204 | -2.44e+64 / -2.27e+64 / -2.69e+65 |
| A | small clumps | 4.70 / 1.30 | -1.702 | -2.463 | -3.82e+64 / -3.55e+64 / -4.22e+65 |
| A | big clumps | 5.61 / 0.39 | -2.609 | -11.669 | -5.85e+64 / -5.45e+64 / -6.47e+65 |
| A | phase separation, row interface | 5.83 / 0.17 | -2.826 | -17.571 | -6.34e+64 / -5.90e+64 / -7.01e+65 |
| A | phase separation, column interface | 5.82 / 0.18 | -2.818 | -19.585 | -6.32e+64 / -5.89e+64 / -6.99e+65 |
| A | MC anneal, nearest-neighbour links [A] | 5.83 / 0.17 | -2.826 | -17.571 | -6.34e+64 / -5.90e+64 / -7.01e+65 |
| A | MC anneal, all pairs (Ewald) [A] | 5.82 / 0.18 | -2.818 | -19.585 | -6.32e+64 / -5.89e+64 / -6.99e+65 |
| B | random 50/50 | 2.98 / 3.02 | -18.479 | +13.734 | -4.14e+65 / -3.86e+65 / -4.58e+66 |
| B | alternating rows (= best alternation) | 2.00 / 4.00 | -6.471 | +21.714 | -1.45e+65 / -1.35e+65 / -1.60e+66 |
| B | double-row stripes | 4.09 / 1.91 | -29.708 | +9.031 | -6.66e+65 / -6.20e+65 / -7.36e+66 |
| B | small clumps | 4.70 / 1.30 | -39.498 | -16.986 | -8.86e+65 / -8.25e+65 / -9.79e+66 |
| B | big clumps | 5.61 / 0.39 | -50.588 | -129.527 | -1.13e+66 / -1.06e+66 / -1.25e+67 |
| B | phase separation, row interface | 5.83 / 0.17 | -53.245 | -201.671 | -1.19e+66 / -1.11e+66 / -1.32e+67 |
| B | phase separation, column interface | 5.82 / 0.18 | -53.149 | -226.296 | -1.19e+66 / -1.11e+66 / -1.32e+67 |
| B | MC anneal, nearest-neighbour links [B] | 5.83 / 0.17 | -53.245 | -201.671 | -1.19e+66 / -1.11e+66 / -1.32e+67 |
| B | MC anneal, all pairs (Ewald) [B] | 5.82 / 0.18 | -53.149 | -226.296 | -1.19e+66 / -1.11e+66 / -1.32e+67 |

All-pairs size scaling (case A, per node, G M²/a): random: N=528: +0.09, N=2024: +0.11, N=4620: -0.02; alternating rows: N=528: +0.70, N=2024: +0.70, N=4620: +0.70; phase separation (best interface): N=528: -8.78, N=2024: -19.58, N=4620: -30.40

## Compensated layout (voids of many sizes/contrasts; clusters hold each void's deficit in a wall at R–1.3R)

Nodes 1868 (397 voids, 1471 clusters, 3.71 clusters per void), 85 Gpc periodic 2D box, equivalent lattice spacing a = 2.11 Gpc.
Relation: N = 0.92 + M_def / 2.10e+50 kg (effective cluster mass close to the expected mean 2.21e50 kg; the +0.9 is the partly-filled last cluster). Log-log slope 0.67. All voids fully compensated.

| layout | U/node all pairs (J) | U/node all pairs (G M_c²/a) | U/node Delaunay links (J) | median triangle ÷ lattice | triangle-area cv | median link ÷ a | mean degree (void / cluster) | links c–c / c–v / v–v | clusters at r/R (mean) |
|---|---|---|---|---|---|---|---|---|---|
| placed (random angles) | +1.19e+65 | +5.95 | +9.75e+64 | 0.74 | 0.87 | 1.04 | 6.0 (6.5 / 5.9) | 3166 / 2304 / 134 | 1.15 |
| annealed (constrained ground state) | -1.72e+65 | -8.55 | -1.84e+65 | 0.18 | 1.98 | 0.51 | 6.0 (8.2 / 5.4) | 3012 / 1936 / 656 | 1.21 |

Compensated all-pairs energy uses the minimum-image 1/r sum (each void + wall is mass-neutral, so the sum converges); masses: ρ = 2.886e-27 kg/m³ from the web's median void, clusters' excess mass δ_c ρ (4/3)π r_c³.
