# FCC crystal ladder vs foam s2 web (3D Delaunay)

Lattice: FCC clusters, cubic cell 171 Mpc; voids in tetrahedral holes. 'web ratio' = voids in 1/12 of tetra holes (1 void per 6 clusters, like foam s2); 'fluorite' = every tetra hole (2 voids per cluster).
Masses: cluster = web median 1.40e50 kg, void = web median deficit −8.39e50 kg. Web: whole modeled ball, stats inside 28 Gpc (hull slivers dropped). Lattices: 14-cell box (2.39 Gpc), stats inside a 2-cell margin; polycrystal 18-cell box, 10 grains.
S_hist = entropy of the histogram of log10(V/⟨V⟩), 0.1-dex bins (disorder: delta spikes → low). S_vol/lnN = size evenness (equal cells → 1, so the PERFECT lattice scores HIGHEST on this one). cv = std/mean of tet volume. u = U per node ÷ (G M_cl² / own mean spacing).

| structure | nodes | tets | S_hist (disorder) | S_vol | evenness S_vol/lnN | cv vol | mean degree | mean link / spacing | u (dimensionless U/N) | U/N (J) | links c–c / c–v / v–v |
|---|---|---|---|---|---|---|---|---|---|---|---|
| web ratio: perfect | 4648 | 25609 | 1.245 | 10.047 | 0.9898 | 0.417 | 15.34 | 1.210 | +1.93 | 7.96e+65 | 0.72 / 0.26 / 0.02 |
| web ratio: jitter 5% | 4645 | 24960 | 1.855 | 10.015 | 0.9892 | 0.429 | 15.40 | 1.208 | +1.91 | 7.90e+65 | 0.72 / 0.26 / 0.02 |
| web ratio: jitter 15% | 4623 | 24011 | 2.069 | 9.992 | 0.9906 | 0.406 | 14.87 | 1.205 | +1.35 | 5.56e+65 | 0.76 / 0.23 / 0.01 |
| web ratio: jitter 30% | 4665 | 24131 | 2.286 | 9.969 | 0.9879 | 0.487 | 14.87 | 1.225 | +1.09 | 4.49e+65 | 0.77 / 0.22 / 0.01 |
| web ratio: defects 3% | 4664 | 25568 | 1.546 | 10.033 | 0.9886 | 0.458 | 15.32 | 1.215 | +1.64 | 6.77e+65 | 0.73 / 0.25 / 0.02 |
| web ratio: defects 10% | 4650 | 25306 | 1.985 | 9.999 | 0.9862 | 0.521 | 15.28 | 1.228 | +1.56 | 6.46e+65 | 0.73 / 0.25 / 0.02 |
| web ratio: polycrystal | 13076 | 76289 | 2.270 | 11.086 | 0.9861 | 0.549 | 15.42 | 1.227 | +1.48 | 6.18e+65 | 0.73 / 0.25 / 0.02 |
| web ratio: random | 4686 | 25751 | 2.806 | 9.837 | 0.9686 | 0.865 | 15.55 | 1.269 | -0.11 | -4.51e+64 | 0.72 / 0.26 / 0.02 |
| web BEFORE | 13516 | 77711 | 3.064 | 10.696 | 0.9498 | 1.314 | 15.06 | 1.259 | +4.93 | 1.10e+65 | 0.67 / 0.32 / 0.01 |
| web AFTER | 13516 | 77792 | 3.061 | 10.714 | 0.9514 | 1.272 | 15.08 | 1.258 | +4.72 | 1.05e+65 | 0.67 / 0.32 / 0.01 |
| fluorite: perfect | 12000 | 62376 | 0.663 | 10.975 | 0.9940 | 0.385 | 14.66 | 1.260 | -92.73 | -5.26e+67 | 0.00 / 0.38 / 0.62 |
| fluorite: jitter 5% | 12002 | 64221 | 1.479 | 10.986 | 0.9924 | 0.412 | 14.65 | 1.255 | -93.48 | -5.30e+67 | 0.00 / 0.37 / 0.63 |
| fluorite: jitter 15% | 12011 | 66013 | 2.003 | 11.014 | 0.9925 | 0.416 | 14.93 | 1.259 | -89.52 | -5.07e+67 | 0.01 / 0.41 / 0.58 |
| fluorite: jitter 30% | 11978 | 69926 | 2.260 | 11.037 | 0.9894 | 0.489 | 15.51 | 1.270 | -80.02 | -4.53e+67 | 0.06 / 0.45 / 0.49 |
| fluorite: defects 3% | 12075 | 64065 | 0.864 | 10.995 | 0.9934 | 0.405 | 14.79 | 1.258 | -89.42 | -5.08e+67 | 0.02 / 0.39 / 0.59 |
| fluorite: defects 10% | 11896 | 65129 | 1.102 | 11.000 | 0.9924 | 0.436 | 15.05 | 1.257 | -85.17 | -4.81e+67 | 0.05 / 0.41 / 0.53 |
| fluorite: polycrystal | 32842 | 186705 | 2.187 | 11.997 | 0.9884 | 0.541 | 14.68 | 1.261 | -92.30 | -5.23e+67 | 0.04 / 0.40 / 0.56 |
| fluorite: random | 11967 | 69591 | 2.796 | 10.841 | 0.9722 | 0.846 | 15.48 | 1.272 | -88.67 | -5.02e+67 | 0.11 / 0.45 / 0.44 |
| fluorite: FCC clusters only (no voids) | 864 | 3768 | -0.000 | 8.234 | 1.0000 | 0.000 | 14.02 | 1.189 | -5.02 | -1.97e+66 | 1.00 / 0.00 / 0.00 |

- web BEFORE: ladder position on disorder entropy (0 = perfect, 1 = random) = 1.166

- web AFTER: ladder position on disorder entropy (0 = perfect, 1 = random) = 1.163

Notes: FCC = maximal-order (minimum configurational entropy) reference. Web lengths/volumes normalized by their own spacing (1.89 Gpc; web voids are 1.5 Gpc-radius lumps) vs 75–110 Mpc for the lattice: shapes comparable, absolute scales are not.
Link energy depends on which Delaunay links exist, so it is not monotonic along the fluorite ladder (void–void links dominate there). Toy estimate.
