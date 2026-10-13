# Compensated random web: move + resize with the entropy rule (10 Oct)
Code: compweb_gen.py (generator), minedit_comp.py (fitter), comp_figs.py (figures). Region and kernel are the same as the minedit_ent runs (ball of 30 Gpc around us including the core cap, point-mass tidal kernel, same residual and floor).

| seed | bodies (voids + clusters) | start | end | bodies edited | units (group / cluster) | stop |
|---|---|---|---|---|---|---|
| 3 (best) | 9,062 (1,645 + 7,417) | 789× | **0.998×** | 1,099 (12.1%) | 324 (191 / 133) | reached 1× |
| 2 | 9,028 (1,656 + 7,372) | 1,077× | 1.18× | 1,927 (21.3%) | 1,046 (483 / 563) | plateau |
| 1 | 9,146 (1,674 + 7,472) | 1,260× | 7.40× | 2,538 (27.7%) | 1,122 | plateau |

Best run, seed 3:
- **Moves and resizes.** 1,067 bodies moved (178 voids, 889 clusters) and 930 resized. Median move is 0.6 Gpc, 90th percentile 2.3 Gpc, largest 6.5 Gpc. The median mass factor is 1.02.
- **Leftover by scale.** m1–12 goes 640× → 0.02×. m13+ goes 461× → 1.00×, so all of what is left is fine-scale.
- **Symmetric share.** 67% → 35%.
- **Entropy, before → after.** No component went down.

| component | before | after |
|---|---|---|
| radial | −0.0456 | −0.0441 |
| ring | −0.0054 | −0.0050 |
| angular | −0.0110 | −0.0107 |
| link_cc | 3.0009 | 3.0037 |
| link_cv | 2.3176 | 2.3449 |
| link_vv | 2.0223 | 2.0666 |
| ncount | 2.4804 | 2.5002 |

- **Ring density.** Bodies centred 0–3 Gpc from the ring: 151 → 160 in 2,530 Gpc³, which is 0.060 → 0.063 per Gpc³ (the web average is 0.080). Bodies whose surface is within 3 Gpc: 217 → 216. No ring formed and nothing piled up.
- **Compensation.** It stayed exact: the largest net mass of any group divided by its deficit is 4e-16.
- **Mesh (slab).** Triangle median 1.45 → 1.38 Gpc² over the whole slab and 1.84 → 1.59 Gpc² in the 0–3 Gpc band. Connections per node stay at median 6.0 (mean 5.93). Link length median 2.09 → 2.05 Gpc.
