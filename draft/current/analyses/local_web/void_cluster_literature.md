# How voids and clusters are connected in the literature

## Short answer

In the standard picture, voids and clusters are not two separate groups of objects. They are the two ends of one structure, the **cosmic web**.
Voids are the empty cells, and clusters are the knots where the cell walls meet. Matter that leaves the voids is what builds the walls,
then the filaments, then the clusters. So the mass missing from voids and the extra mass in clusters are the same matter, moved.
The pattern also grows step by step: today's web is the early, faint pattern made sharper. Each stage depends on the one before it.

## The main links

1. **Same foam, opposite ends.** Space divides into cells (voids). Cells share walls, walls meet along filaments, and filaments meet at
   points. Rich clusters sit at those meeting points (Icke & van de Weygaert 1987; van de Weygaert & Icke 1989; Bond, Kofman & Pogosyan 1996).
2. **The mass balances.** Each void is a dip in density, surrounded by a ring of extra mass. Averaged over a large enough region,
   the dip and the ring cancel. That is why the books balance to zero in our packing (Hamaus, Sutter & Wandelt 2014; Clampitt & Jain 2015).
3. **Voids feed clusters over time.** Voids expand faster than the average, pushing matter outward. It collects into walls, slides
   along filaments, and piles up at the clusters (Zel'dovich 1970; Icke 1984; Gurbatov, Saichev & Shandarin 1989; Cautun et al. 2014).
   Whether a void survives depends on its surroundings: voids next to a collapsing region get crushed (Sheth & van de Weygaert 2004).
4. **We see it in the sky and in motions.** Stacked voids show a rim of galaxies around them (S-type "shell" voids vs R-type
   "rising" voids; Ceccarelli et al. 2013). Locally, our motion is a pull toward the Shapley concentration plus a push away from a big void,
   the "Dipole Repeller", and the two contribute about equally (Hoffman, Pomarède, Tully & Courtois 2017).
5. **Regular spacing (contested).** Einasto et al. (1997) reported superclusters and voids in a rough 3-D grid with ~120 Mpc spacing.
   This is the clearest "voids and clusters as a pattern" claim, but it has not become the standard view.

## Key citations (all checked)

| Paper | What it shows about the void–cluster link |
|-|-|
| Icke & van de Weygaert 1987, A&A 184, 16 ("Fragmenting the universe I") and **van de Weygaert & Icke 1989, A&A 213, 1** ("II: Voronoi vertices as Abell clusters") | The Voronoi foam picture. Voids are the cells, and the vertices match rich (Abell) clusters in both clustering strength and slope. This is the model our foam runs use. |
| **Bond, Kofman & Pogosyan 1996, Nature 380, 603** ("How filaments of galaxies are woven into the cosmic web") | Clusters are the nodes, filaments are bridges between them, and voids fill the space in between. The web is already present, faintly, in the early density pattern. Later growth only sharpens it. |
| **Zel'dovich 1970, A&A 5, 84**; Gurbatov, Saichev & Shandarin 1989, MNRAS 236, 385 (adhesion model); Icke 1984, MNRAS 206, 1P ("Voids and filaments") | Growth story: matter collapses first into sheets, then filaments, then knots. Underdense regions expand and push their matter into the surrounding walls. |
| **Sheth & van de Weygaert 2004, MNRAS 350, 517** ("A hierarchy of voids: much ado about nothing") | Voids form and merge the way clusters do, but in mirror image. Small voids inside a collapsing overdense region get squeezed out ("void-in-cloud"). Typical size today is ~15–20 h⁻¹ Mpc. |
| **Cautun et al. 2014, MNRAS 441, 2923** ("Evolution of the cosmic web") | In simulations, matter flows from voids into walls, then along filaments into cluster nodes. Today voids fill ~77% of the volume but hold ~15% of the mass. Nodes fill <0.1% of the volume but hold ~11% of the mass. Filaments hold ~50% of the mass in ~6% of the volume. |
| **Hamaus, Sutter & Wandelt 2014, PRL 112, 251302** ("Universal density profile for cosmic voids") | Every void has an empty centre and an overdense "compensation wall" just outside its edge. Small voids are overcompensated and large voids are undercompensated. |
| **Clampitt & Jain 2015, MNRAS 454, 3357**; Amendola, Frieman & Waga 1999, MNRAS 309, 465 | Voids lens light like negative mass. Clampitt & Jain measured this with SDSS voids. The deficit reaches ~3 void radii before it is fully compensated, and no sharp wall shows up. |
| **Hoffman, Pomarède, Tully & Courtois 2017, Nature Astronomy 1, 0036** ("The dipole repeller") | Observed motions: clusters pull and voids push. The local flow is about half Shapley attractor and half void repeller. |
| Ceccarelli et al. 2013, MNRAS (arXiv:1306.5798, "Clues on void evolution I") | In SDSS, voids come in two kinds: shell-type voids are surrounded by overdensity and are collapsing (void-in-cloud), while rising-type voids are expanding (void-in-void). |
| Einasto et al. 1997, Nature 385, 139 | Claimed a ~120 Mpc quasi-regular supercluster–void network. Contested, so mention only briefly. |

Wider reviews: van de Weygaert & Platen 2011 (IJMP Conf. Ser. 1, 41, "Cosmic voids: structure, dynamics and galaxies");
Libeskind et al. 2018 (MNRAS 473, 1195, "Tracing the cosmic web"), which compares 12 web finders.

## How this sits with our runs

- The idea behind our foam runs (Voronoi cells, clusters at the vertices, net mass zero) is the standard one: Icke/van de Weygaert, then Bond et al., then Cautun et al.
- Our shares differ from Cautun's: voids 78% vol / 23% mass, clusters 29% mass, filaments 40%, walls 8%. Cautun's are 77% / 15%, nodes 11%, filaments 50%, walls 24%.
  The volume share matches, but we put more mass in clusters and less in walls. Libeskind et al. 2018 show that these splits vary a lot from one method to another.
  Our total cluster mass is on the high side of that range.
- The "t+1 depends on t" point: in the literature the web is **not** a random scatter of pieces. Its layout is inherited from the early density pattern
  (Bond et al. 1996) and grows from it (Zel'dovich, adhesion). Our random foam ignored that inheritance, which fits with the finding that a random foam
  does not cancel the leftover the way an arranged pattern does.
