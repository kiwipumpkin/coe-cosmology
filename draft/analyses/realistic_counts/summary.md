# How many real clusters and voids would it take to reach 1×?

Chris: "How many of those do you need to place and where? Symmetrically and asymmetrically, inside, outside or both to reach CMB."

## Method
The effect is linear in mass, so each 1× solution (point-mass fits and smooth minimal fields) reduces to two numbers:
- **M+**, the extra mass on the dense side.
- **M−**, the missing mass on the empty side.

Dividing these by the excess or deficit of one realistic object gives the count. Every count means **extra objects beyond normal**: a rearrangement of the universe's normal content.

**Object values used** (literature-approximate; ρ = 2.9e-27 kg/m³ = 4.3e10 M☉/Mpc³):

|object|size, contrast|excess or deficit per object|
|-|-|-|
|rich cluster|1e15 M☉, R ~2 Mpc, δ ~200|+2.0e45 kg|
|supercluster (Shapley-like)|R 50 Mpc, δ +1.2|+5.4e46 kg (2.7e16 M☉ excess, ~5e16 M☉ total)|
|void|R 28.6 Mpc (20 h⁻¹), δ −0.4|−3.3e45 kg (1.7e15 M☉, Clampitt & Jain scale)|
|supervoid (Eridanus-like)|R 200 Mpc, δ −0.2|−5.7e47 kg (2.9e17 M☉)|

**What each region holds anyway** (`../local_web/volumes_counts.md`), as cluster mass and void deficit:
- Visible sphere: 3.0e53 kg and 5.7e53 kg, or 1.5e8 cluster-equivalents and 1.7e8 void-equivalents (foam count: 6.6e8 clusters, 9.8e7 voids).
- Outside, the 14.3–30 Gpc tidal region: 1.1e9 cluster-equivalents and 1.2e9 void-equivalents.

**New runs:**
- **Smooth fields** (`field_multi.py`, `run_fields.py`): minimum-norm smooth fields at 2 Gpc smoothing.
  - Inside: 0.5–13.8 Gpc on a 0.3 Gpc grid.
  - Outside: 14.8–30 Gpc on a 0.6 Gpc grid, CoE interior excluded.
  - Both: the two regions together.
  - Net mass is zero, the target is exactly 1×, and the kernel and residual are the same as before.
- **"No mirror" for fields**, tested two ways:
  - **One-sided:** the field is allowed only on one side of the mirror plane, and nowhere within 10° of it.
  - **Staggered:** the dense side and the empty side sit at clearly different distances, so nothing mirrors anything.
- **"No mirror" for points inside:** the box rule of `../solution_space/asym` (20°, 15%, 25%, 10° off the plane) applied to inside fits (`run_inside_asym.py`). Classes 7, 8 and 10 masses, 6 + 2 and 4 + 4, with 10 starts each.

## Counts per scenario (`objects_needed_by_scenario.png`)
Best 1× case per scenario. Clusters and voids are the extra objects needed. Superclusters and supervoids are the alternative counts if only those are used.

|scenario|best 1× case|rich clusters|or superclusters|voids|or supervoids|share of what's there anyway|peak δ|
|-|-|-|-|-|-|-|-|
|**inside, mirror allowed**|smooth field|**1.8e6**|6.5e4|**1.0e6**|6.1e3|1.2% / 0.6%|±5.5%|
|inside, no mirror|points: best is 6 + 2 at **3.2×** (no 1×). One-sided field: **can't reach 1×**. Staggered field: 1×, but δ reaches **−1.59 (unphysical)**|1.7e7|6.2e5|9.9e6|5.8e4|11% / 6%|δ < −1|
|outside, mirror allowed|points, 6 masses + 5 voids (0.58×)|3.4e7|1.3e6|2.0e7|1.2e5|3% / 1.6%|core δ +100, voids −0.7|
| |smooth field, for comparison|3.9e7|1.4e6|2.3e7|1.4e5|3.5% / 1.9%|+26% / −21%|
|outside, no mirror|points, 7 masses + 5 voids (0.73×)|2.9e7|1.1e6|1.7e7|1.0e5|2.6% / 1.4%|core δ +100, voids −0.7|
| |one-sided field, for comparison|8e10|–|5e10|–|~7,000%|δ ~800 (impossible)|
|both, mirror allowed|the inside field is the best. The joint minimum-norm field uses more mass (3.0e6 / 1.8e6) because it spreads δ² over a larger volume|1.8e6|6.5e4|1.0e6|6.1e3|0.13% / 0.07% of the larger region|±5.5%|
|**both, no mirror**|staggered field: dense side inside near the ring, empty side outside (or the reverse: 3.0e52 kg)|**1.4e7**|5.1e5|**8.1e6**|4.8e4|1.1% / 0.6%|+13% / −17%|

## Where they go (best case, inside with mirror allowed: `map_where_objects_go.png`)
Eight lobes, 10–11 Gpc from us, close to the ring plane, each about 400–700 Gpc³:
- **Over-dense** at ℓ ≈ 15°, 128°, 195° and 263°. Each needs about **+4.3e5 extra rich clusters** (ℓ 195°: 2.9e5), or about 1–2e4 superclusters.
- **Under-dense** at ℓ ≈ 75°, 162°, 228° and 315°. Each needs about **2–2.7e5 extra voids**, or about 1–2e3 supervoids.
- These are the same lobes as the earlier smooth-field result. The empty lobes sit at the outside four-mass longitudes, and the dense lobes sit at the inside point-mass longitudes.

## Can it physically be built?
**The smooth inside field: packing yes, density no (by observed standards).**
- **Dense lobes:** a lobe naturally holds about 7e6 cluster-equivalents. Adding 4e5 is about **+6%** more cluster mass there. That's easy to fit, since clusters fill well under 1% of the volume.
- **Empty lobes:** 2–3e5 extra voids take **4–5% more** of the lobe's volume. That fits too, since voids already fill 50–78%.
- **But the density is far too high.** It implies lobe-wide means of **±1.4–2%** and peaks of **±5.5%**, coherent over 400–700 Gpc³ and arranged in a fixed m = 3–5 pattern.
  - ΛCDM expects rms 1.5e-4 to 8e-4 at 2 Gpc smoothing, so the peaks are about 70–370× typical.
  - Observed large-scale homogeneity is at or below about 1%, so the lobe means are above that and the peaks well above.
- **Distance.** In standard cosmology, 10–11 Gpc comoving is roughly z ≈ 10–30 (approximate), where linear fluctuations were 10–30× smaller still.

**Point-mass fits: not realistically buildable.**
- **Masses:** one fitted mass of 2e52 kg equals 1e7 rich clusters. Packed into a 1 Gpc sphere it means **δ ≈ 56**. In a 2 Gpc sphere, δ ≈ 7. Even in a 3 Gpc sphere it is δ ≈ 2, far beyond any supercluster complex.
- **Voids:** the fitted voids are δ = −0.7 over 3.5–4 Gpc. Built from real voids (δ −0.4) they would need **175% of the volume**; from supervoids (δ −0.2), 350%. The emptiest a region full of real voids can get is about −0.2 to −0.3.
- **Real share of mass.** The counts look like a small share of what's there (1.6–3.5% outside, 10–18% inside), but each fitted object concentrates that share into a single Gpc-scale lump.

**Asymmetric cases.**
- Outside points are buildable on the same terms as the symmetric ones.
- Inside without mirror images is not buildable: δ < −1 is impossible.
- The staggered "both" field is physically allowed (δ between −0.17 and +0.13), but its 13–17% peaks over Gpc scales are even further from observed homogeneity.

## Symmetric vs asymmetric
- **Allowing mirror images gives the fewest objects.** The residual itself is mirror-symmetric about ℓ 15°/195°, so the cheapest solution is symmetric too.
- **Banning mirror copies costs about 8×** in objects (1.8e6 → 1.4e7 clusters).
- **The cheapest asymmetric layout splits the job across the ring:** the dense side just inside it, the empty side just beyond it.
- **Inside alone,** the ban stops point fits at 3.2×, and smooth fields reach 1× only with impossible (δ < −1) densities.
- **The strictest version (everything on one side of the mirror plane)** is impossible inside, and outside needs δ ~ 800. The pattern needs structure on both sides of the CoE axis.
- **Outside points barely care:** the no-mirror 7 + 5 needs 2.9e7 clusters against 3.4e7 for the mirror-allowed 6 + 5, because the asymmetric fit reuses near-mirror positions with different masses.

## Caveats
- **Object values** are typical literature numbers, not fitted. Counts scale inversely with them.
- **Counts are extra objects beyond normal:** a net rearrangement, not total population.
- **Minimum-norm isn't minimum-mass.** The fields minimize ∫δ², not mass, so mass-minimal layouts could use somewhat fewer objects. That's why the joint "both" field used more mass than inside alone.
- **Outside field limits:** limited to 14.8–30 Gpc on a 0.6 Gpc grid. Inside no-mirror points: 10 starts per class only.
- **Missing physics still applies:** time-varying potentials, lensing and the potential difference between ring and us are left out. They matter most for the inside lobes, which sit right on the line of sight.

## Files
- **Charts:** `objects_needed_by_scenario.png`, `map_where_objects_go.png`
- **Data:** `counts.json`, `fields.json`, `lobes_inside.json`, `inside_asym.json`, `field_*.npz`
- **Code:** `field_multi.py`, `run_fields.py`, `run_stagger.py`, `run_inside_stagger.py`, `run_inside_asym.py`, `counts.py`, `lobes.py`, `plots.py`
