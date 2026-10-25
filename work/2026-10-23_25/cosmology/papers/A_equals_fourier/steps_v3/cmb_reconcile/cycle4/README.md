# Cycle 4: two ways out for the CoE core + shell model

## Setup (Chris's parameters unchanged unless flagged)
- Chris's Newtonian integrator (prof_core.integrate, same start and Ht iteration). Omega_s = 0.049 shell with R_in = D - 1.028 R and R_out = D + 1.475 R, R = 14.3 Gpc. <H> = 71.6 sets the unit conversion. A = |dV|/separation is compared linearly with Migkas A (never H^2). H0 = 70 is never an input.
- Cluster side: for each profile and D, the interior mass, CoE longitude l0, and sine amplitude and phase are re-solved exactly on the slice's dipole and quadrupole terms (same equations as termslice). R2 is scored against the rebuilt Migkas slice profile. Slices: all clusters (0.31 Gpc) and z>0.07 (0.47 Gpc).
- CMB side: Legendre a_l about the CoE axis on the 14.3 Gpc sphere around us, using today's field.
  - Static redshift: (Phi_emit - Phi_obs)/c^2.
  - Doppler: -n.(V_emit - V_obs)/c, from the integrator's velocities (linear convention, as in cycle 1).
  - Pass means |a_l| <= 1e-5 for l >= 2. The dipole is judged by the observer velocity needed to absorb it (|a1| c).
- D grid: 14.7 (the new lower bound), 16, 20, 22.3, 30, 49.17, 94 Gpc (way out B: 14.7 to 200).

## Two exact statements that frame both ways out
- **Static term.** The whole 14.3 Gpc sphere (D-R <= r <= D+R) lies outside the interior (r < D-R) and inside the uniform shell. So by the shell theorem the interior profile cannot matter. Only M_eff = (interior mass) - (cavity mass 0.049 r_in^3) does, and a2 = -5.84 M_eff/Mcrit. A uniform shell adds only l <= 1.
- **Doppler term.** For any flow centred on the CoE, V = V(r) r-hat, the l >= 2 part of n.V on the sphere vanishes only if V/r is constant over D-R..D+R (cycle 1 no-go). That same condition removes the shear at Earth. So no interior profile can null the Doppler term while keeping the cluster pattern; the numbers below measure how far off each profile is.

## Way out A: change the mass profile inside D - R (solveA.py, solveA.json, cycle4_A_profiles.png)
- **Cluster fit is kept by almost every profile.** After re-solving, R2 is 0.86-0.875 (all clusters) and 0.73-0.744 (z>0.07). The quad axis sits at 8.3 deg (all) and 15.5 deg (z>0.07), the same as the headline. A few solutions land on the l0 + 180 branch (CoE direction reversed, same axis); this is flagged.
- **Only the compact part drives the fit.** In Chris's integrator the flow is sheared at early times, when all the matter sits near the CoE. Extended mass is nearly dead weight for the cluster pattern, but it still counts in M_eff.
- **Headline compact core.**
  - Interior mass 0.0019 Mcrit at every D.
  - Static a2 is -0.011 at D = 14.7 and +0.087 at D = 49; the sign flips because the cavity (0.0169) outweighs the core at D = 49.
  - Doppler a2 is -0.039 at D = 14.7 and -0.38 at D >= 21.
- **(1) Density rising toward the CoE (Chris's intuition), rho ~ r^-2.**
  - Keeps the fit (R2 0.867-0.874).
  - Needs 1.1x to 20x the mass: 0.0021 (D 14.7), 0.033 (D 49), 0.039 (D 94) Mcrit.
  - Static a2 -0.013 / -0.093 / -0.055. Doppler a2 -0.046 / -0.37 / -0.38.
  - The r^-1 version is worse: 0.15 Mcrit and static -0.76 at D = 49.
  - No gain on either CMB term.
- **(2) Thin Sedov-Taylor interior (mass piled just inside D - R).**
  - Keeps the fit only for D <= 30 (R2 0.86-0.87).
  - At D = 49 it cannot reach the quadrupole amplitude (best R2 0.51). At D = 94 it fits only with the CoE axis turned 90 deg, a change of Chris's geometry.
  - Mass 0.004 (D 14.7), 0.074 (D 20), 0.31 (D 30) Mcrit.
  - Static a2 -0.024 / -0.42 / -1.8, the worst of all.
  - Doppler a2 -0.041 / -0.17 / -0.20, the lowest at D 20-30. It is the only genuinely scale-dependent profile: the anisotropy at 14.3 Gpc is 0.38-0.45 of the cluster-scale value (about 1.0 for the others at D >= 20). That is a 2.3x reduction against the roughly 17,000x still needed.
- **(3) Compact remnant inside a thin Sedov interior.**
  - Keeps the fit at every D (R2 0.863-0.875). The remnant carries the same approximately 0.0019-0.0020 Mcrit whatever its share (50%: total 0.0039; 5%: total 0.040).
  - The extended mass is free, so it can be sized to fill the cavity mass (variant 3c). The static term is then exactly zero at every D >= 22.3, with R2 0.869 (all) / 0.740 (z>0.07) at D = 49.
  - Doppler a2 stays at -0.37 to -0.41.
- **Which profile does what.**
  - Best static reduction: (3c), exactly zero.
  - Best Doppler reduction: D near 14.7 with any profile (-0.039 to -0.054, about 4,000x over).
  - Among the three at fixed D: (2) cuts the Doppler term by at most 2.3x, but only where it still fits.
- **Dipole.** Absorbing it needs an observer velocity far above the CMB dipole's 370 km/s:
  - (3c) at D = 49: a1 -0.021, about 6,300 km/s.
  - Compact core at D = 22.3: -0.0095, about 2,900 km/s.
  - Compact core at D = 49: -0.31 (static-dominated), about 93,000 km/s.
  - D = 14.7: -0.22, about 66,000 km/s.
- **Verdict A: no profile passes both.** The static roughly 1% can be removed (3c, or the compact core at D = 22.3). The Doppler term, 0.04-0.4, cannot be removed by any profile inside our radius.

## Way out B: terms that treat the CMB differently from clusters (partB.py, partB.json, cycle4_B_terms.png)
- The headline compact core + shell throughout. The core mass is re-solved at each D so the all-cluster quadrupole stays at 0.084.
- **Static redshift (the "~1%" term).**
  - Core alone: a2 = -0.011 at every D.
  - The shell cavity adds a positive term. The net vanishes, all l together, at D = 22.28 Gpc, where core mass = cavity mass (0.00193 Mcrit). No parameter changes except D, which the clusters do not constrain. The tolerance is +/- 3 Mpc for |a2| < 1e-5.
  - Flag: this D differs from the 49.17 Gpc headline.
- **Adiabatic 1/3 factor.** Same sign, one third the size. It reduces the term but does not cancel it.
- **Doppler from emitter motion.** a2 = -0.39 for D >= 21 and -0.039 at D = 14.7, with a3 positive and a4 negative. Its sign is fixed by the observed A pattern (maximum along the CoE axis). No other term has the same angular shape, so nothing cancels it.
  - Flag: HR/c = 3.4, so the linear v/c is outside its range. A fractional (anisotropic-expansion) estimate gives 0.1-0.8, which does not change the verdict.
- **ISW / Rees-Sciama.** Zero in Chris's static-table picture: the potential is static in the CoE frame, so only static + Doppler terms exist. With a growing-mode history (open, Omega = 0.049): SW factor 9.5 (cycle 2: 9.6), ISW -5.9, net 3.6x the static term linear (cycle 2 exact GR: 6.1x). The ISW has the opposite sign and removes about 60% of the SW, not all.
- **SW + ISW + Doppler together.**
  - Cancelling the Doppler quadrupole needs the static term multiplied by +4.4 at D = 49 (+6.3 at D 40, +3.0 at D 70). For D < 22.3 the factor must be negative, which no history gives.
  - The growing-mode range 3.6-6.1 does cancel l = 2 at D about 40-62 Gpc. It leaves a3 = 4.5e-2 and a4 = 1.6e-2 at D = 49 (4,500x over). The smallest residuals on the grid are about 1e-2 (at D 94-200), still 1,000x over.
  - It is also inconsistent: the growing mode needs Mc/Mcrit about 0.06-0.08 (cycle 2), not 0.002.
- **Lensing.** A uniform source stays uniform (surface brightness is conserved), so lensing cannot remove a quadrupole. It only remaps the existing pattern, at about the 1% (convergence) level.
- **Ionised-shell scattering (kSZ).** A fully ionised shell at today's density over the 14.3 Gpc path gives tau = 0.0073. That shifts the Doppler quadrupole by about +1.4e-3 (0.4% of it), with the same sign structure, so it cannot cancel.
- **Moving the emission sphere closer.** To pass, the Doppler term would need emission within 0.37 Mpc, and the core's static term within 0.42 Gpc. That is the dead late-thermalisation case.
- **Already closed, not redone:** plasma temperature and release-timing patterns (wrong sign), the Tolman / static-potential route (tuning), and the far core in exact GR (GM/c^2 > D).
- **Verdict B.**
  - What can cancel: the static roughly 1% quadrupole, exactly at D = 22.28 +/- 0.003 Gpc (cavity), and l = 2 alone via growing-mode SW + ISW against the Doppler term at D about 40-60.
  - What cannot: the Doppler term of the emitters' anisotropic motion (0.04-0.4), and the l >= 3 residuals of any partial cancellation.
  - Nothing passes the CMB while keeping the cluster fit.

## Caveats
- Everything is Newtonian with today's field on a sphere of comoving-like radius R (ballpark and partial GR). The linear Doppler uses v/c at v > c.
- The cluster-scale pattern is only as good as the 2D longitude-only model (R2 as in termslice). Mass fits use the <H> = 71.6 conversion.
- At D <= 22.3 some matter on the 14.3 Gpc sphere has turned around (minimum radial velocity < 0), so the Doppler numbers there include infalling matter.
- The static-term numbers are exact given M_eff (shell theorem). An earlier numeric version smeared the cavity edge on the radial grid and was replaced by the exact formula (old/ keeps the first run).

## Files
- lib4.py: profile builder, integrator wrapper, CMB diagnostics.
- solveA.py -> solveA.json, solveA.log (way out A).
- partB.py -> partB.json, partB.log (way out B).
- make_figs.py -> cycle4_A_profiles.png, cycle4_B_terms.png.
