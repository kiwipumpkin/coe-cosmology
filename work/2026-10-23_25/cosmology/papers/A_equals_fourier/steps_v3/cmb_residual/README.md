# CMB residual multipoles in the CoE toy cosmology

Run: `python3 cmb_residual.py` (about 5 s). It writes cmb_residual.json, cmb_residual_legendre.csv, cmb_residual_profiles.csv, cmb_residual.png and cmb_residual_dscan.png.

- Parameters come from ../dscan/dscan.json (headline points: free D=93.72, omb D=49.17, rel/LTB D=32.94 Gpc) and from ../notebook_results.json anchor <H> (71.79 / 71.61 / 71.63). Nothing is refitted.
- theta is the angle from the CoE direction. a_l are Legendre coefficients of dT/T(theta). Residual means l>=2.
- l>=2 only comes from the 1/r potential terms: the core Mc, plus the empty interior of the shell (hole mass = Os rho_c(70) 4/3 pi R_in^3). "core" runs use Mc only. "core+hole" runs use Meff = Mc - Mhole.
- Components:
  - SW = Phi_e/c^2 (full Phi, no 1/3; with a Laplace-exterior field and Phi ~ 1/a the Newtonian-gauge delta is 0, so there is no recombination-timing term).
  - ISW = 2/c^2 int dPhi/dt.
  - Doppler = n.(u_o - u_e)/c, using the linear peculiar velocity in the decaying field. The velocity is zero before a_i, the start of pure V = k r flow.
- The null-ray shape comes from a coasting (Om=0) or decelerated (Friedmann, Os) background. It is normalised so that the comoving radius is R = 14.3 Gpc at a = 1/1101.
- SW_frozen and Doppler_today are charitable lower bounds. They take the potential or velocity pattern at today's amplitude and add no growth back to recombination.
- The linear calculation is formally invalid because GMc/(c^2 D) is about 0.1 to 2.5 today, and about 1100 times larger at recombination. The frozen bounds already exceed 1e-5 by 3 to 5 orders of magnitude.
