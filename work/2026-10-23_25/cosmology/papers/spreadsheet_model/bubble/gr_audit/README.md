# Independent audit of bubble/gr_ltb (other worker's LTB) vs bubble/universe_thin (Newtonian), 2026-09-26
Files: newton_audit.py/.json (own RK4 replica of the fixed-profile Newtonian model plus ablations), ltb_local.py (own minimal LTB:
uniform bang, R(r,t0)=r, local H_r/H_perp), quad_fit.py/.json (best possible fit of an axisymmetric local shear to Migkas),
scan_mc.py/.json (coarse Mc x D x Os scan), refit.py/.json/.log (Mc solved for the optimal shear at each D, Os).
Checks: newton replica matches P.run af to 8e-4 (free) and 6e-5 (baryons); the LTB replica reproduces GR-A (H_perp 89.86, H_r 30.29, age 9.03 Gyr, 2E=598).
Key results: see the final report (the parent agent's message). The Newtonian orientation comes from the core's point-mass tidal field acting on test
particles that spend 98% (free) / 94% (baryons) of their time in the empty interior of a mass shell that is FIXED in space. Remove the core and it
flips; remove the shell and it stays. In comoving dust (LTB, and identically Newtonian comoving dust) the sign follows mean enclosed density vs local density.
With l0 refit, the local-shear LTB family reaches R2 ~0.80 (with sine) / ~0.53 (geometry only) for either sign; D and Mc are degenerate (Mc ~ 0.03-0.08 Mcrit(D) at Os=0.05).
