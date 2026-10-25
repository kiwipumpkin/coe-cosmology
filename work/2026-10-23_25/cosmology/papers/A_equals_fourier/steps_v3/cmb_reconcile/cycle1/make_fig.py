"""One figure: (left) l=2 sky amplitude vs sphere radius s for the static potential piece (today's core field) and the matter-velocity
pattern the A fit needs; (right) best-case |a2| / 1e-5 for each cycle-1 candidate (omb fit unless noted). Reads estimates.json, c5_kinematic_fit.json."""
import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
E = json.load(open("estimates.json")); K = json.load(open("c5_kinematic_fit.json"))
c = 299792.458; beta = E["data"]["beta_needed"]
fits = {"free": (0.002651, 71.79, "C0"), "omb": (0.001742, 71.61, "C1"), "rel": (0.0837, 71.63, "C2")}
s = np.geomspace(0.05, 14.3, 200)
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
for f, (mcc, H, col) in fits.items():
    ax[0].loglog(s, mcc / 2 * (H * s * 1e3 / c)**2 * (70 / H)**2, col, label=f"static potential piece, {f} fit")
ax[0].loglog(s, 2 * beta / 3 * (71.7 * s * 1e3 / c), "k", lw=2, label="matter velocity pattern needed by A (all fits)")
ax[0].axhline(1e-5, color="r", ls="--", label="CMB residual bound 1e-5")
for x, t in ((0.5, "Migkas bulk-flow depth"), (14.3, "emission sphere R")):
    ax[0].axvline(x, color="gray", ls=":"); ax[0].text(x * 1.05, 2e-6, t, rotation=90, fontsize=8, color="gray")
ax[0].set_xlabel("sphere radius around Earth (Gpc)"); ax[0].set_ylabel("|quadrupole| of dT/T or v/c"); ax[0].legend(fontsize=7.5); ax[0].set_ylim(1e-7, 3)
ax[0].set_title("Static redshift piece vs velocity effect", fontsize=10)
o = "omb"; P = E["priority_static_ISW_ratio"][o]; C2 = E["c2_SW_ISW_cancellation"][o]
bars = [("core as fitted, SW+ISW (constant Mc)", abs(P["constant_Mc_Phi_1_over_a"]["net_a2"])),
        ("static potential only (today's field)", abs(P["static_Phi_emit_minus_obs_today_field"]["a2"])),
        ("1 velocity-only, today's pattern at R", abs(E["c1_velocity_only"][o]["doppler_a2_today_pattern"])),
        ("1 velocity-only, ballistic back to emission", abs(E["c1_velocity_only"][o]["doppler_a2_if_ballistic_at_emission"])),
        ("2 Mc grows as a^0.635 (octupole left)", abs(C2["residual_a3"])),
        ("2 potential halves just after emission", abs(C2["step_half_after_emission"]["net_a2"])),
        ("2 flow start time tuned (octupole left)", abs(C2["a3_at_l2_zero"])),
        ("3 fixed-T front, adiabatic one-third", abs(E["c3_recombination_front"]["frozen_a2_with_1_3"][o])),
        ("4 core switched on at a=0.99", abs(E["c4_late_thermalisation_or_core"][o]["late_core_ISW_a2"]["a_f=0.99"])),
        ("7 rotation, transverse Doppler at R", E["c7_rotation"]["transverse_v2_over_c2_at_R"]),
        ("local flow L=1.25 Gpc (growing mode)", E["local_confined_flow"]["Om=0.049"]["rows"]["sig=0.05|L=1.25Gpc"]["late_ISW_a2"]),
        ("local flow L=0.5 Gpc (growing mode)", E["local_confined_flow"]["Om=0.049"]["rows"]["sig=0.05|L=0.5Gpc"]["late_ISW_a2"])]
y = np.arange(len(bars)); v = np.array([b[1] for b in bars]) / 1e-5
ax[1].barh(y, v, color=["C3" if x > 1 else "C2" for x in v]); ax[1].set_xscale("log"); ax[1].axvline(1, color="k", ls="--")
ax[1].set_yticks(y); ax[1].set_yticklabels([b[0] for b in bars], fontsize=8); ax[1].invert_yaxis()
ax[1].set_xlabel("best-case |residual| / 1e-5  (omb fit)")
ax[1].set_title(f"Candidates (5 and 6 pass the CMB but keep only R² ≈ {K['dipole_only_and_harmonic_refs']['dipole']:.2f}–0.34)", fontsize=10)
plt.tight_layout(); plt.savefig("cycle1_summary.png", dpi=130)
