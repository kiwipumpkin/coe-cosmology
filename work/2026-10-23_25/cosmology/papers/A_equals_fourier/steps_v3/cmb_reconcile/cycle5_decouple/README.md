# Cycle 5: CMB–Hubble decoupling from the same Φ

Chris's ask, computed exactly: from **one** Newtonian Φ of the compact-core + shell model, compare the matter (Migkas) anisotropy to the CMB anisotropy term by term.

## Setup
- Model: gravity only, Λ = 0, Ω_s = 0.049 shell (R_in = D − 1.028 R, R_out = D + 1.475 R), R = 14.3 Gpc fixed (R_see track left to cycle4/dtt_rec_track). Compact point core; Mc re-solved per D so cluster quadrupole a₂ = 0.084 at R_s = 0.31 Gpc. ⟨H⟩ = 71.6. Reuses `../cycle4/lib4.py` / termslice integrator.
- D = 14.7 (shell minimum), **22.28** (static null from cycle4/partB: M_core = M_cavity), 49.
- Angle θ = 0 toward CoE, 180 away.

## Same Φ → two integrals

1. **Matter (clusters).**  
   Δv = −∫ ∇Φ dt (the integrator).  
   ΔH/H = v_grav · n̂ / r ≡ af/⟨af⟩ − 1 at R_s = 0.31 Gpc.

2. **CMB (R = 14.3 Gpc sphere).**  
   ΔT/T = endpoint Φ/c² + (2/c²)∫ Φ̇ dt + Doppler (−n · Δv/c)  
   (+ lensing not needed: remaps only, ~1%).

**Nearly-static Φ limit** (Chris's table / this Newtonian run): substantial ∇Φ, Φ̇ ≈ 0.  
Matter accumulates velocity over cosmic time → finite ΔH/H.  
ISW ∝ Φ̇ ≈ 0. Endpoint redshift ∝ M_eff = M_interior − M_cavity.

## Key bullets (|a₂| unless noted)

### D = 14.7 Gpc (Mc/Mcrit = 0.00187, M_eff = +0.00187)
- |ΔH/H| a₂ = **0.084** (ptp 0.171) — Migkas match by construction
- |ΔT/T|_static (endpoint) a₂ = **0.011**
- |ΔT/T|_ISW (static Φ) a₂ = **0**
- |ΔT/T|_Doppler a₂ = **0.039** (ptp 0.616; turnaround near θ ≲ 25°)
- |ΔT/T|_total (static-Φ) a₂ = **0.050** — **3,900×** over CMB limit 1e−5

### D = 22.28 Gpc — static null (Mc/Mcrit = 0.00193, M_eff ≈ 0)
- |ΔH/H| a₂ = **0.084** (ptp 0.170) — **still matches Migkas**
- |ΔT/T|_static a₂ ≈ **0** (−4.7e−10; core = cavity)
- |ΔT/T|_ISW (static Φ) a₂ = **0**; growing-mode SW/ISW also ∝ M_eff → **0**
- |ΔT/T|_Doppler a₂ = **0.411** (ptp 3.03; turnaround spike)
- |ΔT/T|_total ≈ Doppler = **0.411** — **41,000×** over 1e−5

### D = 49 Gpc (Mc/Mcrit = 0.00196, M_eff = −0.0149)
- |ΔH/H| a₂ = **0.084** (ptp 0.169)
- |ΔT/T|_static a₂ = **0.087**
- |ΔT/T|_ISW (static Φ) a₂ = **0**
- |ΔT/T|_Doppler a₂ = **0.385**
- |ΔT/T|_total (static-Φ) a₂ = **0.298** — **30,000×** over 1e−5

### Growing-mode contrast (Om = 0.049 open; not Chris's static table)
- Multipliers on today's static tidal pattern (cycle4/partB): SW ×9.52, ISW ×(−5.91), net ×3.61.
- At D = 22.28: M_eff ≈ 0 ⇒ SW and ISW still vanish; Doppler alone remains a₂ = −0.411.
- At D = 49: growing net + Doppler leaves |a₂| ≈ 0.072 — still ≫ 1e−5; higher multipoles from Doppler are not cancelled.

## Verdict on decoupling
**Yes — they decouple.** The same Φ that builds Migkas-scale |ΔH/H| a₂ = 0.084 does **not** force a small CMB.

- In the nearly-static limit Chris prefers: ISW is killed, and at D = 22.28 the endpoint is killed too.
- **Doppler from emitter motion survives and dominates** (|a₂| ≈ 0.04–0.41, thousands to tens of thousands times the CMB quadrupole limit).
- Static-null D therefore answers Chris's question directly: **Doppler is still large while ΔH/H still matches Migkas.** Nulling Φ on the last-scattering sphere does not erase the velocity field that Φ has already integrated into the matter.

Physical reason: matter responds to ∫∇Φ dt (cumulative); photons in a static potential see only the endpoint Φ/c² (and Φ̇). Killing the endpoint leaves the accumulated velocity — which is exactly what makes the clusters work.

## Files
- `compute_decouple.py` → `decouple.json`, `decouple.log`
- `plot_decouple.py` → `decouple_vs_theta.png`, `decouple_a2_bars.png`, `decouple_growing_contrast.png`
- Reuses `../cycle4/lib4.py`

## Caveats
- Linear Doppler uses v/c with HR/c ≈ 3.4 (outside linear range); fractional estimate |a₂|/(HR/c) still ≫ 1e−5 (cycle4/dtt).
- D ≤ 22.28 has turnaround (vrmin < 0) near the CoE axis; quadrupole is the robust metric.
- R = 14.3 fixed here; the recombination-track R_see ≈ 4.25 Gpc (cycle4/dtt_rec_track) lowers Doppler |a₂| to ~0.11 but does not change the verdict.
