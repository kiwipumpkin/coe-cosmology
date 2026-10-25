# Heuristics: CoE red universe + clump/void web, light-path CMB (low ℓ 2–30)

1. **Observations first.** Planck's low-ℓ D_ℓ is about 200–1500 μK², so ΔT/T rms is about 1e-5. Any mechanism has to land at Φ/c² ~ 3e-5 (SW) or at an equivalent integrated effect. Check this before building anything.
2. **Look for the pattern, not the decimals.** The shape tests are a flat log-slope over ℓ 4–30, no runaway at ℓ 2–3, and the level relative to Planck's median. χ² and r² are labels; the call is Chris's.
3. **Scalings with body mass M, distance d and path length L** (Φ ~ GM/(c² r)):
   - Sachs–Wolfe at the source: ΔT/T = Φ_src/3. It is linear in M and only sees mass near the last-scattering shell (~14 Gpc).
   - Linear ISW in a matter-only (EdS) universe: Φ̇ = 0, so it is **zero** at linear order. Any no-Λ time-varying-potential signal is second order.
   - Rees–Sciama / moving lens: Φ̇ = −v·∇Φ, with v = −(2/3aH)∇Φ, so ΔT/T = (2/c³)∫(2/3aH)|∇Φ|² dχ. This is **quadratic** in M (∝E²) for web×web and **linear** in M for the web×CoE-gradient cross term. It grows with L but is weighted by 1/(aH), so the late (near) path dominates.
   - Lensing: α = −2∫(χs−χ)/χs ∇⊥Φ dχ, roughly 4GM/(c²b) per body. It only **remaps** existing anisotropy, δT ≈ α·∇T, so it cannot create power where the primary sky has none. At low ℓ (pixel ~2°) it's a 10⁻³–10⁻² relative effect.
   - CoE gradient: the large-scale potential gradient is smooth, so its own power sits at ℓ ≤ 3. It matters mainly through the cross term with the web (an amplifier: RS_wc ∝ E·|∇Φ_CoE|).
4. **Order-of-magnitude check before trusting a run.** For each component, print the sky rms: SW ~ E·8e-5, RS_ww ~ E²·3e-8, RS_wc ~ E·5e-7, RS_cc ~ 2.6e-6 (fixed), α_web ~ E·3e-4 rad, α_CoE ~ 3.6e-3 rad. The dominant term sets the required E. Here SW dominates, so required E ≈ (Planck rms)/(SW rms per unit E).
5. **The 1/d² trap.** A single body near a sampling sphere dominates a static read. Path integrals smooth this out, so ladder curves should scale smoothly with E. A jump means a single-body artifact.
6. **Decompose.** Run the components separately (SW only, +RS, +lens) and see which one moves the required strength. A no-Λ RS term enters as E² and so becomes negligible at small E.
7. **Same standard both ways.** Use the same Planck band, χ² definition, r² and bracket rule for ΛCDM, smooth red and the web.

## Added after los_flow (web riding the red LTB flow)
8. **Find out where the light cone sits in Lagrangian (today-label) space before choosing a grid.** In the red LTB (famP poly build, frac 20000, u = 4), the past light cone back to 1+z = 1100 only reaches Lagrangian labels within 0.228 Gpc of us: x = 0.9946–1.0054, φ ≤ 5.4e-3. That is inside the 0.5 Gpc empty zone, so in this picture the web acts only through its long-range potential.
9. **A field seen from inside an empty region is smooth (near-harmonic).** On the sky that means a constant, a gradient and a tidal term, so the power piles up at ℓ = 2–4 and the spectrum falls steeply: log-slope −0.25 (growth) or −0.16 (fixed mass) over ℓ 4–30. Spectral flatness needs structure at a range of distances along the path.
10. **How the potential's amplitude evolves matters more than where the bodies are.** Φ = A(t)Φ0(X), and A(t) depends on the assumption:
    - Fixed-mass bodies: A = 1/a, about 1100 at recombination.
    - Linear growth in the local red background (local Ω ≈ 2e4 at our shell, closed): A = D/a, about 0.20 at recombination, rising to 1 today.
    The ISW-like term ≈ 2ΔA·Φ0, with ΔA = A_today − A_rec. It exceeds Sachs–Wolfe (A_rec·Φ0/3) in both cases.
11. **Differential CoE stretching is small near us.** S(t,x)/(x S(t,1)) differs from 1 by at most 0.35% across x = 0.65–1.35 at recombination, so the anisotropic-stretch term is about 1e-4 of the main ISW-like term. Local q = a_r/a_t differs from 1 by under 3e-6.
12. **Check the direction convention (ψ from the CoE versus from the anti-CoE)** between sm_red (run.py) and the geodesic tracer before combining them.
