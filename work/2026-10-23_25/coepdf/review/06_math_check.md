# Appendix M — Mathematics Check

Method: each equation re-derived by hand; numbers recomputed in `mathcheck.py` (this folder) and spot-checked against box files. "Verified" = re-derived and/or reproduced to the quoted precision. "Error" = wrong as written. "Unverified" = not checkable from the paper or box in this review. No code was re-run except the small checks in `mathcheck.py`; the run9 rerun result is the earlier rerun's log (`/workspace/coepdf/run_out/run9/run9.log`), not re-executed here.

| # | Item (location) | Status | Note |
|---|---|---|---|
| M1 | Volumes, §2.2 (V_R=268.1, V_{D+R}=1.008e6, 0.0266 %, 1.618e6, 0.0166 %) | Verified | 268.08; 1.00800e6; 0.02660 %; 1.61815e6; 0.01657 %. 72.83 ≠ 1.25×58.2=72.75 (paper already notes). |
| M2 | Eq. 1 ρ_r/ρ_m0 = a_eq/a⁴; example a=1e-4 → 2.94e12 | Verified | Standard. |
| M3 | Eq. 2 closed form of ∫a′da′/√(a_eq+a′) | Verified | Closed form = quadrature to 1e-11 (a=0.5). Limits a^{1/2}, a^{2/3} correct. |
| M4 | t(a_eq)=40.8 kyr, t(1/1090)=297.5 kyr with t(1)=13.8 Gyr | Verified (arithmetic) | 40.79 kyr, 297.47 kyr. |
| M5 | Implied H₀ of Eq. 2 normalisation | **Error (inconsistency)** | Forcing t₀=13.8 Gyr with no Λ gives H* ≈ 47.2 km/s/Mpc. The CatWISE lock uses H₀=70, Ω_m=0.3 (EdS age at H₀=70 is 9.31 Gyr). The thermal history and the density normalisation use incompatible H₀. a_eq=1/3400 is itself a ΛCDM value (depends on Ω_m h² including CDM). |
| M6 | Recombination radius R = 4 (8, 14) Gpc | **Error / unsupported** | R is a free input, not derived. Matter-only comoving distance to z=1090: 8.31 Gpc (H₀=70) or 12.37 Gpc (H₀=47.2). χ*=4 Gpc would need H₀≈145. The final CMB fit (R≈4, R_c=3.57 Gpc) is not tied to the model's own kinematics. |
| M7 | Volume growth n/t; 2/13.8 = 0.145 Gyr⁻¹ | Verified | |
| M8 | CoE kinematics: r(s,n̂)=√(D²+s²+2Ds cosψ), first-order dipole in cosψ | Verified | Correct; note expansion velocity at D=58.2 Gpc with H₀=70 is ≈13.6c — Newtonian kinematics invalid there (paper concedes "superluminal"). |
| M9 | Eq. 3 A/A_all functional | Unverified | f and p(s) never specified (paper self-flags). The forward map cannot be reproduced from the text. |
| M10 | Eq. 4 δ=0.031 = 2×CatWISE dipole | Arithmetic verified; **physics error** | 2×0.01554 = 0.03108 (rounded to 0.031: L becomes 920.3 vs 922.65 Gpc, negligible). But a number-count dipole is mostly kinematic (∝[2+x(1+α)]β) plus a clustering term integrated over the quasar redshift range (z≈1–2, i.e. few Gpc), not a density difference at D±14.3 Gpc. Equating it to [ρ(D−d)−ρ(D+d)]/ρ(D) is an unjustified identification. |
| M11 | Eq. 5 L=922.65 Gpc; ρ₀=0.3195; ρ′/ρ=−0.001084; d=3 → 193.56 | Verified | 922.650; 0.31953; −0.0010838; 193.556. Matches catwise-density-levelset-fit.json. |
| M12 | Law table (power, cored, Gauss, soft, stretched, η-edge): δ=0.031 and ρ₀/ρ_c | Verified | All reproduce δ=0.0310 at the quoted L (to the 4-5 s.f. given); ρ₀ values 0.3201, 0.3098, 0.3096, 0.3784, 0.3847, 1.0319 reproduced. |
| M13 | "all laws ρ′/ρ ≈ −0.00106 to −0.00109 Gpc⁻¹" | Verified | At D: −0.00108 (exp, power), −0.00109 (cored, Gauss), −0.00107 (soft, stretched), −0.00106 (η-edge). |
| M14 | Null χ² = 45.14 over 71 bins | Verified | 45.137 using mean of asymmetric errors (fig6_digitized.csv). With low/high/side errors: 60.5 / 40.7 / 43.4 — the choice matters and is not stated. Note Eq. 9's side rule is not what produced 45.14. |
| M15 | χ²/dof and significance of 45.14 → 7.67–11.6 | **Error in interpretation** | Treating bins as independent: P(χ²≥45.14 \| 71) = 0.993, i.e. the null itself is "too good". Best fits give χ²/dof ≈ 0.11 (7.67/68), 0.16 (10.79/67), 0.17 (11.6/69); P(χ²≤7.67\|69) ≈ 2e-21. Bins are 90° sliding cones at 5° steps — each cluster enters ~18 bins; effective independent bins ~4. χ² values cannot be read as fit quality or significance. No significance is quoted in the paper, but "χ² falls from 45 to 10–12" is presented as the strongest evidence. |
| M16 | A_F baseline (AF_smooth) | Verified (it is the axisymmetric CoE model, not a Migkas fit) | 3-harmonic Fourier series reproduces AF_smooth to 1.5e-5 rms; a 3-harmonic fit to Migkas A differs from AF by 0.055 rms. **But**: data vs AF χ² = 33.2 (71 bins, mean errors); a free 3-harmonic (7-parameter) Fourier fit reaches χ² = 3.17. Generic low-order harmonics beat every CoE law (7.67–11.6). |
| M17 | m(r) modulation (Fig. d5): slope −9.8e-4 at r=42, m(56.3)=0.987, R_b=2.011 | Verified | c₀/42 = −9.76e-4; 0.9874; 2.0107. **Inconsistency**: web uses D=42 Gpc, cluster fit D=58.2 (refits 35–60 Gpc). |
| M18 | Eq. 6–7 harmonics, D_ℓ; example 1.30e4 μK²; θ≈180°/ℓ → 0.82° at ℓ=220 | Verified | 13005 μK². |
| M19 | Eq. 8 cosmic variance values | Verified | 0.632, 0.309, 0.181, 0.0673, 0.0316. |
| M20 | Eq. 9 asymmetric-error χ² | Verified as stated (convention ok) | Side-switching makes χ² non-smooth in parameters; fine as comparison only (paper says so). |
| M21 | Eq. 10 Pearson r² | Verified | Correct; r² comparison 0.976 vs 0.983 is uninformative (paper acknowledges). |
| M22 | Eq. 11 Press–Schechter, bias b=1+(ν²−1)/δ_c | Verified (standard) | But uses ΛCDM σ₈=0.8, Ω_m=0.31 with Eisenstein–Hu CDM transfer: dark matter enters the "no dark matter" model. |
| M23 | Void deficit R_v=30 Mpc → 4.0e15 M☉ | Verified | 4.05e15. |
| M24 | Eq. 13 tidal signal s_R and K_G scaling | Unverified | K_G=7.7e51/0.0039 physical origin unknown (self-flagged). s_R is a radial *acceleration* difference, not a temperature; the map is rescaled by ε fitted per R, so its amplitude is free. |
| M25 | Eq. 14 ladder values for universe #3 (ε, χ²) | Verified vs guesses.json | 29.0 %/441.6, 25.0/411.5, 25.0/177.2, 21.6/69.6, 52.5/85.0, 45.3/47.2. **Confirmed inconsistency**: R=8 lowl_over=2.994 > 2.5 yet "PASS"; also universe #0 R=8 (2.91, PASS). R=14 for #3 is FAIL (slope 0.060>0.05) — paper does not state this. |
| M26 | run9 shear step (272,175 bodies; pattern/noise 149.9) | **Not reproduced** | Earlier rerun log: 272,638 bodies, 48.9 vs repo 149.9. Cause unknown. Shear step is not reported in the paper at all, yet is the first fit step of the method. |
| M27 | Eq. 15 Poisson; 4πGρ̄ = 2.278e9; Φ/c² = 2.5e-6 | Verified | 2.2777e9; 2.534e-6. |
| M28 | SW Φ/3c², ISW (2/c²)∫Φ̇dt | Verified (standard linear GR, but inserted into a Newtonian superluminal flow) | |
| M29 | ⟨(k̂·û)²⟩=1/3 | Verified | |
| M30 | Eq. 17 exact C_ℓ = (2/π)∫k²dk P_Φ\|Δ_ℓ\|² | Verified | Prefactor re-derived from the stated Fourier convention. |
| M31 | Eq. 18 Limber | Verified (form) | Standard; W′ term dropped in Limber part (not stated). |
| M32 | R_c/w = 10.6; source 1.28w above centre | Verified | 10.645; 1.277. |
| M33 | Eq. 19 example P(k=100) = 2.2e-9 Gpc³ | **Error (minor)** | With the damping factor stated (e^{−0.089}) it is 2.06e-9; 2.2(5)e-9 is the value *without* damping. |
| M34 | Eq. 20 Δℓ = 302.6, 282; χ_eff 3.8; 57.1/114.2/199.9 Mpc | Verified | 302.57, 282.27, 3.796, 57.1/114.2/199.9. |
| M35 | Residual mapping r₂=41.7 Mpc, k_B=251 Gpc⁻¹ | Verified | 41.75, 251.3. |
| M36 | Fit table (3787/0.976/137; 4131/0.966/49.1) and parameters | Verified vs acoustic_half_fixed.json | All match. Note: 2499 unbinned points used (npts), not the binned ℓ≥30 file the Data section implies. |
| M37 | χ² 3787 vs ΛCDM 2558 | Verified; interpretation added | Δχ²=1229 with 12 vs ~6 parameters; χ²/dof 1.52 vs 1.03 — CoE decisively disfavoured (paper says "short of ΛCDM", understates). The 3787 is the **ensemble mean** spectrum; the seed-42 one-sky realisation gives 4601 (2–2500) and 148 (2–30), but Fig. 18 caption calls the red curve a one-sky realisation with χ²=3787. |
| M38 | Planck peak positions | **Error (internal inconsistency)** | §Acoustic fit: Planck 224, 535, 821, 1122, 1412; Table 2: 228, 544, 823, 1121, 1417. Same data, two lists. |
| M39 | "half" variant parameters | Verified; flag | s_w=0.0500 sits on its lower bound (boundary solution) — not disclosed. |
| M40 | Δχ²=2.3, 6.2, 11.8 for 2 parameters | Verified | 2.296, 6.180, 11.829. |
| M41 | Harmonic decomposition A: 25.54→15.30→12.06→11.46; B 13.68→8.67→10.04→10.24; phases C₂=109.0° ↔ peaks 54.5/234.5, C₃=248.2° ↔ 82.7° | Phases verified (109/2=54.5, 248.2/3=82.7) ; χ² values unverified | Scripts dumpM.py/recon_n3.py/cum2.py not re-run. n=1 peak at "ℓ=102.8°" while C₁ phase 102.8° — consistent. Notation clash: ℓ used for Galactic longitude and multipole. |
| M42 | 11/12 sign agreement | Arithmetic: binomial P(≥11/12 \| p=0.5) = 0.0032 — **but not a valid significance** | Spots are not independent (antipodal pairs: dipole hot/cold, bulk flow toward/away), orientation (15° lift, 45° rotation) and ±25° smoothing chosen after the fact, A/B assignment is the author's, n=1–3 removal chosen. No null simulation. |
| M43 | CMB isotropy: "95.9 % cancellation" of a several-percent geometric multipole | **Error in sufficiency** | Remaining residual = 4.1 % of 1–5 % ≈ 4e-4–2e-3, i.e. 40–200× the observed ~1e-5 low-ℓ anisotropy. Literature bounds (Alnes & Amarzguioui 2006, PRD 74 103520; Foreman, Moss, Zibin & Scott 2010, PRD 82 103532) put an observer in a Gpc-scale LTB structure within ~tens of Mpc of the centre (≲ ~1 % of the void radius) from the CMB dipole/quadrupole — here the observer is 58 Gpc off-centre. [exact literature numbers quoted from memory — verify before citing] |
| M44 | Power-spectrum stage omits tidal shear | Verified (code) | fwd.py/acoustic_model.py: SW+ISW+moving lens only; moving lens ≥99.99 %. Lensing/shear not in the TT fit; the shear (run9) step is separate and not reproduced (M26). |
| M45 | H/H_ref = 1/√(A/A_all) | Verified | Follows from A ∝ H₀⁻² at fixed flux. |
