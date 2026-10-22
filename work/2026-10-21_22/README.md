# 2026-10-21/22: stochastic field, ℓ≈220 hump, acoustic fit, residual harmonics

- `stochfield/`: Stage A Gaussian stochastic field fit to ℓ 2–30 (χ² 17.6 vs ΛCDM 20.1, Planck errors); axial m=0 analysis; combined low-ℓ + hump spectra (combo*.png).
- `rings220/`: HEALPix rebuild of the ring-hump field normalised to D_ℓ peak 5000 μK².
- Acoustic fit `acoustic_*.py`: P_δ = P0 kⁿ[1+A cos(k r_s+φ)]e^{-(k/k_D)²}; r_s=81.9 Mpc, shell 3.944 Gpc. χ² ℓ 2–2500: 11,772 (ΛCDM 2,558); Pearson r² 0.871 (ΛCDM 0.983). Ladders R=2,4,8,14 and 4,5,6 Gpc.
- Residual (Planck−CoE)/σ in `resid_free_W.npy`: damped sine Pearson r² 0.80 (Δℓ 593, L≈1019).
  Damped harmonic series 1–9 reaches Pearson r² 0.837 (fundamental Δℓ≈593, 34.6σ, L≈989); undamped plateaus at 0.551. Figures `resid_fourier*.png`, `resid_damped_1to*.png`, `resid_undamped_1to9.png`.
