# R² maximise note (local-density constraint)

**Date:** 2026-09-18 (America/Chicago)

## Constraint (unchanged)

ρ(D) = Ω_m ρ_crit = **0.3 ρ_crit** (verified = 0.3000).
κ is **derived** from M_enc — **not** forced to 0.15.
For decreasing ρ with this local match, κ > Ω_m/2 = 0.15 always.

## Baseline to beat

Prior balanced (shared |Amp| opp. sign): **R²_A ≈ R²_B ≈ 0.685**, κ≈0.231, D/L≈5.9, f_u≈0.94, AF=`hratio`.

Chris sheet ballpark ~**0.74 / 0.73**.

## Three optima (this scan)

### 1) Best A alone (independent Amp_A, φ_A)

| | |
|---|---|
| R²_A | **0.7413** |
| Profile | `twoscale_fu0.97` / AF=`invsq` |
| D/L, f_u | 9.1000, 0.9719801679556928 |
| α, d/D | 0.464, 0.240 |
| κ / ρ_mean/ρ_c | **0.4438** / 0.8875 |
| Amp_A, φ_A | -0.2200, 226.37° |

### 2) Best B alone (independent Amp_B, φ_B)

| | |
|---|---|
| R²_B | **0.7868** |
| Profile | `plummer_n2` / AF=`hratio` |
| D/L | 3.3201 |
| α, d/D | 0.327, 0.190 |
| κ | **2.7052** (high — Plummer pathology for B-solo) |
| Amp_B, φ_B | -0.0200, 240.33° |

Preferable B-solo with κ<1: R²_B=**0.7831**, κ=0.7012, `exponential`, AF=`invsq`.

### 3a) Best pair — max min(R²_A, R²_B)

| | |
|---|---|
| R²_A / R²_B | **0.7413 / 0.7473** |
| R²_min | **0.7413** (beats 0.685) |
| Mode | `indep_amps_shared_phi` (indep Amps, shared φ) |
| Profile | `twoscale_fu0.97` / AF=`invsq` |
| D/L, f_u | 9.1000, 0.9719801679556928 |
| α, d/D | 0.464, 0.240 |
| κ | **0.4438** (not 0.15) |
| Amp_A, Amp_B, φ | -0.2200, -0.1400, 226.37° |

### 3b) Best pair — max (R²_A+R²_B)/2

| | |
|---|---|
| R²_A / R²_B | **0.7285 / 0.7822** |
| R²_mean | **0.7554** |
| Mode | `indep_amps_shared_phi` |
| Profile | `exponential` / AF=`invsq` |
| D/L | 4.4000 |
| α, d/D | 0.327, 0.240 |
| κ | **0.7012** |
| Amp_A, Amp_B, φ | -0.2200, -0.0900, 218.39° |

## Honest ceiling vs Chris ~0.74

- **Joint min ≈ 0.741** under ρ(D)=0.3ρ_c — matches Chris ballpark numerically, but with **κ≈0.44** and **independent Amp_A ≠ −Amp_B**.
- Strict **shared-|Amp| opposite-sign** (sheet-style) only reaches **R²_min ≈ 0.687** (κ≈0.331) — barely above the 0.685 baseline.
- Near-κ band [0.20, 0.28]: best joint min ≈ **0.734** (still with indep Amps).
- Solo B absolute max (0.787) is a high-κ Plummer; treat as diagnostic, not preferred physics.
- Do **not** fake κ=0.15.

## What moved the needle

Scanning **α**, **d/D**, AF=`invsq`, and especially **independent Amp_A/Amp_B** with shared φ. Profile stayed two-scale / exp under the density constraint.

## Regenerable

```bash
/workspace/cosmology/.venv/bin/python \
  /workspace/cosmology/density_local_match/r2_maximize/run_r2_maximize.py
```
