# Migkas A/A_all: sigma vs flat null (A/A_all=1)

**Date:** 2026-10-06 (America/Chicago)

**Data:** Migkas Fig. 6 via `bubble_core` (`MA`, `MLO`, `MHI`); N=71 bins; sigma_i=(lo+hi)/2; WLS chi2 = sum w_i (y_i - m_i)^2.

## Null

Flat A/A_all=1: chi2_null=45.1366, nu=71, R2_unw=-0.0038.

## Clarification (CatWISE vs all-free)

The CatWISE-constrained Migkas runs **fixed** L = d / arcsinh(0.031/2) from 2 sinh(d/L)=0.031. Free params were only {D, mc, l0}. The R^2 ≈ 0.71 there is **not** an all-free best fit.

**Re-run with ALL free** {D, L, mc, l0} (pure exp, Ob=0, no sine): R2_unw=0.7589, chi2=8.4753, D=27.95 Gpc, L=6.073 Gpc, mc=0.1895, l0=18.6 deg — confirms prior free-L R^2 ≈ 0.76 is the best pure-exp geometry fit.

**Core density:** when L is free, (rho'/rho)_D = -1/L is set by the Migkas gradient fit; rho(0)=rho(D) exp(D/L) then follows from that L and is **not** independently constrained by the local gradient alone.

## Table vs null

| model | free params | n_free | chi2 | nu | R2 unw. | Delta chi2 | sigma_nested | sqrt(Delta chi2) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| null_flat_A_Aall_eq_1 | — | 0 | 45.137 | 71 | -0.0038 | 0.000 | — | 0.00 |
| ptshell_WLS_geometry_only | D,mcc,kin,kout,l0 | 5 | 18.747 | 66 | 0.5456 | 26.390 | 3.79 | 5.14 |
| catwise_fixed_L_d14.3 | D,mc,l0 | 3 | 10.572 | 68 | 0.7160 | 34.564 | 5.12 | 5.88 |
| catwise_fixed_L_d3 | D,mc,l0 | 3 | 10.489 | 68 | 0.7184 | 34.647 | 5.13 | 5.89 |
| pure_exp_ALL_free_DL_mc_l0 | D,L,mc,l0 | 4 | 8.475 | 67 | 0.7589 | 36.661 | 5.06 | 6.05 |

**Method:** Delta chi2 = chi2_null - chi2; p = P(chi2_{Delta nu} > Delta chi2) with Delta nu = n_free; sigma_nested = norm.isf(p) (one-sided Gaussian equivalent). sqrt(Delta chi2) equals nested sigma only for Delta nu = 1.

## Density plot labels

`catwise-exp-density.png`: checked — legend outside right; rho(0), D±d, and Earth-slope callouts do not overlap.

## Files

| Path | Content |
|---|---|
| `catwise-sigma-null.json` | machine-readable table + free re-run |
| `catwise-sigma-null.md` | this note |
| `catwise-sigma-null.py` | script |
| `catwise-exp-migkas-fit.json` | CatWISE fixed-L fits (source) |
| `exponential-migkas-fit.json` | prior free-L WLS |
