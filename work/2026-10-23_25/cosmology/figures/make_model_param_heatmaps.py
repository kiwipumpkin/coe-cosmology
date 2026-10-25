#!/usr/bin/env python3
"""R^2 heat maps of the PHYSICAL Sedov-Taylor + Gauss enclosed-mass model
against the 71 digitised Migkas et al. (2020) Fig. 6 points A/A_all.

Model = the "ideal drag" first-principles model used for
map_compare/amp_mean_param_noise/ideal_fit_M_d_D_L_migkas.npz (2026-09-19 run):
  * Sedov-Taylor homologous background, R ∝ t^{2/5}, H = 2/(5t), T = 2/(5H).
  * Two-scale density rho = rho_core0 e^{-r/L} + rho_u, f_u = 0.94, local match
    rho(D) = 0.3 rho_crit (density_local_match.params_twoscale_local_match),
    masses scaled by mass_scale  ->  kappa = mass_scale * kappa_base(D/L).
  * Peculiar radial velocity from du/dt + H u = -G M(<r)/r^2 integrated along
    homologous paths from t_i = 1e-4 T to T (explicit Euler, 800 steps),
    tabulated on r in [0.2, 2.0] D (50 radii) and interpolated.
  * H_obs(n) = H + [(u(r)-u(D)) . n]/d, sources at depth d on a sphere about
    the observer; A ∝ H_obs/H, mean-renormalised, Δℓ = 90° boxcar, renormalised.
  * NO Amp×sin template, NO 2ℓ correction.
Two R^2 definitions are produced:
  * PRIMARY ("direct"): R^2 = 1 - SS_res/SS_tot of the mean-renormalised model
    against the data with NO rescaling -- the definition behind the paper's
    physical-model value R^2 = 0.307. Amplitude (hence kappa) matters.
  * "_corr2": squared Pearson correlation (= R^2 after the best linear rescale
    a*model+b), the objective of the original ideal fit. Because u is exactly
    linear in the mass scale, corr^2 is exactly independent of kappa.
The d/D grid stops at 0.80 so that all source radii stay inside the tabulated
u(r) range [0.2, 2.0] D (no clamped extrapolation).
The only change from the original code: the time integration is vectorised over
radii (identical arithmetic), and the axis may be tilted out of the Galactic
plane (mu = cos b0 cos(ℓ-ℓ0) along the plane; the boxcar is applied in ℓ after
evaluating the model along the plane, so b0=0 reproduces the original exactly).
"""
import json
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/workspace/cosmology/density_local_match")
from run_density_local_match import params_twoscale_local_match, m_enc_at_r, H, G, D  # noqa: E402

OUT = "/workspace/cosmology/figures"
CSV = "/workspace/cosmology/migkas/fig6_digitized.csv"
FU, ETA, NRAD, NT = 0.94, 1e-4, 50, 800
T = 2.0 / (5.0 * H)
plt.rcParams.update({"font.size": 11})

dat = np.genfromtxt(CSV, delimiter=",", names=True)
ell_m, A_m = dat["l_deg"].astype(float), dat["A_over_Aall"].astype(float)

# best physical fit (ideal_fit_M_d_D_L_migkas.npz)
best = np.load("/workspace/cosmology/map_compare/amp_mean_param_noise/ideal_fit_M_d_D_L_migkas.npz")
MS0, EPS0, X0, L0DEG = float(best["mass_scale"]), float(best["eps"]), float(best["x"]), float(best["phi0_deg"])
B0DEG = 0.0

_ucache = {}
def u_base(x):
    """Radial peculiar velocity table for mass_scale = 1 (u is linear in mass)."""
    key = round(float(x), 6)
    if key in _ucache:
        return _ucache[key]
    p = params_twoscale_local_match(float(x), FU)
    r_grid = np.linspace(0.2 * D, 2.0 * D, NRAD)
    chi = r_grid / T ** 0.4
    t = np.linspace(ETA * T, T, NT)
    dt = t[1] - t[0]
    u = np.zeros(NRAD)
    for ti in t[1:]:
        Hi = 2.0 / (5.0 * ti)
        ri = chi * ti ** 0.4
        gi = G * m_enc_at_r(ri, p) / np.maximum(ri, 1e-30) ** 2
        u = u + dt * (-Hi * u - gi)
    _ucache[key] = (r_grid, u, float(p["kappa"]))
    return _ucache[key]

def kappa_of(ms, x):
    return ms * u_base(x)[2]

def boxcar_dl90(y):
    n = len(y); half = int(round(45.0 / (360.0 / n)))
    k = np.ones(2 * half + 1) / (2 * half + 1)
    return np.convolve(np.pad(y, (half, half), mode="wrap"), k, mode="valid")

def model_A(ms, eps, x, l0, b0, n_phi=181):
    r_grid, ub, _ = u_base(x)
    u = lambda r: ms * np.interp(r, r_grid, ub)
    l_deg = np.linspace(0.0, 360.0, n_phi, endpoint=False)
    # direction along the Galactic plane, expressed in the model frame whose
    # +x axis is the anti-CoE direction: mu = cos b0 cos(l - l0)
    mu = np.cos(np.radians(b0)) * np.cos(np.radians(l_deg - l0))
    s = np.sqrt(np.clip(1 - mu ** 2, 0, None))
    n_hat = np.stack([mu, s, np.zeros_like(mu)], axis=-1)   # axisymmetric: azimuth irrelevant
    d = eps * D
    r_vec = np.array([D, 0, 0]) + d * n_hat
    r = np.linalg.norm(r_vec, axis=-1)
    u_r = u(r)[:, None] * r_vec / r[:, None]
    u_D = u(D) * np.array([1.0, 0, 0])
    Hobs = H + np.sum((u_r - u_D) * n_hat, axis=-1) / d
    AF = Hobs / H; AF = AF / AF.mean(); AF = boxcar_dl90(AF); AF = AF / AF.mean()
    return np.interp(ell_m % 360.0, l_deg, AF, period=360.0)

MODE = "corr2"
SST = float(((A_m - A_m.mean()) ** 2).sum())
def r2(ms, eps, x, l0=L0DEG, b0=B0DEG, mode=None):
    m = model_A(ms, eps, x, l0, b0)
    if (mode or MODE) == "corr2":
        return float(np.corrcoef(m, A_m)[0, 1] ** 2)
    return float(1 - ((A_m - m) ** 2).sum() / SST)

R2_CHECK = r2(MS0, EPS0, X0, mode="corr2")
K0 = kappa_of(MS0, X0)
print(f"reproduce best fit: R2={R2_CHECK:.4f} (stored {float(best['R2']):.4f}), kappa={K0:.4f} (stored {float(best['kappa']):.4f})")

def edges(v):
    dv = np.diff(v).mean(); return np.concatenate([v - dv / 2, [v[-1] + dv / 2]])

FOOT = ("Neighbouring Migkas $\\Delta\\ell=90^\\circ$ sectors overlap (71 points are not independent); "
        "$R^2$ differences are indicative.")

def heat(xv, yv, Z, xlab, ylab, sub, name, marks=(), hlines=(), vlines=()):
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    pm = ax.pcolormesh(edges(xv), edges(yv), Z, cmap="magma", shading="flat", vmin=0.0,
                       edgecolors="face", linewidth=0, rasterized=True)
    iy, ix = np.unravel_index(np.nanargmax(Z), Z.shape)
    for y0, lab in hlines:
        ax.axhline(y0, color="#17becf", ls=":", lw=1.6, label=lab)
    for x0, lab in vlines:
        ax.axvline(x0, color="#17becf", ls=":", lw=1.6, label=lab)
    for (mx, my, mk, lab) in marks:
        ax.plot(mx, my, marker=mk, ms=10, color="#17becf", mec="k", mew=0.7, ls="none", label=lab)
    ax.plot(xv[ix], yv[iy], marker="*", ms=16, color="w", mec="k", mew=0.8, ls="none",
            label=f"best cell: ({xv[ix]:.3g}, {yv[iy]:.3g}), $R^2={Z[iy, ix]:.3f}$")
    ax.set_xlabel(xlab); ax.set_ylabel(ylab)
    ax.set_title(sub, fontsize=8.8)
    cb = fig.colorbar(pm, ax=ax, extend="min" if np.nanmin(Z) < 0 else "neither")
    cb.set_label(r"$R^2_A$ vs Migkas $A/A_{\mathrm{all}}$" + (" (after linear rescale)" if MODE == "corr2" else " (no rescale)"))
    ax.legend(loc="best", fontsize=8, framealpha=0.9)
    fig.text(0.01, 0.01, FOOT, fontsize=7, color="#444")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/heatmap_R2_{name}{'_corr2' if MODE == 'corr2' else ''}.{ext}", dpi=300)
    plt.close(fig)
    return {"best": {"x": float(xv[ix]), "y": float(yv[iy]), "R2": float(Z[iy, ix])}}

res = {"model": __doc__, "R2_definition": "squared Pearson correlation (= R^2 after best linear rescale)",
       "best_fit_reference": {"mass_scale": MS0, "d_over_D": EPS0, "D_over_L": X0, "l0_deg": L0DEG,
                              "b0_deg": B0DEG, "kappa": K0, "R2_recomputed": R2_CHECK,
                              "R2_stored": float(best["R2"])}}

eps_g = np.round(np.arange(0.02, 0.8001, 0.02), 3)
x_g = np.round(np.arange(3.0, 10.001, 0.2), 3)
k_g = np.round(np.arange(0.05, 1.0001, 0.025), 4)
fixed = lambda s: s

ALL = {}
for MODE in ("direct", "corr2"):
    R = {}
    # 1. d/D vs D/L  (kappa held at best value, so mass_scale varies with D/L)
    Z1 = np.array([[r2(K0 / u_base(x)[2], e, x) for x in x_g] for e in eps_g])
    R["dD_vs_DL"] = heat(x_g, eps_g, Z1, r"$D/L$", r"$d/D$",
        rf"$\kappa={K0:.3f}$, anti-CoE $(\ell_0,b_0)=({L0DEG:.1f}^\circ,0^\circ)$, $f_u=0.94$, drag, $\eta=10^{{-4}}$",
        "dD_vs_DL", marks=[(X0, EPS0, "o", "best physical fit")])
    R["dD_vs_DL"].update({"x_grid": x_g.tolist(), "y_grid": eps_g.tolist(), "Z": Z1.tolist(),
                            "fixed": {"kappa": K0, "l0": L0DEG, "b0": 0.0}})

    # 2. d/D vs kappa (D/L held)
    Z2 = np.array([[r2(k / u_base(X0)[2], e, X0) for k in k_g] for e in eps_g])
    R["dD_vs_kappa"] = heat(k_g, eps_g, Z2, r"$\kappa = GM(<D)/(D^3H^2)$", r"$d/D$",
        rf"$D/L={X0:.2f}$, anti-CoE $(\ell_0,b_0)=({L0DEG:.1f}^\circ,0^\circ)$, $f_u=0.94$, drag, $\eta=10^{{-4}}$",
        "dD_vs_kappa", marks=[(K0, EPS0, "o", "best physical fit")],
        vlines=[(0.15, r"$\kappa=\Omega_m/2=0.15$")])
    R["dD_vs_kappa"].update({"x_grid": k_g.tolist(), "y_grid": eps_g.tolist(), "Z": Z2.tolist(),
                               "fixed": {"D_over_L": X0, "l0": L0DEG, "b0": 0.0}})

    # 3. axis orientation
    l0_g = np.arange(0.0, 360.0, 5.0)
    b0_g = np.arange(-60.0, 60.001, 5.0)
    Z3 = np.array([[r2(MS0, EPS0, X0, l0, b0) for l0 in l0_g] for b0 in b0_g])
    R["axis_l0_b0"] = heat(l0_g, b0_g, Z3, r"anti-CoE longitude $\ell_0$ [deg]", r"anti-CoE latitude (tilt) $b_0$ [deg]",
        rf"$d/D={EPS0:.3f}$, $D/L={X0:.2f}$, $\kappa={K0:.3f}$; model evaluated along the Galactic plane",
        "axis_l0_b0", marks=[(0.0, 15.0, "^", r"paper pose: anti-CoE $(0^\circ,+15^\circ)$")])
    R["axis_l0_b0"].update({"x_grid": l0_g.tolist(), "y_grid": b0_g.tolist(), "Z": Z3.tolist(),
                              "fixed": {"d_over_D": EPS0, "D_over_L": X0, "kappa": K0},
                              "R2_at_paper_pose_0_15": r2(MS0, EPS0, X0, 0.0, 15.0)})

    # 4. D/L vs kappa (d/D held)
    Z4 = np.array([[r2(k / u_base(x)[2], EPS0, x) for x in x_g] for k in k_g])
    R["DL_vs_kappa"] = heat(x_g, k_g, Z4, r"$D/L$", r"$\kappa = GM(<D)/(D^3H^2)$",
        rf"$d/D={EPS0:.3f}$, anti-CoE $(\ell_0,b_0)=({L0DEG:.1f}^\circ,0^\circ)$, $f_u=0.94$, drag, $\eta=10^{{-4}}$",
        "DL_vs_kappa", marks=[(X0, K0, "o", "best physical fit")],
        hlines=[(0.15, r"$\kappa=\Omega_m/2=0.15$")])
    R["DL_vs_kappa"].update({"x_grid": x_g.tolist(), "y_grid": k_g.tolist(), "Z": Z4.tolist(),
                               "fixed": {"d_over_D": EPS0, "l0": L0DEG, "b0": 0.0}})

    ALL[MODE] = R
MODE = "direct"
res["maps"] = ALL
# extra diagnostics
D_ = ALL["direct"]
ik = int(np.argmin(abs(k_g - 0.15)))
Z2d, Z4d = np.array(D_["dD_vs_kappa"]["Z"]), np.array(D_["DL_vs_kappa"]["Z"])
res["diagnostics"] = {
    "kappa_base_of_DL": {str(x): u_base(x)[2] for x in x_g},
    "direct_R2_at_best_fit_params": r2(MS0, EPS0, X0, mode="direct"),
    "direct_R2_at_kappa0.15_best_eps_x": r2(0.15 / u_base(X0)[2], EPS0, X0, mode="direct"),
    "direct_max_R2_along_kappa0.15_dD_map": {"R2": float(Z2d[:, ik].max()), "d_over_D": float(eps_g[np.argmax(Z2d[:, ik])])},
    "direct_max_R2_along_kappa0.15_DL_map": {"R2": float(Z4d[ik, :].max()), "D_over_L": float(x_g[np.argmax(Z4d[ik, :])])},
    "corr2_kappa_spread_dD_map": float(np.ptp(np.array(ALL["corr2"]["dD_vs_kappa"]["Z"]), axis=1).max()),
}
json.dump(res, open(f"{OUT}/heatmaps_R2_model_params.json", "w"), indent=1)
for m in ALL:
    for k in ("dD_vs_DL", "dD_vs_kappa", "axis_l0_b0", "DL_vs_kappa"):
        print(m, k, ALL[m][k]["best"], ALL[m][k].get("R2_at_paper_pose_0_15", ""))
print({k: v for k, v in res["diagnostics"].items() if k != "kappa_base_of_DL"})
