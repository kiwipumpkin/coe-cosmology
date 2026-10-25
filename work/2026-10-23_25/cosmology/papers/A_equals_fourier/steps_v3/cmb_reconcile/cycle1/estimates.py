"""Cycle 1 quick numbers for all reconciliation candidates. Reuses ../../cmb_residual/cmb_residual.py (unchanged) for geometry,
fit parameters (dscan.json headline points) and the model <H> anchors. Writes estimates.json.
Conventions: theta from CoE direction, a_l = Legendre coefficients of dT/T; 'residual' = l>=2; pass = |a_l| <~ 1e-5.
Coasting (Om=0) null-ray shape, R = 14.3 Gpc at a_e = 1/1101, as in the prior step."""
import sys, os, json, numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "..", "cmb_residual"))
import cmb_residual as CR
G, c, GPC, MPC = CR.G, CR.c, CR.GPC, CR.MPC; R = CR.R_GPC * GPC; A_E = CR.A_E; MU = CR.MU; leg = CR.legendre
PASS = 1e-5; OUT = {}
FITS = {k: CR.point(k, v) for k, v in CR.HEAD.items()}
HMOD = {k: CR.anchor[k]["H"] for k in FITS}

# ---------- data-driven shear from the Migkas Fourier form (model-independent requirement) ----------
AF = {0: (0.9965, 0.0), 1: (0.0296, 0.0434), 2: (0.0636, 0.0366), 3: (-0.0299, -0.0114)}
C1, C2, C3 = (np.hypot(*AF[m]) for m in (1, 2, 3))
def a2_series(b):
    B = 2 * b + b * b; m0 = 1 + B / 4 - 3 * B * B / 64; return (B / 4 - B * B / 16) / m0
BETA = brentq(lambda b: a2_series(b) - C2 / AF[0][0], 1e-4, 1.0)      # h'/h at Earth needed for the 2l amplitude
OUT["data"] = dict(C1=C1, C2=C2, C3=C3, beta_needed=BETA, note="beta = dln h/dln q at Earth (radial vs transverse rate: H_r/H_perp = 1+beta)")

def grid(ch):                                   # comoving points on sphere radius ch around Earth; CoE at +z, distance D
    return ch
def h_point(q, beta):                           # point-core radial rate profile normalised h(1)=1, h'(1)=beta, h''(1)=-4beta
    return 1 + beta / 3 * (1 - q**-3)

def doppler_profile(p, H, beta, s):
    """l-pattern of n.(u_o - u_e)/c for emitters on a sphere of comoving radius s, flow V = k x h(|x|/D) about the CoE (linear velocities)."""
    D = p["D"] * GPC; k = H * 1e3 / MPC; mu = MU
    xz = s * mu - D; xp = s * np.sqrt(1 - mu**2); q = np.hypot(xz, xp) / D     # emitter position relative to CoE, theta from CoE dir.
    n_dot_x = s - D * mu                                                           # n.(x_e - x_CoE) with n = (sin, 0, cos)... (n.x_e = s, n.x_CoE = D mu)
    u_e = k * (h_point(q, beta) - 1) * n_dot_x
    return -u_e / c                                                                # u_o = 0 because h(1)=1 (observer's own motion only gives a dipole anyway)

def G_run(p, H, Gfun, s_emit=R, a_emit=A_E, na=6000):
    """SW + ISW for Phi(x,a) = Phi0_today(x) * G(a), G(1)=1 (Phi0 = -G Mc/|x - x_CoE| with today's Mc). Coasting ray shape, emission at s_emit, a_emit."""
    D = p["D"] * GPC; GM = G * p["Mc"]; a = np.geomspace(a_emit, 1.0, na)
    chi = s_emit * np.log(a) / np.log(a_emit)                                      # coasting: chi proportional to ln(1/a)
    m = MU[:, None]; X = np.sqrt(D**2 + chi[None, :]**2 - 2 * D * chi[None, :] * m)
    P0 = -GM / X / c**2; Ga = Gfun(a)
    sw = P0[:, 0] * Ga[0]; isw = 2 * np.trapezoid(P0, Ga[None, :], axis=1)
    return leg(sw), leg(isw), leg(sw + isw), leg(P0[:, 0])

# ---------- PRIORITY (Chris): static piece, time-varying piece, ratio to matter-velocity effect ----------
pri = {}
for f, p in FITS.items():
    H = HMOD[f]; e = {}
    sw_static = leg((-G * p["Mc"] / np.sqrt((p["D"] * GPC)**2 + R**2 - 2 * p["D"] * GPC * R * MU) + G * p["Mc"] / (p["D"] * GPC)) / c**2)
    e["static_Phi_emit_minus_obs_today_field"] = dict(a2=sw_static[2], a3=sw_static[3], a4=sw_static[4],
        analytic_a2=-p["mcc"] * (H * 1e3 / MPC * R / c)**2 / 2 * (70 / H)**2)
    sw, isw, tot, _ = G_run(p, H, lambda a: 1 / a)                                 # constant Mc: Phi ~ 1/a
    e["constant_Mc_Phi_1_over_a"] = dict(SW_a2=sw[2], ISW_a2=isw[2], net_a2=tot[2], net_a3=tot[3], ISW_over_SW=isw[2] / sw[2])
    sw1, isw1, tot1, _ = G_run(p, H, lambda a: np.ones_like(a))                    # static potential (Mc grows as a): ISW = 0
    e["static_potential_all_time"] = dict(SW_a2=sw1[2], ISW_a2=isw1[2], net_a2=tot1[2])
    dop_fit = leg(doppler_profile(p, H, BETA, R))                                  # matter velocity pattern the A fit needs, at the emission sphere
    s_mig = 0.5 * GPC; dop_mig = leg(doppler_profile(p, H, BETA, s_mig))
    stat_mig = leg((-G * p["Mc"] / np.sqrt((p["D"] * GPC)**2 + s_mig**2 - 2 * p["D"] * GPC * s_mig * MU)) / c**2)
    e["velocity_effect_fit_beta"] = dict(at_R_a2=dop_fit[2], at_R_a3=dop_fit[3], at_500Mpc_a2=dop_mig[2])
    e["ratio_static_over_velocity_at_R"] = sw_static[2] / dop_fit[2]
    e["ratio_net_constMc_over_velocity_at_R"] = tot[2] / dop_fit[2]
    e["ratio_static_over_velocity_at_500Mpc"] = stat_mig[2] / dop_mig[2]
    pri[f] = e
OUT["priority_static_ISW_ratio"] = pri

# ---------- candidate 1: velocity-only ----------
c1 = {}
for f, p in FITS.items():
    d = pri[f]["velocity_effect_fit_beta"]["at_R_a2"]
    c1[f] = dict(doppler_a2_today_pattern=d, doppler_a2_if_ballistic_at_emission=d / A_E,  # peculiar v ~ 1/a with no force (coasting: x1101)
                 suppression_needed=PASS / abs(d), sigmaH_today_allowed_ballistic=BETA / 3 * PASS / abs(d / A_E))
OUT["c1_velocity_only"] = c1

# local confined flow (gravity-only, growing mode, Om_b = 0.049, Lambda = 0): potential and late-ISW estimate
def growth_open(Om):
    # Lambda=0 dust: D'' + (3/a + H'/H) D' - 1.5 Om/(a^5 E^2) D = 0 in a; return f(1)
    E2 = lambda a: Om / a**3 + (1 - Om) / a**2
    def rhs(a, y):
        D, dD = y; dE2 = -3 * Om / a**4 - 2 * (1 - Om) / a**3
        return [dD, -(3 / a + dE2 / (2 * E2(a))) * dD + 1.5 * Om / (a**5 * E2(a)) * D]
    s = solve_ivp(rhs, [1e-4, 1], [1e-4, 1.0], rtol=1e-9); D1, dD1 = s.y[0, -1], s.y[1, -1]; return dD1 / D1
loc = {}
for Om in (0.049, 1.0):
    f = growth_open(Om); rows = {}
    for sig in (0.036, 0.05):
        for Lg in (0.3, 0.5, 0.7, 1.0, 1.25):
            x = HMOD["omb"] * 1e3 * Lg * 1e3 / (c)                                 # H L / c
            Phi2 = 1.5 * Om / f * sig * (x / np.pi)**2
            a2 = 2 * abs(1 - f) * x * Phi2 if Om < 1 else 2 * 0.0 * x * Phi2     # Om=1: Phi constant, late ISW = 0 (SW not involved: structure is local)
            rows[f"sig={sig}|L={Lg}Gpc"] = dict(HL_c=x, Phi2=Phi2, delta_needed=3 * sig / f, late_ISW_a2=a2)
        Lmax = brentq(lambda Lg: 2 * abs(1 - f) * (HMOD["omb"] * Lg * 1e6 / c) * 1.5 * Om / f * sig * (HMOD["omb"] * Lg * 1e6 / c / np.pi)**2 - PASS, 0.01, 20) if Om < 1 else None
        rows[f"sig={sig}|Lmax_Gpc_for_1e-5"] = Lmax
    loc[f"Om={Om}"] = dict(f_today=f, rows=rows)
OUT["local_confined_flow"] = loc

# ---------- candidate 2: SW-ISW cancellation by core-mass history ----------
c2 = {}
for f, p in FITS.items():
    H = HMOD[f]
    a2p = lambda pp: G_run(p, H, lambda a: a**(pp - 1))[2][2]
    pz = brentq(a2p, 0.0, 1.0); sw, isw, tot, _ = G_run(p, H, lambda a: a**(pz - 1))
    dp = 1e-4; slope = (a2p(pz + dp) - a2p(pz - dp)) / (2 * dp)
    # step: Phi = 2*Phi_today just before a_s, Phi_today after, a_s just after emission (1.02 a_e); then static
    step = lambda a, r=0.5, a_s=1.02 * A_E: np.where(a < a_s, 1 / r, 1.0)
    swS, iswS, totS, _ = G_run(p, H, step)
    ratio_tol = brentq(lambda r: G_run(p, H, lambda a: step(a, r))[2][2] - PASS, 0.3, 0.5) if G_run(p, H, lambda a: step(a, 0.3))[2][2] * G_run(p, H, lambda a: step(a, 0.5))[2][2] < 0 else None
    # a_i scan: constant Mc, SW+ISW+Doppler (prior linear machinery); zero crossings of l=2 and l=3
    ais = np.geomspace(1e-6, A_E * 0.999, 60); A2 = []; A3 = []
    for ai in ais:
        r = CR.rays(p, H * 1e3 / MPC, p["Mc"], "coast", ai, MU); L = leg(r["total"]); A2.append(L[2]); A3.append(L[3])
    A2, A3 = np.array(A2), np.array(A3)
    z2 = [float(ais[i]) for i in range(len(ais) - 1) if A2[i] * A2[i + 1] < 0]; z3 = [float(ais[i]) for i in range(len(ais) - 1) if A3[i] * A3[i + 1] < 0]
    a3_at_z2 = float(np.interp(np.log(z2[0]), np.log(ais), A3)) if z2 else None
    c2[f] = dict(p_Mc_prop_a_p=pz, SW_a2=sw[2], ISW_a2=isw[2], residual_a3=tot[3], residual_a4=tot[4], da2_dp=slope, p_tolerance=PASS / abs(slope),
                 step_half_after_emission=dict(net_a2=totS[2], net_a3=totS[3], net_a4=totS[4], SW_a2=swS[2]), step_ratio_for_a2_1em5=ratio_tol,
                 ai_scan_zero_l2=z2, ai_scan_zero_l3=z3, a3_at_l2_zero=a3_at_z2)
OUT["c2_SW_ISW_cancellation"] = c2

# ---------- candidate 3: fixed-temperature front / Tolman ----------
H0s = HMOD["omb"] * 1e3 / MPC
c3 = dict(external_Laplace_potential_intrinsic_term=0.0, SW_factor=1.0, adiabatic_selfsourced_SW_factor=1 / 3,
          frozen_a2_with_1_3={f: pri[f]["static_Phi_emit_minus_obs_today_field"]["a2"] / 3 for f in FITS})
for f, p in FITS.items():   # Tolman equilibrium needs sound crossing of the potential scale D before a_e (coasting: r_s = c/(sqrt3 H0) ln(a_e/a_i))
    for N in (1, 10):
        c3[f"{f}|ln(a_e/a_i)_for_{N}_crossings_of_D"] = N * p["D"] * GPC * np.sqrt(3) * H0s / c
OUT["c3_recombination_front"] = c3

# ---------- candidate 4: late last scattering / late core ----------
c4 = {}
for f, p in FITS.items():
    H = HMOD[f]
    def a2_ls(z):
        a_ls = 1 / (1 + z); s = R * np.log(1 + z) / np.log(1 / A_E); return G_run(p, H, lambda a: 1 / a, s_emit=s, a_emit=a_ls)[2][2]
    zmax = brentq(lambda z: abs(a2_ls(z)) - PASS, 1e-3, 1100)
    K = p["mcc"] * (1 / 0.01 - 1) / 2                                              # fit's velocity-shear parameter (tau_i = 0.01)
    def a2_late(af):                                                               # core switched on at a_f, mass rescaled to keep K fixed
        scale = (1 / 0.01 - 1) / (1 / af - 1); return scale * G_run(p, H, lambda a: np.where(a < af, 0.0, 1 / a))[2][2]
    rows = {f"a_f={af}": a2_late(af) for af in (0.01, 0.1, 0.5, 0.9, 0.99)}
    c4[f] = dict(z_ls_max_for_1em5=zmax, K_fit=K, late_core_ISW_a2=rows)
OUT["c4_late_thermalisation_or_core"] = c4

# ---------- candidate 7: rotation (transverse velocity, invisible to redshift) ----------
OUT["c7_rotation"] = dict(omega_over_k_needed=np.sqrt(4 * C2), transverse_v2_over_c2_at_R=(np.sqrt(4 * C2) * HMOD["omb"] * 1e3 * R / MPC / c)**2,
                          transverse_v2_over_c2_at_500Mpc=(np.sqrt(4 * C2) * HMOD["omb"] * 500 / 299792.458)**2,
                          omega_over_H_at_emission_coasting=np.sqrt(4 * C2) / A_E)
# tilt / bulk flow numbers
OUT["c6_tilt"] = dict(v_bulk_kms_if_A_prop_H_zeff0p1=C1 * 299792.458 * 0.1, v_bulk_kms_if_A_prop_H2=C1 / 2 * 299792.458 * 0.1,
                      sigmaH_allowed_superhorizon=PASS / (HMOD["omb"] * 1e3 * R / MPC / c)**2)
json.dump(OUT, open(os.path.join(HERE, "estimates.json"), "w"), indent=1, default=float)
def pr(o, ind=0):
    for k, v in o.items():
        if isinstance(v, dict): print(" " * ind + k); pr(v, ind + 2)
        else: print(" " * ind + f"{k}: {v:.4g}" if isinstance(v, (float, np.floating)) else " " * ind + f"{k}: {v}")
pr(OUT)
