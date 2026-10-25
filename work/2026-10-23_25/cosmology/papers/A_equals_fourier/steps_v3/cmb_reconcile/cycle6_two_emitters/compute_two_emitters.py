#!/usr/bin/env python3
"""
Cycle 6: celestial sphere between two emitters — parallel plates OR concentric shells.

Rest-frame blackbody emission at T: specific intensity along each ray = B_ν(T_emitter).
Brightness temperature T_obs = T_emitter for the surface the ray hits (no intervening opacity).

Geometry decides WHICH surface is seen; if both surfaces share the same rest-frame T,
the entire sky is T → perfectly isotropic. Anisotropy only from unequal T or Doppler.
"""
from __future__ import annotations
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- geometry params (dimensionless; only ratios matter for rest-frame) ----
D_PLATE = 1.0          # half-separation of plates (±d); observer z in (-d, d)
R_IN = 1.0
R_OUT = 3.0
D_SHELL = 1.8          # R_in < D < R_out
T0 = 1.0               # rest-frame temperature units
BETA_OBS = 0.01        # optional observer motion (plates) β = v/c along +z
BETA_REC = 0.02        # optional shell recession scale: v = β_rec * (r/R_out) * c  (Hubble-like)


def plates_T_obs(theta_deg, z_obs=0.0, T_plus=T0, T_minus=T0, beta_obs=0.0):
    """
    Infinite parallel plates at z = ±d. Observer at z = z_obs, |z_obs| < d.
    theta = 0 looks toward +z (upper plate), theta = 180 toward -z (lower plate).
    Every ray with cosθ > 0 hits +plate; cosθ < 0 hits -plate; equator is measure-zero.
    Rest-frame: T_obs = T_plate. With observer velocity +β ĉ_z: dipole on top.
    """
    th = np.asarray(theta_deg, float)
    mu = np.cos(np.deg2rad(th))  # +1 toward +plate
    T_rest = np.where(mu >= 0.0, T_plus, T_minus)
    # Special-rel Doppler for observer moving at β along +z through isotropic field
    # of local T_rest (piecewise). Photon direction in lab: n̂; T_obs = T_rest / [γ(1 - β·n̂)]
    # For β along +z: β·n̂ = β mu. Approaching upper plate (mu>0) → blueshift if β>0? 
    # Observer velocity +z means moving toward upper plate → blueshift for +plate rays.
    # Standard CMB dipole formula: T'(n) = T / [γ(1 - β·n̂)] where n̂ is photon arrival direction
    # (= look direction). If β points toward a hot spot, that direction is blueshifted.
    if abs(beta_obs) < 1e-15:
        return T_rest
    g = 1.0 / np.sqrt(1.0 - beta_obs**2)
    return T_rest / (g * (1.0 - beta_obs * mu))


def shell_hit(theta_deg, D=D_SHELL, R_in=R_IN, R_out=R_OUT):
    """
    Observer at (0,0,D) between concentric shells. theta=0 toward CoE (−ẑ),
    theta=180 away (+ẑ). Returns which shell is hit ('in'|'out'), λ, and
    unit radial direction of the hit point from CoE.
    """
    th = np.asarray(theta_deg, float)
    phi = np.deg2rad(th)  # angle from toward-CoE
    # n̂ = look direction: toward CoE is −ẑ
    nz = -np.cos(phi)
    nx = np.sin(phi)
    # r_obs · n̂ = D * nz
    rdotn = D * nz

    def intersect(R):
        disc = rdotn**2 - (D**2 - R**2)
        out = {
            "hit": np.zeros(th.shape, dtype=bool),
            "lam": np.full(th.shape, np.nan),
            "rx": np.full(th.shape, np.nan),
            "rz": np.full(th.shape, np.nan),
        }
        ok = disc >= -1e-14
        disc = np.maximum(disc, 0.0)
        sqrt_d = np.sqrt(disc)
        # two roots; take smallest λ > 0
        lam1 = -rdotn - sqrt_d
        lam2 = -rdotn + sqrt_d
        # choose positive root(s)
        cand = np.full(th.shape, np.nan)
        for L in (lam1, lam2):
            pos = (L > 1e-12) & ok
            use = pos & (np.isnan(cand) | (L < cand))
            cand = np.where(use, L, cand)
        out["hit"] = ~np.isnan(cand)
        out["lam"] = cand
        # hit position
        hx = nx * cand
        hz = D + nz * cand
        rn = np.sqrt(hx**2 + hz**2)
        out["rx"] = np.where(out["hit"], hx / rn, np.nan)
        out["rz"] = np.where(out["hit"], hz / rn, np.nan)
        return out

    inn = intersect(R_in)
    # Prefer inner shell when it is hit (closer toward CoE disk)
    which = np.where(inn["hit"], "in", "out")
    # For rays that miss inner, they must hit outer (observer inside outer shell)
    outr = intersect(R_out)
    # Build hit kinematics
    hit_R = np.where(inn["hit"], R_in, R_out)
    lam = np.where(inn["hit"], inn["lam"], outr["lam"])
    rx = np.where(inn["hit"], inn["rx"], outr["rx"])
    rz = np.where(inn["hit"], inn["rz"], outr["rz"])
    assert np.all(outr["hit"] | inn["hit"]), "some rays miss both shells"
    return {
        "which": which,
        "hit_R": hit_R,
        "lam": lam,
        "rx": rx,
        "rz": rz,
        "alpha_deg": float(np.rad2deg(np.arcsin(R_in / D))),  # angular radius of inner disk
        "nx": nx,
        "nz": nz,
    }


def shells_T_obs(theta_deg, D=D_SHELL, R_in=R_IN, R_out=R_OUT,
                 T_in=T0, T_out=T0, beta_rec=0.0, beta_obs_radial=0.0):
    """
    Rest-frame: T_obs = T_in or T_out by which shell is hit.
    Optional Hubble-like Doppler: emitter vel = β_rec * c * (r_hit / R_out) * ê_r
    (radial from CoE); observer vel = β_obs_radial * c * ê_z (ê_z = CoE→observer).
    Relative emitter velocity in CoE frame: v_em; observer: v_obs.
    Photon look direction n̂; first-order: ΔT/T = − n̂ · (v_em − v_obs)/c
    (CMB dipole convention: motion of observer toward emitter → blueshift).
    Equivalently T_obs ≈ T_rest * (1 − n̂·v_em/c) / (1 − n̂·v_obs/c) ≈ T_rest*(1 − n̂·(v_em−v_obs)/c).
    """
    th = np.asarray(theta_deg, float)
    geo = shell_hit(th, D=D, R_in=R_in, R_out=R_out)
    T_rest = np.where(geo["which"] == "in", T_in, T_out)
    if abs(beta_rec) < 1e-15 and abs(beta_obs_radial) < 1e-15:
        return T_rest, geo

    # emitter β vector (in units of c): β_rec * (R_hit/R_out) * ê_r
    bem_scale = beta_rec * (geo["hit_R"] / R_out)
    bem_x = bem_scale * geo["rx"]
    bem_z = bem_scale * geo["rz"]
    bobs_x = 0.0
    bobs_z = beta_obs_radial
    # n̂ · (β_em − β_obs)
    ndb = geo["nx"] * (bem_x - bobs_x) + geo["nz"] * (bem_z - bobs_z)
    # Use exact collinear Doppler along each ray for small β: T' = T * sqrt((1-β_los)/(1+β_los))
    # with β_los = n̂·(β_em − β_obs)  (positive β_los = emitter receding along look dir → redshift)
    # For |β|≪1: T' ≈ T (1 − β_los)
    T = T_rest * (1.0 - ndb)
    return T, geo


def sky_map_shells(nside_lon=720, nside_lat=360, **kw):
    """lon 0..360, lat −90..90; θ_from_CoE = 90 − lat for axis along z… 
    We put CoE direction at lon=0, lat=0 for great-circle friendliness? 
    Better: put toward-CoE at north pole for Mollweide (standard dipole style).
    θ = 0 toward CoE = +z polar; Mollweide uses lon, lat with lat = 90 − θ.
    """
    lon = np.linspace(-180, 180, nside_lon, endpoint=False)
    lat = np.linspace(-90, 90, nside_lat)
    Lon, Lat = np.meshgrid(lon, lat)
    # θ from CoE axis: toward CoE at north pole ⇒ θ = 90° − lat
    theta = 90.0 - Lat
    # azimuthal symmetry: T depends only on θ
    T1d, geo = shells_T_obs(np.linspace(0, 180, 181), **kw)
    # interpolate
    th_grid = np.linspace(0, 180, 181)
    Tflat = np.interp(theta.ravel(), th_grid, T1d).reshape(theta.shape)
    return Lon, Lat, Tflat, geo


def main():
    th360 = np.linspace(0, 360, 721)
    th180 = np.linspace(0, 180, 361)

    # --- plates ---
    # map θ∈[0,360] great circle: 0 = +z, 180 = −z, 360 = +z
    # For plates, use polar angle from +z in [0,180] then mirror for full circle
    def plate_circle(z_obs, T_plus, T_minus, beta):
        # great circle through poles: θ_from_+z goes 0→360, with mu = cos(θ) for θ in 0..180
        # and for 180..360, mu = cos(360−θ) wait: standard polar: at φ=0 plane,
        # direction angle α from +z: α∈[0,360) with mu = cos(α) only if α in [0,180];
        # for α in [180,360], cos(α) still works (cos 270=0, cos 360=1). Yes cos is fine on [0,360].
        return plates_T_obs(th360, z_obs=z_obs, T_plus=T_plus, T_minus=T_minus, beta_obs=beta)

    plates = {
        "theta_deg": th360.tolist(),
        "mid_sameT": plate_circle(0.0, T0, T0, 0.0).tolist(),
        "offset_sameT": plate_circle(0.5 * D_PLATE, T0, T0, 0.0).tolist(),
        "mid_diffT": plate_circle(0.0, 1.1 * T0, 0.9 * T0, 0.0).tolist(),
        "mid_sameT_doppler": plate_circle(0.0, T0, T0, BETA_OBS).tolist(),
        "offset_diffT": plate_circle(0.5 * D_PLATE, 1.1 * T0, 0.9 * T0, 0.0).tolist(),
        "params": {
            "d": D_PLATE, "T0": T0, "beta_obs": BETA_OBS,
            "T_plus_diff": 1.1 * T0, "T_minus_diff": 0.9 * T0,
            "z_offset": 0.5 * D_PLATE,
        },
    }

    # --- shells ---
    def shell_circle(T_in, T_out, beta_rec, beta_obs_radial, D=D_SHELL):
        # full 0..360: θ_from_CoE goes 0→180→0 along a great circle (forward then back other side)
        # Use |α| folded: for α in [0,180] use α; for α in [180,360] use 360−α (other meridian, same θ)
        th_from_coe = np.where(th360 <= 180, th360, 360.0 - th360)
        T, geo = shells_T_obs(th_from_coe, D=D, T_in=T_in, T_out=T_out,
                              beta_rec=beta_rec, beta_obs_radial=beta_obs_radial)
        return T, geo

    T_same, geo0 = shell_circle(T0, T0, 0.0, 0.0)
    T_diff, geo1 = shell_circle(1.1 * T0, 0.9 * T0, 0.0, 0.0)
    T_dop, geo2 = shell_circle(T0, T0, BETA_REC, 0.0)  # Hubble-like on emitters; observer comoving at D ⇒ set beta_obs = beta_rec*(D/R_out)
    # Proper comoving observer in same Hubble flow:
    T_dop_comoving, geo3 = shell_circle(T0, T0, BETA_REC, BETA_REC * (D_SHELL / R_OUT))
    # Observer at rest w.r.t. CoE (not comoving):
    T_dop_coe_rest, geo4 = shell_circle(T0, T0, BETA_REC, 0.0)

    # also centred observer (D→0 limit not allowed; use D very close to centre — skip;
    # instead D such that still between: report alpha)
    shells = {
        "theta_deg": th360.tolist(),
        "alpha_deg": geo0["alpha_deg"],
        "D": D_SHELL, "R_in": R_IN, "R_out": R_OUT,
        "sameT": T_same.tolist(),
        "diffT": T_diff.tolist(),
        "doppler_coe_rest_obs": T_dop_coe_rest.tolist(),
        "doppler_comoving_obs": T_dop_comoving.tolist(),
        "which_0_180": shells_T_obs(th180)[1]["which"].tolist(),
        "theta_0_180": th180.tolist(),
        "params": {
            "T_in_diff": 1.1 * T0, "T_out_diff": 0.9 * T0,
            "beta_rec": BETA_REC,
            "beta_obs_comoving": BETA_REC * (D_SHELL / R_OUT),
        },
    }

    # summary stats
    def aniso(y):
        y = np.asarray(y, float)
        return {
            "mean": float(y.mean()),
            "ptp": float(y.max() - y.min()),
            "std": float(y.std()),
            "is_isotropic": bool(np.ptp(y) < 1e-12 * max(abs(y.mean()), 1.0)),
        }

    summary = {
        "plates_mid_sameT": aniso(plates["mid_sameT"]),
        "plates_offset_sameT": aniso(plates["offset_sameT"]),
        "plates_mid_diffT": aniso(plates["mid_diffT"]),
        "plates_mid_sameT_doppler": aniso(plates["mid_sameT_doppler"]),
        "shells_sameT": aniso(shells["sameT"]),
        "shells_diffT": aniso(shells["diffT"]),
        "shells_doppler_coe_rest_obs": aniso(shells["doppler_coe_rest_obs"]),
        "shells_doppler_comoving_obs": aniso(shells["doppler_comoving_obs"]),
        "shell_inner_angular_radius_deg": geo0["alpha_deg"],
        "verdict": (
            "Rest-frame equal T: sky is PERFECTLY ISOTROPIC for both infinite plates "
            "(any offset) and concentric shells (any off-centre D). "
            "Anisotropy appears only when (i) the two emitters have different T, or "
            "(ii) there is relative motion (observer Doppler / recession Doppler)."
        ),
    }

    out = {"plates": plates, "shells": shells, "summary": summary}
    path = os.path.join(HERE, "two_emitters.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print("wrote", path)
    print(json.dumps(summary, indent=2))
    return out


if __name__ == "__main__":
    main()
