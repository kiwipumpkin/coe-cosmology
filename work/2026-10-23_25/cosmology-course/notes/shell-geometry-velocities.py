"""Export absolute velocities from the Newtonian shell-geometry (ptshell) integrator
for the same two R runs as smooth-symmetric-A.py.

UNIT CONVERSION (README)
========================
Integrator units (bubble_core / prof_core / integ_core):
  - length unit = D  (Earth-CoE distance today)
  - time unit   = t0 (age today; dimensionless time tau = t/t0 ends at tau=1)
  - velocity V is therefore dimensionless:  V = v_phys * t0 / D

Physical conversion (matches newt_omb.py describe()):
  Ht  = |V_E| / |X_E|          (dimensionless; ~ H_E * t0)
  age = Ht * TH0               (Gyr; TH0 = 1/H0 in Gyr)
  t0  = age * 3.156e16         (seconds)  = Ht / H0
  D_m = D_Gpc * GPC * 1e3      (metres; GPC=3.0857e22 is 1 Mpc in m)

  v_phys [m/s] = V_dimless * D_m / t0
               = V_dimless * D_m * H0 / Ht
  v_kms        = v_phys / 1000
  beta         = v_phys / c      (c = 2.99792458e8 m/s)

CAVEAT: for cosmologically large D the Newtonian |V_E| ~ H0*D can be
superluminal (beta > 1). That is a known feature of this Newtonian model
(see universe_thin / newt_omb VO_over_c), not a conversion bug. Use beta
only as the absolute scale for an SR temperature curve construction; do
not interpret beta>1 as a physical Lorentz factor of the Newtonian run.
"""
import json, os, sys, csv
import numpy as np

NOTES = os.path.dirname(os.path.abspath(__file__))
BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"
sys.path.insert(0, BUB)
cwd = os.getcwd()
os.chdir(BUB)
import ptshell  # noqa: F401  -- patches P.tables / P.run for "ptshell"
import prof_core as P
import bubble_core as BC
from bubble_core import geometry, H0, GPC, TH0
os.chdir(cwd)

# Fine 0.25 deg sky grid (psi from anti-CoE; psi=180 points at CoE)
FINE = np.arange(0, 360, 0.25)
P.TH1 = FINE
BC.TH1 = FINE

# Fixed parameters from migkas-geometry-fit.json wls_fit (must match smooth-symmetric-A.py)
fit = json.load(open(os.path.join(NOTES, "migkas-geometry-fit.json")))["wls_fit"]
D = float(fit["D"])          # Gpc
mcc = float(fit["mcc"])
Rin = float(fit["Rin"])      # Gpc
Rout = float(fit["Rout"])    # Gpc
OS, LC, MC, OB = 0.049, 0.005, 0.0, 0.0

C_MS = 2.99792458e8          # m/s
C_KMS = C_MS / 1000.0
GM = GPC * 1e3               # 1 Gpc in metres
YR_S = 3.156e16              # seconds per Gyr (codebase convention)
TAUI, NSTEPS = 0.01, 200

R_VALUES = (3.0, 14.30048287260589)


def run_one(R):
    """Run ptshell at observation distance R (Gpc); return absolute velocities."""
    DR = D / R
    a = dict(DR=DR, Os=OS, mcc=mcc, rin=Rin / D, rout=Rout / D, Lc=LC, mc=MC, Ob=OB)
    os.chdir(BUB)
    rr = P.run("ptshell", a, 0.0, 0.0)

    # Re-integrate once with the converged (mu, mub, Rb) to recover full V vectors.
    # P.run does not return X,V; same integrator / tables / geometry as inside run.
    q, rho = geometry(FINE, DR)
    Rb = rr["Rb"]
    _S, MT = P.tables("ptshell", a, Rb)
    _X, V, _VR = P.integrate(q, rr["mu"], P.RG, MT, rr["mub"], Rb, TAUI, NSTEPS)
    os.chdir(cwd)

    # Sanity: af from captured V should match rr["af"]
    vrel = V[1:] - V[0]
    af_chk = np.hypot(vrel[:, 0], vrel[:, 1]) / rho
    Ht = float(rr["Ht"])
    age = float(rr["age"])
    Dm = D * GM                              # metres
    # Conversion factor: v_phys [m/s] = V_dimless * factor_ms
    # factor_ms = Dm / t0 = Dm * H0 / Ht
    t0_s = age * YR_S
    factor_ms = Dm / t0_s
    factor_ms_alt = Dm * H0 / Ht
    assert abs(factor_ms - factor_ms_alt) / factor_ms < 1e-12

    VE = V[0].copy()                         # dimensionless 2-vector
    Vs = V[1:].copy()                        # (N, 2)
    VE_kms = VE * factor_ms / 1000.0
    Vs_kms = Vs * factor_ms / 1000.0
    beta_E = VE * factor_ms / C_MS
    beta_s = Vs * factor_ms / C_MS

    af = np.asarray(rr["af"], dtype=float)
    A_Aall = af / af.mean()
    af_err = float(np.max(np.abs(af_chk - af) / np.maximum(af, 1e-30)))

    VE_kms_mag = float(np.hypot(VE_kms[0], VE_kms[1]))
    Vs_kms_mag = np.hypot(Vs_kms[:, 0], Vs_kms[:, 1])
    beta_E_mag = float(np.hypot(beta_E[0], beta_E[1]))
    beta_s_mag = np.hypot(beta_s[:, 0], beta_s[:, 1])

    meta = dict(
        R_Gpc=R, D_Gpc=D, DR=DR, Rin_Gpc=Rin, Rout_Gpc=Rout,
        mcc=mcc, Os=OS, Lc=LC, mc=MC, Ob=OB,
        rin=Rin / D, rout=Rout / D,
        Ht=Ht, age_Gyr=age, vrmin=float(rr["vrmin"]),
        rE=float(rr["rE"]), vE_dimless=float(rr["vE"]),
        mu=float(rr["mu"]), mub=float(rr["mub"]), Rb=float(rr["Rb"]),
        D_m=Dm, t0_s=t0_s, H0_SI=float(H0), TH0_Gyr=float(TH0),
        conversion="v_phys[m/s] = V_dimless * D_m / t0_s = V_dimless * D_m * H0 / Ht",
        factor_ms=float(factor_ms),
        factor_kms=float(factor_ms / 1000.0),
        c_ms=C_MS,
        VE_dimless=[float(VE[0]), float(VE[1])],
        VE_kms=[float(VE_kms[0]), float(VE_kms[1])],
        VE_kms_mag=VE_kms_mag,
        beta_E=[float(beta_E[0]), float(beta_E[1])],
        beta_E_mag=beta_E_mag,
        mean_Vs_kms=float(Vs_kms_mag.mean()),
        mean_beta_s=float(beta_s_mag.mean()),
        max_beta_s=float(beta_s_mag.max()),
        mean_af=float(af.mean()),
        A_Aall_max=float(A_Aall.max()),
        A_Aall_min=float(A_Aall.min()),
        af_reintegrate_max_rel_err=af_err,
        note_beta=(
            "Newtonian |beta| can exceed 1 for large D (VO ~ H0*D); "
            "documented in newt_omb VO_over_c."
        ),
    )
    arrays = dict(
        psi=FINE.copy(),
        Vx_s=Vs_kms[:, 0], Vy_s=Vs_kms[:, 1], Vs_mag=Vs_kms_mag,
        Vx_E=float(VE_kms[0]), Vy_E=float(VE_kms[1]), VE_mag=VE_kms_mag,
        af=af, A_Aall=A_Aall,
        beta_sx=beta_s[:, 0], beta_sy=beta_s[:, 1], beta_s_mag=beta_s_mag,
        beta_Ex=float(beta_E[0]), beta_Ey=float(beta_E[1]), beta_E_mag=beta_E_mag,
        Vx_s_dimless=Vs[:, 0], Vy_s_dimless=Vs[:, 1],
        Vx_E_dimless=float(VE[0]), Vy_E_dimless=float(VE[1]),
    )
    return meta, arrays


def extrema_report(A, psi):
    lines = []
    n = len(A)
    for k in range(n):
        a0, a1, a2 = A[k - 1], A[k], A[(k + 1) % n]
        if a1 > a0 and a1 >= a2:
            lines.append(("max", float(psi[k]), float(a1)))
        if a1 < a0 and a1 <= a2:
            lines.append(("min", float(psi[k]), float(a1)))
    return lines


def main():
    os.chdir(NOTES)
    summary = dict(
        conversion_formula=(
            "v_phys[m/s] = V_dimless * D_m / t0_s "
            "with D_m = D_Gpc * GPC * 1e3, t0_s = age_Gyr * 3.156e16, "
            "equivalently V_dimless * D_m * H0 / Ht; beta = v_phys / c"
        ),
        constants=dict(
            H0_km_s_Mpc=70.0, H0_SI=float(H0), GPC_Mpc_m=float(GPC),
            Gpc_m=float(GM), c_ms=C_MS, YR_S=YR_S, TH0_Gyr=float(TH0),
        ),
        fixed=dict(
            D_Gpc=D, mcc=mcc, Rin_Gpc=Rin, Rout_Gpc=Rout,
            Os=OS, Lc=LC, mc=MC, Ob=OB,
        ),
        runs={},
    )
    npz = {}
    print("UNIT CONVERSION: v[m/s] = V_dimless * D_m * H0 / Ht  (= V_dimless * D_m / t0)")
    print(f"  D = {D:.6f} Gpc,  H0 = 70 km/s/Mpc,  c = {C_KMS:.3f} km/s")
    print()

    for R in R_VALUES:
        meta, arr = run_one(R)
        tag = f"R{R:.1f}"
        key = f"R{R:.1f}"
        summary["runs"][key] = meta

        npz[f"{tag}_psi"] = arr["psi"]
        npz[f"{tag}_Vx_s"] = arr["Vx_s"]
        npz[f"{tag}_Vy_s"] = arr["Vy_s"]
        npz[f"{tag}_Vs"] = arr["Vs_mag"]
        npz[f"{tag}_Vx_E"] = np.array([arr["Vx_E"]])
        npz[f"{tag}_Vy_E"] = np.array([arr["Vy_E"]])
        npz[f"{tag}_VE"] = np.array([arr["VE_mag"]])
        npz[f"{tag}_af"] = arr["af"]
        npz[f"{tag}_A_Aall"] = arr["A_Aall"]
        npz[f"{tag}_beta_sx"] = arr["beta_sx"]
        npz[f"{tag}_beta_sy"] = arr["beta_sy"]
        npz[f"{tag}_beta_s"] = arr["beta_s_mag"]
        npz[f"{tag}_beta_Ex"] = np.array([arr["beta_Ex"]])
        npz[f"{tag}_beta_Ey"] = np.array([arr["beta_Ey"]])
        npz[f"{tag}_beta_E"] = np.array([arr["beta_E_mag"]])
        npz[f"{tag}_Vx_s_dimless"] = arr["Vx_s_dimless"]
        npz[f"{tag}_Vy_s_dimless"] = arr["Vy_s_dimless"]

        csv_name = f"shell-geometry-velocities-R{R:.1f}.csv"
        sidecar = {
            "V_E_kms": meta["VE_kms"],
            "V_E_kms_mag": meta["VE_kms_mag"],
            "beta_E": meta["beta_E"],
            "beta_E_mag": meta["beta_E_mag"],
            "fixed_parameters": summary["fixed"],
            "R_Gpc": R, "DR": meta["DR"], "Ht": meta["Ht"],
            "age_Gyr": meta["age_Gyr"], "vrmin": meta["vrmin"],
            "conversion": meta["conversion"],
            "factor_kms": meta["factor_kms"],
            "factor_ms": meta["factor_ms"],
            "caveat": meta["note_beta"],
        }
        with open(csv_name, "w", newline="") as f:
            f.write("# shell-geometry absolute velocities; psi from anti-CoE (psi=180 at CoE)\n")
            f.write("# " + json.dumps(sidecar) + "\n")
            w = csv.writer(f)
            w.writerow([
                "psi_deg", "Vx_s_kms", "Vy_s_kms", "Vs_kms",
                "beta_sx", "beta_sy", "beta_s", "af", "A_over_Aall",
            ])
            for j, p in enumerate(arr["psi"]):
                w.writerow([
                    f"{p:.2f}",
                    f"{arr['Vx_s'][j]:.8e}", f"{arr['Vy_s'][j]:.8e}",
                    f"{arr['Vs_mag'][j]:.8e}",
                    f"{arr['beta_sx'][j]:.8e}", f"{arr['beta_sy'][j]:.8e}",
                    f"{arr['beta_s_mag'][j]:.8e}",
                    f"{arr['af'][j]:.8e}", f"{arr['A_Aall'][j]:.8e}",
                ])
        with open(csv_name.replace(".csv", ".meta.json"), "w") as f:
            json.dump(sidecar, f, indent=2)

        print(f"R = {R:.6f} Gpc   DR = {meta['DR']:.6f}")
        print(f"  |V_E| = {meta['VE_kms_mag']:.6e} km/s   |beta_E| = {meta['beta_E_mag']:.6f}")
        print(f"  max|beta_s| = {meta['max_beta_s']:.6f}   mean|V_s| = {meta['mean_Vs_kms']:.6e} km/s")
        print(f"  mean a_f = {meta['mean_af']:.6e}   Ht = {meta['Ht']:.6f}   age = {meta['age_Gyr']:.4f} Gyr")
        print(f"  A/A_all max = {meta['A_Aall_max']:.4f}  min = {meta['A_Aall_min']:.4f}")
        print(f"  factor = {meta['factor_kms']:.6e} km/s per dimless V")
        print(f"  af reintegrate max rel err = {meta['af_reintegrate_max_rel_err']:.2e}")
        for kind, ang, val in extrema_report(arr["A_Aall"], arr["psi"]):
            print(f"  {kind} at psi={ang:.2f}  A/A_all={val:.4f}")
        print()

    np.savez(os.path.join(NOTES, "shell-geometry-velocities.npz"), **npz)
    with open(os.path.join(NOTES, "shell-geometry-velocities-summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("Wrote:")
    print(" ", os.path.join(NOTES, "shell-geometry-velocities.py"))
    print(" ", os.path.join(NOTES, "shell-geometry-velocities.npz"))
    print(" ", os.path.join(NOTES, "shell-geometry-velocities-summary.json"))
    for R in R_VALUES:
        print(" ", os.path.join(NOTES, f"shell-geometry-velocities-R{R:.1f}.csv"))
        print(" ", os.path.join(NOTES, f"shell-geometry-velocities-R{R:.1f}.meta.json"))


if __name__ == "__main__":
    main()
