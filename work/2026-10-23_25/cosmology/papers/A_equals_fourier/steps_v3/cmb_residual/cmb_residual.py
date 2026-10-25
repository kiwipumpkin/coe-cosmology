"""CMB temperature anisotropy of the CoE toy cosmology (core + uniform Sedov-Taylor shell, Lambda=0, no DM).
Linear (weak-field) Sachs-Wolfe / Doppler / ISW split along null rays, axisymmetric about the CoE axis.
theta = angle between the line of sight and the direction to the CoE (theta=0 looks at the CoE).

Model ingredients (taken from files, not refitted):
  - D, Mc [kg], Os, R_in, R_out from steps_v3/dscan/dscan.json (headline points: free D~94, omb D~49, rel/LTB D~33).
  - <H> (model's own sky-averaged H) from steps_v3/notebook_results.json["anchor"] -> background expansion rate today.
  - R = 14.3 Gpc (CMB emission sphere, comoving, centred on us), emission at scale a_e = 1/1101 (T=3000 K).
Potential: inside the shell the uniform shell part is quadratic in r (-> only l<=1 on the sky, absorbed by the
homologous V = k r background, Taylor 1950 / Sedov 1959). l>=2 comes only from the 1/r pieces:
  Phi_pert = -G Meff / r_phys,  Meff = Mc (core only)  or  Mc - Mhole (core + the empty interior r<R_in of the shell).
In the homologous background Phi_pert at fixed comoving x decays as 1/a (Mc, Mhole conserved) -> ISW.
"""
import json, os, numpy as np
from numpy.polynomial.legendre import leggauss, legval
from scipy.integrate import cumulative_trapezoid
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); S3 = os.path.dirname(HERE)
G = 6.67430e-11; c = 2.99792458e8; GPC = 3.0856775814913673e25; MPC = GPC / 1e3; T0 = 2.7255
R_GPC = 14.3; A_E = 1 / 1101.0
H70 = 70e3 / MPC; RHOC70 = 3 * H70**2 / (8 * np.pi * G)          # dscan shell-mass convention (Os x rho_c(70))

dscan = json.load(open(os.path.join(S3, "dscan", "dscan.json")))
anchor = json.load(open(os.path.join(S3, "notebook_results.json")))["anchor"]
HEAD = {"free": 6.553635010401956, "omb": 3.4386644382849716, "rel": 2.3034233550510663}   # D/R of the headline fits
LABEL = {"free": "Newtonian, Os free (D~94 Gpc)", "omb": "Newtonian, Os=0.049 (D~49 Gpc)", "rel": "relativistic/LTB, Os=0.049 (D~33 Gpc)"}

def point(fit, DR):
    for p in dscan["points"][fit]:
        if abs(p["DR"] - DR) < 1e-6: return p
    raise KeyError((fit, DR))

def background(Om, H0, a):
    """H(a) for Lambda=0 dust (Om=0 -> coasting a ~ t)."""
    return H0 * np.sqrt(Om * a**-3 + (1 - Om) * a**-2)

def rays(p, H0, meff, null="coast", a_i=1e-4, mu=None, na=6000):
    """Return dict of Delta T / T components on directions mu = cos(theta)."""
    D = p["D"] * GPC; R = R_GPC * GPC
    Om = 0.0 if null == "coast" else p["Os"]
    a = np.geomspace(A_E, 1.0, na)                                   # along the ray, a_e -> 1
    Hs = background(Om, H0, a)
    # null-ray shape chi(a) from the background, normalised so chi(a_e) = R (user spec: emission sphere radius R at a = 1/1101)
    shape = cumulative_trapezoid((1 / (a**2 * Hs))[::-1], a[::-1], initial=0.0)[::-1] * -1   # int_a^1 da/(a^2 H)
    chi = R * shape / shape[0]
    # conformal-time integral for peculiar velocities (physical dynamics, real H): eta(a) = int_{a_i}^a da/(a^2 H)
    ag = np.geomspace(min(a_i, A_E) * 0.999, 1.0, 20000); eg = cumulative_trapezoid(1 / (ag**2 * background(Om, H0, ag)), ag, initial=0.0)
    ei = np.interp(a_i, ag, eg)
    deta = lambda aa: np.maximum(np.interp(aa, ag, eg) - ei, 0.0)   # zero before the homologous (pure V=kr) start a_i
    GM = G * meff
    out = {}
    X = lambda m, ch: np.sqrt(D**2 + ch**2 - 2 * D * ch * m)        # comoving distance from CoE (theta from CoE direction)
    Phi0 = lambda m, ch: -GM / X(m, ch)                               # today's-amplitude potential at comoving point
    m = mu[:, None]
    # SW (emission; full Phi: Laplace-exterior pattern with Phi ~ 1/a has zero Newtonian-gauge density contrast -> no intrinsic term)
    sw = Phi0(mu, R) / A_E / c**2
    sw_o = Phi0(np.array([1.0]), 0.0)[0] / c**2                       # observer potential (monopole)
    # ISW: 2/c^2 int dPhi/dt dt = 2/c^2 int Phi0(x(a)) d(1/a)
    P = Phi0(m, chi[None, :]); inv = 1 / a
    isw = 2 / c**2 * np.trapezoid(P, inv[None, :], axis=1)            # integrate in the direction of increasing t (a_e -> 1)
    # Doppler: u = -grad_x Phi0 * deta(a)/a (linear growth in the decaying external field); dT/T = n.(u_o - u_e)/c
    def n_dot_grad(mm, ch):                                           # n . grad Phi0 at comoving point D z + ch n, n=(sin, 0, -cos)
        x = X(mm, ch); return GM * (ch - D * mm) / x**3
    u_e = -n_dot_grad(mu, R) * deta(A_E) / A_E
    u_o = -n_dot_grad(mu, 0.0) * deta(1.0)
    dop = (u_o - u_e) / c
    out = dict(SW=sw - sw_o, ISW=isw, Doppler=dop)
    out["total"] = out["SW"] + out["ISW"] + out["Doppler"]
    # charitable Doppler: emitters on the R sphere carrying TODAY's model peculiar-velocity pattern (the pattern the A(angle) fits see)
    out["Doppler_today"] = (u_o + n_dot_grad(mu, R) * deta(1.0)) / c
    out["SW_frozen"] = (Phi0(mu, R) - Phi0(np.array([1.0]), 0.0)[0]) / c**2   # charitable: potential frozen at today's amplitude, no growth/ISW
    return out

def tidal_dH(p, H0, meff, null, a_i):
    """linear-theory H_radial - H_perp (w.r.t. CoE axis) today at Earth from the tidal tensor of Meff: 3 G M / D^3 * eta(1)"""
    Om = 0.0 if null == "coast" else p["Os"]; ag = np.geomspace(a_i, 1.0, 20000)
    eta = np.trapezoid(1 / (ag**2 * background(Om, H0, ag)), ag)
    return float(3 * G * meff / (p["D"] * GPC)**3 * eta * MPC / 1e3)

NG = 256; MU, W = leggauss(NG)
def legendre(f, lmax=6):
    return np.array([(2 * l + 1) / 2 * np.sum(W * f * legval(MU, [0] * l + [1])) for l in range(lmax + 1)])

def analyse(fit, DR, meff_mode, null, a_i):
    p = point(fit, DR); H0 = anchor[fit]["H"] * 1e3 / MPC
    Mhole = 4 * np.pi / 3 * p["Os"] * RHOC70 * (p["R_in"] * GPC)**3
    meff = p["Mc"] if meff_mode == "core" else p["Mc"] - Mhole
    comp = rays(p, H0, meff, null, a_i, MU)
    res = {k: legendre(v) for k, v in comp.items()}
    eps = G * meff / (c**2 * p["D"] * GPC)
    return p, H0, Mhole, meff, eps, comp, res

def fmt(v): return float(v)

if __name__ == "__main__":
    out = {"meta": {"R_Gpc": R_GPC, "a_e": A_E, "T0_K": T0, "theta": "angle from CoE direction; a_l = Legendre coeffs of dT/T(theta)",
                    "H_source": "notebook_results.json anchor <H>", "params_source": "dscan/dscan.json points",
                    "variants": "null=coast (Om=0 ray shape & dynamics) or decel (Friedmann with Os); a_i = start of pure V=kr flow (peculiar velocity zero before)"},
           "headline": {}, "dscan": {}}
    rows = []
    for fit, DR in HEAD.items():
        out["headline"][fit] = {}
        for mm in ("core", "core+hole"):
            for null in ("coast", "decel"):
                for a_i in (1e-4, 1e-2):
                    p, H0, Mh, meff, eps, comp, res = analyse(fit, DR, mm, null, a_i)
                    key = f"{mm}|{null}|ai={a_i:g}"
                    d = {"D_Gpc": p["D"], "Mc_kg": p["Mc"], "Mhole_kg": Mh, "Meff_kg": meff, "Os": p["Os"], "R_in": p["R_in"], "R_out": p["R_out"],
                         "H_model": anchor[fit]["H"], "VO_c_file": p["VO_c"], "R2_file": p["R2"], "eps_GMeff_c2D": eps, "R_over_D": R_GPC / p["D"],
                         "H_file_r_perp": [p.get("H_r"), p.get("H_perp")], "linear_tidal_dH_r_minus_perp_kms_Mpc": tidal_dH(p, H0, meff, null, a_i)}
                    for k, v in res.items():
                        d[k] = {f"a{l}": fmt(v[l]) for l in range(7)}
                        d[k + "_uK"] = {f"a{l}": fmt(v[l] * T0 * 1e6) for l in range(5)}
                    out["headline"][fit][key] = d
                    for k, v in res.items():
                        rows.append([fit, mm, null, a_i, k] + [v[l] for l in range(7)])
    # D scan (coasting null shape, a_i = 1e-4), all dscan points
    for fit, pts in dscan["points"].items():
        out["dscan"][fit] = []
        for p in pts:
            e = {"DR": p["DR"], "D": p["D"], "mcc": p["mcc"], "Mc": p["Mc"], "R2": p["R2"]}
            for mm in ("core", "core+hole"):
                _, H0, Mh, meff, eps, comp, res = analyse(fit, p["DR"], mm, "coast", 1e-4)
                e[mm] = {"eps": eps, "frozen_a2": fmt(res["SW_frozen"][2]), "frozen_a3": fmt(res["SW_frozen"][3]), "frozen_a4": fmt(res["SW_frozen"][4]),
                         "total_a2": fmt(res["total"][2]), "total_a3": fmt(res["total"][3]), "total_a4": fmt(res["total"][4]),
                         "Dop_a2": fmt(res["Doppler"][2]), "Dop_today_a2": fmt(res["Doppler_today"][2]), "tidal_pred_a2": eps * (R_GPC / p["D"])**2}
            out["dscan"][fit].append(e)
    json.dump(out, open(os.path.join(HERE, "cmb_residual.json"), "w"), indent=1)
    with open(os.path.join(HERE, "cmb_residual_legendre.csv"), "w") as f:
        f.write("fit,mass,null_shape,a_i,component," + ",".join(f"a{l}_dT/T" for l in range(7)) + "\n")
        for r in rows: f.write(",".join(str(x) if not isinstance(x, float) else f"{x:.6e}" for x in r) + "\n")
    # theta profiles (headline, core-only and core+hole, coasting, a_i=1e-4)
    th = np.linspace(0, np.pi, 361); mu = np.cos(th)
    with open(os.path.join(HERE, "cmb_residual_profiles.csv"), "w") as f:
        f.write("fit,mass,theta_deg,SW,ISW,Doppler,total,residual_l>=2,SW_frozen,SW_frozen_residual_l>=2\n")
        prof = {}
        for fit, DR in HEAD.items():
            for mm in ("core", "core+hole"):
                p = point(fit, DR); H0 = anchor[fit]["H"] * 1e3 / MPC
                Mh = 4 * np.pi / 3 * p["Os"] * RHOC70 * (p["R_in"] * GPC)**3
                meff = p["Mc"] if mm == "core" else p["Mc"] - Mh
                comp = rays(p, H0, meff, "coast", 1e-4, mu)
                L = legendre_tot = legendre(rays(p, H0, meff, "coast", 1e-4, MU)["total"])
                Lf = legendre(rays(p, H0, meff, "coast", 1e-4, MU)["SW_frozen"])
                resid = comp["total"] - L[0] - L[1] * mu
                residf = comp["SW_frozen"] - Lf[0] - Lf[1] * mu
                prof[(fit, mm)] = (comp, resid, residf)
                for i in range(len(th)):
                    f.write(f"{fit},{mm},{np.degrees(th[i]):.1f},{comp['SW'][i]:.6e},{comp['ISW'][i]:.6e},{comp['Doppler'][i]:.6e},{comp['total'][i]:.6e},{resid[i]:.6e},{comp['SW_frozen'][i]:.6e},{residf[i]:.6e}\n")
    # figure
    fig, axs = plt.subplots(3, 3, figsize=(13, 10.5))
    deg = np.degrees(th)
    for j, fit in enumerate(HEAD):
        comp, resid, residf = prof[(fit, "core")]
        ax = axs[0, j]
        for k, col in (("SW", "C0"), ("ISW", "C1"), ("Doppler", "C2"), ("total", "k")):
            ax.plot(deg, comp[k], color=col, lw=1.4 if k == "total" else 1.1, label=k)
        ax.set_yscale("symlog", linthresh=1e-3); ax.set_title(LABEL[fit] + "\ncore only, linear decaying-potential calc", fontsize=9)
        ax.set_ylabel("ΔT/T"); ax.legend(fontsize=7, frameon=False)
        ax = axs[1, j]
        ax.plot(deg, resid, "k", lw=1.4, label="total, monopole+dipole removed (core only)")
        comp2, resid2, residf2 = prof[(fit, "core+hole")]
        ax.plot(deg, resid2, "m--", lw=1.1, label="same, core + shell hole")
        ax.set_yscale("symlog", linthresh=1e-3); ax.set_ylabel("residual ΔT/T (ℓ≥2)"); ax.legend(fontsize=7, frameon=False)
        ax.set_title("residual after removing ℓ=0,1", fontsize=9)
        ax = axs[2, j]
        ax.plot(deg, residf, "C3", lw=1.4, label="frozen-potential SW only, core (charitable bound)")
        ax.plot(deg, residf2, "C4--", lw=1.1, label="frozen-potential SW only, core + hole")
        ax.axhspan(-1e-5, 1e-5, color="g", alpha=0.35, label="±1e-5 target (±27 μK)")
        ax.set_yscale("symlog", linthresh=1e-6); ax.set_ylabel("residual ΔT/T (ℓ≥2)"); ax.legend(fontsize=7, frameon=False)
        ax.set_xlabel("θ from CoE direction [deg]"); ax.set_title("charitable lower bound vs target", fontsize=9)
        for a_ in axs[:, j]: a_.set_xlim(0, 180); a_.set_xticks(range(0, 181, 30)); a_.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "cmb_residual.png"), dpi=130)
    # D-scan figure
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for fit, col in (("free", "C3"), ("omb", "C2"), ("rel", "C1")):
        E = out["dscan"][fit]; Dv = [e["D"] for e in E]
        ax[0].loglog(Dv, [abs(e["core"]["frozen_a2"]) for e in E], "o-", color=col, label=f"{fit}: |a2| core, frozen SW")
        ax[0].loglog(Dv, [abs(e["core+hole"]["frozen_a2"]) for e in E], "s--", color=col, mfc="none", label=f"{fit}: |a2| core+hole, frozen SW")
        ax[1].loglog(Dv, [abs(e["core"]["total_a2"]) for e in E], "o-", color=col, label=f"{fit}: |a2| core, linear total")
    for a_ in ax: a_.axhline(1e-5, color="g", lw=1.2, label="1e-5 target"); a_.set_xlabel("D [Gpc]"); a_.grid(alpha=0.3, which="both")
    ax[0].set_ylabel("|a2| (ΔT/T)"); ax[0].legend(fontsize=6.5, frameon=False); ax[1].legend(fontsize=7, frameon=False)
    ax[0].set_title("quadrupole vs D along the dscan best fits (charitable)", fontsize=9); ax[1].set_title("quadrupole vs D, linear decaying-potential calc", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "cmb_residual_dscan.png"), dpi=130)
    # print summary
    for fit in HEAD:
        for key, d in out["headline"][fit].items():
            if not key.endswith("ai=0.0001") and "decel" not in key: pass
            print(fit, key, "eps=%.3g R/D=%.3f" % (d["eps_GMeff_c2D"], d["R_over_D"]))
            print("   linear dH_r-dH_perp [km/s/Mpc] =", d["linear_tidal_dH_r_minus_perp_kms_Mpc"], "file H_r,H_perp", d["H_file_r_perp"])
            for k in ("SW", "ISW", "Doppler", "total", "SW_frozen", "Doppler_today"):
                print("   %-9s" % k, " ".join("a%d=%+.3e" % (l, d[k]["a%d" % l]) for l in range(5)))
