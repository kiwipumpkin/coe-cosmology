"""CMB ΔT/T along Chris's recombination track (variable R_see, NOT fixed 14.3 Gpc).

Track (coasting homologous + continuous T∝1/a cooling to T_rec=3000 K):
  f = t_rec/t_now = D_rec/D_now = T0/T_rec
  Δt = t_now − t_rec
  R_see = c × Δt
  D_now − D_rec = V_shell × Δt  with V_shell = D_now/t_now  (lock: R_see/(D_now−D_rec)=c/V_shell)

t_now = age from the Newtonian integrator (Ht×TH0). Gravity only, Λ=0, Ω_s=0.049 shell
with R_in=D−KIN·R_see, R_out=D+KOUT·R_see. Compact core Mc re-solved for cluster a2≈0.084.
Same Doppler convention as cycle4: ΔT/T = −n·(V_emit−V_obs)/c; kill = −(Doppler−⟨Doppler⟩).
"""
import json, os, sys
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import lib4 as L

T = L.T
P = L.P
HM = L.HM
CKM = L.CKM
TH1 = L.TH1
KIN, KOUT = float(T.KIN), float(T.KOUT)
OS = 0.049
RG = L.RG
TH0 = float(P.TH0)           # Gyr age scale
C_GPC_GYR = CKM * 1e3 / 3.0857e22 * 3.15576e16 / 1e9 * 1e9  # messy; use exact below
# c [km/s] * 1 Gyr [s] / (1 Gpc [km]) with GPC_km = 3.0857e22
GPC_KM = 3.0857e22
GYR_S = 3.15576e16
C_GPC_PER_GYR = CKM * GYR_S / GPC_KM   # ≈ 0.3066 Gpc/Gyr

T_REC = 3000.0
T0 = 2.725
F_REC = T0 / T_REC          # t_rec/t_now = D_rec/D_now
RS_CL = 0.31
TARGET_A2 = 0.084
PASS = 1e-5

# D_now grid: from near shell-min upward, including headline-like and static-null-ish
D_GRID = [4.5, 5.0, 6.0, 8.0, 10.0, 12.0, 14.7, 18.0, 22.28, 30.0, 40.0, 49.17, 70.0, 94.0]


def build_var(D, mcc, R_see, kind="point", a=None, Os=OS):
    """Density table with shell edges set by R_see (not fixed RV)."""
    a = a or {}
    rc = 1.0 - R_see / D
    rin = 1.0 - KIN * R_see / D
    rout = 1.0 + KOUT * R_see / D
    s_in = L.shape_mass(kind, max(rc, 1e-4), a)
    s = mcc * s_in + np.where((RG > rin) & (RG < rout), Os * 3 / (4 * np.pi), 0.0)
    m = cumulative_trapezoid(4 * np.pi * RG**2 * s, RG, initial=0.0)
    m1 = np.interp(1.0, RG, m)
    return s / m1, m / m1, m1, rin, rout


def run_var(D, mcc, R_see, Rs_list, kind="point", a=None, taui=0.01, nsteps=200, iters=3):
    a = a or {}
    S, MT, m1, rin, rout = build_var(D, mcc, R_see, kind, a)
    DR0 = D / min(Rs_list)
    Rb = max(P.rb_of(DR0, 1.5), 1.02 * rout)
    out = []
    Ht = 0.95
    qs = [P.geometry(TH1, D / Rs) for Rs in Rs_list]
    for _ in range(iters):
        mu = m1 * Ht**2 / 2
        X, V, _ = P.integrate(qs[0][0], mu, RG, MT, 0.0, Rb, taui, nsteps)
        Ht = np.hypot(*V[0]) / np.hypot(*X[0])
    for (q, rho) in qs:
        X, V, VR = P.integrate(q, mu, RG, MT, 0.0, Rb, taui, nsteps)
        vrel = V[1:] - V[0]
        af = np.hypot(vrel[:, 0], vrel[:, 1]) / rho
        out.append(dict(af=af, X=X, V=V, vrmin=float(VR.min())))
    age = float(Ht * TH0)
    return out, dict(mu=mu, MT=MT, S=S, Ht=float(Ht), Rb=Rb, m1=m1, rin=rin, rout=rout, age=age)


def meff_var(D, mcc, R_see, kind="point", Os=OS):
    rin = 1.0 - KIN * R_see / D
    Mh = Os * rin**3
    Min = mcc
    return Min - Mh, Mh, rin


def static_pattern_var(D, Me, R_see, H=HM):
    th = np.radians(TH1)
    rs = R_see / D
    rr = np.hypot(1 + rs * np.cos(th), rs * np.sin(th))
    return -0.5 * Me * (H * D * 1e3 / CKM)**2 * (1 / rr - 1)


def track_from_D(D, age_gyr):
    """Given today's D and age, return recombination-track kinematics."""
    t_now = age_gyr
    t_rec = F_REC * t_now
    dt = t_now - t_rec
    R_see = C_GPC_PER_GYR * dt
    D_rec = F_REC * D
    V_over_c = D / (C_GPC_PER_GYR * t_now)   # homologous V_shell/c = D/(c t)
    R_lock = (1.0 / V_over_c) * (D - D_rec)  # (c/V)*(D−D_rec); should ≈ R_see
    D_min = KIN * R_see
    return dict(
        t_now=t_now, t_rec=t_rec, dt=dt, R_see=R_see, D_rec=D_rec,
        V_shell_over_c=V_over_c, R_see_lock=R_lock, D_min=D_min,
        below_shell_min=bool(D < D_min - 1e-9),
        f_rec=F_REC, frac_visible=R_see / (R_see + D),
    )


def mcc_for(D, R_see, target=TARGET_A2, Rs=RS_CL):
    def f(lm):
        af = run_var(D, 10**lm, R_see, [Rs])[0][0]["af"]
        return np.log(T.dq(af)["a2"] / target)
    return 10 ** brentq(f, -3.5, -2.0, xtol=1e-6)


def curves(D, mcc, R_see):
    res, info = run_var(D, mcc, R_see, [RS_CL, R_see])
    cl, far = res
    Hconv = HM / cl["af"].mean()
    vunit = Hconv * D * 1e3
    X, V = far["X"], far["V"]
    d = X[1:] - X[0]
    n = d / np.hypot(d[:, 0], d[:, 1])[:, None]
    kin = -np.sum(n * (V[1:] - V[0]), axis=1) * vunit / CKM
    Me, Mh, rin = meff_var(D, mcc, R_see)
    stat = static_pattern_var(D, Me, R_see)
    tot = kin + stat
    HRC = HM * R_see * 1e3 / CKM
    kin_frac = kin / HRC
    i_toward = np.arange(180, -1, -1)
    th = np.arange(0.0, 181.0, 1.0)
    return dict(
        theta_toward_deg=th.tolist(),
        dtt_static=stat[i_toward].tolist(),
        dtt_doppler=kin[i_toward].tolist(),
        dtt_total=tot[i_toward].tolist(),
        dtt_doppler_frac=kin_frac[i_toward].tolist(),
        kin_code=kin.tolist(),
        stat_code=stat.tolist(),
        Meff=float(Me), Mhole=float(Mh), rin=float(rin),
        Hconv=float(Hconv), HRC=float(HRC),
        vrmin=float(min(cl["vrmin"], far["vrmin"])),
        Ht=float(info["Ht"]), age=float(info["age"]),
        cl_a2=float(T.dq(cl["af"])["a2"]),
    )


def summarize(dop, stat, tot, dop_frac, kin_code, stat_code, HRC):
    kin_leg = L.legendre(np.array(kin_code))
    stat_leg = L.legendre(np.array(stat_code))
    tot_leg = L.legendre(np.array(kin_code) + np.array(stat_code))
    frac_leg = L.legendre(np.array(kin_code) / HRC)
    dop = np.asarray(dop); stat = np.asarray(stat); tot = np.asarray(tot)
    kill = -(dop - dop.mean())
    return dict(
        doppler_ptp=float(dop.max() - dop.min()),
        doppler_peak_abs=float(np.max(np.abs(dop - dop.mean()))),
        doppler_a2=float(kin_leg[2]),
        static_a2=float(stat_leg[2]),
        total_a2=float(tot_leg[2]),
        frac_a2=float(frac_leg[2]),
        kill_a2=float(-kin_leg[2]),
        kill_peak_abs=float(np.max(np.abs(kill))),
        kill_energy_pct_peak=float(100.0 * np.max(np.abs(kill))),
        kill_energy_pct_min=float(100.0 * kill.min()),
        kill_energy_pct_max=float(100.0 * kill.max()),
        kin_legendre=kin_leg,
        stat_legendre=stat_leg,
        tot_legendre=tot_leg,
        passes_cmb=bool(abs(kin_leg[2]) <= PASS and abs(stat_leg[2]) <= PASS),
        passes_cmb_kill_needed=bool(abs(kin_leg[2]) <= PASS),
    )


def bootstrap_age(D, R_see_guess, mcc_guess=0.002):
    """One integrator pass to get age at this D with approximate R_see."""
    _, info = run_var(D, mcc_guess, R_see_guess, [RS_CL])
    return info["age"], info["Ht"]


def main():
    # First-pass R_see from TH0 (refined per point)
    R_see_seed = C_GPC_PER_GYR * TH0 * (1.0 - F_REC)
    rows = []
    for D in D_GRID:
        print("=== D_now=%.4f Gpc ===" % D, flush=True)
        # Iterate age ↔ R_see once or twice (weak dependence)
        R_see = R_see_seed
        age, Ht = bootstrap_age(D, R_see)
        tr = track_from_D(D, age)
        R_see = tr["R_see"]
        # re-age with updated R_see
        age, Ht = bootstrap_age(D, R_see)
        tr = track_from_D(D, age)
        R_see = tr["R_see"]
        print("  age=%.4f Gyr  R_see=%.4f  D_rec=%.6f  D_min=%.4f  below_min=%s"
              % (age, R_see, tr["D_rec"], tr["D_min"], tr["below_shell_min"]), flush=True)
        try:
            mcc = mcc_for(D, R_see)
        except Exception as e:
            print("  mcc solve FAILED:", e, flush=True)
            rows.append(dict(D=D, error=str(e), **tr))
            continue
        print("  mcc=%.6f  running curves..." % mcc, flush=True)
        c = curves(D, mcc, R_see)
        # final track with age from full solve
        tr = track_from_D(D, c["age"])
        # if R_see shifted, note it (do not re-solve unless large)
        R_see_final = tr["R_see"]
        s = summarize(
            c["dtt_doppler"], c["dtt_static"], c["dtt_total"], c["dtt_doppler_frac"],
            c["kin_code"], c["stat_code"], c["HRC"],
        )
        row = dict(
            D=D, D_label=("%.2f" % D).rstrip("0").rstrip("."),
            mcc=mcc,
            R_see=R_see_final,
            R_see_used=R_see,
            D_rec=tr["D_rec"],
            t_now=tr["t_now"], t_rec=tr["t_rec"], dt=tr["dt"],
            V_shell_over_c=tr["V_shell_over_c"],
            R_see_lock=tr["R_see_lock"],
            D_min=tr["D_min"],
            below_shell_min=tr["below_shell_min"],
            frac_visible=tr["frac_visible"],
            f_rec=F_REC,
            Meff=c["Meff"], Mhole=c["Mhole"], rin=c["rin"],
            Hconv=c["Hconv"], HRC=c["HRC"],
            vrmin=c["vrmin"], Ht=c["Ht"], age=c["age"],
            cl_a2=c["cl_a2"],
            theta_toward_deg=c["theta_toward_deg"],
            dtt_static=c["dtt_static"],
            dtt_doppler=c["dtt_doppler"],
            dtt_total=c["dtt_total"],
            dtt_doppler_frac=c["dtt_doppler_frac"],
            **s,
        )
        rows.append(row)
        print(
            "  kill_a2=%+.4e  peak|dE/E|=%.3f%%  cl_a2=%.4f  Meff=%+.3e  vrmin=%.3f  pass=%s"
            % (s["kill_a2"], s["kill_energy_pct_peak"], c["cl_a2"], c["Meff"], c["vrmin"], s["passes_cmb"]),
            flush=True,
        )

    out = dict(
        note=(
            "Chris recombination track: continuous T∝1/a cooling, T_rec=3000 K, T0=2.725 K; "
            "coasting homologous lock R_see=c·Δt=(c/V_shell)·(D_now−D_rec); "
            "Ω_s=0.049 shell with KIN/KOUT rule on R_see; compact core; Newtonian today's field; "
            "ΔT/T_Doppler=−n·(V_emit−V_obs)/c; kill=−(Doppler−mean)."
        ),
        formulas=dict(
            f_rec="T0/T_rec = t_rec/t_now = D_rec/D_now",
            T0=T0, T_rec=T_REC, f_rec_value=F_REC,
            R_see="c * (t_now - t_rec)",
            D_rec="f_rec * D_now",
            V_shell="D_now / t_now",
            lock="R_see / (D_now - D_rec) = c / V_shell",
            D_min="KIN * R_see",
            shell="R_in = D - KIN*R_see, R_out = D + KOUT*R_see",
            KIN=KIN, KOUT=KOUT, Os=OS,
            c_Gpc_per_Gyr=C_GPC_PER_GYR,
            TH0=TH0,
            expansion="coasting homologous (R∝t); Sedov T∝t^{-6/5} would change f but keep R_see≈c·t_now",
        ),
        H_mean=HM,
        cluster_target_a2=TARGET_A2,
        cmb_pass=PASS,
        rows=rows,
    )
    path = os.path.join(HERE, "track.json")
    json.dump(out, open(path, "w"), default=float)
    print("wrote", path)


if __name__ == "__main__":
    main()
