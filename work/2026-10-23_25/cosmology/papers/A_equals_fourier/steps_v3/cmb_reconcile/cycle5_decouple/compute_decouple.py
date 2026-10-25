"""CMB–Hubble decoupling from ONE Φ for Chris's CoE compact-core + shell.

Matter:  Δv = −∫ ∇Φ dt  →  ΔH/H = v_grav·n̂/r  at cluster scale (Rs=0.31 Gpc).
CMB:    ΔT/T = endpoint Φ/c² + (2/c²)∫ Φ̇ dt + Doppler (−n·Δv/c).

Nearly-static Φ limit (Chris's table / Newtonian integrator): Φ̇ ≈ 0 → ISW ≈ 0;
matter still accumulates velocity from ∇Φ; endpoint redshift ∝ M_eff (null at D≈22.28).

Reuses cycle4/lib4.py (same integrator, shell rule, linear Doppler convention).
Mc re-solved per D so cluster a2 = 0.084 (Migkas all-slice quad). D = 14.7, 22.28, 49.
"""
import json, os, sys
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "cycle4"))
import lib4 as L

T = L.T
RV = L.RV
HM = L.HM
CKM = L.CKM
TH1 = L.TH1
HRC = HM * RV * 1e3 / CKM
RS_CL = 0.31
TARGET_A2 = 0.084
PASS_CMB = 1e-5

PARTB = os.path.join(HERE, "..", "cycle4", "partB.json")
_pb = json.load(open(PARTB))
D_NULL = float(_pb["static_null"]["D"])
MULT = _pb["multipliers"]  # growing-mode multipliers on today's static tidal pattern

DS = [14.7, D_NULL, 49.0]
D_LABELS = {14.7: "14.7", D_NULL: "22.28", 49.0: "49"}


def mcc_for(D, target=TARGET_A2, Rs=RS_CL):
    f = lambda lm: np.log(T.dq(L.run(D, 10**lm, "point", {}, [Rs])[0][0]["af"])["a2"] / target)
    return 10 ** brentq(f, -3.5, -2.3, xtol=1e-6)


def toward_index():
    """Code θ=0 = away from CoE → plot θ=0 toward CoE."""
    return np.arange(180, -1, -1)


def demean(y):
    y = np.asarray(y, float)
    return y - y.mean()


def run_one(D, mcc):
    """Same Φ → cluster ΔH/H and CMB ΔT/T components."""
    res, info = L.run(D, mcc, "point", {}, [RS_CL, RV])
    cl, far = res
    Hconv = HM / cl["af"].mean()
    vunit = Hconv * D * 1e3  # km/s per code velocity

    # --- Matter: ΔH/H from relative velocity / separation (cluster scale) ---
    af = cl["af"]
    dHH = af / af.mean() - 1.0  # ΔH/H(θ) code-ordered
    dqc = T.dq(af)

    # --- CMB: Doppler + endpoint static from same run / same Meff ---
    X, V = far["X"], far["V"]
    d = X[1:] - X[0]
    n = d / np.hypot(d[:, 0], d[:, 1])[:, None]
    dop = -np.sum(n * (V[1:] - V[0]), axis=1) * vunit / CKM
    Me = L.meff(D, mcc, "point")
    stat = L.static_pattern(D, Me)  # endpoint (Φ_e − Φ_o)/c², today's field

    # Nearly-static Φ: Φ̇ ≈ 0 → ISW ≈ 0. Growing-mode estimate scales with today's static.
    isw_static = np.zeros_like(stat)  # Chris's preferred limit
    isw_growing = MULT["growing_ISW"] * stat  # open Om=0.049 history (partB)
    sw_growing = MULT["growing_SW_full"] * stat
    net_growing = MULT["growing_net_linear"] * stat

    it = toward_index()
    th = np.arange(0.0, 181.0, 1.0)

    dHH_t = dHH[it]
    dop_t = dop[it]
    stat_t = stat[it]
    isw_g_t = isw_growing[it]
    tot_static_lim = dop + stat  # ISW=0
    tot_growing = dop + net_growing

    # Legendre on code-ordered 0..180 (lib4.legendre)
    leg_dHH = L.legendre(dHH)
    leg_dop = L.legendre(dop)
    leg_stat = L.legendre(stat)
    leg_isw_g = L.legendre(isw_growing)
    leg_tot_s = L.legendre(tot_static_lim)
    leg_tot_g = L.legendre(tot_growing)
    leg_net_g = L.legendre(net_growing)

    rin = 1 - T.KIN * RV / D
    Mh = 0.049 * rin**3

    def peak_abs(y):
        y = demean(y)
        return float(np.max(np.abs(y)))

    return dict(
        D=float(D),
        D_label=D_LABELS[D],
        mcc=float(mcc),
        Mhole=float(Mh),
        Meff=float(Me),
        Hconv=float(Hconv),
        Ht=float(info["Ht"]),
        vrmin=float(min(cl["vrmin"], far["vrmin"])),
        HRC=float(HRC),
        # matter
        cluster_a2=float(dqc["a2"]),
        cluster_a1=float(dqc["a1"]),
        dHH_a2=float(leg_dHH[2]),
        dHH_ptp=float(dHH_t.max() - dHH_t.min()),
        dHH_peak_abs=peak_abs(dHH_t),
        dHH_min=float(dHH_t.min()),
        dHH_max=float(dHH_t.max()),
        # CMB nearly-static
        static_a2=float(leg_stat[2]),
        static_ptp=float(stat_t.max() - stat_t.min()),
        static_peak_abs=peak_abs(stat_t),
        isw_static_a2=0.0,  # Φ̇=0
        doppler_a2=float(leg_dop[2]),
        doppler_ptp=float(dop_t.max() - dop_t.min()),
        doppler_peak_abs=peak_abs(dop_t),
        total_static_lim_a2=float(leg_tot_s[2]),
        total_static_lim_ptp=float((dop_t + stat_t).max() - (dop_t + stat_t).min()),
        total_static_lim_peak_abs=peak_abs(dop_t + stat_t),
        # growing-mode contrast (same Meff, not Chris's static table)
        isw_growing_a2=float(leg_isw_g[2]),
        sw_growing_a2=float(L.legendre(sw_growing)[2]),
        net_growing_a2=float(leg_net_g[2]),
        total_growing_a2=float(leg_tot_g[2]),
        # ratios / verdict helpers
        ratio_doppler_a2_over_dHH_a2=float(abs(leg_dop[2]) / max(abs(leg_dHH[2]), 1e-30)),
        ratio_doppler_a2_over_CMB_pass=float(abs(leg_dop[2]) / PASS_CMB),
        ratio_dHH_a2_over_Migkas=float(dqc["a2"] / TARGET_A2),
        # curves toward-CoE
        theta_toward_deg=th.tolist(),
        dHH=dHH_t.tolist(),
        dtt_static=stat_t.tolist(),
        dtt_isw_static=isw_static[it].tolist(),
        dtt_isw_growing=isw_g_t.tolist(),
        dtt_doppler=dop_t.tolist(),
        dtt_total_static_lim=(dop_t + stat_t).tolist(),
        dtt_total_growing=(dop[it] + net_growing[it]).tolist(),
        # code arrays for Legendre checks
        legendre=dict(
            dHH=leg_dHH,
            doppler=leg_dop,
            static=leg_stat,
            isw_growing=leg_isw_g,
            total_static_lim=leg_tot_s,
            total_growing=leg_tot_g,
        ),
    )


def main():
    rows = []
    for D in DS:
        print("solving mcc at D=%.4f ..." % D, flush=True)
        mcc = mcc_for(D)
        print("  mcc=%.6f  integrating..." % mcc, flush=True)
        row = run_one(D, mcc)
        rows.append(row)
        print(
            "  D=%s  Meff=%+.3e  |ΔH/H|a2=%.4f (ptp=%.3f)  "
            "stat_a2=%+.3e  ISW_static=0  dop_a2=%+.4f  tot_a2=%+.4f  "
            "dop/CMB_pass=%.0fx  vrmin=%.3f"
            % (
                row["D_label"],
                row["Meff"],
                row["cluster_a2"],
                row["dHH_ptp"],
                row["static_a2"],
                row["doppler_a2"],
                row["total_static_lim_a2"],
                row["ratio_doppler_a2_over_CMB_pass"],
                row["vrmin"],
            ),
            flush=True,
        )

    out = dict(
        note=(
            "CMB–Hubble decoupling from one Φ. Compact core (point) + Ω_s=0.049 shell. "
            "Matter: ΔH/H = af/<af>−1 at Rs=0.31 Gpc (v_grav·n̂/r from −∫∇Φ dt). "
            "CMB: endpoint static (Φ_e−Φ_o)/c² + ISW (0 if Φ̇=0; growing-mode estimate otherwise) "
            "+ Doppler −n·(V_e−V_o)/c. Mc re-solved so cluster a2=0.084. ⟨H⟩=71.6. R=14.3 Gpc. "
            "θ=0 toward CoE."
        ),
        R_Gpc=float(RV),
        Rs_cluster_Gpc=RS_CL,
        H_mean=HM,
        HRC=float(HRC),
        target_cluster_a2=TARGET_A2,
        CMB_pass=PASS_CMB,
        D_static_null=D_NULL,
        multipliers=MULT,
        shell=dict(Omega_s=0.049, kin=float(T.KIN), kout=float(T.KOUT)),
        rows=rows,
    )
    path = os.path.join(HERE, "decouple.json")
    json.dump(out, open(path, "w"), default=float)
    print("wrote", path)


if __name__ == "__main__":
    main()
