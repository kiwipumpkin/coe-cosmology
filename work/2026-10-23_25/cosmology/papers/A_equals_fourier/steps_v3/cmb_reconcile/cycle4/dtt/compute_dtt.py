"""ΔT/T vs angle on the R=14.3 Gpc CMB emission sphere for Chris's compact-core + shell model.
Reuses cycle4/lib4 (Newtonian integrator, same linear Doppler convention). Headline compact core,
Mc re-solved per D so the all-cluster quadrupole stays at 0.084. D = 14.7, 22.28 (static null), 49."""
import json, os, sys
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import lib4 as L

T = L.T
RV = L.RV
HM = L.HM
CKM = L.CKM
TH1 = L.TH1  # deg; code θ=0 is AWAY from CoE
HRC = HM * RV * 1e3 / CKM  # HR/c ≈ 3.415

# Prefer these three D (user). Static null taken from cycle4/partB.json when available.
PARTB = os.path.join(HERE, "..", "partB.json")
if os.path.isfile(PARTB):
    _pb = json.load(open(PARTB))
    D_NULL = float(_pb["static_null"]["D"])  # ≈ 22.276
else:
    D_NULL = 22.28

DS = [14.7, D_NULL, 49.0]
D_LABELS = {14.7: "14.7", D_NULL: "22.28", 49.0: "49"}


def mcc_for(D, target=0.084, Rs=0.31):
    f = lambda lm: np.log(T.dq(L.run(D, 10**lm, "point", {}, [Rs])[0][0]["af"])["a2"] / target)
    return 10 ** brentq(f, -3.5, -2.3, xtol=1e-6)


def curves(D, mcc):
    """Full θ arrays on the R sphere. Code θ=0 = away from CoE; we also return toward-CoE indexing."""
    res, info = L.run(D, mcc, "point", {}, [0.31, RV])
    cl, far = res
    Hconv = HM / cl["af"].mean()
    vunit = Hconv * D * 1e3
    X, V = far["X"], far["V"]
    d = X[1:] - X[0]
    n = d / np.hypot(d[:, 0], d[:, 1])[:, None]
    kin = -np.sum(n * (V[1:] - V[0]), axis=1) * vunit / CKM  # linear Doppler, cycle4 convention
    Me = L.meff(D, mcc, "point")
    stat = L.static_pattern(D, Me)
    tot = kin + stat
    # Fractional anisotropic-expansion estimate: divide by HR/c (linear v/c outside range).
    kin_frac = kin / HRC
    # Code θ=0 = away from CoE; plot θ=0 = toward CoE → reverse indices 180..0
    i_toward = np.arange(180, -1, -1)
    th = np.arange(0.0, 181.0, 1.0)
    kin_t = kin[i_toward]
    stat_t = stat[i_toward]
    tot_t = tot[i_toward]
    kin_frac_t = kin_frac[i_toward]
    return dict(
        theta_toward_deg=th.tolist(),
        dtt_static=stat_t.tolist(),
        dtt_doppler=kin_t.tolist(),
        dtt_total=tot_t.tolist(),
        dtt_doppler_frac=kin_frac_t.tolist(),
        # also keep code-ordered full 360 for Legendre (same as lib4.cmb)
        kin_code=kin.tolist(),
        stat_code=stat.tolist(),
        Meff=float(Me),
        Hconv=float(Hconv),
        vrmin=float(min(cl["vrmin"], far["vrmin"])),
        Ht=float(info["Ht"]),
    )


def summarize(th, stat, dop, tot, dop_frac, kin_code, stat_code):
    kin_leg = L.legendre(np.array(kin_code))
    stat_leg = L.legendre(np.array(stat_code))
    tot_leg = L.legendre(np.array(kin_code) + np.array(stat_code))
    frac_leg = L.legendre(np.array(kin_code) / HRC)
    dop = np.array(dop)
    stat = np.array(stat)
    tot = np.array(tot)
    return dict(
        doppler_ptp=float(dop.max() - dop.min()),
        doppler_peak_abs=float(np.max(np.abs(dop - dop.mean()))),
        doppler_min=float(dop.min()),
        doppler_max=float(dop.max()),
        doppler_mean=float(dop.mean()),
        static_ptp=float(stat.max() - stat.min()),
        static_min=float(stat.min()),
        static_max=float(stat.max()),
        total_ptp=float(tot.max() - tot.min()),
        kin_legendre=kin_leg,
        stat_legendre=stat_leg,
        tot_legendre=tot_leg,
        frac_legendre=frac_leg,
        doppler_a2=float(kin_leg[2]),
        static_a2=float(stat_leg[2]),
        total_a2=float(tot_leg[2]),
        frac_a2=float(frac_leg[2]),
        frac_ptp=float(np.max(dop_frac) - np.min(dop_frac)),
    )


def main():
    rows = []
    for D in DS:
        print("solving mcc at D=%.4f ..." % D, flush=True)
        mcc = mcc_for(D)
        print("  mcc=%.6f  running integrator..." % mcc, flush=True)
        c = curves(D, mcc)
        s = summarize(
            c["theta_toward_deg"],
            c["dtt_static"],
            c["dtt_doppler"],
            c["dtt_total"],
            c["dtt_doppler_frac"],
            c["kin_code"],
            c["stat_code"],
        )
        rin = 1 - T.KIN * RV / D
        Mh = 0.049 * rin**3
        row = dict(
            D=D,
            D_label=D_LABELS[D],
            mcc=mcc,
            Mhole=float(Mh),
            Meff=c["Meff"],
            HRC=float(HRC),
            Hconv=c["Hconv"],
            vrmin=c["vrmin"],
            Ht=c["Ht"],
            theta_toward_deg=c["theta_toward_deg"],
            dtt_static=c["dtt_static"],
            dtt_doppler=c["dtt_doppler"],
            dtt_total=c["dtt_total"],
            dtt_doppler_frac=c["dtt_doppler_frac"],
            **s,
        )
        rows.append(row)
        print(
            "  D=%s  Meff=%+.4e  dop_ptp=%.4f  dop_a2=%+.4e  stat_a2=%+.4e  frac_a2=%+.4e  vrmin=%.3f"
            % (D_LABELS[D], c["Meff"], s["doppler_ptp"], s["doppler_a2"], s["static_a2"], s["frac_a2"], c["vrmin"]),
            flush=True,
        )

    out = dict(
        note="Compact core (point) + Omega_s=0.049 shell. Angle θ=0 toward CoE, θ=180 away. "
        "ΔT/T_Doppler = -n·(V_emit-V_obs)/c linear (cycle4). Static = (Φ_emit-Φ_obs)/c². "
        "Fractional = linear / (HR/c). Mc re-solved so cluster a2=0.084 at Rs=0.31 Gpc. <H>=71.6.",
        R_Gpc=float(RV),
        H_mean=HM,
        HRC=float(HRC),
        shell=dict(Omega_s=0.049, kin=float(T.KIN), kout=float(T.KOUT)),
        D_static_null=D_NULL,
        rows=rows,
    )
    path = os.path.join(HERE, "dtt.json")
    # strip huge code arrays already represented; keep plot arrays
    json.dump(out, open(path, "w"), default=float)
    print("wrote", path)


if __name__ == "__main__":
    main()
