#!/usr/bin/env python3
"""Time derivatives of the curves on notes/cross-a.png.

The generator cross-a.py is not on disk. Functions are the ones recorded as
matching that PNG in notes/universe-storyboard.json (and checked against the
PNG pixels): E=0 shell, R scaled so R(61.2 Gyr)=14.3 Gpc, not R=ct.
"""
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

OUT_PNG = Path("/workspace/cosmology-course/notes/cross-a-deriv.png")
OUT_JSON = Path("/workspace/cosmology-course/notes/cross-a-deriv.json")

# Settled geometry used by the cross-a chart.
T0 = 0.550  # Gyr from the singularity to the crossing
T_TODAY = 61.2  # Gyr after the crossing; a=1 here
T_TODAY_S = T_TODAY + T0  # 61.75 Gyr from the singularity
D0 = 2.5  # Gpc at the crossing
T_KILL = 1.10  # Gyr after the crossing; D_rec hits zero
R_TODAY = 14.3  # Gpc; R = R_TODAY * t / T_TODAY
C = 0.3066  # Gpc/Gyr
A_CROSS = 1.0 / 1000.0
P = math.log(1.0 / A_CROSS) / math.log(T_TODAY_S / T0)


def shell(t):
    """Analytic curves and d/dt. t is time after the crossing (Gyr)."""
    t = np.asarray(t, dtype=float)
    ts = t + T0
    D = D0 * (ts / T0) ** (2.0 / 3.0)
    v = (2.0 / 3.0) * D / ts  # Gpc/Gyr; equals sqrt(2GM/D) for this normalisation
    a = (ts / T_TODAY_S) ** P
    da = P * a / ts
    R = R_TODAY * t / T_TODAY
    dR = R_TODAY / T_TODAY
    dDa = v * a + D * da
    dRa = dR * a + R * da
    d_minus = (v - dR) * a + (D - R) * da
    d_plus = (v + dR) * a + (D + R) * da
    # Purple curve on cross-a.png is v/c (right axis), not v in Gpc/Gyr.
    d_vc = (-v / (3.0 * ts)) / C
    out = {
        "t": t,
        "D": D,
        "v": v,
        "v_over_c": v / C,
        "a": a,
        "R": R,
        "Da": D * a,
        "(D-R)a": (D - R) * a,
        "(D+R)a": (D + R) * a,
        "Ra": R * a,
        "d(Da)/dt": dDa,
        "d((D-R)a)/dt": d_minus,
        "d((D+R)a)/dt": d_plus,
        "d(Ra)/dt": dRa,
        "d(v/c)/dt": d_vc,
        "da/dt": da,
    }
    u = 1.0 - t / T_KILL
    alive = t < T_KILL
    Drec = np.full_like(t, np.nan)
    dDrec = np.full_like(t, np.nan)
    dDrec_a = np.full_like(t, np.nan)
    Drec_a = np.full_like(t, np.nan)
    uu = np.clip(u, 1e-30, None)
    Drec[alive] = D0 * uu[alive] ** (1.0 / 3.0)
    dDrec[alive] = D0 * (1.0 / 3.0) * uu[alive] ** (-2.0 / 3.0) * (-1.0 / T_KILL)
    Drec_a[alive] = Drec[alive] * a[alive]
    dDrec_a[alive] = dDrec[alive] * a[alive] + Drec[alive] * da[alive]
    out["D_rec"] = Drec
    out["D_rec a"] = Drec_a
    out["d(D_rec)/dt"] = dDrec
    out["d(D_rec a)/dt"] = dDrec_a
    return out


def main():
    # Log grid over the chart window, plus a fine linear grid up to the kill
    # so the cube-root derivative is analytic and not a finite-difference spike.
    t_log = np.logspace(-3, math.log10(T_TODAY * 1.02), 1200)
    t_rec = np.linspace(1e-3, T_KILL * (1.0 - 1e-9), 8000)
    s_log = shell(t_log)
    s_rec = shell(t_rec)

    keys_finite = [
        "d(Da)/dt",
        "d((D+R)a)/dt",
        "d((D-R)a)/dt",
        "d(Ra)/dt",
        "d(v/c)/dt",
        "da/dt",
    ]
    finite = np.concatenate([s_log[k] for k in keys_finite])
    ymin = float(np.min(finite))
    ymax = float(np.max(finite))
    pad = 0.12 * (ymax - ymin)
    ymin -= pad
    ymax += pad

    # The cube-root derivative is finite and small until the last ~0.001 Gyr,
    # then runs to -infinity. On a log time axis that runaway is one vertical
    # pixel column. Stop at the last point still above a finite floor so the
    # trace does not become that spike, and do not draw anything at t>=1.10.
    Y_FLOOR = -0.25  # Gpc/Gyr; still finite, before the divergence
    yrec = s_rec["d(D_rec a)/dt"]
    crossed = np.where(yrec <= Y_FLOOR)[0]
    if len(crossed):
        last = int(crossed[0]) - 1
    else:
        last = len(t_rec) - 1
    last = max(last, 1)
    t_rec_plot = t_rec[: last + 1]
    y_rec_plot = yrec[: last + 1]

    plt.rcParams.update({
        "font.size": 14,
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
    })
    fig, ax = plt.subplots(figsize=(14.5, 7.2), layout="constrained")

    # Same colours and linestyles as cross-a.png.
    ax.plot(t_log, s_log["d(Da)/dt"], color="tab:red", lw=2.2, label=r"$d(Da)/dt$ (Gpc/Gyr)")
    ax.plot(t_rec_plot, y_rec_plot, color="tab:blue", lw=2.2, label=r"$d(D_{\mathrm{rec}}a)/dt$ (Gpc/Gyr)")
    ax.plot(t_log, s_log["d((D+R)a)/dt"], color="tab:green", lw=2.0, ls="--", label=r"$d((D+R)a)/dt$ (Gpc/Gyr)")
    ax.plot(t_log, s_log["d((D-R)a)/dt"], color="tab:green", lw=2.0, ls=":", label=r"$d((D-R)a)/dt$ (Gpc/Gyr)")
    ax.plot(t_log, s_log["d(Ra)/dt"], color="tab:green", lw=2.2, ls="-", label=r"$d(Ra)/dt$ (Gpc/Gyr)")
    ax.plot(t_log, s_log["d(v/c)/dt"], color="tab:purple", lw=2.2, label=r"$d(v/c)/dt$ (1/Gyr)")
    ax.plot(t_log, s_log["da/dt"], color="black", lw=2.2, label=r"$da/dt$ (1/Gyr)")

    ax.axhline(0.0, color="0.5", lw=0.8, zorder=0)
    ax.axvline(T_KILL, color="0.45", lw=1.0, zorder=0)
    ax.axvline(T_TODAY, color="0.2", lw=1.0, ls="--", zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 74)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("time after the crossing (Gyr)")
    ax.set_ylabel("time derivative")
    ax.set_title("Time derivative of each curve")
    ax.legend(loc="upper left", framealpha=0.95)
    ax.annotate(
        r"$d(D_{\mathrm{rec}}a)/dt$ cut at $-0.25$ Gpc/Gyr; diverges as $t\to 1.10$ Gyr",
        xy=(0.99, 0.02),
        xycoords="axes fraction",
        ha="right",
        va="bottom",
        fontsize=14,
        color="tab:blue",
    )
    fig.savefig(OUT_PNG, dpi=140)
    plt.close(fig)

    # Numerical gradient as a check, not as the plotted curves.
    t_check = np.logspace(-3, math.log10(T_TODAY), 5000)
    sc = shell(t_check)
    grad_err = {}
    for name, qty in (
        ("Da", "d(Da)/dt"),
        ("(D+R)a", "d((D+R)a)/dt"),
        ("(D-R)a", "d((D-R)a)/dt"),
        ("Ra", "d(Ra)/dt"),
        ("v/c", "d(v/c)/dt"),
        ("a", "da/dt"),
    ):
        g = np.gradient(sc[name if name != "v/c" else "v_over_c"], t_check)
        analytic = sc[qty]
        # ignore the first/last few samples of the one-sided gradient
        sl = slice(20, -20)
        rel = np.max(np.abs(g[sl] - analytic[sl]) / np.maximum(np.abs(analytic[sl]), 1e-8))
        grad_err[qty] = float(rel)

    t_rec_check = np.linspace(0.05, 1.05, 4000)
    sr = shell(t_rec_check)
    g = np.gradient(sr["D_rec a"], t_rec_check)
    sl = slice(20, -20)
    grad_err["d(D_rec a)/dt"] = float(
        np.max(np.abs(g[sl] - sr["d(D_rec a)/dt"][sl]) / np.maximum(np.abs(sr["d(D_rec a)/dt"][sl]), 1e-8))
    )

    samples = {}
    for t in (0.1, 1.0, 5.0, 13.8, 61.2):
        row = shell(np.array([t]))
        entry = {}
        for k in (
            "Da",
            "D_rec a",
            "(D-R)a",
            "(D+R)a",
            "Ra",
            "v_over_c",
            "a",
            "d(Da)/dt",
            "d(D_rec a)/dt",
            "d((D-R)a)/dt",
            "d((D+R)a)/dt",
            "d(Ra)/dt",
            "d(v/c)/dt",
            "da/dt",
        ):
            val = float(row[k][0])
            entry[k] = None if math.isnan(val) else val
        samples[str(t)] = entry

    # Where the plotted recombination derivative leaves the frame.
    t_cut = float(t_rec_plot[-1])
    y_cut = float(y_rec_plot[-1])

    sidecar = {
        "png": str(OUT_PNG),
        "source_chart": "/workspace/cosmology-course/notes/cross-a.png",
        "script_found": False,
        "note": (
            "cross-a.py is not on disk. Curves are the analytic functions "
            "recorded as matching cross-a.png in universe-storyboard.json, "
            "and checked against the PNG: Ra at 61.2 Gyr sits on 14.3 Gpc, not on c*61.2."
        ),
        "R_law": "R = 14.3 * t / 61.2 Gpc, with t the time after the crossing in Gyr. Not R = c t.",
        "R_at_61.2_Gpc": R_TODAY,
        "c_t_at_61.2_Gpc": C * T_TODAY,
        "a_of_t": (
            "a = (t_sing / 61.75)^p, t_sing = t + 0.550, "
            "p = ln(1000)/ln(61.75/0.550). a(0)=1/1000, a(61.2 Gyr)=1."
        ),
        "p": P,
        "quantities_differentiated": {
            "Da": "Gpc; derivative Gpc/Gyr",
            "D_rec a": "Gpc; derivative Gpc/Gyr; only for t < 1.10",
            "(D-R)a": "Gpc; derivative Gpc/Gyr",
            "(D+R)a": "Gpc; derivative Gpc/Gyr",
            "Ra": "Gpc; derivative Gpc/Gyr",
            "v/c": "dimensionless purple curve on the chart; derivative 1/Gyr",
            "a": "dimensionless; derivative 1/Gyr",
        },
        "formulas": {
            "t_sing": "t + 0.550",
            "D": "2.5 * (t_sing/0.550)^(2/3)",
            "v": "(2/3) * D / t_sing",
            "v_over_c": "v / 0.3066",
            "a": "(t_sing/61.75)^p",
            "R": "14.3 * t / 61.2",
            "D_rec": "2.5 * (1 - t/1.10)^(1/3) for t < 1.10; absent afterwards. No floor segment.",
            "d(Da)/dt": "v*a + D*da/dt",
            "da/dt": "p * a / t_sing",
            "d(Ra)/dt": "(14.3/61.2)*a + R*da/dt",
            "d((D-R)a)/dt": "(v - dR/dt)*a + (D-R)*da/dt",
            "d((D+R)a)/dt": "(v + dR/dt)*a + (D+R)*da/dt",
            "d(v/c)/dt": "-(v/(3*t_sing))/c",
            "d(D_rec)/dt": "2.5 * (1/3) * (1-t/1.10)^(-2/3) * (-1/1.10)",
            "d(D_rec a)/dt": "d(D_rec)/dt * a + D_rec * da/dt",
        },
        "y_window": {"ymin": ymin, "ymax": ymax, "reason": "span of the non-divergent derivatives, with 12% padding"},
        "D_rec_cutoff": {
            "t_last_plotted_Gyr": t_cut,
            "value_last_plotted": y_cut,
            "how": (
                "Analytic derivative drawn only while it is greater than -0.25 Gpc/Gyr. "
                "Past that it runs to -infinity as t approaches 1.10 Gyr; that spike is not drawn. "
                "Nothing is drawn at or after 1.10 Gyr, and there is no horizontal floor."
            ),
            "floor_Gpc_per_Gyr": -0.25,
        },
        "gradient_check_max_rel_err": grad_err,
        "samples_Gyr": samples,
    }
    OUT_JSON.write_text(json.dumps(sidecar, indent=2) + "\n")
    print("wrote", OUT_PNG, OUT_PNG.stat().st_size)
    print("p", P)
    print("ylim", ymin, ymax)
    print("D_rec cut", t_cut, y_cut)
    print("grad", json.dumps(grad_err, indent=2))
    print("samples", json.dumps(samples, indent=2))


if __name__ == "__main__":
    main()
