#!/usr/bin/env python3
"""
Extrema of the model A/A_all(l) curve vs each input variable (one at a time).

Model: imported unchanged from
  /workspace/figs_letter_exact_python/makers/make_one_at_a_time_levelsets.py
  (model_A_curve -> af_an_alpha + model_AB_fixed, Amp_A=-0.06, phi_A=0.35 rad,
   f_u=0.94, af_mode="hratio", 90-deg wrap boxcar), i.e. exactly the curves drawn
   in preprint Figs A10-A13.  No fitting of any kind is done here.

Hold values match the level-set script, including its d/D override:
  D/L, alpha, kappa sweeps : D/L=5.9, d/D=0.1, alpha=1.0, kappa=0.230588
  d/D sweep (Fig A13)      : same but alpha=1.3 (code: hold["alpha"]=1.3)

Method: the model is evaluated on its native 361-point grid (as in the figures),
interpolated with a periodic cubic spline onto a 0.1-deg l grid, local maxima
and minima found on that grid (circular neighbours, so 0/360 wraps), then each
extremum refined to the spline's dA/dl=0 root inside its 0.1-deg bracket.
Extrema are tracked as continuous branches outward from the baseline value;
a branch ends where a max and min meet and annihilate (the turning point is
refined by bisection on the swept variable).
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

NOTES = Path("/workspace/cosmology-course/notes")
SRC = Path("/workspace/figs_letter_exact_python/makers/make_one_at_a_time_levelsets.py")
_spec = importlib.util.spec_from_file_location("levelsets", SRC)
lv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lv)

N_SWEEP = 200
DL_GRID = 0.1  # deg

HOLD = lv.baseline_hold()
SWEEPS = [  # key, label, lo, hi, baseline, file tag, held overrides
    ("x", r"$D/L$", 3.0, 10.0, 5.9, "DL", {}),
    ("alpha", r"$\alpha$", 0.2, 2.0, 1.0, "alpha", {}),
    ("kappa", r"$\kappa$", 0.15, 0.70, lv.KAPPA0, "kappa", {}),
    ("eps", r"$d/D$", 0.06, 0.74, 0.10, "dD", {"alpha": 1.3}),  # as Fig A13 code
]
PLAIN = {"x": "D/L", "alpha": "alpha", "kappa": "kappa", "eps": "d/D"}
# colours by family (order of baseline angle): max-low, max-high, min-low, min-high
COLS = {"max": ["#d62728", "#ff7f0e", "#8c564b", "#e377c2"],
        "min": ["#1f77b4", "#2ca02c", "#17becf", "#9467bd"]}


def curve(kw):
    return lv.model_A_curve(x=kw["x"], f_u=kw["f_u"], alpha=kw["alpha"],
                            eps=kw["eps"], kappa=kw["kappa"])


def extrema(kw):
    """Return list of (kind, l_deg, value) for every local max/min, 0<=l<360."""
    lf, ma = curve(kw)
    cs = CubicSpline(np.r_[lf, 360.0], np.r_[ma, ma[0]], bc_type="periodic")
    d1 = cs.derivative()
    g = np.arange(0.0, 360.0, DL_GRID)
    y = cs(g)
    yp, yn = np.roll(y, 1), np.roll(y, -1)
    out = []
    for kind, mask in (("max", (y > yp) & (y >= yn)), ("min", (y < yp) & (y <= yn))):
        for i in np.where(mask)[0]:
            a, b = g[i] - DL_GRID, g[i] + DL_GRID
            try:
                l0 = brentq(lambda t: d1(t % 360.0), a, b, xtol=1e-6)
            except ValueError:
                l0 = g[i]
            l0 %= 360.0
            out.append((kind, float(l0), float(cs(l0))))
    return out


def counts(kw):
    e = extrema(kw)
    return sum(k == "max" for k, *_ in e), sum(k == "min" for k, *_ in e)


def cdist(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def track(key, lo, hi, base, over):
    hold = dict(HOLD); hold.update(over)
    vals = np.unique(np.r_[np.linspace(lo, hi, N_SWEEP), base])
    ib = int(np.argmin(abs(vals - base)))
    ex = []
    for v in vals:
        kw = dict(hold); kw[key] = float(v)
        ex.append(extrema(kw))
    # branches seeded at baseline
    branches = []
    for kind in ("max", "min"):
        seeds = sorted([e for e in ex[ib] if e[0] == kind], key=lambda e: e[1])
        for j, e in enumerate(seeds):
            branches.append(dict(kind=kind, fam=j, l0=e[1],
                                 l=np.full(len(vals), np.nan),
                                 A=np.full(len(vals), np.nan)))
            branches[-1]["l"][ib] = e[1]; branches[-1]["A"][ib] = e[2]
    ends = []
    for step in (+1, -1):
        alive = {id(b): True for b in branches}
        i = ib
        while 0 <= i + step < len(vals):
            i2 = i + step
            pool = list(ex[i2])
            for b in branches:
                if not alive[id(b)] or np.isnan(b["l"][i]):
                    continue
                prev = b["l"][i]
                cand = [(cdist(prev, e[1]), n) for n, e in enumerate(pool)
                        if e is not None and e[0] == b["kind"]]
                cand = [c for c in cand if c[0] < 30.0]
                if not cand:
                    alive[id(b)] = False
                    ends.append(dict(branch=b, i_last=i, i_gone=i2))
                    continue
                _, n = min(cand)
                e = pool[n]; pool[n] = None
                # unwrap so the branch is continuous through 0/360
                l_new = prev + ((e[1] - prev + 180.0) % 360.0 - 180.0)
                b["l"][i2] = l_new; b["A"][i2] = e[2]
            for e in pool:
                if e is not None:
                    raise RuntimeError(f"new extremum appeared at {key}={vals[i2]}: {e}")
            i = i2
    # refine annihilation points (max+min ending together)
    merges = []
    used = set()
    for e1 in ends:
        for e2 in ends:
            if e1 is e2 or id(e1) in used or id(e2) in used:
                continue
            if e1["branch"]["kind"] == "max" and e2["branch"]["kind"] == "min" \
                    and e1["i_gone"] == e2["i_gone"]:
                used |= {id(e1), id(e2)}
                va, vb = vals[e1["i_last"]], vals[e1["i_gone"]]
                n_ok = counts({**hold, key: float(va)})
                for _ in range(40):
                    vm = 0.5 * (va + vb)
                    if counts({**hold, key: float(vm)}) == n_ok:
                        va = vm
                    else:
                        vb = vm
                exs = extrema({**hold, key: float(va)})
                lm = [e for e in exs if e[0] == "max"
                      and cdist(e[1], e1["branch"]["l"][e1["i_last"]] % 360) < 30]
                ln = [e for e in exs if e[0] == "min"
                      and cdist(e[1], e2["branch"]["l"][e2["i_last"]] % 360) < 30]
                lmax = min(lm, key=lambda e: cdist(e[1], ln[0][1]))
                lmin = min(ln, key=lambda e: cdist(e[1], lmax[1]))
                lmerge = (lmax[1] + 0.5 * ((lmin[1] - lmax[1] + 180) % 360 - 180)) % 360
                merges.append(dict(var=key, value=float(va), l_deg=float(lmerge),
                                   A=float(0.5 * (lmax[2] + lmin[2])),
                                   max_branch=f"max near {e1['branch']['l0']:.0f} deg",
                                   min_branch=f"min near {e2['branch']['l0']:.0f} deg",
                                   max_last=(float(lmax[1]), float(lmax[2])),
                                   min_last=(float(lmin[1]), float(lmin[2]))))
    for b in branches:
        b["name"] = f"{b['kind']} near {b['l0']:.0f}°"
        b["col"] = COLS[b["kind"]][b["fam"]]
    return dict(key=key, vals=vals, ib=ib, hold=hold, branches=branches,
                merges=merges, ex=ex)


def plot_panels(ax_l, ax_a, T, label, base, big=True, legend_handles=None):
    fs = 15 if big else 14
    for b in T["branches"]:
        ls = "-" if b["kind"] == "max" else "--"
        ax_l.plot(T["vals"], b["l"], ls=ls, color=b["col"], lw=2.6, label=b["name"])
        ax_a.plot(T["vals"], b["A"], ls=ls, color=b["col"], lw=2.6)
    for m in T["merges"]:
        lab = "max–min merge (" + PLAIN[T["key"]].replace("alpha", "α").replace("kappa", "κ") + f" = {m['value']:.3g})"
        ax_l.plot(m["value"], m["l_deg"], "o", mfc="none", mec="k", mew=2, ms=11,
                  label=lab, zorder=6)
        ax_a.plot(m["value"], m["A"], "o", mfc="none", mec="k", mew=2, ms=11, zorder=6)
    for ax in (ax_l, ax_a):
        ax.axvline(base, color="k", ls=":", lw=1.8, label="baseline" if ax is ax_l else None)
        ax.grid(alpha=0.35)
        ax.tick_params(labelsize=fs - 1)
        ax.set_xlim(T["vals"][0], T["vals"][-1])
    ax_a.axhline(1.0, color="#888888", ls="-", lw=0.9)
    allL = np.concatenate([b["l"] for b in T["branches"]])
    lo = 0.0 if np.nanmin(allL) >= 0 else 30 * np.floor((np.nanmin(allL) - 5) / 30); hi = 30 * np.ceil((np.nanmax(allL) + 5) / 30)
    ax_l.set_ylim(lo, hi)
    ax_l.set_yticks(np.arange(lo, hi + 1, 30 if big else 60))
    ax_l.tick_params(labelbottom=False)
    ax_a.set_xlabel(label, fontsize=fs + 1)


def held_str(T):
    h = T["hold"]; parts = []
    for k, f in (("x", "D/L={:g}"), ("eps", "d/D={:g}"), ("alpha", "α={:g}"), ("kappa", "κ={:.6g}")):
        if k != T["key"]:
            parts.append(f.format(h[k]))
    return ", ".join(parts)


def fig_single(T, label, base, tag):
    fig, (al, aa) = plt.subplots(2, 1, figsize=(9.5, 9.6), sharex=True,
                                 gridspec_kw=dict(hspace=0.08))
    plot_panels(al, aa, T, label, base)
    al.set_ylabel(r"$\ell$ of extremum [deg]", fontsize=16)
    aa.set_ylabel(r"$A/A_{\mathrm{all}}$ at extremum", fontsize=16)
    fig.suptitle(f"Extrema vs {PLAIN[T['key']].replace('alpha','α').replace('kappa','κ')}"
                 f" (others at baseline)\nheld: {held_str(T)}; "
                 "Amp$_A$=−0.06, $\\phi_A$=0.35 rad", fontsize=16, y=0.995)
    h, l = al.get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.915), ncol=2,
               fontsize=14, frameon=False)
    fig.subplots_adjust(top=0.79, bottom=0.08, left=0.13, right=0.97)
    p = NOTES / f"extrema-vs-{tag}.png"
    fig.savefig(p, dpi=200); plt.close(fig)
    return p


def fig_overview(Ts):
    fig, axs = plt.subplots(2, 4, figsize=(16, 8.4), sharex="col",
                            gridspec_kw=dict(hspace=0.08, wspace=0.28))
    for c, (T, (key, label, lo, hi, base, tag, over)) in enumerate(zip(Ts, SWEEPS)):
        plot_panels(axs[0, c], axs[1, c], T, label, base, big=False)
        ttl = f"vary {label}" + (r" ($\alpha$=1.3 held)" if key == "eps" else "")
        axs[0, c].set_title(ttl, fontsize=15)
    axs[0, 0].set_ylabel(r"$\ell$ of extremum [deg]", fontsize=15)
    axs[1, 0].set_ylabel(r"$A/A_{\mathrm{all}}$ at extremum", fontsize=15)
    # shared y range for the value row so columns compare directly
    ymin = min(np.nanmin(b["A"]) for T in Ts for b in T["branches"])
    ymax = max(np.nanmax(b["A"]) for T in Ts for b in T["branches"])
    for a in axs[1]:
        a.set_ylim(ymin - 0.03, ymax + 0.03)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], color=COLS["max"][0], lw=2.6, label="max near 10°"),
          Line2D([], [], color=COLS["max"][1], lw=2.6, label="max near 170°"),
          Line2D([], [], color=COLS["min"][0], lw=2.6, ls="--", label="min near 95°"),
          Line2D([], [], color=COLS["min"][1], lw=2.6, ls="--", label="min near 266°"),
          Line2D([], [], color="k", lw=1.8, ls=":", label="baseline"),
          Line2D([], [], marker="o", mfc="none", mec="k", mew=2, ms=11, ls="none",
                 label="max–min merge")]
    fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=6,
               fontsize=14, frameon=False, columnspacing=1.2, handlelength=2.2)
    fig.subplots_adjust(top=0.86, bottom=0.10, left=0.07, right=0.99)
    p = NOTES / "extrema-vs-inputs-overview.png"
    fig.savefig(p, dpi=200); plt.close(fig)
    return p


def summary_rows(T):
    rows = []
    vals = T["vals"]
    for tagname, i in (("low end", 0), ("baseline", T["ib"]), ("high end", len(vals) - 1)):
        for b in T["branches"]:
            rows.append(dict(var=PLAIN[T["key"]], point=tagname, value=float(vals[i]),
                             branch=b["name"].replace("°", " deg"),
                             l_deg=None if np.isnan(b["l"][i]) else round(float(b["l"][i] % 360), 3),
                             A_over_Aall=None if np.isnan(b["A"][i]) else round(float(b["A"][i]), 5)))
    return rows


def main():
    plt.rcParams.update({"mathtext.fontset": "dejavusans", "font.size": 14})
    Ts, out, allrows = [], {}, []
    for key, label, lo, hi, base, tag, over in SWEEPS:
        T = track(key, lo, hi, base, over)
        Ts.append(T)
        fig_single(T, label, base, tag)
        rows = summary_rows(T); allrows += rows
        # sanity: raw 361-pt plotted curve at baseline
        lf, ma = curve({**T["hold"], key: base})
        raw = {"argmax_raw_deg": [], "argmin_raw_deg": []}
        n = len(ma)
        for i in range(n):
            if ma[i] > ma[i - 1] and ma[i] >= ma[(i + 1) % n]:
                raw["argmax_raw_deg"].append(round(float(lf[i]), 2))
            if ma[i] < ma[i - 1] and ma[i] <= ma[(i + 1) % n]:
                raw["argmin_raw_deg"].append(round(float(lf[i]), 2))
        out[PLAIN[key]] = dict(
            range=[lo, hi], baseline=base, n_sweep=len(T["vals"]),
            hold={k: T["hold"][k] for k in ("x", "eps", "alpha", "kappa")},
            points=rows, merges=T["merges"],
            baseline_check_on_plotted_361pt_curve=raw)
        print(PLAIN[key], "merges:", T["merges"])
        print("  raw-curve baseline extrema:", raw)
    (NOTES / "extrema-vs-inputs.json").write_text(json.dumps(out, indent=2))
    with open(NOTES / "extrema-vs-inputs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(allrows[0].keys()))
        w.writeheader(); w.writerows(allrows)
    fig_overview(Ts)
    for r in allrows:
        print(r)


if __name__ == "__main__":
    main()
