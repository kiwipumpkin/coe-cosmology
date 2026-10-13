"""Simple charts + summary for sequential scenarios A and B."""
import os, json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BASE0, BASE4 = 237.26, 88.81

def load(name):
    with open(os.path.join(HERE, name)) as f:
        return json.load(f)

def chart_a(a):
    # 1) local greedy trajectory (best seed)
    traj = a["local"]["traj"]
    ns = [t["n_units"] for t in traj]
    left = [t["leftover"] for t in traj]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(ns, left, "o-", color="#1f77b4", lw=2, ms=6, label="local leftover")
    ax.axhline(BASE0, color="0.4", ls="--", lw=1, label=f"baseline {BASE0:.0f}×")
    ax.axhline(BASE4, color="#d62728", ls="--", lw=1, label=f"four clumps {BASE4:.0f}×")
    ax.set_xlabel("local units placed (one at a time)")
    ax.set_ylabel("leftover / floor")
    ax.set_title("Scenario A — fixing the near residual")
    ax.legend(loc="best", frameon=False)
    ax.set_ylim(0, max(BASE0 * 1.05, max(left) * 1.05))
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "a_local_trajectory.png"), dpi=140)
    plt.close()

    # 2) outer configs: median leftover bars + helped/hurt
    keys = ["empty", "random", "foam", "lattice"]
    labels = ["empty\n(control)", "random\npack", "foam-like", "lattice"]
    frozen = a["local"]["leftover"]
    meds, helped, hurt = [], [], []
    for k in keys:
        o = a["outer"][k]
        if k == "empty":
            meds.append(o["left"]); helped.append(0); hurt.append(0)
        else:
            meds.append(o["median"]); helped.append(o["helped"]); hurt.append(o["hurt"])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axes[0]
    colors = ["#8c8c8c", "#1f77b4", "#2ca02c", "#ff7f0e"]
    ax.bar(labels, meds, color=colors, width=0.65)
    ax.axhline(frozen, color="black", ls=":", lw=1.5, label=f"frozen local {frozen:.0f}×")
    ax.axhline(BASE0, color="0.5", ls="--", lw=1)
    ax.axhline(BASE4, color="#d62728", ls="--", lw=1)
    ax.set_ylabel("median leftover / floor")
    ax.set_title("Outer shell configs (local frozen)")
    ax.legend(loc="best", frameon=False, fontsize=8)

    ax = axes[1]
    x = np.arange(len(keys) - 1)
    h = helped[1:]; u = hurt[1:]
    ax.bar(x - 0.18, h, 0.35, color="#2ca02c", label="helped (vs frozen)")
    ax.bar(x + 0.18, u, 0.35, color="#d62728", label="hurt (vs frozen)")
    ax.set_xticks(x); ax.set_xticklabels(labels[1:])
    ax.set_ylabel("layouts")
    ax.set_title("Helped vs hurt vs frozen local")
    ax.legend(loc="best", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "a_outer_effect.png"), dpi=140)
    plt.close()

def chart_b(b):
    # trajectory lines — one panel per rule, mean of seeds, simple
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = {"void_then_rim": "#1f77b4", "alternate": "#2ca02c",
              "anti_coe_bias": "#ff7f0e", "random_scatter": "#7f7f7f"}
    labels = {"void_then_rim": "void then rim", "alternate": "alternate V/C",
              "anti_coe_bias": "anti-CoE bias", "random_scatter": "random scatter"}
    for name, rule in b["rules"].items():
        # plot each seed lightly, then median envelope via longest traj resampled by n_obj
        for run in rule["runs"]:
            t = run["traj"]
            xs = [p.get("n_units", p["n_obj"]) for p in t]
            ys = [p["leftover"] for p in t]
            ax.plot(xs, ys, color=colors[name], alpha=0.25, lw=1)
        # pick best seed as bold
        finals = [r["leftover"] for r in rule["runs"]]
        ib = int(np.argmin(finals))
        t = rule["runs"][ib]["traj"]
        xs = [p.get("n_units", p["n_obj"]) for p in t]
        ys = [p["leftover"] for p in t]
        ax.plot(xs, ys, color=colors[name], lw=2.5, label=labels[name])
    ax.axhline(BASE0, color="0.4", ls="--", lw=1, label=f"baseline {BASE0:.0f}×")
    ax.axhline(BASE4, color="#d62728", ls="--", lw=1, label=f"four clumps {BASE4:.0f}×")
    ax.set_xlabel("units placed (or objects for alternate)")
    ax.set_ylabel("leftover / floor")
    ax.set_title("Scenario B — residual as objects are added")
    ax.legend(loc="best", frameon=False, fontsize=8)
    # keep linear, avoid huge scale if random blows up
    ymax = max(BASE0 * 1.2, BASE4 * 2)
    for rule in b["rules"].values():
        for run in rule["runs"]:
            ymax = max(ymax, min(max(p["leftover"] for p in run["traj"]), BASE0 * 3))
    ax.set_ylim(0, ymax)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "b_growth_trajectory.png"), dpi=140)
    plt.close()

    # final helped/hurt style bar
    fig, ax = plt.subplots(figsize=(7, 4))
    names = list(b["rules"].keys())
    labs = [labels[n] for n in names]
    meds = [b["rules"][n]["median_final"] for n in names]
    bests = [b["rules"][n]["best_final"] for n in names]
    x = np.arange(len(names))
    ax.bar(x - 0.18, meds, 0.35, color="#1f77b4", label="median final")
    ax.bar(x + 0.18, bests, 0.35, color="#2ca02c", label="best final")
    ax.axhline(BASE0, color="0.4", ls="--", lw=1)
    ax.axhline(BASE4, color="#d62728", ls="--", lw=1)
    ax.set_xticks(x); ax.set_xticklabels(labs, rotation=15, ha="right")
    ax.set_ylabel("leftover / floor")
    ax.set_title("Scenario B — final leftover by growth rule")
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "b_final_bars.png"), dpi=140)
    plt.close()

def write_md(a, b):
    frozen = a["local"]["leftover"]
    g = a["geometry"]
    lines = []
    lines.append("# Sequential non-overlapping placement — Scenarios A & B\n")
    lines.append("Objects placed **one at a time** with hard exclusion (centers ≥ R_i + R_j apart; ")
    lines.append("void–rim contact inside a unit allowed). Compensated units: void holds −0.546 M_cell, ")
    lines.append("rim clusters hold +0.546 M_cell (wall+filament folded into clusters so net mass = 0).\n")
    lines.append(f"\nBaseline leftover (SR + central mass): **{a['baseline']['BASE0']:.1f}×** floor. ")
    lines.append(f"Four fitted clumps: **{a['baseline']['BASE4']:.1f}×**.\n")
    lines.append(f"\n**Effective sizes:** local unit cell s = {g['S_local']} Gpc → void R ≈ {g['R_v_local']:.2f} Gpc, ")
    lines.append(f"cluster R ≈ {g['R_c_local']:.2f} Gpc. Outer pack s = {g['S_outer']} Gpc → void R ≈ {g['R_v_outer']:.2f} Gpc. ")
    lines.append("(Real-web voids are ~0.03 Gpc; coarse-graining keeps packing fractions and mass balance runnable.)\n")

    lines.append("\n## Scenario A — fix near dipole/residual, then fill the rest\n")
    lines.append(f"\nLocal region 14.3–30 Gpc, greedy non-overlapping units (seed {a['local']['seed']}): ")
    lines.append(f"**{a['local']['info']['n_units']} units**, leftover **{frozen:.1f}×** ")
    lines.append(f"(best-scaled {a['local']['best_scaled']:.1f}×). ")
    lines.append(f"Overlap check: min margin {a['local']['info']['min_margin_gpc']:.3f} Gpc, ok={a['local']['info']['no_overlap']}.\n")
    lines.append("\nThen the outer shell (30–42.9 Gpc from us, CoE interior still empty) was filled:\n\n")
    lines.append("|outer config|median leftover|Δ vs frozen|helped / hurt|median outer units|\n|-|-|-|-|-|\n")
    lines.append(f"|empty (control)|{frozen:.1f}×|0|—|0|\n")
    for k in ("random", "foam", "lattice"):
        o = a["outer"][k]; d = a["outer_vs_frozen"][k]["median_delta"]
        lines.append(f"|{k}|{o['median']:.1f}×|{d:+.1f}×|{o['helped']}/{o['hurt']}|{o['median_outer_units']:.0f}|\n")
    lines.append("\n**Does the rest of the shell matter once the near residual is set?** ")
    deltas = [abs(a["outer_vs_frozen"][k]["median_delta"]) for k in ("random", "foam", "lattice")]
    lines.append(f"Median move from outer fills: {np.mean(deltas):.1f}× (range {min(deltas):.1f}–{max(deltas):.1f}×). ")
    if np.mean(deltas) < 15:
        lines.append("Once the local arrangement is frozen, outer packing changes our leftover only a little — ")
        lines.append("the near region dominates our view of the CMB residual.\n")
    else:
        lines.append("Outer packing still moves the needle noticeably; the far shell is not negligible.\n")

    lines.append("\n## Scenario B — sequential growth from empty\n")
    lines.append(f"\nCell s = {b['s']} Gpc (void R ≈ {b['R_v']:.2f} Gpc). Five seeds per rule.\n\n")
    lines.append("|growth rule|median final|best final|helped / hurt vs 237×|\n|-|-|-|-|\n")
    labs = {"void_then_rim": "void then rim", "alternate": "alternate void/cluster",
            "anti_coe_bias": "anti-CoE bias", "random_scatter": "random scatter"}
    for k, lab in labs.items():
        r = b["rules"][k]
        lines.append(f"|{lab}|{r['median_final']:.1f}×|{r['best_final']:.1f}×|{r['helped_vs_base0']}/{r['hurt_vs_base0']}|\n")
    lines.append("\nGreedy unit growth (void-then-rim, anti-CoE) can beat random scatter. ")
    lines.append("None of the automatic growth rules reliably reach the four-clump **89×** class — ")
    lines.append("that still needs an arranged fit, not just sequential packing with the right mix.\n")
    lines.append("\n## Charts\n")
    lines.append("- `a_local_trajectory.png` — leftover as local units are placed\n")
    lines.append("- `a_outer_effect.png` — outer configs vs frozen local; helped/hurt bars\n")
    lines.append("- `b_growth_trajectory.png` — leftover vs objects added\n")
    lines.append("- `b_final_bars.png` — final leftover by rule\n")
    with open(os.path.join(HERE, "summary.md"), "w") as f:
        f.write("".join(lines))
    print("wrote summary.md and PNGs")

if __name__ == "__main__":
    a = load("scenario_a.json"); b = load("scenario_b.json")
    chart_a(a); chart_b(b); write_md(a, b)
