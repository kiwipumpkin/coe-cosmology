"""Scenario B: sequential growth from empty — add one void OR one cluster at a time,
each placement conditioned on previous ones (non-overlapping). Track leftover trajectory.
Growth rules compared:
  1. void_then_rim  — place a void, then its rim clusters, repeat (unit growth)
  2. alternate      — void, cluster, void, cluster, ... (clusters prefer junctions near voids)
  3. anti_coe_bias  — prefer placing voids toward anti-CoE / clusters toward CoE
  4. random_scatter — independent random voids/clusters (still non-overlap) as control
"""
import os, sys, time, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from seq_core import *

OUT = HERE
S = 3.0                 # intermediate cell size for growth
N_RIM = 6
MAX_UNITS = 40          # for unit-based rules
MAX_OBJ = 200           # for object-by-object rules
TRACK_EVERY = 1


def _mass_for_void(center, s):
    g = cell_geometry(s)
    return g["f_void"] * float(rho_gpc(center)) * g["V"] * GPC ** 3, g["R_v"]

def _mass_for_clu(center, s, n_share=6):
    g = cell_geometry(s)
    # one cluster's share of excess
    return g["f_excess"] * float(rho_gpc(center)) * g["V"] * GPC ** 3 / n_share, g["R_c"]


def place_if_ok(center, mass, rad, Xcat, radcat, region="full"):
    if center is None:
        return False
    if region == "full":
        if not in_any_shell(center[None])[0]:
            return False
        if mass < 0 and np.linalg.norm(center) < RMIN + rad:
            return False
    if overlaps(center[None], np.array([rad]), Xcat, radcat, gap=0.0):
        return False
    return True


def growth_void_then_rim(rng, s=S, max_units=MAX_UNITS, region="full"):
    """Add one compensated unit at a time (void then rim clusters as a block, but
    recorded object-by-object on the trajectory after each full unit)."""
    g = cell_geometry(s)
    gap = GAP_FRAC * g["R_v"]
    Xcat = np.zeros((0, 3)); Mcat = np.zeros(0); radcat = np.zeros(0)
    w = np.zeros_like(R)
    traj = [dict(n_obj=0, n_units=0, leftover=float(BASE0))]
    for u in range(max_units):
        # try several candidates, pick the one that helps most (gentle greedy + random mix)
        best = None
        for _ in range(40):
            c = propose_center(rng, region if region != "full" else "full", r_void=g["R_v"])
            # propose_center regions: local|outer|full — "full" works
            if c is None:
                continue
            Xu, Mu, ru, ku, _ = make_unit(c, s, float(rho_gpc(c)), n_rim=N_RIM, gap=gap)
            if not in_any_shell(Xu).all():
                continue
            if overlaps(Xu, ru, Xcat, radcat, gap=0.0):
                continue
            dw = pattern(Xu, Mu)
            left = leftover(w + dw)
            if best is None or left < best[0]:
                best = (left, Xu, Mu, ru, dw)
        if best is None:
            break
        # accept even if mild hurt early on (exploration); stop if big hurt after a few
        left, Xu, Mu, ru, dw = best
        if u >= 3 and left > traj[-1]["leftover"] + 5:
            # reject and try pure random accept once
            accepted = False
            for _ in range(60):
                c = propose_center(rng, "full", r_void=g["R_v"])
                if c is None:
                    continue
                Xu, Mu, ru, ku, _ = make_unit(c, s, float(rho_gpc(c)), n_rim=N_RIM, gap=gap)
                if not in_any_shell(Xu).all() or overlaps(Xu, ru, Xcat, radcat, gap=0.0):
                    continue
                dw = pattern(Xu, Mu)
                left = leftover(w + dw)
                Xcat = np.vstack([Xcat, Xu]) if len(Xcat) else Xu.copy()
                Mcat = np.concatenate([Mcat, Mu]); radcat = np.concatenate([radcat, ru])
                w = w + dw
                traj.append(dict(n_obj=len(Mcat), n_units=u + 1, leftover=float(left)))
                accepted = True
                break
            if not accepted:
                break
            continue
        Xcat = np.vstack([Xcat, Xu]) if len(Xcat) else Xu.copy()
        Mcat = np.concatenate([Mcat, Mu]); radcat = np.concatenate([radcat, ru])
        w = w + dw
        traj.append(dict(n_obj=len(Mcat), n_units=u + 1, leftover=float(left)))
    return Xcat, Mcat, radcat, w, traj


def growth_alternate(rng, s=S, max_obj=MAX_OBJ):
    """Alternate void / cluster. Clusters biased toward existing void rims."""
    g = cell_geometry(s)
    Xcat = np.zeros((0, 3)); Mcat = np.zeros(0); radcat = np.zeros(0)
    w = np.zeros_like(R)
    traj = [dict(n_obj=0, leftover=float(BASE0))]
    voids = []  # centers of placed voids
    next_void = True
    stuck = 0
    while len(Mcat) < max_obj and stuck < 100:
        if next_void:
            c = propose_center(rng, "full", r_void=g["R_v"])
            if c is None:
                stuck += 1; continue
            mass, rad = _mass_for_void(c, s)
            if not place_if_ok(c, mass, rad, Xcat, radcat):
                stuck += 1; continue
            dw = pattern(c[None], np.array([mass]))
            Xcat = np.vstack([Xcat, c]) if len(Xcat) else c[None].copy()
            Mcat = np.concatenate([Mcat, [mass]]); radcat = np.concatenate([radcat, [rad]])
            voids.append(c); w = w + dw
            traj.append(dict(n_obj=len(Mcat), leftover=float(leftover(w)), kind="void"))
            next_void = False; stuck = 0
        else:
            # place cluster near a random existing void (junction preference) or random
            placed = False
            for _ in range(30):
                if voids and rng.random() < 0.75:
                    v = voids[rng.integers(len(voids))]
                    u = rng.normal(size=3); u /= np.linalg.norm(u)
                    c = v + u * (g["R_v"] + g["R_c"] + GAP_FRAC * g["R_v"])
                else:
                    c = propose_center(rng, "full", r_void=g["R_c"])
                    if c is None:
                        continue
                mass, rad = _mass_for_clu(c, s)
                if not place_if_ok(c, mass, rad, Xcat, radcat):
                    continue
                if not in_any_shell(c[None])[0]:
                    continue
                dw = pattern(c[None], np.array([mass]))
                Xcat = np.vstack([Xcat, c]); Mcat = np.concatenate([Mcat, [mass]]); radcat = np.concatenate([radcat, [rad]])
                w = w + dw
                traj.append(dict(n_obj=len(Mcat), leftover=float(leftover(w)), kind="clu"))
                placed = True; break
            if not placed:
                stuck += 1
            else:
                stuck = 0
            next_void = True
        # rebalance mass every so often: scale last few clusters — skip, keep raw
    # final mass rebalance note: net may not be zero for alternate; report it
    return Xcat, Mcat, radcat, w, traj


def growth_anti_coe(rng, s=S, max_units=MAX_UNITS):
    """Bias: voids preferred on far side (anti-CoE), clusters on near side (toward CoE).
    Uses compensated units but candidate score includes directional bias + residual help."""
    g = cell_geometry(s)
    gap = GAP_FRAC * g["R_v"]
    ax = C_G / np.linalg.norm(C_G)   # us → CoE
    Xcat = np.zeros((0, 3)); Mcat = np.zeros(0); radcat = np.zeros(0)
    w = np.zeros_like(R)
    traj = [dict(n_obj=0, n_units=0, leftover=float(BASE0))]
    for u in range(max_units):
        best = None
        for _ in range(50):
            c = propose_center(rng, "full", r_void=g["R_v"])
            if c is None:
                continue
            Xu, Mu, ru, ku, _ = make_unit(c, s, float(rho_gpc(c)), n_rim=N_RIM, gap=gap)
            if not in_any_shell(Xu).all() or overlaps(Xu, ru, Xcat, radcat, gap=0.0):
                continue
            dw = pattern(Xu, Mu)
            left = leftover(w + dw)
            # bias: void center with positive projection on -ax (away from CoE) preferred
            bias = float(np.dot(c / np.linalg.norm(c), -ax))  # +1 = anti-CoE
            score = left - 8.0 * bias   # lower is better; reward anti-CoE voids
            if best is None or score < best[0]:
                best = (score, left, Xu, Mu, ru, dw)
        if best is None:
            break
        _, left, Xu, Mu, ru, dw = best
        if u >= 4 and left > traj[-1]["leftover"] + 15:
            break
        Xcat = np.vstack([Xcat, Xu]) if len(Xcat) else Xu.copy()
        Mcat = np.concatenate([Mcat, Mu]); radcat = np.concatenate([radcat, ru])
        w = w + dw
        traj.append(dict(n_obj=len(Mcat), n_units=u + 1, leftover=float(left)))
    return Xcat, Mcat, radcat, w, traj


def growth_random_scatter(rng, s=S, max_units=MAX_UNITS):
    """Control: sequential random units, first-fit (no residual greed)."""
    g = cell_geometry(s)
    gap = GAP_FRAC * g["R_v"]
    Xcat = np.zeros((0, 3)); Mcat = np.zeros(0); radcat = np.zeros(0)
    w = np.zeros_like(R)
    traj = [dict(n_obj=0, n_units=0, leftover=float(BASE0))]
    tries = 0
    u = 0
    while u < max_units and tries < max_units * 40:
        tries += 1
        c = propose_center(rng, "full", r_void=g["R_v"])
        if c is None:
            continue
        Xu, Mu, ru, ku, _ = make_unit(c, s, float(rho_gpc(c)), n_rim=N_RIM, gap=gap)
        if not in_any_shell(Xu).all() or overlaps(Xu, ru, Xcat, radcat, gap=0.0):
            continue
        dw = pattern(Xu, Mu)
        Xcat = np.vstack([Xcat, Xu]) if len(Xcat) else Xu.copy()
        Mcat = np.concatenate([Mcat, Mu]); radcat = np.concatenate([radcat, ru])
        w = w + dw
        u += 1
        traj.append(dict(n_obj=len(Mcat), n_units=u, leftover=float(leftover(w))))
    return Xcat, Mcat, radcat, w, traj


RULES = {
    "void_then_rim": growth_void_then_rim,
    "alternate": growth_alternate,
    "anti_coe_bias": growth_anti_coe,
    "random_scatter": growth_random_scatter,
}


def main():
    t0 = time.time()
    print("=== Scenario B: sequential growth trajectories ===", flush=True)
    g = cell_geometry(S)
    print(f"s={S} R_v={g['R_v']:.3f} Gpc R_c={g['R_c']:.3f} BASE0={BASE0:.2f} BASE4={BASE4:.2f}", flush=True)
    out = dict(s=S, R_v=g["R_v"], R_c=g["R_c"], BASE0=BASE0, BASE4=BASE4, rules={})
    # warm pattern
    _ = pattern(np.array([[20.0, 0.0, 0.0]]), np.array([1e50]))

    for name, fn in RULES.items():
        print(f"\nRule: {name}", flush=True)
        runs = []
        for seed in (1, 2, 3, 5, 8):
            rng = np.random.default_rng(seed)
            X, M, rad, w, traj = fn(rng)
            left = float(leftover(w))
            sc, bl = best_scale(w)
            # overlap verify
            if name in ("void_then_rim", "anti_coe_bias", "random_scatter"):
                n_u = traj[-1].get("n_units", 0)
                pairs = unit_ignore_pairs(n_u, N_RIM) if n_u else set()
            else:
                pairs = set()
            margin, ok = verify_no_overlap(X, rad, ignore_pairs=pairs) if len(X) else (np.inf, True)
            runs.append(dict(seed=seed, leftover=left, best_scaled=float(bl), scale=float(sc),
                             n_obj=len(M), net_mass=float(M.sum()) if len(M) else 0.0,
                             min_margin_gpc=float(margin), no_overlap=bool(ok), traj=traj))
            print(f"  seed {seed}: left={left:.2f}× objs={len(M)} margin={margin:.4f} ok={ok}", flush=True)
        finals = np.array([r["leftover"] for r in runs])
        out["rules"][name] = dict(
            runs=runs,
            median_final=float(np.median(finals)),
            best_final=float(finals.min()),
            worst_final=float(finals.max()),
            helped_vs_base0=int(np.sum(finals < BASE0)),
            hurt_vs_base0=int(np.sum(finals > BASE0)),
        )
        # save best trajectory layout
        ib = int(np.argmin(finals))
        np.savez(os.path.join(OUT, f"b_{name}_best.npz"),
                 leftover=finals[ib], seed=runs[ib]["seed"])

    out["elapsed_s"] = time.time() - t0
    with open(os.path.join(OUT, "scenario_b.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nDone Scenario B in {out['elapsed_s']:.0f}s", flush=True)
    for k, v in out["rules"].items():
        print(f"  {k}: med={v['median_final']:.1f} best={v['best_final']:.1f} helped={v['helped_vs_base0']}/5", flush=True)


if __name__ == "__main__":
    main()
