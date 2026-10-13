"""Scenario A: fit non-overlapping void+cluster units in the local tidal region to cancel
the CMB residual, FREEZE that arrangement, then fill the outer shell (30–42.9 Gpc) with
empty / random / foam-like / lattice configs and measure how much our leftover moves.
"""
import os, sys, time, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from seq_core import *

OUT = HERE
S_LOCAL = 4.0       # coarse units that can move the needle (void deficit ~3e51 kg)
S_OUTER = 2.5       # outer packing cell size
N_RIM = 6
N_GREEDY = 12       # max local units
N_CAND = 80         # candidates per greedy step
N_OUTER_REAL = 24   # realizations per outer config (except empty/lattice variants)


def greedy_local(rng, s=S_LOCAL, n_units=N_GREEDY, n_cand=N_CAND, n_rim=N_RIM, target=None):
    """Place units one at a time; each step pick the non-overlapping candidate that
    most reduces leftover RMS. Returns frozen local layout + pattern + trajectory."""
    if target is None:
        target = R
    g = cell_geometry(s)
    gap = GAP_FRAC * g["R_v"]
    Xcat = np.zeros((0, 3)); Mcat = np.zeros(0); radcat = np.zeros(0); kinds = np.zeros(0, dtype=object)
    w = np.zeros_like(target)
    traj = [dict(n_units=0, leftover=float(leftover(w, target)), best_scaled_left=float(leftover(w, target)), scale=0.0)]
    for step in range(n_units):
        best = None  # (left, Xu, Mu, ru, ku, dw)
        for _ in range(n_cand):
            c = propose_center(rng, "local", r_void=g["R_v"])
            if c is None:
                continue
            rho = float(rho_gpc(c))
            Xu, Mu, ru, ku, _ = make_unit(c, s, rho, n_rim=n_rim, gap=gap)
            if not in_local(Xu).all():
                continue
            if np.linalg.norm(Xu[0]) < RMIN + g["R_v"]:
                continue
            if overlaps(Xu, ru, Xcat, radcat, gap=0.0):
                continue
            dw = pattern(Xu, Mu)
            left = leftover(w + dw, target)
            if best is None or left < best[0]:
                best = (left, Xu, Mu, ru, ku, dw)
        if best is None:
            break
        # accept only if it helps (or first unit even if mild hurt — skip if hurts)
        prev = traj[-1]["leftover"]
        if best[0] >= prev - 0.05 and step > 0:
            # try a few more candidates with larger search before giving up
            improved = False
            for _ in range(n_cand):
                c = propose_center(rng, "local", r_void=g["R_v"])
                if c is None:
                    continue
                rho = float(rho_gpc(c))
                Xu, Mu, ru, ku, _ = make_unit(c, s, rho, n_rim=n_rim, gap=gap)
                if not in_local(Xu).all() or np.linalg.norm(Xu[0]) < RMIN + g["R_v"]:
                    continue
                if overlaps(Xu, ru, Xcat, radcat, gap=0.0):
                    continue
                dw = pattern(Xu, Mu)
                left = leftover(w + dw, target)
                if left < prev - 0.05:
                    best = (left, Xu, Mu, ru, ku, dw); improved = True; break
            if not improved:
                break
        left, Xu, Mu, ru, ku, dw = best
        Xcat = np.vstack([Xcat, Xu]) if len(Xcat) else Xu.copy()
        Mcat = np.concatenate([Mcat, Mu]); radcat = np.concatenate([radcat, ru]); kinds = np.concatenate([kinds, ku])
        w = w + dw
        sc, bl = best_scale(w, target)
        traj.append(dict(n_units=step + 1, leftover=float(left), best_scaled_left=float(bl), scale=float(sc),
                         n_obj=len(Mcat), void_center=Xu[0].tolist()))
        print(f"  greedy step {step+1}: left={left:.2f}×  best-scaled={bl:.2f}× (a={sc:.3f})  objs={len(Mcat)}", flush=True)
    # overlap check (ignore void-cluster contact inside units)
    n_u = len(traj) - 1
    pairs = unit_ignore_pairs(n_u, n_rim)
    margin, ok = verify_no_overlap(Xcat, radcat, gap=0.0, ignore_pairs=pairs)
    info = dict(s=s, R_v=g["R_v"], R_c=g["R_c"], n_units=n_u, n_obj=len(Mcat),
                leftover=float(leftover(w, target)), min_margin_gpc=margin, no_overlap=ok,
                net_mass=float(Mcat.sum()) if len(Mcat) else 0.0)
    return Xcat, Mcat, radcat, kinds, w, traj, info


def outer_empty(local):
    Xl, Ml, rl, kl = local
    return Xl.copy(), Ml.copy(), rl.copy(), kl.copy(), dict(n_units=0, tag="empty")


def run_outer_realizations(local, w_local, tag, pack_fn, nreal, seed0, s=S_OUTER):
    Xl, Ml, rl, kl = local
    existing = (Xl, Ml, rl, kl)
    rows = []
    t0 = time.time()
    for i in range(nreal):
        rng = np.random.default_rng([seed0, i])
        X, M, rad, kinds, info = pack_fn(rng, s, "outer", existing=existing)
        # pattern of NEW objects only
        n0 = len(Ml)
        if len(M) > n0:
            w_new = pattern(X[n0:], M[n0:])
        else:
            w_new = np.zeros_like(w_local)
        w = w_local + w_new
        left = leftover(w)
        sc, bl = best_scale(w)
        helped = left < leftover(w_local)
        rows.append(dict(left=float(left), best_scaled=float(bl), scale=float(sc),
                         helped=bool(helped), n_outer_units=int(info.get("n_units", 0)),
                         n_obj=int(len(M)), rms_w_outer=float(rms(w_new) / FLOOR)))
        if (i + 1) % 6 == 0 or i == nreal - 1:
            L = np.array([r["left"] for r in rows])
            print(f"  {tag} {i+1}/{nreal} med={np.median(L):.2f} helped={np.mean(L<leftover(w_local)):.2f} t={time.time()-t0:.0f}s", flush=True)
    return rows


def main():
    t0 = time.time()
    print("=== Scenario A: local dipole/residual fix, then outer fills ===", flush=True)
    print(f"BASE0={BASE0:.2f} BASE4={BASE4:.2f} S_local={S_LOCAL} S_outer={S_OUTER}", flush=True)
    gL = cell_geometry(S_LOCAL); gO = cell_geometry(S_OUTER)
    print(f"Local unit: R_v={gL['R_v']:.3f} Gpc R_c={gL['R_c']:.3f} Gpc", flush=True)
    print(f"Outer unit: R_v={gO['R_v']:.3f} Gpc R_c={gO['R_c']:.3f} Gpc", flush=True)

    # --- local greedy (several seeds, keep best leftover) ---
    best_pack = None
    all_local = []
    for seed in (1, 2, 3, 7, 11):
        print(f"\nLocal greedy seed={seed}", flush=True)
        rng = np.random.default_rng(seed)
        X, M, rad, kinds, w, traj, info = greedy_local(rng)
        all_local.append(dict(seed=seed, traj=traj, info=info))
        print(f"  -> units={info['n_units']} left={info['leftover']:.2f}× margin={info['min_margin_gpc']:.4f} ok={info['no_overlap']}", flush=True)
        if best_pack is None or info["leftover"] < best_pack["info"]["leftover"]:
            best_pack = dict(X=X, M=M, rad=rad, kinds=kinds, w=w, traj=traj, info=info, seed=seed)

    local = (best_pack["X"], best_pack["M"], best_pack["rad"], best_pack["kinds"])
    w_local = best_pack["w"]
    left0 = leftover(w_local)
    sc0, bl0 = best_scale(w_local)
    print(f"\nFROZEN local (seed {best_pack['seed']}): left={left0:.2f}× best-scaled={bl0:.2f}× (a={sc0:.3f})", flush=True)
    print(f"  units={best_pack['info']['n_units']} objs={best_pack['info']['n_obj']} net_mass={best_pack['info']['net_mass']:.3e}", flush=True)

    # save local snapshot
    np.savez(os.path.join(OUT, "local_frozen.npz"),
             X=best_pack["X"], M=best_pack["M"], rad=best_pack["rad"],
             kinds=best_pack["kinds"], w=w_local, seed=best_pack["seed"])

    # --- outer configs ---
    results = dict(
        geometry=dict(local_gpc=[RMIN, LOCAL_RMAX], outer_gpc=[LOCAL_RMAX, S_MAX_G],
                      S_local=S_LOCAL, S_outer=S_OUTER, R_v_local=gL["R_v"], R_v_outer=gO["R_v"],
                      R_c_local=gL["R_c"], R_c_outer=gO["R_c"], n_rim=N_RIM,
                      note="Compensated units: void δ share −0.546 M_cell + rim clusters holding +0.546 M_cell (wall+fil folded into clusters). Hard non-overlap except void–rim contact inside a unit."),
        baseline=dict(BASE0=BASE0, BASE4=BASE4, FLOOR=FLOOR),
        local=dict(seed=best_pack["seed"], leftover=left0, best_scaled=bl0, scale=sc0,
                   info=best_pack["info"], traj=best_pack["traj"], all_seeds=all_local),
        outer={},
    )

    # empty control
    results["outer"]["empty"] = dict(
        left=left0, best_scaled=bl0, scale=sc0, n_outer_units=0,
        note="frozen local only; outer empty",
    )
    print(f"\nOuter empty (control): left={left0:.2f}×", flush=True)

    # random packing
    print("\nOuter random non-overlap:", flush=True)
    def _rand(rng, s, region, existing=None):
        return pack_random(rng, s, region, n_rim=N_RIM, existing=existing, n_tries_factor=30)
    rows = run_outer_realizations(local, w_local, "random", _rand, N_OUTER_REAL, seed0=100)
    L = np.array([r["left"] for r in rows])
    results["outer"]["random"] = dict(
        nreal=len(rows), median=float(np.median(L)), mean=float(L.mean()),
        p10=float(np.percentile(L, 10)), p90=float(np.percentile(L, 90)),
        min=float(L.min()), max=float(L.max()),
        helped=int(np.sum(L < left0)), hurt=int(np.sum(L > left0)),
        helped_vs_base0=int(np.sum(L < BASE0)),
        median_outer_units=float(np.median([r["n_outer_units"] for r in rows])),
        rows=rows,
    )

    # foam-like
    print("\nOuter foam-like sequential:", flush=True)
    def _foam(rng, s, region, existing=None):
        return pack_foam_like(rng, s, region, n_rim=N_RIM, existing=existing)
    rows = run_outer_realizations(local, w_local, "foam", _foam, N_OUTER_REAL, seed0=200)
    L = np.array([r["left"] for r in rows])
    results["outer"]["foam"] = dict(
        nreal=len(rows), median=float(np.median(L)), mean=float(L.mean()),
        p10=float(np.percentile(L, 10)), p90=float(np.percentile(L, 90)),
        min=float(L.min()), max=float(L.max()),
        helped=int(np.sum(L < left0)), hurt=int(np.sum(L > left0)),
        helped_vs_base0=int(np.sum(L < BASE0)),
        median_outer_units=float(np.median([r["n_outer_units"] for r in rows])),
        rows=rows,
    )

    # lattice (several phases)
    print("\nOuter lattice:", flush=True)
    rows = []
    for i, seed in enumerate(range(300, 300 + N_OUTER_REAL)):
        rng = np.random.default_rng(seed)
        X, M, rad, kinds, info = pack_lattice(S_OUTER, "outer", n_rim=N_RIM, existing=local, rng=rng)
        n0 = len(local[1])
        w_new = pattern(X[n0:], M[n0:]) if len(M) > n0 else np.zeros_like(w_local)
        w = w_local + w_new
        left = leftover(w); sc, bl = best_scale(w)
        rows.append(dict(left=float(left), best_scaled=float(bl), scale=float(sc),
                         helped=bool(left < left0), n_outer_units=int(info["n_units"]),
                         n_obj=int(len(M)), rms_w_outer=float(rms(w_new) / FLOOR)))
        if (i + 1) % 6 == 0:
            print(f"  lattice {i+1}/{N_OUTER_REAL} left={left:.2f} units={info['n_units']}", flush=True)
    L = np.array([r["left"] for r in rows])
    results["outer"]["lattice"] = dict(
        nreal=len(rows), median=float(np.median(L)), mean=float(L.mean()),
        p10=float(np.percentile(L, 10)), p90=float(np.percentile(L, 90)),
        min=float(L.min()), max=float(L.max()),
        helped=int(np.sum(L < left0)), hurt=int(np.sum(L > left0)),
        helped_vs_base0=int(np.sum(L < BASE0)),
        median_outer_units=float(np.median([r["n_outer_units"] for r in rows])),
        rows=rows,
    )

    results["elapsed_s"] = time.time() - t0
    # delta summary
    results["outer_vs_frozen"] = {
        k: dict(median_delta=float(results["outer"][k]["median"] - left0),
                median=results["outer"][k]["median"],
                helped=results["outer"][k]["helped"],
                hurt=results["outer"][k]["hurt"])
        for k in ("random", "foam", "lattice")
    }
    results["outer_vs_frozen"]["empty"] = dict(median_delta=0.0, median=left0, helped=0, hurt=0)

    with open(os.path.join(OUT, "scenario_a.json"), "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nDone Scenario A in {results['elapsed_s']:.0f}s → scenario_a.json", flush=True)
    print("Outer vs frozen:", json.dumps(results["outer_vs_frozen"], indent=1), flush=True)
    return results


if __name__ == "__main__":
    main()
