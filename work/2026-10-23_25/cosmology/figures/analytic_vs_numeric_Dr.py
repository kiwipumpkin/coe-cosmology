#!/usr/bin/env python3
"""Analytic (small r/D tidal expansion) vs NUMERIC (ideal-drag ST + time-integrated Gauss gravity,
exp core + floor, 90-deg sectors) CoE model for Migkas+2020 Fig. 6 A(l), at D/r = 0.9, 1.01, 1.3, 2, 3, 4.
Sources on the vision shell: d = r, eps = r/D = 1/(D/r).  Pure physics: no rotation / 2l phase transforms.

NUMERIC: r2_vs_Dr_density_families.py model (family a_exp_floor, kappa free, non-strict), params from
  R2_vs_Dr_1to10.json (D/r = 2, 3, 4; stored grid r in [0.02, 2] D, 80 pts) and from
  R2_numeric_Dr_0p9_1p01_1p3.json (fresh fits by fit_numeric_Dr_0p9_1p01_1p3.py; r grid [0.001, 2] D, 600 pts).
ANALYTIC: A/A_all = N[1 + K(3mu^2-1) + K eps(9/2 mu - 15/2 mu^3)], mu = n.D_hat (D_hat from CoE to us, same
  sign convention as the numeric code; this is the exact point-mass tidal series to O(eps)),
  mu = cos b cos b0 cos(l-l0) + sin b sin b0, averaged over 90-deg longitude sectors, |b|>20 cut, cos b weight
  (closed-form latitude moments as in papers/A_equals_fourier/solve_v2.py).
  (a) matched: same l0, b0 as numeric; N fixed by mean(A)=1 (numeric curve is mean-renormalised);
      K chosen so the analytic m=2 coefficient along cos 2(l-l0) equals the numeric curve's.
  (b) free: N, K, l0, b0 least-squares fit to the 71 data points at fixed eps.
  (c) diagnostic (table + thin dotted line): same expansion but evaluated in the numeric code's own sky
      geometry (single great circle mu = cos b0 cos(l-l0), same 181-pt 93.5-deg boxcar + renormalisation),
      m=2 matched the same way -> isolates the expansion error from the latitude-averaging difference.
  (d) reference (table only): the numeric H_obs pipeline fed a pure point-mass field u = -c/r^2 (c = 0.02, linear
      regime) at the same l0, b0, vs (c): pure truncation error of the O(eps) series, no profile effects.
      (Verified: -> 0 as eps -> 0, e.g. 0.001 at D/r = 10, b0 = 0, which also confirms the sign convention.)"""
import json, sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import least_squares
sys.path.insert(0, "/workspace/cosmology/figures")
import r2_vs_Dr_density_families as M
for f in ("a_exp_floor", "d_pure_exp"): M.FAM[f]["hi"][0] = 300.0
OUT = "/workspace/cosmology/figures/analytic_vs_numeric_Dr"
DRS = [0.9, 1.01, 1.3, 2.0, 3.0, 4.0]
d2r = np.pi / 180
dat = np.genfromtxt(M.CSV, delimiter=",", names=True)
l_d, A_d, elo, ehi = dat["l_deg"], dat["A_over_Aall"], dat["A_err_low"], dat["A_err_high"]
r2 = lambda p: float(1 - np.sum((A_d - p) ** 2) / np.sum((A_d - A_d.mean()) ** 2))

# ---------------- numeric params ----------------
stored = [r for r in json.load(open("/workspace/cosmology/figures/R2_vs_Dr_1to10.json"))["results"]
          if r["family"] == "a_exp_floor" and not r["strict_local_density"]]
new = json.load(open("/workspace/cosmology/figures/R2_numeric_Dr_0p9_1p01_1p3.json"))["results"]
def num_entry(Dr):
    for r in new:
        if abs(r["D_over_r"] - Dr) < 1e-9: return r, (0.001, 600), "fresh fit (R2_numeric_Dr_0p9_1p01_1p3.json)"
    for r in stored:
        if abs(r["D_over_r"] - Dr) < 1e-9: return r, (0.02, 80), "R2_vs_Dr_1to10.json"
    raise KeyError(Dr)

# ---------------- analytic, latitude-averaged closed form ----------------
bc, hp = 20 * d2r, np.pi / 2
Fi = lambda f: f(hp) - f(bc); nrm = Fi(np.sin)
c1 = Fi(lambda b: b / 2 + np.sin(2 * b) / 4) / nrm
c2 = Fi(lambda b: np.sin(b) - np.sin(b) ** 3 / 3) / nrm
c3 = Fi(lambda b: 3 * b / 8 + np.sin(2 * b) / 4 + np.sin(4 * b) / 32) / nrm
s2 = Fi(lambda b: np.sin(b) ** 3 / 3) / nrm
cs2 = Fi(lambda b: b / 8 - np.sin(4 * b) / 32) / nrm
W = {m: np.sin(m * np.pi / 4) / (m * np.pi / 4) for m in (1, 2, 3)}
def coefs(b0):
    a, s0 = np.cos(b0 * d2r), np.sin(b0 * d2r)
    return dict(q0=1.5 * a * a * c2 + 3 * s0 * s0 * s2 - 1, q2=1.5 * a * a * c2 * W[2],
                o1=a * W[1] * (4.5 * c1 - 45 / 8 * a * a * c3 - 22.5 * s0 * s0 * cs2), o3=-15 / 8 * a ** 3 * c3 * W[3])
def A_lat(l, N, K, l0, b0, eps):
    t = coefs(b0); p = (l - l0) * d2r
    return N * (1 + K * (t["q0"] + t["q2"] * np.cos(2 * p)) + K * eps * (t["o1"] * np.cos(p) + t["o3"] * np.cos(3 * p)))
def brute(l, N, K, l0, b0, eps):   # direct 2-D sector average, for a sanity check of the closed form
    bg = np.radians(np.arange(-89.875, 90, 0.25)); bg = bg[np.abs(bg) > bc]
    Lo = np.radians(np.arange(-44.875, 45, 0.25)); B, Lg = np.meshgrid(bg, Lo, indexing="ij"); w = np.cos(B)
    out = []
    for lc in np.atleast_1d(l):
        mu = np.cos(B) * np.cos(b0 * d2r) * np.cos((lc - l0) * d2r + Lg) + np.sin(B) * np.sin(b0 * d2r)
        out.append(np.sum(w * N * (1 + K * (3 * mu ** 2 - 1) + K * eps * (4.5 * mu - 7.5 * mu ** 3))) / w.sum())
    return np.array(out)
_chk = (np.array([0., 47., 200.]), 1.0, 0.07, 30.0, 25.0, 0.8)
assert np.allclose(A_lat(*_chk), brute(*_chk), atol=2e-5), (A_lat(*_chk), brute(*_chk))

# ---------------- analytic in the numeric code's ring geometry (diagnostic) ----------------
def A_ring(K, l0, b0, eps):
    mu = np.cos(b0 * d2r) * np.cos((M.LG - l0) * d2r)
    y = 1 + K * (3 * mu ** 2 - 1) + K * eps * (4.5 * mu - 7.5 * mu ** 3)
    y = y / y.mean(); y = M.boxcar_dl90(y); return y / y.mean()

def c2_along(y, l0):   # m=2 coefficient along cos 2(l-l0) of a curve on M.LG
    return float(2 * np.mean(y * np.cos(2 * (M.LG - l0) * d2r)))

rows, curves = [], {}
for Dr in DRS:
    eps = 1 / Dr
    e, grid, src = num_entry(Dr)
    M.RG_LO, M.NRAD = grid
    AF, kap, kb = M.evaluate("a_exp_floor", eps, np.array(e["x"]))
    l0, b0 = e["x"][1] % 360, e["x"][2]
    num_d = M.at_data(AF); R2n = r2(num_d)
    # (a) matched K
    t = coefs(b0); c2n = c2_along(AF, l0)
    K = c2n / (t["q2"] - c2n * t["q0"]); N = 1 / (1 + K * t["q0"])
    an_fine = A_lat(M.LG, N, K, l0, b0, eps); R2a = r2(A_lat(l_d, N, K, l0, b0, eps))
    assert abs(c2_along(an_fine, l0) - c2n) < 1e-9
    # (c) ring geometry, matched K (1-D root find; m2 of ring curve is monotone in K for small K)
    from scipy.optimize import brentq
    g = lambda k: c2_along(A_ring(k, l0, b0, eps), l0) - c2n
    Kr = brentq(g, -0.3, 0.3) if g(-0.3) * g(0.3) < 0 else np.nan
    ring = A_ring(Kr, l0, b0, eps); R2r = r2(M.at_data(ring))
    rgp = np.linspace(0.001, 2.2, 4000); AFp = M.model_fine(1.0, eps, l0, b0, (rgp, -0.02 / rgp ** 2)); c2p = c2_along(AFp, l0)
    try:
        Kp = brentq(lambda k: c2_along(A_ring(k, l0, b0, eps), l0) - c2p, -0.5, 0.5)
        frac_pm = float(np.sqrt(np.mean((A_ring(Kp, l0, b0, eps) - AFp) ** 2)) / np.std(AFp))
    except ValueError:
        frac_pm = float("nan")
    # (b) free fit
    best = None
    for l0s in np.arange(0, 360, 30):
        for b0s in (-60, -30, 0, 30, 60):
            for Ks in (0.05, -0.05):
                res = least_squares(lambda p: A_lat(l_d, p[0], p[1], p[2], p[3], eps) - A_d, [1, Ks, l0s, b0s],
                                    bounds=([0.5, -2, -720, -90], [1.5, 2, 720, 90]))
                if best is None or res.cost < best.cost: best = res
    Nf, Kf, l0f, b0f = best.x; l0f %= 360
    free_fine = A_lat(M.LG, Nf, Kf, l0f, b0f, eps); R2f = r2(A_lat(l_d, Nf, Kf, l0f, b0f, eps))
    aniso = np.sqrt(np.mean((AF - AF.mean()) ** 2))
    rms_an = np.sqrt(np.mean((an_fine - AF) ** 2)); rms_ring = np.sqrt(np.mean((ring - AF) ** 2))
    # harmonic content of numeric vs analytic (m=1,3 along cos m(l-l0))
    cm = lambda y, m: float(2 * np.mean(y * np.cos(m * (M.LG - l0) * d2r)))
    row = dict(D_over_r=Dr, eps=eps, numeric_source=src, numeric_r_grid=grid,
               kappa=float(kap), D_over_L=e["density"]["D_over_L"], f_u=e["density"]["f_u"], l0_deg=l0, b0_deg=b0,
               R2_numeric=R2n, R2_numeric_stored=e["R2_direct"], R2_analytic_matched=R2a, R2_analytic_free=R2f,
               R2_analytic_ring_matched=R2r,
               K_matched=float(K), N_matched=float(N), K_ring_matched=float(Kr),
               free=dict(N=float(Nf), K=float(Kf), l0_deg=float(l0f), b0_deg=float(b0f)),
               rms_numeric_anisotropy=float(aniso),
               abs_rms_diff_analytic_vs_numeric=float(rms_an), frac_rms_diff_analytic_vs_numeric=float(rms_an / aniso),
               frac_rms_diff_ring_analytic_vs_numeric=float(rms_ring / aniso),
               frac_rms_pointmass_reference_ring=frac_pm,
               numeric_c1_c2_c3=[cm(AF, 1), cm(AF, 2), cm(AF, 3)],
               analytic_matched_c1_c2_c3=[cm(an_fine, 1), cm(an_fine, 2), cm(an_fine, 3)],
               ring_analytic_c1_c2_c3=[cm(ring, 1), cm(ring, 2), cm(ring, 3)])
    rows.append(row); curves[Dr] = (AF.copy(), an_fine, free_fine, ring)

# ---------------- table ----------------
hdr = f"{'D/r':>5} {'eps':>6} {'R2_num':>7} {'R2_anM':>7} {'R2_anF':>7} {'R2_ring':>7} {'fracRMS':>8} {'fracRMS_ring':>12} {'K_match':>8} {'K_free':>7} {'fracRMS_ptmass':>14}"
print(hdr)
for r in rows:
    print(f"{r['D_over_r']:5.2f} {r['eps']:6.3f} {r['R2_numeric']:7.4f} {r['R2_analytic_matched']:7.4f} {r['R2_analytic_free']:7.4f} "
          f"{r['R2_analytic_ring_matched']:7.4f} {r['frac_rms_diff_analytic_vs_numeric']:8.3f} {r['frac_rms_diff_ring_analytic_vs_numeric']:12.3f} "
          f"{r['K_matched']:8.4f} {r['free']['K']:7.4f} {r['frac_rms_pointmass_reference_ring']:14.3f}")
for r in rows:
    print(f"D/r={r['D_over_r']}: l0={r['l0_deg']:.2f} b0={r['b0_deg']:.2f} kappa={r['kappa']:.4f} D/L={r['D_over_L']:.3g} f_u={r['f_u']:.3g} | "
          f"c1,c2,c3 num={np.round(r['numeric_c1_c2_c3'],4)} anM={np.round(r['analytic_matched_c1_c2_c3'],4)} ring={np.round(r['ring_analytic_c1_c2_c3'],4)} | free={ {k: round(v,3) for k,v in r['free'].items()} }")

# ---------------- figure ----------------
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "legend.fontsize": 9.5})
fig, axs = plt.subplots(2, 3, figsize=(16, 9.6), sharex=True, sharey=True)
for ax, r in zip(axs.flat, rows):
    Dr = r["D_over_r"]; AF, an, fr, ring = curves[Dr]
    Lp = np.append(M.LG, 360.0); cyc = lambda y: np.append(y, y[0])
    ax.errorbar(l_d, A_d, yerr=[elo, ehi], fmt="o", ms=3.5, color="0.35", ecolor="0.7", elinewidth=0.9, capsize=0,
                zorder=1, label="Migkas+20 Fig. 6")
    ax.plot(Lp, cyc(AF), "-", color="C0", lw=2.4, label=f"numeric  $R^2$={r['R2_numeric']:.3f}", zorder=3)
    ax.plot(Lp, cyc(an), "-", color="C3", lw=1.8, label=f"analytic, m=2-matched  $R^2$={r['R2_analytic_matched']:.3f}", zorder=4)
    ax.plot(Lp, cyc(fr), "--", color="C2", lw=1.8, label=f"analytic, free fit  $R^2$={r['R2_analytic_free']:.3f}", zorder=4)
    ax.plot(Lp, cyc(ring), ":", color="k", lw=1.3, label=f"analytic, ring geom., matched  $R^2$={r['R2_analytic_ring_matched']:.3f}", zorder=5)
    ax.set_title(f"D/r = {Dr:g}   (r/D = {r['eps']:.3f})")
    ax.text(0.02, 0.03, f"κ={r['kappa']:.3f}, D/L={r['D_over_L']:.3g}, f$_u$={r['f_u']:.2g}\n"
                        f"ℓ$_0$={r['l0_deg']:.1f}°, b$_0$={r['b0_deg']:.1f}°, K$_{{match}}$={r['K_matched']:.4f}\n"
                        f"frac. RMS(analytic−numeric)={r['frac_rms_diff_analytic_vs_numeric']:.2f}",
            transform=ax.transAxes, fontsize=9, va="bottom", bbox=dict(fc="white", ec="0.8", alpha=0.85))
    ax.legend(loc="upper right", framealpha=0.85, handlelength=2.2)
    ax.grid(alpha=0.25); ax.set_xlim(0, 360); ax.set_xticks(range(0, 361, 60))
for ax in axs[1]: ax.set_xlabel("Galactic longitude ℓ [deg]")
for ax in axs[:, 0]: ax.set_ylabel("A / A$_{all}$")
lo = min(np.min(A_d - elo), min(np.min(c) for v in curves.values() for c in v)); hi = max(np.max(A_d + ehi), max(np.max(c) for v in curves.values() for c in v))
axs[0, 0].set_ylim(lo - 0.02, hi + 0.2)
fig.suptitle("Sedov–Taylor CoE: analytic tidal expansion vs numeric ideal-drag model, sources on vision shell (d = r)\n"
             "Pure physics: no rotation, no 2ℓ phase transform.  Same κ, ℓ$_0$, b$_0$ in both; analytic K set so its m=2 amplitude equals the numeric curve's",
             fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.955))
fig.savefig(OUT + ".png", dpi=150); fig.savefig(OUT + ".pdf")
json.dump(dict(doc=__doc__, frac_rms_definition="RMS_l(A_analytic - A_numeric) / RMS_l(A_numeric - mean), over 0-360 on the 181-pt model grid",
               lat_moments=dict(c1=c1, c2=c2, c3=c3, s2=s2, cs2=cs2), rows=rows,
               curves={str(k): dict(l_deg=list(M.LG), numeric=list(v[0]), analytic_matched=list(v[1]), analytic_free=list(v[2]),
                                    analytic_ring_matched=list(v[3])) for k, v in curves.items()}),
          open(OUT + ".json", "w"), indent=1, default=float)
print("saved", OUT + ".png/.pdf/.json")
