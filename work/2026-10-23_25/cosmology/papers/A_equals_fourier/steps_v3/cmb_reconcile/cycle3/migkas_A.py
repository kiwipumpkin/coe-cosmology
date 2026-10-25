"""Rebuild the Migkas 2020 Fig. 6 profile (best-fit L_X-T normalisation A in Delta l = 90 deg longitude bins, all latitudes, over A_all)
from the per-cluster catalogue, for the full sample and for redshift slices. Relation: log10(L44/E(z)) = log10 A + B log10(T/4 keV),
weighted least squares with intrinsic scatter tuned to reduced chi2 = 1 (as Migkas 2020 Sect. 4), E(z) flat LCDM Om=0.3.
Bootstrap errors on A(l) and on the Fourier harmonic amplitudes."""
import numpy as np
from migkas_data import load
C = load()
Ez = lambda z: np.sqrt(0.3 * (1 + z)**3 + 0.7)
Y = np.log10(C["LX"] / 100 / Ez(C["z"])); X = np.log10(C["T"] / 4.0)
sY = 0.4343 * C["e_LX"] / 100; sX = 0.4343 * 0.5 * (C["E_T"] + C["e_T"]) / C["T"]
def fitAB(idx, B=None):
    y, x, sy, sx = Y[idx], X[idx], sY[idx], sX[idx]
    Bf = 2.0 if B is None else B; sint = 0.2
    for _ in range(30):
        w = 1 / (sy**2 + Bf**2 * sx**2 + sint**2)
        if B is None:
            M = np.c_[np.ones_like(x), x]; coef = np.linalg.solve(M.T @ (M * w[:, None]), M.T @ (w * y)); a, Bf = coef
        else: a = np.sum(w * (y - Bf * x)) / np.sum(w)
        r = y - a - Bf * x; chi = np.sum(r**2 / (sy**2 + Bf**2 * sx**2 + sint**2)) / max(len(y) - 2, 1)
        # adjust intrinsic scatter toward reduced chi2 = 1
        tot = np.var(r) - np.mean(sy**2 + Bf**2 * sx**2); sint = np.sqrt(max(tot, 1e-4))
    return 10**a, Bf
LC = np.arange(0, 360, 5.0)
def profile(sel, B=None, halfwidth=45.0, minn=5):
    Aall, Ball = fitAB(sel, B)
    out = np.full(len(LC), np.nan)
    for i, lc in enumerate(LC):
        d = np.abs((C["GLON"] - lc + 180) % 360 - 180); idx = sel & (d <= halfwidth)
        if idx.sum() >= minn: out[i] = fitAB(idx, B)[0] / Aall
    return out
def harm(y, lc=LC):
    t = np.radians(lc); ok = np.isfinite(y); y = y[ok]; t = t[ok]
    X = np.c_[np.ones_like(t)] ; cols = [np.ones_like(t)]
    for m in (1, 2, 3): cols += [np.cos(m * t), np.sin(m * t)]
    M = np.c_[tuple(cols)]; c, *_ = np.linalg.lstsq(M, y, rcond=None)
    return dict(C0=c[0], amp1=np.hypot(c[1], c[2]), amp2=np.hypot(c[3], c[4]), amp3=np.hypot(c[5], c[6]),
                phase2_deg=(np.degrees(np.arctan2(c[4], c[3])) / 2) % 180, coef=c.tolist())
def bootstrap(sel, n=300, seed=0, B=None):
    rng = np.random.default_rng(seed); idx0 = np.where(sel)[0]; profs = []; hs = []
    for _ in range(n):
        s = np.zeros(len(sel), bool); pick = rng.choice(idx0, len(idx0), replace=True)
        # bootstrap with multiplicity: duplicate by jittering index list -> emulate by weights: simple approach, resample unique set
        cnt = np.bincount(pick, minlength=len(sel)); s = cnt > 0
        p = profile_w(cnt, B); profs.append(p); hs.append(harm(p))
    P = np.array(profs)
    return np.nanstd(P, axis=0), {k: float(np.std([h[k] for h in hs])) for k in ("amp1", "amp2", "amp3")}
def fitAB_w(idx, cnt, B=None):
    y, x, sy, sx, m = Y[idx], X[idx], sY[idx], sX[idx], cnt[idx]
    Bf = 2.0 if B is None else B; sint = 0.2
    for _ in range(30):
        w = m / (sy**2 + Bf**2 * sx**2 + sint**2)
        if B is None:
            M = np.c_[np.ones_like(x), x]; a, Bf = np.linalg.solve(M.T @ (M * w[:, None]), M.T @ (w * y))
        else: a = np.sum(w * (y - Bf * x)) / np.sum(w)
        r = y - a - Bf * x; tot = np.average(r**2, weights=m) - np.average(sy**2 + Bf**2 * sx**2, weights=m); sint = np.sqrt(max(tot, 1e-4))
    return 10**a, Bf
def profile_w(cnt, B=None, halfwidth=45.0, minn=5):
    sel = cnt > 0; Aall, _ = fitAB_w(sel, cnt, B); out = np.full(len(LC), np.nan)
    for i, lc in enumerate(LC):
        d = np.abs((C["GLON"] - lc + 180) % 360 - 180); idx = sel & (d <= halfwidth)
        if cnt[idx].sum() >= minn: out[i] = fitAB_w(idx, cnt, B)[0] / Aall
    return out
