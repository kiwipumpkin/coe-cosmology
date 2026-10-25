"""Pantheon+ supernova test of the CoE / LTB geometry (Lambda = 0)."""
import sys, os, numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from numpy.polynomial.legendre import leggauss, legvander
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../cmb_redshift"))
import ltb, rays, models
C_KMS = 299792.458
REQ2GAL = np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
                    [0.4941094279, -0.4448296300, 0.7469822445],
                    [-0.8676661490, -0.1980763734, 0.4559837762]])

def radec2gal(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    v = np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
    g = REQ2GAL @ v
    return np.degrees(np.arctan2(g[1], g[0])) % 360, np.degrees(np.arcsin(np.clip(g[2], -1, 1)))

def unit(l, b):
    l, b = np.radians(l), np.radians(b)
    return np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])

def load(zmin=0.023, zcol="zHD", zmax=10.0):
    d = np.genfromtxt(os.path.join(HERE, "data/PantheonPlusSH0ES.dat"), names=True, dtype=None, encoding=None)
    cov = np.load(os.path.join(HERE, "data/cov.npy"))
    z = d[zcol]
    m = (z > zmin) & (z < zmax)
    l, b = radec2gal(d["RA"][m], d["DEC"][m])
    C = cov[np.ix_(m, m)]
    return dict(z=z[m], zhel=d["zHEL"][m], mb=d["m_b_corr"][m], l=l, b=b, n=unit(l, b).T, C=C, Cinv=np.linalg.inv(C),
                N=int(m.sum()), cid=d["CID"][m], survey=d["IDSURVEY"][m])

def chi2_marg(S, mu_model):
    """chi2 with the additive offset (absolute magnitude + 5 log10 H0) marginalised analytically (profile = marginal up to const)."""
    r = S["mb"] - mu_model; Ci = S["Cinv"]
    a = r @ Ci @ r; b = np.sum(Ci @ r); c = np.sum(Ci)
    return a - b * b / c, b / c

# ---------- reference FRW models (c/H0 = 1 units; offset absorbs H0)
def mu_frw(z, Om, OL):
    Ok = 1 - Om - OL
    zz = np.linspace(0, z.max() * 1.001, 4000)
    Ez = np.sqrt(Om * (1 + zz)**3 + Ok * (1 + zz)**2 + OL)
    chi = np.interp(z, zz, cumulative_trapezoid(1 / Ez, zz, initial=0))
    if abs(Ok) < 1e-8: dm = chi
    elif Ok > 0: dm = np.sinh(np.sqrt(Ok) * chi) / np.sqrt(Ok)
    else: dm = np.sin(np.sqrt(-Ok) * chi) / np.sqrt(-Ok)
    return 5 * np.log10((1 + z) * dm)

# ---------- LTB distance tables
def dl_table(p, t0, ro, Npsi=12, zmax=2.4, nt=260):
    """D_L(z) [Gpc] on Gauss-Legendre nodes in cos(psi) (psi = angle from the CoE direction)."""
    mu, _ = leggauss(Npsi); psi = np.arccos(mu)
    tq = t0 * np.geomspace(0.9995, 0.05, nt)
    zg = np.geomspace(1e-3, zmax, 300); tab = np.zeros((Npsi, zg.size))
    for k, ps in enumerate(psi):
        z, DA = rays.area_distance(p, t0, ro, ps, tq)
        ok = np.isfinite(z) & np.isfinite(DA)
        z, DA = z[ok], DA[ok]
        if z.max() < zmax: raise RuntimeError("zmax not reached: %.3f" % z.max())
        tab[k] = np.exp(np.interp(np.log(zg), np.log(z), np.log((1 + z)**2 * DA)))
    Vinv = np.linalg.inv(legvander(mu, Npsi - 1))
    return dict(mu=mu, zg=zg, tab=tab, Vinv=Vinv, t0=t0, p=p, ro=ro)

def dl_nodes(T, z):
    """(N_sn, Npsi) D_L at each SN redshift for every node direction."""
    lz = np.log(z); lzg = np.log(T["zg"])
    return np.exp(np.array([np.interp(lz, lzg, np.log(row)) for row in T["tab"]]).T)

def mu_ltb(T, DLn, nhat, lc, bc):
    """distance modulus (up to constant) for SN directions nhat given CoE direction (lc, bc)."""
    cpsi = nhat @ unit(lc, bc)
    W = legvander(cpsi, len(T["mu"]) - 1) @ T["Vinv"]
    return 5 * np.log10(np.sum(W * DLn, axis=1))

def ltb_model(D, L, MD, Om):
    """Om = 0 and pure lump -> Chris's case-B construction; otherwise open background + lump.
    MD = lump mass inside D in geometric Gpc."""
    if Om <= 0:
        return models.model_B(M=MD, D=D, L=L)
    dM = MD / ltb.Fexp(D / L)[0]
    return models.model_C(dM, Om=Om, D=D, L=L)

def fib_sphere(n):
    i = np.arange(n) + 0.5; th = np.arccos(1 - 2 * i / n); ph = np.pi * (1 + 5**0.5) * i
    return np.degrees(ph) % 360, 90 - np.degrees(th)

def best_direction(S, T, DLn, ngrid=1200, fix=None):
    from scipy.optimize import minimize
    if fix is not None:
        return fix[0], fix[1], chi2_marg(S, mu_ltb(T, DLn, S["n"], *fix))[0]
    L, B = fib_sphere(ngrid)
    c = np.array([chi2_marg(S, mu_ltb(T, DLn, S["n"], l, b))[0] for l, b in zip(L, B)])
    i = np.argmin(c)
    f = lambda v: chi2_marg(S, mu_ltb(T, DLn, S["n"], v[0], np.clip(v[1], -90, 90)))[0]
    r = minimize(f, [L[i], B[i]], method="Nelder-Mead", options=dict(xatol=0.05, fatol=1e-4))
    return r.x[0] % 360, float(np.clip(r.x[1], -90, 90)), r.fun

H0 = models.H0
def geom_from(D, LD, mc, Om):
    """D [Gpc], L/D, lump mass inside D in units of the critical mass inside D, background Omega."""
    MD = mc * H0**2 * D**3 / 2
    p, t0, ro = ltb_model(D, LD * D, MD, Om)
    return p, t0, ro

# ---------- directional statistics
DIRS = dict(CMB=(264.02, 48.25), Migkas=(280.0, -20.0), CoE=(193.17, 0.0))

def sub(S, m):
    C = S["C"][np.ix_(m, m)]
    return dict(z=S["z"][m], mb=S["mb"][m], n=S["n"][m], l=S["l"][m], b=S["b"][m], C=C, Cinv=np.linalg.inv(C), N=int(m.sum()))

def gls(X, y, Ci):
    A = X.T @ Ci @ X; beta = np.linalg.solve(A, X.T @ Ci @ y); cov = np.linalg.inv(A)
    r = y - X @ beta
    return beta, cov, r @ Ci @ r

def residual_dipole(S, mu_model):
    """fit r = a + d.n (GLS, full covariance); returns amplitude, direction, errors and delta chi2 vs no dipole."""
    y = S["mb"] - mu_model; one = np.ones(S["N"])
    b0, _, c0 = gls(one[:, None], y, S["Cinv"])
    X = np.column_stack([one, S["n"]]); b1, cv, c1 = gls(X, y, S["Cinv"])
    d = b1[1:]; A = np.linalg.norm(d); cd = cv[1:, 1:]
    sigA = np.sqrt(d @ cd @ d) / A
    l = np.degrees(np.arctan2(d[1], d[0])) % 360; b = np.degrees(np.arcsin(d[2] / A))
    # direction uncertainty from Monte Carlo of the 3-vector
    rng = np.random.default_rng(1); smp = rng.multivariate_normal(d, cd, 4000)
    ang = np.degrees(np.arccos(np.clip(smp @ d / (np.linalg.norm(smp, axis=1) * A), -1, 1)))
    return dict(A=A, sigA=sigA, l=l, b=b, dang68=float(np.percentile(ang, 68)), dchi2=c0 - c1, vec=d.tolist())

def cosmo_mu(z, q, j):
    return 5 * np.log10(z * (1 + z) * (1 + 0.5 * (1 - q) * z - (1 - q - 3 * q**2 + j) * z**2 / 6) / (1 + z))

def colin_fit(S, lc, bc, Sz=0.0262, zmax=0.8):
    """q0 = qm + qd (n.c) exp(-z/S), j0 free, offset marginalised (Colin et al. 2019 form)."""
    from scipy.optimize import minimize
    m = S["z"] < zmax; T = sub(S, m); cth = T["n"] @ unit(lc, bc); F = np.exp(-T["z"] / Sz)
    f = lambda v: chi2_marg(T, cosmo_mu(T["z"], v[0] + v[1] * cth * F, v[2]))[0]
    r1 = minimize(f, [-0.5, 0.0, 1.0], method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-6, maxiter=4000))
    f0 = lambda v: f([v[0], 0.0, v[1]])
    r0 = minimize(f0, [-0.5, 1.0], method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-6))
    # 1-sigma on qd from profile curvature
    h = 0.5; g = lambda qd: minimize(lambda v: f([v[0], qd, v[1]]), [r1.x[0], r1.x[2]], method="Nelder-Mead").fun
    curv = (g(r1.x[1] + h) + g(r1.x[1] - h) - 2 * r1.fun) / h**2
    return dict(qm=r1.x[0], qd=r1.x[1], j0=r1.x[2], sig_qd=float(np.sqrt(2 / curv)) if curv > 0 else np.nan, dchi2=r0.fun - r1.fun, N=T["N"])

def hemi_q0(S, lc, bc, zmax=0.8):
    from scipy.optimize import minimize_scalar
    out = []
    for sgn in (1, -1):
        m = (S["z"] < zmax) & (sgn * (S["n"] @ unit(lc, bc)) > 0); T = sub(S, m)
        f = lambda q: chi2_marg(T, cosmo_mu(T["z"], q, 1.0))[0]
        r = minimize_scalar(f, bounds=(-3, 3), method="bounded")
        h = 0.02; curv = (f(r.x + h) + f(r.x - h) - 2 * r.fun) / h**2
        out.append((r.x, np.sqrt(2 / curv), T["N"]))
    return out

# ---------- Migkas A(l) with A proportional to H0^2 = D_L^-2 at fixed z
def migkas_data():
    import csv
    rows = list(csv.DictReader(open(os.path.join(HERE, "../../migkas/fig6_digitized.csv"))))
    l = np.array([float(r["l_deg"]) for r in rows]); A = np.array([float(r["A_over_Aall"]) for r in rows])
    e = 0.5 * (np.array([float(r["A_err_low"]) for r in rows]) + np.array([float(r["A_err_high"]) for r in rows]))
    return l, A, e

def migkas_model(T, lc, bc, z=0.1, ell=None):
    ell = np.arange(0, 360, 1.0) if ell is None else ell
    nh = unit(ell, np.zeros_like(ell)).T
    DLn = dl_nodes(T, np.full(len(ell), z))
    dl = 10 ** (mu_ltb(T, DLn, nh, lc, bc) / 5)
    A = dl ** -2
    return ell, A / A.mean()

def migkas_chi2(T, lc, bc, z=0.1):
    l, A, e = migkas_data()
    _, Am = migkas_model(T, lc, bc, z, ell=l)
    w = 1 / e**2; k = np.sum(w * A * Am) / np.sum(w * Am * Am)      # free normalisation
    return float(np.sum(w * (A - k * Am) ** 2)), float(np.corrcoef(A, Am)[0, 1])
