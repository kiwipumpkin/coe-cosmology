"""Cycle 4 helpers: Chris's Newtonian integrator (prof_core.integrate, unchanged) with a general interior profile inside D-R plus the
Omega_s = 0.049 shell (headline rule R_in = D - 1.028 R, R_out = D + 1.475 R). Returns cluster-scale a_f and CMB-scale diagnostics."""
import sys, os, numpy as np
sys.path.insert(0, "../termslice"); import tcore as T          # sets up Chris's modules (ptshell, prof_core)
P = T.P
from scipy.integrate import cumulative_trapezoid
RG = P.RG; RV = T.RV; TH1 = T.TH1; CKM = 299792.458; HM = 71.6
ST53 = P.ST["st53"]
def inner_shape(kind, r, rc, a):
    L = 0.005
    if kind == "point":   return np.exp(-r / L)
    if kind == "pl":      return (r * r + L * L) ** (-a["n"] / 2) * (r <= rc)           # rising toward CoE (n>0), uniform n=0
    if kind == "st":      xi, G = ST53; return np.interp(r / rc, xi, G, left=0.0, right=0.0)
    if kind == "hollow":  return np.exp(-0.5 * ((r - a["rh"] * rc) / (0.03 * rc)) ** 2)
    raise ValueError(kind)
def shape_mass(kind, rc, a):
    s = inner_shape(kind, RG, rc, a); m = cumulative_trapezoid(4 * np.pi * RG**2 * s, RG, initial=0.0); return s / m[-1]   # unit total mass
def build(D, mcc, kind, a, Os=0.049):
    """density table (units: mass/Mcrit per unit volume*(3/4pi) as ptshell), normalised to M(<D)=1; mcc = interior mass / Mcrit(D)."""
    rc = 1 - RV / D; rin = 1 - T.KIN * RV / D; rout = 1 + T.KOUT * RV / D
    if kind == "remnant":        # compact point core carrying fraction f of the interior mass + thin Sedov interior
        s_in = a["f"] * shape_mass("point", rc, a) + (1 - a["f"]) * shape_mass("st", rc, a)
    elif kind.startswith("comp_"):   # compensated: point core mcc + extended interior (st / pl0 / pl1) carrying M_hole - mcc, so M(<D-R) = hole mass
        Mh = Os * rin**3; ext = kind[5:]
        se = shape_mass("st", rc, a) if ext == "st" else shape_mass("pl", rc, dict(n=float(ext[2:])))
        s_in = shape_mass("point", rc, a) + max(Mh - mcc, 0.0) / mcc * se
    else: s_in = shape_mass(kind, rc, a)
    s = mcc * s_in + np.where((RG > rin) & (RG < rout), Os * 3 / (4 * np.pi), 0.0)
    m = cumulative_trapezoid(4 * np.pi * RG**2 * s, RG, initial=0.0); m1 = np.interp(1.0, RG, m)
    return s / m1, m / m1, m1, rin, rout
def run(D, mcc, kind, a, Rs_list, taui=0.01, nsteps=200, iters=3):
    """a_f(theta) for each sampling radius in Rs_list (Gpc), plus X, V (for CMB). Same integrator, start and Ht iteration as prof_core.run."""
    S, MT, m1, rin, rout = build(D, mcc, kind, a)
    DR0 = D / min(Rs_list); Rb = max(P.rb_of(DR0, 1.5), 1.02 * rout)
    out = []; Ht = 0.95
    qs = [P.geometry(TH1, D / Rs) for Rs in Rs_list]
    for it in range(iters):
        mu = m1 * Ht**2 / 2
        X, V, _ = P.integrate(qs[0][0], mu, RG, MT, 0.0, Rb, taui, nsteps)
        Ht = np.hypot(*V[0]) / np.hypot(*X[0])
    for (q, rho) in qs:
        X, V, VR = P.integrate(q, mu, RG, MT, 0.0, Rb, taui, nsteps)
        vrel = V[1:] - V[0]; af = np.hypot(vrel[:, 0], vrel[:, 1]) / rho
        out.append(dict(af=af, X=X, V=V, vrmin=float(VR.min())))
    return out, dict(mu=mu, MT=MT, S=S, Ht=Ht, Rb=Rb, m1=m1, rin=rin, rout=rout)
def meff(D, mcc, kind, Os=0.049):
    """interior mass minus the shell's cavity mass (units Mcrit(D)); by the shell theorem this alone sets the static field on the R sphere,
    because the whole sphere (D-R <= r <= D+R) lies outside the interior (r < D-R) and inside the uniform shell."""
    rin = 1 - T.KIN * RV / D; Mh = Os * rin**3
    Min = max(mcc, Mh) if kind.startswith("comp_") else mcc
    return Min - Mh
def static_pattern(D, Me, H=HM):
    """(Phi_emit - Phi_obs)/c^2 on the R sphere for a point mass Me*Mcrit(D) at the CoE (uniform shell adds only l<=1): GM/(c^2 D) = Me (HD/c)^2 / 2"""
    th = np.radians(TH1); rs = RV / D; rr = np.hypot(1 + rs * np.cos(th), rs * np.sin(th))
    return -0.5 * Me * (H * D * 1e3 / CKM)**2 * (1 / rr - 1)
def legendre(f_theta, lmax=4):
    """Legendre coefficients about the CoE axis of an axisymmetric pattern sampled at TH1 (deg; theta from the away-from-CoE direction)."""
    th = np.radians(TH1[:181]); mu = np.cos(th); f = f_theta[:181]
    o = np.argsort(mu); mu, f = mu[o], f[o]
    from numpy.polynomial import legendre as L
    return [float((2 * l + 1) / 2 * np.trapezoid(f * L.legval(mu, [0] * l + [1]), mu)) for l in range(lmax + 1)]
def cmb(D, mcc, kind, a, Rs_cl=0.31):
    """cluster-scale pattern and CMB-scale (R = 14.3 Gpc sphere) diagnostics, today's field.
    kinematic: -n.(V_emit - V_obs)/c (linear, same convention as cycle 1); static: (Phi_emit - Phi_obs)/c^2; both Legendre a_l."""
    res, info = run(D, mcc, kind, a, [Rs_cl, RV])
    cl, far = res
    Hconv = HM / cl["af"].mean()               # code 1/time -> km/s/Mpc, set by the cluster-scale mean (<H> = 71.6)
    Dmpc = D * 1e3; vunit = Hconv * Dmpc        # km/s per code velocity
    X, V = far["X"], far["V"]; d = X[1:] - X[0]; n = d / np.hypot(d[:, 0], d[:, 1])[:, None]
    kin = -np.sum(n * (V[1:] - V[0]), axis=1) * vunit / CKM
    stat = static_pattern(D, meff(D, mcc, kind))
    stat_num_check = None
    dqc = T.dq(cl["af"]); dqf = T.dq(far["af"])
    return dict(kind=kind, a=a, D=D, mcc=mcc, cl_a2=dqc["a2"], cl_p2=dqc["p2"], far_a2=dqf["a2"], ratio=dqf["a2"] / max(dqc["a2"], 1e-30),
                kin=legendre(kin), stat=legendre(stat), vrmin=min(cl["vrmin"], far["vrmin"]), Ht=info["Ht"], Hconv=Hconv)
