"""CMB velocity/acceleration cancellation test (Chris's procedure, sections 1-13).

Target: G(ell) = curve/<curve>, curve = full-precision WLS geometry-only A/A_all (1-deg, 0..359).
Units c = 1, d1 = 1. theta = ell - phi.

Geometry / sign conventions
  Origin O (preferred centre, central mass) at (0,0). Observer initially at (d1,0) = (1,0).
  Source at distance q from the observer in direction theta:
      x_e0 = (1 + q cos theta, q sin theta),  r_e0 = sqrt(1 + q^2 + 2 q cos theta)
  so theta = 0 points away from O (anti-centre), theta = 180 deg toward O.
  Velocities are radial from O, outward positive:
      beta_o = v_o * x_hat,   beta_e = v_e * x_e / |x_e|.
  Run 1: v = k r with beta1 = k d1/c, no gravity, positions unscaled (s = 1).
  Runs 2-3: common expansion r = s r0 (r_o = s, r_e = s r_e0), current speeds from
      v^2 = beta1^2 r0^2 + 2 mu (1/r - 1/r0).
  Photon intercept: d = x_o - x_e (source -> observer at emission), observer moves at v_o along +x;
      (v_o^2 - 1) t^2 + 2 d_x v_o t + |d|^2 = 0, positive root; n = (d + v_o t x_hat)/t
  (n is the photon propagation direction, source -> observer, |n| = 1).
  D_motion = gamma_o (1 - beta_o.n) / [gamma_e (1 - beta_e.n)]
  D_grav   = sqrt[(1 - 2 mu/r_e)/(1 - 2 mu/r_o)]
  S_model  = D_motion * D_grav;  R = G S/<G S> - 1;  objective RMS(R).
"""
import json
import numpy as np
from scipy.optimize import differential_evolution, minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NOTES = "/workspace/cosmology-course/notes/"
ELL = np.arange(360.0)
full = np.loadtxt(NOTES + "migkas-best-fit-curve-0-360-full.csv", delimiter=",", skiprows=1)[:, 1]
rnd = np.loadtxt(NOTES + "migkas-best-fit-curve-0-360.csv", delimiter=",", skiprows=1)[:, 1]
G_FULL = full / full.mean()
G_RND = rnd / rnd.mean()
BAD = 1.0


def S_model(beta1, q, mu, s, phi, ell=ELL):
    """Return S(ell) or None if a physical guard fails."""
    th = np.radians(ell - phi)
    xe0 = 1.0 + q * np.cos(th)
    ye0 = q * np.sin(th)
    re0 = np.sqrt(xe0 ** 2 + ye0 ** 2)
    if np.any(re0 < 1e-6):
        return None
    ro, re = s * 1.0, s * re0
    vo2 = beta1 ** 2 + 2 * mu * (1.0 / ro - 1.0)
    ve2 = beta1 ** 2 * re0 ** 2 + 2 * mu * (1.0 / re - 1.0 / re0)
    if vo2 <= 0 or vo2 >= 1 or np.any(ve2 <= 0) or np.any(ve2 >= 1):
        return None
    if 2 * mu / ro >= 1 or np.any(2 * mu / re >= 1):
        return None
    vo, ve = np.sqrt(vo2), np.sqrt(ve2)
    # positions at emission (common s)
    xe, ye = s * xe0, s * ye0
    dx, dy = ro - xe, -ye
    a = vo2 - 1.0
    b = 2 * dx * vo
    c = dx ** 2 + dy ** 2
    disc = b * b - 4 * a * c
    tf = (-b - np.sqrt(disc)) / (2 * a)          # positive root (a < 0, c > 0)
    nx, ny = (dx + vo * tf) / tf, dy / tf
    bex, bey = ve * xe / re, ve * ye / re
    go, ge = 1 / np.sqrt(1 - vo2), 1 / np.sqrt(1 - ve2)
    Dm = go * (1 - vo * nx) / (ge * (1 - (bex * nx + bey * ny)))
    Dg = np.sqrt((1 - 2 * mu / re) / (1 - 2 * mu / ro))
    return Dm * Dg, Dm, Dg


def resid(G, S):
    p = G * S
    return p / p.mean() - 1.0


def rms(x):
    return float(np.sqrt(np.mean(x ** 2)))


def obj(p, G, kind, bounds=None):
    if bounds is not None and any((v < lo) or (v > hi) for v, (lo, hi) in zip(p, bounds)):
        return BAD
    if kind == 1:
        beta1, q, phi = p
        out = S_model(beta1, q, 0.0, 1.0, phi)
    else:
        beta1, q, mu, s, phi = p
        out = S_model(beta1, q, mu, s, phi)
    if out is None:
        return BAD
    return rms(resid(G, out[0]))


def fit(G, kind, bounds, seeds=(1, 2, 3, 4), x0=None, de_x0=None):
    best = None
    if x0 is not None:  # local polish only, from a given start (used to check the note's run-1 minimum)
        r2 = minimize(obj, x0, args=(G, kind, bounds), method="Nelder-Mead",
                      options=dict(xatol=1e-12, fatol=1e-16, maxiter=40000, maxfev=40000, adaptive=True))
        print(f"  local from {x0}: rms={r2.fun:.6e} x={np.round(r2.x, 6).tolist()}", flush=True)
        return r2.x, r2.fun
    for sd in seeds:
        r = differential_evolution(obj, bounds, args=(G, kind), seed=sd, popsize=40, maxiter=3000,
                                   tol=1e-12, atol=0, mutation=(0.5, 1.0), recombination=0.9,
                                   polish=False, init="sobol", x0=de_x0)
        r2 = minimize(obj, r.x, args=(G, kind, bounds), method="Nelder-Mead",
                      options=dict(xatol=1e-12, fatol=1e-16, maxiter=40000, maxfev=40000, adaptive=True))
        x, f = (r2.x, r2.fun) if r2.fun < r.fun else (r.x, r.fun)
        print(f"  kind {kind} seed {sd}: rms={f:.6e} x={np.round(x, 6).tolist()}", flush=True)
        if best is None or f < best[1]:
            best = (x, f)
    return best


def fourier(R, M=12):
    a0 = float(R.mean())
    rad = np.radians(ELL)
    rows, cum = [], []
    var = float(np.mean((R - a0) ** 2))
    left = R - a0
    for m in range(1, M + 1):
        am = float(2 * np.mean(R * np.cos(m * rad)))
        bm = float(2 * np.mean(R * np.sin(m * rad)))
        Am = float(np.hypot(am, bm))
        ph = float((np.degrees(np.arctan2(bm, am)) / m) % (360.0 / m))
        rows.append(dict(m=m, a=am, b=bm, A=Am, phi_deg=ph, period_deg=360.0 / m,
                         var_frac=float(0.5 * Am ** 2 / var)))
        left = left - am * np.cos(m * rad) - bm * np.sin(m * rad)
        cum.append(dict(m_max=m, leftover_rms=rms(left)))
    return dict(a0=a0, rms=rms(R), modes=rows, cumulative=cum)


def pack(name, G, kind, x, f):
    if kind == 1:
        beta1, q, phi = x
        mu, s = 0.0, 1.0
    else:
        beta1, q, mu, s, phi = x
    S, Dm, Dg = S_model(beta1, q, mu, s, phi)
    R = resid(G, S)
    th = np.radians(ELL - phi)
    re0 = np.sqrt(1 + q * q + 2 * q * np.cos(th))
    vo = np.sqrt(beta1 ** 2 + 2 * mu * (1 / s - 1))
    ve = np.sqrt(beta1 ** 2 * re0 ** 2 + 2 * mu * (1 / (s * re0) - 1 / re0))
    return dict(name=name, beta1=float(beta1), q=float(q), mu=float(mu), s=float(s),
                phi_deg=float(phi % 360), rms=rms(R), max_abs_R=float(np.abs(R).max()),
                v_o_now=float(vo), v_e_now_min=float(ve.min()), v_e_now_max=float(ve.max()),
                two_mu_over_r_max=float(max(2 * mu / s, (2 * mu / (s * re0)).max())),
                fourier=fourier(R)), S, R


def main():
    out = {}
    curves = {}
    Sreq = (1 / G_FULL) / np.mean(1 / G_FULL)
    out["target"] = dict(source="migkas-best-fit-curve-0-360-full.csv (12 sig. digits, regenerated from wls_fit)",
                         G_min=float(G_FULL.min()), G_max=float(G_FULL.max()),
                         max_abs_round_check=float(np.abs(np.round(full, 4) - rnd).max()),
                         rms_G_round_minus_full=rms(G_RND - G_FULL))

    print("Run 1", flush=True)
    b1 = [(1e-4, 0.999), (1e-4, 5.0), (0.0, 360.0)]
    x, f = fit(G_FULL, 1, b1)
    r, S, R = pack("run1_v=kr", G_FULL, 1, x, f)
    out["run1"] = r; curves["run1"] = (S, R)
    o = S_model(0.667, 0.143, 0.0, 1.0, 15.16)
    Rn = resid(G_FULL, o[0])
    out["run1_note_values_evaluated_phi15.16"] = dict(rms=rms(Rn), max_abs_R=float(np.abs(Rn).max()))
    o = S_model(0.667, 0.143, 0.0, 1.0, 195.16)
    Rn = resid(G_FULL, o[0])
    out["run1_note_values_evaluated_phi195.16"] = dict(rms=rms(Rn), max_abs_R=float(np.abs(Rn).max()))
    xl, fl = fit(G_FULL, 1, b1, x0=np.array([0.667, 0.143, 15.2]))
    out["run1_local_from_note_start"] = pack("run1_local_from_note_start", G_FULL, 1, xl, fl)[0]

    b2 = [(1e-4, 0.999), (1e-4, 3.0), (0.0, 2.0), (1.0, 3.0), (0.0, 360.0)]
    print("Run 2", flush=True)
    x2, f = fit(G_FULL, 2, b2)
    out["run2_note_values"] = dict(beta1=0.5865, q=0.2304, mu=0.09327, s=3.0, phi_deg=15.16, rms=2.37e-3, max=3.80e-3)
    r, S, R = pack("run2_central_mass_s<=3", G_FULL, 2, x2, f)
    out["run2"] = r; curves["run2"] = (S, R)

    # mu < beta1^2 s / (2(s-1)) < 0.5 is needed for v_o^2 > 0 once s > 1, so mu <= 0.5 loses nothing.
    # The feasible set is thin in this box, so the DE population is seeded with the run-2 optimum.
    b3 = [(1e-4, 0.999), (1e-4, 3.0), (0.0, 0.5), (1.0, 50.0), (0.0, 360.0)]
    print("Run 3", flush=True)
    x3, f = fit(G_FULL, 2, b3, de_x0=x2)
    r, S, R = pack("run3_central_mass_s<=50", G_FULL, 2, x3, f)
    out["run3"] = r; curves["run3"] = (S, R)

    print("Run 2 on rounded CSV", flush=True)
    xr, f = fit(G_RND, 2, b2, de_x0=x2)
    r, S, R = pack("run2_rounded_csv", G_RND, 2, xr, f)
    out["run2_rounded"] = r

    # note's representative run-2 parameters evaluated (no fit) on our G
    nb = (0.5865, 0.2304, 0.09327, 3.0, 15.16)
    o = S_model(*nb)
    out["note_run2_params_evaluated"] = dict(params=dict(beta1=nb[0], q=nb[1], mu=nb[2], s=nb[3], phi_deg=nb[4]),
                                             rms_full=None if o is None else rms(resid(G_FULL, o[0])),
                                             max_full=None if o is None else float(np.abs(resid(G_FULL, o[0])).max()))
    o = S_model(0.667, 0.143, 0.0, 1.0, out["run1"]["phi_deg"])
    out["note_run1_params_evaluated_at_our_phi"] = dict(rms_full=None if o is None else rms(resid(G_FULL, o[0])))

    out["notes"] = {'conventions': 'O at (0,0), observer at (1,0), source x_e0=(1+q cos th, q sin th); th=0 is anti-centre. Velocities radial outward from O. n = photon direction source->observer from the intercept root. D_motion=g_o(1-b_o.n)/[g_e(1-b_e.n)].', 'run1': "Global DE optimum is the degenerate limit beta1->1, q->0 (all speeds ~0.995c, source 0.0015 d1 from observer). The note's point (0.667,0.143) is not a stationary point under these conventions: evaluated at phi=195.16 it gives RMS 4.20e-2, max 6.42e-2; a local polish from it slides to the beta1=0.999 bound.", 'run3': 's runs to the 50 bound. With the other parameters fixed, s enters only as 1/s (intercept geometry is scale-free), so the s->inf limit is finite; during debugging an out-of-box Nelder-Mead polish reached s~1e12 with RMS 2.1425e-3 (beta1 0.5839, q 0.2037, mu 0.0724). No physical guard failed (2mu/r <= 0.0037 at s=50).', 'fourier_phases': 'All phi_m equal 15.16 deg + k*180/m: every residual mode is a pure cosine about the fitted axis (residual mirror-symmetric about ell=15.16), signs alternating.'}
    json.dump(out, open(NOTES + "cmb-cancel.json", "w"), indent=1)
    np.savez("/tmp/cc/curves.npz", Sreq=Sreq, **{k + "_S": v[0] for k, v in curves.items()},
             **{k + "_R": v[1] for k, v in curves.items()})
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "fourier"} if isinstance(v, dict) else v
                      for k, v in out.items()}, indent=1))


# ---------------------------------------------------------------- plots (python3 cmb-cancel.py --plots)
def sci(x, d=3):
    m, e = f"{x:.{d}e}".split("e")
    return rf"{m}\times10^{{{int(e)}}}"


def params_of(r):
    return r["beta1"], r["q"], r["mu"], r["s"], r["phi_deg"]


def plots():
    plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "axes.labelsize": 16,
                         "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14,
                         "mathtext.fontset": "cm"})
    out = json.load(open(NOTES + "cmb-cancel.json"))
    best_key = min(("run1", "run2", "run3"), key=lambda k: out[k]["rms"])
    show = [best_key] + (["run2"] if best_key != "run2" else [])
    Sreq = (1 / G_FULL) / np.mean(1 / G_FULL)
    cols = {"run1": "C2", "run2": "C0", "run3": "C3"}
    lab = {"run1": "run 1 ($v=kr$)", "run2": r"run 2 ($s\leq3$)", "run3": r"run 3 ($s\leq50$)"}

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 10.5), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.25, 1]))
    a1.plot(ELL, G_FULL, color="k", lw=2.2, label=r"$G(\ell)$")
    a1.plot(ELL, Sreq, color="0.55", lw=5, alpha=0.6, label=r"$S_{\rm req}=(1/G)/\langle1/G\rangle$")
    for k in show:
        S = S_model(*params_of(out[k]))[0]
        a1.plot(ELL, S / S.mean(), color=cols[k], lw=1.8, ls="--",
                label=r"$S_{\rm model}/\langle S\rangle$, " + lab[k])
        R = resid(G_FULL, S)
        a2.plot(ELL, R, color=cols[k], lw=2,
                label=lab[k] + rf", RMS $={sci(out[k]['rms'], 2)}$")
    a1.set_ylabel("normalised factor")
    a1.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=False)
    a2.axhline(0, color="0.6", lw=0.8)
    a2.set_ylabel(r"$R(\ell)=\Delta T/T$")
    a2.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    a2.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=1, frameon=False)
    a2.set_xlim(0, 360)
    a2.set_xticks(np.arange(0, 361, 45))
    a2.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    fig.tight_layout()
    fig.savefig(NOTES + "cmb-cancel-curves.png", dpi=110, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 6.5))
    m = np.arange(1, 13)
    w = 0.27
    for i, k in enumerate(("run1", "run2", "run3")):
        A = [r["A"] for r in out[k]["fourier"]["modes"]]
        ax.bar(m + (i - 1) * w, A, width=w, color=cols[k], label=lab[k])
    ax.set_yscale("log")
    ax.set_xticks(m)
    ax.set_xlabel(r"longitude Fourier mode $m$")
    ax.set_ylabel(r"$A_m$ of residual $R(\ell)$")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(NOTES + "cmb-cancel-fourier.png", dpi=110, bbox_inches="tight")
    plt.close(fig)

    b = out[best_key]; r2 = out["run2"]; r1 = out["run1"]
    lines = [
        (r"$\mathbf{Equations\ used}$  ($c=1,\ d_1=1,\ \theta=\ell-\phi$)", 0),
        (r"$G(\ell)=\dfrac{A/A_{\rm all}(\ell)}{\langle A/A_{\rm all}\rangle}$", 1),
        (r"$S_{\rm req}=\dfrac{1/G}{\langle 1/G\rangle}$", 1),
        (r"$R(\ell)=\dfrac{G\,S_{\rm model}}{\langle G\,S_{\rm model}\rangle}-1$", 1),
        (r"$r_{e0}=\sqrt{1+q^2+2q\cos\theta}$", 1),
        (r"$\mathbf{x}_{e0}=(1+q\cos\theta,\ q\sin\theta)$", 1),
        (r"$r_o=s,\quad r_e=s\,r_{e0}$", 1),
        (r"$v^2=\beta_1^2 r_0^2+2\mu\left(\dfrac{1}{r}-\dfrac{1}{r_0}\right)$", 1),
        (r"$\boldsymbol{\beta}_o=v_o\,\hat{\mathbf{x}},\quad \boldsymbol{\beta}_e=v_e\,\mathbf{x}_e/r_e$", 1),
        (r"$\mathbf{d}=\mathbf{x}_o-\mathbf{x}_e$", 1),
        (r"$(v_o^2-1)\,t_f^2+2d_x v_o\,t_f+|\mathbf{d}|^2=0$", 1),
        (r"$t_f>0\ \mathrm{root}$", 1),
        (r"$\hat{\mathbf{n}}=\dfrac{\mathbf{d}+v_o t_f\,\hat{\mathbf{x}}}{t_f}$", 1),
        (r"$D_{\rm motion}=\dfrac{\gamma_o(1-\boldsymbol{\beta}_o\cdot\hat{\mathbf{n}})}{\gamma_e(1-\boldsymbol{\beta}_e\cdot\hat{\mathbf{n}})}$", 1),
        (r"$D_{\rm grav}=\sqrt{\dfrac{1-2\mu/r_e}{1-2\mu/r_o}}$", 1),
        (r"$S_{\rm model}=D_{\rm motion}\,D_{\rm grav}$", 1),
        (r"$\mathbf{Best\ fits}$", 0),
        (rf"Run 1: $\beta_1={r1['beta1']:.4f},\ q={r1['q']:.5f}$", 1),
        (rf"$\phi={r1['phi_deg']:.2f}^\circ,\ \mathrm{{RMS}}={sci(r1['rms'])}$", 2),
        (r"(degenerate limit: $v\to c$, $q\to 0$)", 2),
        (rf"Run 2: $\beta_{{1,0}}={r2['beta1']:.4f},\ q={r2['q']:.4f}$", 1),
        (rf"$\mu={r2['mu']:.5f},\ s={r2['s']:.3f}\ ({{\rm at\ bound}}),\ \phi={r2['phi_deg']:.2f}^\circ$", 2),
        (rf"$\mathrm{{RMS}}={sci(r2['rms'])},\ \max|R|={sci(r2['max_abs_R'])}$", 2),
    ]
    if best_key == "run3":
        lines += [
            (rf"Run 3: $\beta_{{1,0}}={b['beta1']:.4f},\ q={b['q']:.4f}$", 1),
            (rf"$\mu={b['mu']:.4f},\ s={b['s']:.1f}\ ({{\rm at\ bound}}),\ \phi={b['phi_deg']:.2f}^\circ$", 2),
            (rf"$\mathrm{{RMS}}={sci(b['rms'])},\ \max|R|={sci(b['max_abs_R'])}$", 2),
        ]
    else:
        r3 = out["run3"]
        lines += [
            (rf"Run 3: $\beta_{{1,0}}={r3['beta1']:.4f},\ q={r3['q']:.4f}$", 1),
            (rf"$\mu={r3['mu']:.4f},\ s={r3['s']:.1f}\ ({{\rm at\ bound}}),\ \phi={r3['phi_deg']:.2f}^\circ$", 2),
            (rf"$\mathrm{{RMS}}={sci(r3['rms'])},\ \max|R|={sci(r3['max_abs_R'])}$", 2),
        ]
    tall = [1.75 if "dfrac" in t else (1.2 if "sqrt" in t else 0.9) for t, _ in lines]
    H = sum(tall) * 0.62 + 0.4
    fig = plt.figure(figsize=(9, H), dpi=100)
    y = 1 - 0.2 / H
    for (t, ind), h in zip(lines, tall):
        y -= h * 0.62 / H / 2
        fig.text(0.04 + 0.05 * ind, y, t, fontsize=20, va="center", ha="left", math_fontfamily="cm",
                 family="serif")
        y -= h * 0.62 / H / 2
    fig.savefig(NOTES + "cmb-cancel-eqs.png", dpi=100)
    plt.close(fig)
    print("plots written; best =", best_key)


if __name__ == "__main__":
    import sys
    if "--plots" not in sys.argv:
        main()
    plots()
