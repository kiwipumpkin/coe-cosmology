from tm_core import *
from scipy.optimize import least_squares, brentq
import relfree as RF, rm_core as C
OUT = {}
# ---------- series (O(B^2) even, O(eps) odd) ----------
def series(beta, gam, eps):
    B = 2 * beta + beta**2; Pp = beta * (3 + beta); Q = (1 + beta) * (gam + beta) - 2 * beta
    m0 = 1 + B / 4 - 3 * B**2 / 64; m2 = B / 4 - B**2 / 16
    m1 = eps / 2 * (Pp + 0.75 * (Q - Pp * B / 2) - Q * B / 2 * 10 / 16); m3 = eps / 2 * (0.25 * (Q - Pp * B / 2) - Q * B / 2 * 5 / 16)
    return dict(a1=m1 / m0, a2=m2 / m0, a3=m3 / m0, B=B, P=Pp, Q=Q)
def lin_h(a, Ht, taui=0.01):
    """linear-perturbation expansion-rate profile in the shell (labels q, units D): h ~ 1 - (Ht^2/2)[mcc(1/taui-1)/q^3 + Os(1/2 - 3 rin^2/(2q^2) + rin^3/q^3)]"""
    k = Ht**2 / 2; Kc = k * a["mcc"] * (1 / taui - 1); ks = k * a["Os"] * a["rin"]**3; b = k * 1.5 * a["Os"] * a["rin"]**2; a0 = k * a["Os"] / 2
    h = lambda q: 1 - a0 - (Kc + ks) / q**3 + b / q**2
    h1 = h(1.0); hp = 3 * (Kc + ks) - 2 * b; hpp = -12 * (Kc + ks) + 6 * b
    return dict(beta=hp / h1, gamma=hpp / h1, Kc=Kc, ks=ks, b=b, a0=a0, h=h)
# ---------- data rotated to l0 ----------
def data_rot(l0):
    return {m: rot(*AFc[m], m, l0) for m in range(4)}
OUT["AF_amp"] = {m: list(amp_ph(*AFc[m], m)) for m in (1, 2, 3)}
# ---------- fits ----------
fits = {}
for key in ("free", "omb"):
    a, d = newt_params(key); rr, m1 = newt_run(a); eps = 1 / a["DR"]; af = rr["af"] / rr["af"].mean()
    qs = np.linspace(1 - 1.2 * eps, 1 + 1.2 * eps, 41); h, _ = radial_h(a, rr, qs); cf = np.polyfit(qs - 1, h, 6); h1 = np.polyval(cf, 0)
    beta = np.polyval(np.polyder(cf), 0) / h1; gam = np.polyval(np.polyder(cf, 2), 0) / h1
    L = lin_h(a, rr["Ht"]); Alin = A_radial(L["h"], eps)
    fits[key] = dict(kind="newt", eps=eps, D=a["DR"] * RV, beta=beta, gamma=gam, beta_lin=L["beta"], gamma_lin=L["gamma"], Kc=L["Kc"], ks=L["ks"], b=L["b"],
                     mcc=a["mcc"], Os=a["Os"], rin=a["rin"], Meff=a["mcc"] - a["Os"] * a["rin"]**3, geo=af, Alin=Alin, l0=d["l0"], amp=d["amp"], phi=d["phi"], R2=d.get("R2_mig", d.get("R2")))
o = json.load(open("universe_rel_mcfree.json"))
e = RF.ev(o["Mc_kg"], o["D"], o["kin"], o["kout"], full=True); Ad = np.array(e["Ad"]); geo = C.curve(Ad, 0.0)
fits["rel"] = dict(kind="rel", eps=1.19 / o["D"], D=o["D"], beta=o["H_r"] / o["H_perp"] - 1, gamma=np.nan, mcc=o["mcc"], Os=0.049, rin=o["R_in"] / o["D"],
                   Meff=o["mcc"] - 0.049 * (o["R_in"] / o["D"])**3, geo=geo, l0=o["l0_s"], amp=o["amp"], phi=o["phi"], R2=o["R2_s"])
for k, f in fits.items():
    g = f["geo"]; f["num"] = {m: list(harm(g, m)) for m in range(5)}
    if f["kind"] == "newt":
        f["ser"] = series(f["beta"], f["gamma"], f["eps"]); f["ser_lin"] = series(f["beta_lin"], f["gamma_lin"], f["eps"]); f["lin_num"] = {m: list(harm(f["Alin"], m)) for m in range(4)}
    else:
        b = f["beta"]; f["ser"] = dict(a1=0.0, a2=(b / 2) / (1 + b / 2), a3=0.0)
    full = A_from(shifted(g, f["l0"]), TH1 - f["l0"], f["amp"], f["phi"])
    f["full_num"] = {m: list(harm(full, m)) for m in range(4)}; f["R2F"] = r2(AFfour(TH1), full); f["R2mig_check"] = r2(MA, full[IMIG])
    print(k, "eps %.4f D %.2f beta %.4f gamma %.3f" % (f["eps"], f["D"], f["beta"], f["gamma"]), "lin b/g", f.get("beta_lin"), f.get("gamma_lin"), "Meff", f["Meff"], "R2chk", f["R2mig_check"], f["R2"])
    print("   num a_m", [round(f["num"][m][0], 5) for m in range(5)], " series", {kk: round(v, 5) for kk, v in f["ser"].items() if kk[0] == "a"})
    if f["kind"] == "newt": print("   lin-h num", [round(f["lin_num"][m][0], 5) for m in range(4)], " lin series", {kk: round(v, 5) for kk, v in f["ser_lin"].items() if kk[0] == "a"})
# ---------- term-by-term solve against A_F ----------
def model_harm(p, shape, eps_fixed=None):
    """p = (l0, beta, eps, amp, phi); geometric part = exact radial field with h(q) = 1 + beta (q-1) + gamma/2 (q-1)^2 ... use shape family"""
    l0, beta, eps, amp, phi = p
    if shape == "rel": g = 1 + beta * np.cos(TH)**2; g = g / g.mean()
    else:
        gr = shape * beta
        g = A_radial(lambda q: 1 + beta * (q - 1) + 0.5 * gr * (q - 1)**2, eps)
    A = A_from(shifted(g, l0), TH1 - l0, amp, phi)
    return A
def resid(p, shape):
    A = model_harm(p, shape); out = []
    for m in range(4):
        a, b = harm(A, m); c, s = AFc[m]; out += [a - c, b - s] if m else [a - c]
    return np.array(out)
sol = {}
for name, shape, x0 in (("pointmass", -4.0, [15, 0.16, 0.3, -0.05, 0.9]), ("shape_free", fits["free"]["gamma"] / fits["free"]["beta"], [13, 0.18, 0.15, -0.12, 1.2]),
                        ("shape_omb", fits["omb"]["gamma"] / fits["omb"]["beta"], [13, 0.16, 0.29, -0.06, 0.9]), ("rel", "rel", [13.5, 0.18, 0.0, -0.05, 0.85])):
    lb = [-30, 0.0, 0.0 if shape == "rel" else 0.02, -0.3, -10]; ub = [60, 1.0, 1e-9 if shape == "rel" else 0.95, 0.3, 10]
    if shape == "rel": x0[2] = 5e-10
    r = least_squares(resid, x0, args=(shape,), bounds=(lb, ub))
    A = model_harm(r.x, shape); sol[name] = dict(p=list(r.x), D=RV / r.x[2] if shape != "rel" else None, harm={m: list(harm(A, m)) for m in range(4)}, R2F=r2(AFfour(TH1), A), R2mig=r2(MA, A[IMIG]),
                                                 res=list(resid(r.x, shape)), shape=shape)
    print("SOLVE", name, np.round(r.x, 4), "D", sol[name]["D"], "R2F %.3f R2mig %.3f" % (sol[name]["R2F"], sol[name]["R2mig"]), "res", np.round(r.fun, 4))
# sequential (hand) solution with the series, point-mass shape gamma = -4 beta, and with each fit's gamma/beta
from scipy.optimize import fsolve
def hand(l0, gb):
    Dt = data_rot(l0); C0 = Dt[0][0]; C1, S1 = Dt[1]; C2, S2 = Dt[2]; C3, S3 = Dt[3]
    beta = brentq(lambda b: series(b, gb * b, 0)["a2"] - C2 / C0, 1e-4, 1.0)
    s1 = series(beta, gb * beta, 1.0); a2 = s1["a2"]; g1, g3 = s1["a1"], s1["a3"]      # a1 = eps*g1, a3 = eps*g3 (divided by C0 below)
    def F(x):
        amp, phi, eps = x
        return [C0 * (eps * g1 - amp * np.sin(phi) * (1 + a2 / 2)) - C1, -C0 * amp * np.cos(phi) * (1 - a2 / 2) - S1, C0 * (eps * g3 - amp * a2 / 2 * np.sin(phi)) - C3]
    amp, phi, eps = fsolve(F, [-0.06, 0.9, 0.3])
    if amp > 0: amp, phi = -amp, phi + np.pi
    S3p = -C0 * amp * a2 / 2 * np.cos(phi)
    N = C0 / (1 - amp * eps * g1 * np.sin(phi) / 2)
    return dict(l0=l0, gb=gb, C={m: list(Dt[m]) for m in range(4)}, beta=beta, a2=a2, g1=g1, g3=g3, amp=amp, phi=phi % (2 * np.pi), eps=eps, D=RV / eps, S3_pred=S3p, S3_data=S3, N=N, epsQ=eps * s1["Q"])
OUT["hand"] = {nm: hand(14.96, gb) for nm, gb in (("pointmass", -4.0), ("shape_free", fits["free"]["gamma"] / fits["free"]["beta"]), ("shape_omb", fits["omb"]["gamma"] / fits["omb"]["beta"]))}
for k, v in OUT["hand"].items(): print("HAND", k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "C"})
print("Crot", OUT["hand"]["pointmass"]["C"])
# relativistic: numeric LTB geometric shape scaled, g = 1 + s (g_rel - 1)
grel = fits["rel"]["geo"]
def resid_rel(p):
    l0, s, amp, phi = p; g = 1 + s * (grel - 1); A = A_from(shifted(g, l0), TH1 - l0, amp, phi); out = []
    for m in range(4):
        a, b = harm(A, m); c, s_ = AFc[m]; out += [a - c, b - s_] if m else [a - c]
    return np.array(out), A
r = least_squares(lambda p: resid_rel(p)[0], [13.5, 1.0, -0.05, 0.85]); A = resid_rel(r.x)[1]
sol["rel_numshape"] = dict(p=list(r.x), res=list(r.fun), R2F=r2(AFfour(TH1), A), R2mig=r2(MA, A[IMIG]), harm={m: list(harm(A, m)) for m in range(4)})
print("SOLVE rel_numshape", np.round(r.x, 4), "R2F %.3f R2mig %.3f" % (sol["rel_numshape"]["R2F"], sol["rel_numshape"]["R2mig"]), np.round(r.fun, 4))
OUT["fits"] = {k: {kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv) for kk, vv in f.items() if kk not in ("geo", "Alin")} for k, f in fits.items()}
for k in OUT["fits"]:
    for kk in ("ser", "ser_lin"):
        if kk in OUT["fits"][k]: OUT["fits"][k][kk] = {a: float(b) for a, b in OUT["fits"][k][kk].items()}
OUT["solve"] = sol
json.dump(OUT, open("/workspace/cosmology/papers/A_equals_fourier/thinshell/termmatch.json", "w"), indent=1, default=float)
print("hand", OUT["hand"])
