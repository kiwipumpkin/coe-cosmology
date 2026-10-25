"""Independent audit of the Newtonian fixed-profile ptshell model (universe_thin 'free' and 'baryons').
Own RK4 integrator (numpy), checked against the original P.run. Variants:
 fixed  : mass profile fixed in space (original)
 como   : same profile but COMOVING with the flow (radii scale with tau, enclosed mass per label fixed) = Newtonian analogue of LTB dust
 nocore : fixed, core removed ; noshell : fixed, shell removed
Observables per direction theta (theta=0 = away from CoE):
 af      : |V_P - V_O| / rho_label (what the Newtonian code uses, A ~ af, linear)
 Hrad    : (V_P-V_O).n / |X_P-X_O|  at the FINAL positions (true radial Hubble rate at the galaxy's actual sky angle)
Also Earth's history: fraction of time spent inside the empty interior."""
import sys, json, numpy as np
sys.path.insert(0, ".."); sys.path.insert(0, "../../integrated")
import ptshell as S, prof_core as P
from bubble_core import MA, IMIG, TH1, r2, A_from, best_sine, shifted
from scipy.integrate import cumulative_trapezoid
RG = P.RG
def mtab(d, core=True, shell=True):
    Lc = d["Lc"]; c = np.exp(-RG / Lc); mc_ = cumulative_trapezoid(4*np.pi*RG**2*c, RG, initial=0); mc_ = mc_/np.interp(1, RG, mc_)*d["mcc"]
    ms = d["Os"]*(np.clip(RG, d["rin"], d["rout"])**3 - d["rin"]**3)
    return (mc_ if core else 0*RG) + (ms if shell else 0*RG)          # M(<r)/Mcrit (today's units)
def run(d, mode="fixed", core=True, shell=True, taui=0.01, n=400, iters=4, track=False):
    DR = d["DR"]; rho = 1/DR; th = np.radians(TH1)
    q = np.vstack([[1, 0], np.c_[1 + rho*np.cos(th), rho*np.sin(th)]])
    m = mtab(d, core, shell); Ht = 0.95
    for _ in range(iters):
        k = Ht**2/2
        def acc(x, tau):
            r = np.hypot(x[:, 0], x[:, 1])
            if mode == "fixed": mm = np.interp(r, RG, m)
            else:              # comoving: mass enclosed within r at time tau = profile evaluated at r/tau (labels)
                mm = np.interp(r/tau, RG, m)
            return -k*mm[:, None]*x/r[:, None]**3
        s = np.log(taui); h = -s/n; x = q*taui; v = q.copy(); hist = []
        for i in range(n):
            t1, t2, t4 = np.exp(s), np.exp(s+h/2), np.exp(s+h)
            a1 = acc(x, t1); k1x, k1v = t1*v, t1*a1
            a2 = acc(x+h/2*k1x, t2); k2x, k2v = t2*(v+h/2*k1v), t2*a2
            a3 = acc(x+h/2*k2x, t2); k3x, k3v = t2*(v+h/2*k2v), t2*a3
            a4 = acc(x+h*k3x, t4); k4x, k4v = t4*(v+h*k3v), t4*a4
            x = x + h/6*(k1x+2*k2x+2*k3x+k4x); v = v + h/6*(k1v+2*k2v+2*k3v+k4v); s += h
            if track: hist.append((np.exp(s), np.hypot(*x[0])))
        Ht = np.hypot(*v[0])/np.hypot(*x[0])
    dv = v[1:]-v[0]; dx = x[1:]-x[0]; dist = np.hypot(dx[:, 0], dx[:, 1])
    af = np.hypot(dv[:, 0], dv[:, 1])/rho
    Hrad = np.sum(dv*dx, 1)/dist**2
    ang = np.degrees(np.arctan2(dx[:, 1], dx[:, 0])) % 360           # actual sky angle
    Hrad_grid = np.interp(TH1, np.r_[ang[np.argsort(ang)]-360, ang[np.argsort(ang)], ang[np.argsort(ang)]+360],
                          np.r_[Hrad[np.argsort(ang)]]*1 if False else np.tile(Hrad[np.argsort(ang)], 3))
    out = dict(af=af, Hrad=Hrad_grid, dist=dist/rho, Ht=Ht, rE=np.hypot(*x[0]), vE=np.hypot(*v[0]))
    if track: out["hist"] = np.array(hist)
    return out
def score(y, d, sq=False):
    """y on theta grid (theta = l - l0). Returns R2 no sine (mean-norm), R2 with fitted sine, corr, along/across ratio."""
    y = y**2 if sq else y; y = y/y.mean(); l0 = d["l0"]
    ys = shifted(y, l0); th = TH1 - l0; x = ys[IMIG]
    a, ph = best_sine(x, th[IMIG], MA)
    return dict(R2_nosine=r2(MA, x), R2_sine=r2(MA, A_from(x, th[IMIG], a, ph)), amp=a, corr=np.corrcoef(x, MA)[0, 1],
                along_over_across=float((y[0]+y[180])/(y[90]+y[270])))
if __name__ == "__main__":
    res = {}
    for key, fn in (("free", "../thin_b_free_hug.json"), ("baryons", "../thin_b_0.05_hug.json")):
        d = json.load(open(fn)); R = {}
        orig = P.run("ptshell", dict(d, mc=0, Ob=0), 0, 0)["af"]
        R["orig_code_R2_sine"] = score(orig, d)["R2_sine"]
        base = run(d, track=True)
        R["replication_max_rel_diff_af"] = float(np.max(np.abs(base["af"]/base["af"].mean() - orig/orig.mean())))
        h = base["hist"]; R["earth_frac_time_in_empty_interior"] = float(np.mean(h[:, 1] < d["rin"]))
        R["earth_tau_enter_shell"] = float(h[np.argmax(h[:, 1] >= d["rin"]), 0])
        R["label_vs_final_dist_range"] = [float(base["dist"].min()), float(base["dist"].max())]
        for name, kw in (("fixed", {}), ("fixed_nocore", dict(core=False)), ("fixed_noshell", dict(shell=False)), ("comoving", dict(mode="como"))):
            o = base if name == "fixed" else run(d, **kw)
            R[name] = {"af(linear,code)": score(o["af"], d), "af^2": score(o["af"], d, True), "Hrad^2(true radial, final pos)": score(o["Hrad"], d, True),
                       "Ht": o["Ht"]}
        res[key] = R
    json.dump(res, open("newton_audit.json", "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=lambda v: round(float(v), 4)))
