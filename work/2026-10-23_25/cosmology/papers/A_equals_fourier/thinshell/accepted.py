from tm_core import *
G = 6.674e-11; c = 2.99792458e5; GPCm = 3.0857e25; MPC = 3.0857e19
def rhoc(H): h = H / MPC; return 3 * h * h / (8 * np.pi * G)
RV_m = RV * GPCm; VOL = 4 * np.pi / 3 * RV_m**3
P = dict(H0=67.36, H0e=0.54, Obh2=0.02237, Om=0.3153, Ome=0.0073, age=13.797, agee=0.023)
Ob = P["Obh2"] / (P["H0"] / 100)**2; rc_p = rhoc(P["H0"]); rc70 = rhoc(70.0)
acc = dict(Ob=Ob, rho_b=Ob * rc_p, rhoc_planck=rc_p, rhoc70=rc70, Mb=Ob * rc_p * VOL, Mm=P["Om"] * rc_p * VOL, Om=P["Om"], H0p=P["H0"], H0s=73.04, H0se=1.04,
           vcmb=369.82, vlg=620.0, bulk=395.0, age=P["age"], invH70=1 / (70 / MPC) / 3.156e16)
T = json.load(open("universe_thin.json")); O = json.load(open("universe_newt_omb.json"))["Rs14.3|sine"]
out = dict(accepted=acc, P=P, VOL=VOL)
for key, tk in (("free", "free"), ("omb", "newt_omb")):
    t = T[tk]; a, d = newt_params(key); rr, m1 = newt_run(a); eps = 1 / a["DR"]
    qs = np.linspace(1 - 0.05, 1 + 0.05, 41); h, _ = radial_h(a, rr, qs); cf = np.polyfit(qs - 1, h, 6); h1 = np.polyval(cf, 0)
    hp = np.polyval(np.polyder(cf), 0) / h1; hpp = np.polyval(np.polyder(cf, 2), 0) / h1
    rb = 0.150 / 0.70   # 150 h^-1 Mpc in Gpc, h = 0.70 (model H0)
    eb = rb / t["D"]; bulk_frac = eb**2 / 10 * (hpp + 4 * hp)       # ball average of V minus V_E, / V_E
    VO = t["vE"] * c
    out[key] = dict(Mvis=t["Mvis"], rho_s=t["Os"] * rc70, Os=t["Os"], VO_c=t["vE"], VO=VO, ratio_cmb=VO / 369.82, bulk=abs(bulk_frac) * VO, bulk_sign=np.sign(bulk_frac),
                    Mc=t["Mc"], Mtot=t["Mtot"], Mc_ratio=t["Mc"] / acc["Mb"], Mtot_ratio=t["Mtot"] / acc["Mb"], Mvis_ratio=t["Mvis"] / acc["Mb"], age=t["age"], D=t["D"], hp=hp, hpp=hpp)
    print(key, {k: (float("%.4g" % v) if isinstance(v, float) else v) for k, v in out[key].items()})
print({k: float("%.4g" % v) for k, v in acc.items()})
json.dump(out, open("/workspace/cosmology/papers/A_equals_fourier/thinshell/accepted.json", "w"), indent=1, default=float)
