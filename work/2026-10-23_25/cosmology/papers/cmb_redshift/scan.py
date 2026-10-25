import numpy as np, json, models as m, analyze as an, ltb, warnings, sys
warnings.filterwarnings("ignore")
rows = []
def run(fam, label, p, t0, ro):
    try:
        c = an.cmb(p, t0, ro, N=32); nb_ = an.nearby(p, t0, ro, zlist=(0.05, 0.1, 0.3), N=16)
    except Exception as e:
        print(label, "FAILED", e); return
    four = an.fourier(nb_["ell"], nb_["A"][:, 1])
    rhoc = 3 * m.H0**2 / (8 * np.pi)
    r = dict(family=fam, label=label, t0_Gyr=t0 * ltb.GYR, Om_loc=ltb.density(ro, t0, p) / rhoc,
             A_m1=four[0], A_m2=four[1], A_m3=four[2], dH_c2=nb_["dH_leg"][1][2], dH_c1=nb_["dH_leg"][1][1],
             sky_b=c["sky_frac_b"], c1_a=c["c_a"][1], c2_a=c["c_a"][2], beta_a=c["beta_a"], c2_a_boost=c["cb_a"][2],
             c3_a_boost=c["cb_a"][3], c1_b=c["c_b"][1], c2_b=c["c_b"][2], beta_b=c["beta_b"], c2_b_boost=c["cb_b"][2],
             c3_b_boost=c["cb_b"][3])
    rows.append(r); print(json.dumps(r, default=float)); sys.stdout.flush()
# (i) Chris's family: pure e^{-r/L} lump, Lagrangian labels, Milne energy, uniform bang
for mc in (1.0, 0.1, 0.01):
    for DL in (14.23, 5.0, 2.0):
        L = m.D / DL; Mx = m.MC * mc * m.MCRIT
        p, t0, ro = m.model_B(M=Mx, L=L); run("lump", "M=%.3g Mfit, D/L=%.1f" % (mc, DL), p, t0, ro)
p, t0, ro = m.model_B(M=m.MC * 0.1 * m.H0**2 * (2 * m.D)**3 / 2, D=2 * m.D, L=m.L * 2); run("lump", "D=2Dfit, L=2Lfit, M/Mcrit=0.0024", p, t0, ro)
# (ii) homogeneous open background (Omega=0.3) plus the lump at the fitted D, L
for f in (1.0, 0.3, 0.1, 0.03, 0.01, 0.001):
    p, t0, ro = m.model_C(f * m.M); run("bg+lump", "Om_bg=0.3, dM=%.3g Mfit" % f, p, t0, ro)
for f in (0.1, 0.01):
    p, t0, ro = m.model_C(f * m.M, L=3 * m.L); run("bg+lump", "Om_bg=0.3, dM=%.3g Mfit, L=3Lfit" % f, p, t0, ro)
    p, t0, ro = m.model_C(f * m.M, Om=0.9); run("bg+lump", "Om_bg=0.9, dM=%.3g Mfit" % f, p, t0, ro)
json.dump(rows, open("scan.json", "w"), indent=1, default=float)
