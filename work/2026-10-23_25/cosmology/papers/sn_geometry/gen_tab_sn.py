import numpy as np, json
from scipy.integrate import quad
E = json.load(open("extra.json")); P = json.load(open("post_z0.023.json"))["variants"]; J = json.load(open("joint.json"))
D = json.load(open("direction.json")); F = json.load(open("fit_z0.023.json"))
V = {v["name"]: v for v in P}; e = E["z0.023"]; e1 = E["z0.01"]; N = e["N"]; lnN = np.log(N)
H70 = 977.8 / 70
def age(Om, OL): return H70 * quad(lambda a: 1 / np.sqrt(Om / a + (1 - Om - OL) + OL * a * a), 0, 1)[0]
c0 = e["LCDM"][0]; c01 = e1["LCDM"][0]
rows = [(r"flat $\Lambda$CDM, $\Omega_m=%.3f$" % e["LCDM"][1], 2, c0, c01, age(e["LCDM"][1], 1 - e["LCDM"][1])),
        (r"Milne ($\Omega_m=0$)", 1, e["Milne"], e1["Milne"], age(0, 0)),
        (r"open FRW, baryons only ($\Omega_m=0.05$)", 1, e["FRWb"], e1["FRWb"], age(0.05, 0)),
        (r"Einstein--de Sitter", 1, e["EdS"], e1["EdS"], age(1, 0)),
        (r"CoE lump, no background, axis fixed", 1, e["ChrisCoE"], e1["ChrisCoE"], V["Chris Om0 CoE"]["t0_Gyr"]),
        (r"CoE lump + baryons only ($\Omega_b=0.05$), axis fixed", 1, e["ChrisCoE_b"], e1["ChrisCoE_b"], V["Chris Om0.05 CoE"]["t0_Gyr"]),
        (r"CoE lump + $\Lambda$CDM-equivalent matter ($\Omega=0.3$), axis fixed", 1, e["ChrisCoE_03"], e1["ChrisCoE_03"], V["Chris Om0.3 CoE"]["t0_Gyr"]),
        (r"CoE lump + free $\Omega$ ($%.3f$), axis fixed" % V["Chris Omfree CoE"]["Om"], 2, V["Chris Omfree CoE"]["chi2"], None, V["Chris Omfree CoE"]["t0_Gyr"]),
        (r"CoE lump, no background, axis free", 3, V["Chris Om0 free"]["chi2"], None, V["Chris Om0 free"]["t0_Gyr"]),
        (r"CoE lump + baryons only, axis free", 3, V["Chris Om0.05 free"]["chi2"], None, V["Chris Om0.05 free"]["t0_Gyr"]),
        (r"CoE lump + $\Lambda$CDM-equivalent matter, axis free", 3, V["Chris Om0.3 free"]["chi2"], None, V["Chris Om0.3 free"]["t0_Gyr"]),
        (r"CoE lump + free $\Omega$ ($%.3f$), axis free" % V["Chris Omfree free"]["Om"], 4, V["Chris Omfree free"]["chi2"], None, V["Chris Omfree free"]["t0_Gyr"]),
        (r"free LTB ($D,L,M,\Omega,\ell_0,b_0$)", 7, e["best"]["chi2"], e1["best"]["chi2"], e["best"]["t0_Gyr"]),
        (r"joint SN+cluster geometry, no background", 5, J[0]["sn"], None, J[0]["t0"]),
        (r"joint SN+cluster geometry, baryons only", 5, J[1]["sn"], None, J[1]["t0"])]
f = lambda x, d=1: "--" if x is None else ("%.*f" % (d, x))
L = []
for n, k, c, c1, t in rows:
    L.append("%s & %d & %s & %s & %s & %s & %s & %s & %s \\\\" % (n, k, f(c), f(c - c0), f(c - c0 + 2 * (k - 2)), f(c - c0 + (k - 2) * lnN),
             f(None if c1 is None else c1 - c01), f(t, 2), f(t * 70 / 73, 2)))
open("tab_fits.tex", "w").write("\n".join(L))
# parameters
pr = [("free LTB, $D$ free", F), ]
best = min(F, key=lambda r: r["chi2"]); pure = min([r for r in F if r["mode"] == "pure"], key=lambda r: r["chi2"])
second = min([r for r in F if r["mode"] == "pure" and r["LD"] > 0.1], key=lambda r: r["chi2"])
omf = [r for r in F if r["mode"] == "Omfree"][0]
P2 = [("CoE lump (Chris geometry)", 41.53, 2.919 / 41.53, 0.02395, 0.0, 193.17, 0.0, 1337.82),
      ("free LTB, $D$ free", best["D"], best["LD"], best["mc"], best["Om"], best["l"], best["b"], best["chi2"]),
      ("free LTB, $D$ fixed, compact", 41.53, pure["LD"], pure["mc"], 0.0, pure["l"], pure["b"], pure["chi2"]),
      ("free LTB, $D$ fixed, broad", 41.53, second["LD"], second["mc"], 0.0, second["l"], second["b"], second["chi2"]),
      ("free LTB, $\\Omega$ free", 41.53, omf["LD"], omf["mc"], omf["Om"], omf["l"], omf["b"], omf["chi2"]),
      ("joint SN+cluster, no background", 41.53, J[0]["LD"], J[0]["mc"], 0.0, J[0]["l"], J[0]["b"], J[0]["sn"]),
      ("joint SN+cluster, baryons only", 41.53, J[1]["LD"], J[1]["mc"], 0.05, J[1]["l"], J[1]["b"], J[1]["sn"])]
def ang(l1, b1, l2=193.17, b2=0.0):
    u = lambda l, b: np.array([np.cos(np.radians(b)) * np.cos(np.radians(l)), np.cos(np.radians(b)) * np.sin(np.radians(l)), np.sin(np.radians(b))])
    return np.degrees(np.arccos(np.clip(u(l1, b1) @ u(l2, b2), -1, 1)))
open("tab_par.tex", "w").write("\n".join("%s & %.2f & %.3f & %.4f & %.3f & %.0f & %.0f & %.0f & %.1f \\\\" % (n, d, ld, mc, om, l, b, ang(l, b), c) for n, d, ld, mc, om, l, b, c in P2))
# dipoles
d = D["z0.023"]; R = []
for nm, x in [(r"$\Lambda$CDM, all $z$", d["dip_LCDM"]), (r"$\Lambda$CDM, $z<0.1$", d["dip_LCDM_0_0.1"]), (r"$\Lambda$CDM, $0.1<z<0.3$", d["dip_LCDM_0.1_0.3"]),
              (r"$\Lambda$CDM, $z>0.3$", d["dip_LCDM_0.3_10"]), ("Milne, all $z$", d["dip_Milne"]), ("free LTB, all $z$", e["best_dip"]),
              (r"$\Lambda$CDM, all $z$, $z_{\rm min}=0.01$", D["z0.01"]["dip_LCDM"]), (r"$\Lambda$CDM, $z<0.1$, $z_{\rm min}=0.01$", D["z0.01"]["dip_LCDM_0_0.1"])]:
    R.append("%s & %d & %.4f & %.4f & %.0f & %.0f & %.0f & %.0f & %.1f \\\\" % (nm, x.get("N", N) if "N" in x else (N if "0.01" not in nm else 1590), x["A"], x["sigA"], x["l"], x["b"], x["dang68"], ang(x["l"], x["b"]), x["dchi2"]))
open("tab_dip.tex", "w").write("\n".join(R))
Q = []
for zk, lab in [("z0.023", "0.023"), ("z0.01", "0.01")]:
    for dn in ["CMB", "Migkas", "CoE", "bestdip"]:
        c = D[zk]["colin_" + dn]; h = D[zk].get("hemi_" + dn)
        hs = "%.3f & %.3f & %.3f & %.3f" % (h[0][0], h[0][1], h[1][0], h[1][1]) if h else "-- & -- & -- & --"
        name = {"CMB": "CMB dipole", "Migkas": "Migkas cluster axis", "CoE": "CoE", "bestdip": "best $q_d$ axis"}[dn]
        Q.append("%s & %s & %.3f & %.2f & %.2f & %.1f & %s \\\\" % (lab, name, c["qm"], c["qd"], c["sig_qd"], c["dchi2"], hs))
open("tab_q.tex", "w").write("\n".join(Q))
M = [("isotropic", E["mig_iso_chi2"], None, None),
     ("CoE lump, no background, CoE axis", V["Chris Om0 CoE"]["mig_chi2"], V["Chris Om0 CoE"]["mig_r"], V["Chris Om0 CoE"]["chi2"]),
     ("CoE lump + baryons only, CoE axis", V["Chris Om0.05 CoE"]["mig_chi2"], V["Chris Om0.05 CoE"]["mig_r"], V["Chris Om0.05 CoE"]["chi2"]),
     ("CoE lump + $\\Lambda$CDM-equivalent matter, CoE axis", V["Chris Om0.3 CoE"]["mig_chi2"], V["Chris Om0.3 CoE"]["mig_r"], V["Chris Om0.3 CoE"]["chi2"]),
     ("CoE lump, no background, SN-preferred axis", V["Chris Om0 free"]["mig_chi2"], V["Chris Om0 free"]["mig_r"], V["Chris Om0 free"]["chi2"]),
     ("free LTB (SN best)", e["best_mig"][0], e["best_mig"][1], e["best"]["chi2"]),
     ("joint SN+cluster, no background", J[0]["mig"], J[0]["mig_r"], J[0]["sn"]),
     ("joint SN+cluster, baryons only", J[1]["mig"], J[1]["mig_r"], J[1]["sn"])]
open("tab_mig.tex", "w").write("\n".join("%s & %.1f & %s & %s & %s \\\\" % (n, m, f(r, 2), f(c), f(None if c is None else c - c0)) for n, m, r, c in M))
print(open("tab_fits.tex").read()); print(open("tab_par.tex").read()); print(open("tab_dip.tex").read()); print(open("tab_q.tex").read()); print(open("tab_mig.tex").read())
