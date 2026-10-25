import json, numpy as np
S = json.load(open("res_B.json")); T = json.load(open("tables.json")); scan = json.load(open("scan.json"))
F = json.load(open("check_frw.json")); Cc = json.load(open("check_centre.json")); I = S["info"]
def e(x, n=2):
    if x is None or not np.isfinite(x): return "--"
    if x == 0: return "0"
    s = "%.*e" % (n - 1, x); m, ex = s.split("e"); ex = int(ex)
    return "$%s$" % m if ex == 0 else "$%s\\times10^{%d}$" % (m, ex)
def f(x, n=3): return "--" if not np.isfinite(x) else ("$%.*f$" % (n, x))
out = {}
# multipoles
O = T["obs"]; lin = [r for r in scan if r["label"] == "Om_bg=0.3, dM=0.01 Mfit"][0]; lin9 = [r for r in scan if r["label"] == "Om_bg=0.9, dM=0.01 Mfit"][0]
fa, f9 = T["scale_fac"]
rows = [("Best-fit LTB, fixed $t_{\\rm rec}$, dust frame", "--", S["c_a"][1], S["c_a"][2], S["c_a"][3], S["c_a"][4]),
        ("Best-fit LTB, fixed $t_{\\rm rec}$, boosted", f(S["beta_a"], 3), 0.0, S["cb_a"][2], S["cb_a"][3], S["cb_a"][4]),
        ("$\\Omega=0.3$ background + lump, scaled, boosted", e(fa * lin["beta_a"]), fa * lin["c1_a"], fa * lin["c2_a_boost"], fa * lin["c3_a_boost"], np.nan),
        ("$\\Omega=0.9$ background + lump, scaled, boosted", e(f9 * lin9["beta_a"]), f9 * lin9["c1_a"], f9 * lin9["c2_a_boost"], f9 * lin9["c3_a_boost"], np.nan),
        ("Observed (axisymmetric equivalent)", "$1.23\\times10^{-3}$", O["1"], O["2"], O["3"], O["4"])]
s = ""
for r in rows:
    s += "%s & %s & %s & %s & %s & %s\\\\\n" % (r[0], r[1], e(abs(r[2])) if r[2] != 0 else "0", e(abs(r[3])), e(abs(r[4])), e(abs(r[5])))
out["mult"] = s
# nearby
AC = json.load(open("angle_check.json"))
def row(lab, v, ph, c=None, k=None): return "%s & %s & %s & %s & %.0f$^\\circ$ & %s & %s\\\\\n" % (lab, f(v[0], 4), f(v[1], 4), f(v[2], 4), ph, "--" if c is None else f(c, 2), f(AC[k], 3))
s = row("Data, $A_F$", T["AF4"], T["AF4_phase2"], None, "data")
s += row("Newtonian, AF$\\times$sine", T["newt_A"], T["newt_A_phase2"], T["corr_rows"]["newt_A"], "newtA")
s += row("Newtonian, AF only", T["newt_AF"], T["newt_AF_phase2"], T["corr_rows"]["newt_AF"], "newtAF")
for z in ("0.05", "0.1", "0.3"):
    s += row("LTB, $A\\propto D_L^{-2}$, $z=%s$" % z, T["ltb_A"][z], T["ltb_A_phase2"][z], T["corr_rows"]["ltb_" + z], "ltb_%g" % float(z))
s += row("LTB, $\\delta H_0/H_0$, $z=0.1$", T["ltb_dH"]["0.1"], T["ltb_A_phase2"]["0.1"], T["corr_rows"]["ltb_dH"], "ltb_dH")
out["nearby"] = s
# scan
sel = ["M=1 Mfit, D/L=14.2", "M=1 Mfit, D/L=5.0", "M=1 Mfit, D/L=2.0", "M=0.1 Mfit, D/L=14.2", "M=0.01 Mfit, D/L=14.2", "M=0.01 Mfit, D/L=2.0",
       "D=2Dfit, L=2Lfit, M/Mcrit=0.0024", "Om_bg=0.3, dM=1 Mfit", "Om_bg=0.3, dM=0.1 Mfit", "Om_bg=0.3, dM=0.01 Mfit", "Om_bg=0.3, dM=0.001 Mfit",
       "Om_bg=0.3, dM=0.01 Mfit, L=3Lfit", "Om_bg=0.9, dM=0.01 Mfit"]
nice = {"M=1 Mfit, D/L=14.2": "lump, $M$, $D/L=14.2$ (fit)", "M=1 Mfit, D/L=5.0": "lump, $M$, $D/L=5$", "M=1 Mfit, D/L=2.0": "lump, $M$, $D/L=2$",
        "M=0.1 Mfit, D/L=14.2": "lump, $0.1M$, $D/L=14.2$", "M=0.01 Mfit, D/L=14.2": "lump, $0.01M$, $D/L=14.2$", "M=0.01 Mfit, D/L=2.0": "lump, $0.01M$, $D/L=2$",
        "D=2Dfit, L=2Lfit, M/Mcrit=0.0024": "lump, $2D$, $2L$, $0.1M/M_{\\rm crit}$",
        "Om_bg=0.3, dM=1 Mfit": "$\\Omega=0.3$ + lump $M$", "Om_bg=0.3, dM=0.1 Mfit": "$\\Omega=0.3$ + lump $0.1M$", "Om_bg=0.3, dM=0.01 Mfit": "$\\Omega=0.3$ + lump $0.01M$",
        "Om_bg=0.3, dM=0.001 Mfit": "$\\Omega=0.3$ + lump $0.001M$", "Om_bg=0.3, dM=0.01 Mfit, L=3Lfit": "$\\Omega=0.3$ + lump $0.01M$, $3L$", "Om_bg=0.9, dM=0.01 Mfit": "$\\Omega=0.9$ + lump $0.01M$"}
s = ""
for lab in sel:
    r = [x for x in scan if x["label"] == lab][0]
    s += "%s & %s & %s & %s & %s & %s & %s\\\\\n" % (nice[lab], e(r["Om_loc"]), f(r["sky_b"], 3), e(r["A_m2"]), f(r["beta_a"], 3) if abs(r["beta_a"]) > 1e-3 else e(r["beta_a"]),
                                                   e(abs(r["c2_a_boost"])), f(abs(r["c2_a_boost"]) / r["A_m2"], 2))
    if lab == "D=2Dfit, L=2Lfit, M/Mcrit=0.0024": s += "\\midrule\n"
out["scan"] = s
# model table
out["model"] = ("Case A: areal radii with $1+2E<0$ [Gpc] & 1.67--14.29\\\\\n\\midrule\nCase B: age $t_0$ [Gyr] & %.2f\\\\\nLabel of observer shell $r_{\\rm obs}$ [Gpc] & %.2f\\\\\n$2GM/(c^2D)$ & %.2f\\\\\n"
                "Local density / critical & %s\\\\\n$H_r/H_0$, $H_\\perp/H_0$ at observer & %.4f, %.4f\\\\\nCentral density today / critical & %.2f\\\\\n"
                "$t_{\\rm rec}$, definition (a) [Myr] & %.1f\\\\\nSky fraction with a density-defined LSS & %.3f\\\\\nLargest $\\psi$ reaching $\\rho_{\\rm rec}$ & %.1f$^\\circ$\\\\\n") % (
    I["t0_Gyr"], I["r_obs"], I["twoM_over_R"], e(I["rho_loc_over_crit"]), I["Hr"], I["Hperp"], I["rho_c_centre_today"], S["t_rec_Myr"], S["sky_frac_b"], S["psi_max_b_deg"])
out["valid"] = ("FRW: $\\max|\\Delta z|/(1+z)$, $z=%.0f$ & %s\\\\\nFRW, largest $|c_\\ell|$, $\\ell\\le3$ & %s\\\\\n"
                "FRW, $D_A$ relative error, $z\\le1$ & %s\\\\\nNull drift, FRW / best fit & %s, %s\\\\\n"
                ) % (F["z_rec_true"], e(F["z_rec_rays_max_rel_err"]), e(max(abs(x) for x in F["multipoles"])), e(max(abs(x) for x in F["DA_rel_err"])), e(F["null_err"]), e(S["null_err"]))
for rr, lab in (("0.001", "1"), ("0.01", "10"), ("0.1", "100")):
    out["valid"] += "Centred, $r_{\\rm obs}=%s$ Mpc: $c_1$, $c_2$ & %s, %s\\\\\n" % (lab, e(Cc[rr]["c1"]), e(Cc[rr]["c2"]))
for k, v in out.items(): open("tab_%s.tex" % k, "w").write(v)
print(out["mult"]); print(out["nearby"]); print(out["scan"]); print(out["model"]); print(out["valid"])
