import json, numpy as np
B = json.load(open("bubble_results.json")); AR = json.load(open("analysis_rows.json"))
def e(x):
    m, ex = ("%.2e" % x).split("e"); return r"$%s\times10^{%d}$" % (m, int(ex))
def pct(x): 
    v = 100 * x
    return ("%.3g" % v) if abs(v) < 1000 else "%.0f" % v
def s(x): return ("%.3f" % x).replace("-", "$-$")
def age(d): return "%.1f" % d["age"] + (r"$^\dagger$" if d["age"] < 12 else "")
def tt(d): return r"$\infty$" if not np.isfinite(d["t_turn"]) else "%.1f" % d["t_turn"]
H = [("Previous best, no background", "--", "prev|Ob=0|mean"), ("Previous best + background", "LCDM-equiv.", "prev|Ob=0.3|mean"),
     ("Unconstrained", "LCDM-equiv.", "mean|c=None|a|Ob=0.3|mig_l0"),
     (r"$\le2\%$, no $\ell_0$", "LCDM-equiv.", "mean|c=0.02|a|Ob=0.3|mig"), (r"$\le2\%$", "LCDM-equiv.", "mean|c=0.02|a|Ob=0.3|mig_l0"),
     (r"$\le2\%$, vs $A_F$", "LCDM-equiv.", "mean|c=0.02|a|Ob=0.3|four_l0"), (r"$\le2\%$, $\Omega_b$ free", "free", "mean|c=0.02|a|Ob=free|mig_l0"),
     (r"$\le2\%$, $\Omega_b$ free, bound, age$\ge$12", "free", "mean|c=0.02|b|Ob=free|mig_l0"),
     (r"$\le0.5\%$, no $\ell_0$", "LCDM-equiv.", "mean|c=0.005|a|Ob=0.3|mig"), (r"$\le0.5\%$", "LCDM-equiv.", "mean|c=0.005|a|Ob=0.3|mig_l0"),
     (r"$\le0.5\%$, vs $A_F$", "LCDM-equiv.", "mean|c=0.005|a|Ob=0.3|four_l0"), (r"$\le0.5\%$, $\Omega_b$ free", "free", "mean|c=0.005|a|Ob=free|mig_l0"),
     (r"$\le0.5\%$, $\Omega_b$ free, bound, age$\ge$12", "free", "mean|c=0.005|b|Ob=free|mig_l0"),
     (r"Literal $kR$, $\le2\%$, no $\ell_0$", "LCDM-equiv.", "kR|c=0.02|a|Ob=0.3|mig"), (r"Literal $kR$, $\le2\%$", "LCDM-equiv.", "kR|c=0.02|a|Ob=0.3|mig_l0")]
PRB = json.load(open("prof_results.json"))
extra = [(r"$\le2\%$, no $\ell_0$", "baryons only", PRB["exp|c=0.02|mig|mean|Ob=0.05"]), (r"$\le2\%$", "baryons only", PRB["exp|c=0.02|mig_l0|mean|Ob=0.05"])]
rows = [(a, b, B.get(k)) for a, b, k in H if B.get(k)]
rows = rows[:11] + extra + rows[11:]
t1, t2 = [], []
for lab, mat, d in rows:
    L = d.get("L_Gpc", d["LD"] * d["D_Gpc"])
    t1.append(f'{lab} & {mat} & {d["DR"]:.2f} & {d["LD"]:.3f} & {d["D_Gpc"]:.1f} & {L:.2f} & {d["mc"]:.3g} & {e(d["M_kg"])} & {e(d["M_sun"])} & {d["Ob"]:.2f} & {d["l0"]:.1f} & {s(d["amp"])} & {d["phi"]:.2f}\\\\')
    t2.append(f'{lab} & {mat} & {s(d["R2_mig"])} & {s(d["R2_four"])} & {pct(d["contrast"])} & {pct(d["c34"])} & {d["vesc2"]:.3g} & {age(d)} & {tt(d)} & {"yes" if d["vrmin"] > 0 else "no"}\\\\')
open("../ub_t1.tex", "w").write("\n".join(t1) + "\n"); open("../ub_t2.tex", "w").write("\n".join(t2) + "\n")
t3 = []
for kind, lab, mat, d in AR:
    if not d: continue
    shape = {"exp": "$L$=%.3f" % d.get("LD", 0), "pl": "$L$=%.3f, $n$=%.1f" % (d.get("LD", 0), d.get("n", 0)), "gauss": "$L$=%.3f" % d.get("LD", 0),
             "cored": "$L$=%.3f" % d.get("LD", 0), "shell": "$r_s$=%.2f, $w$=%.3f" % (d.get("rs", 0), d.get("w", 0)),
             "st53": "$R_s$=%.3f" % d.get("rst", 0), "st43": "$R_s$=%.3f" % d.get("rst", 0)}.get(kind)
    if kind == "ptshell": shape = "$L_c$=%.3f, $R$=%.2f--%.2f" % (d["Lc"], d["rin"], d["rout"])
    mcoe = d["mcc"] if kind == "ptshell" else d["mc"]
    om = d["Os"] if kind == "ptshell" else d["Ob"]
    lab = {"power law $(1+r/L)^{-n}$": "power law", "cored $1/(1+r^2/L^2)$": "cored"}.get(lab, lab)
    lab = lab.replace("Sedov--Taylor $\\gamma=5/3$", "ST $\\gamma=5/3$").replace("Sedov--Taylor $\\gamma=4/3$", "ST $\\gamma=4/3$").replace("point + ST shell", "pt+shell")
    mat = {"LCDM-equiv.": "LCDM", "baryons only": "baryons"}.get(mat, mat)
    p3 = lambda x: "0" if abs(100*x) < 1e-3 else pct(x).replace("-", "$-$")
    t3.append(f'{lab} & {mat} & {shape} & {d["DR"]:.2f} & {mcoe:.3g} & {om:.2f} & {d["l0"]:.1f} & {s(d["R2_mig"])} & {s(d["R2_four"])} & {p3(d["contrast"])} & {p3(d["c34"])} & {d["psi"]:.2f} & {d["vesc2"]:.3g} & {age(d)}\\\\')
open("../ub_t3.tex", "w").write("\n".join(t3) + "\n")
print("\n".join(t3))
