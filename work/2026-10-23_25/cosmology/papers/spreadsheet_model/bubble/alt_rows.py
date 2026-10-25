import json
G = 6.674e-11; c = 2.99792458e8; GPCm = 3.0857e25
U = json.load(open("universe_real.json")); PS = json.load(open("ptshell_results.json"))["inside|Os=0.05|mig_l0|c=0.02"]
T = lambda k: json.load(open("thin_%s.json" % k.replace("|", "_"))); C = lambda k: json.load(open("cored_%s.json" % k))
def row(lab, mat, d, shape, Om, delta):
    rsD = 2 * G * d["M_kg"] / c**2 / GPCm / d["D_Gpc"]
    return dict(lab=lab, mat=mat, DR=d["DR"], D=d["D_Gpc"], shape=shape, M=d["M_sun"], Om=Om, l0=d["l0"], R2=d["R2_mig"], R2F=d["R2_four"],
                delta=delta, v2=d["vesc2"], age=d["age"], rsD=rsD, amp=d["amp"], Mc=d.get("Mc_sun"), net=(d["mcc"] - d["Os"] * d["rin"]**3) if "mcc" in d else None,
                hole=d["Os"] * d["rin"]**3 if "Os" in d else None, mcc=d.get("mcc"))
rows = [row("point + ST shell (best real matter)", "baryons", PS, "$L_c$ 0.23; shell 3.7--187", 0.05, "0"),
        row("cored $1/(1+(r/L)^2)$ + background", "LCDM", C("Ob03_c02"), "$L$ %.2f" % (C("Ob03_c02")["LD"] * C("Ob03_c02")["D_Gpc"]), 0.3, "2"),
        row("cored + background", "baryons", C("Ob005_c02"), "$L$ %.1f" % (C("Ob005_c02")["LD"] * C("Ob005_c02")["D_Gpc"]), 0.05, "2"),
        row("cored alone", "none", C("Ob0_c02"), "$L$ %.1f" % (C("Ob0_c02")["LD"] * C("Ob0_c02")["D_Gpc"]), 0.0, "2"),
        row("thin shell + core", "free", T("b|free|hug"), "shell %.0f--%.0f" % (T("b|free|hug")["Rin_Gpc"], T("b|free|hug")["Rout_Gpc"]), T("b|free|hug")["Os"], "0"),
        row("thin shell + core", "baryons", T("b|0.05|hug"), "shell %.0f--%.0f" % (T("b|0.05|hug")["Rin_Gpc"], T("b|0.05|hug")["Rout_Gpc"]), 0.05, "0"),
        row("thin shell only", "free", T("a|free|hug"), "shell %.0f--%.0f" % (T("a|free|hug")["Rin_Gpc"], T("a|free|hug")["Rout_Gpc"]), T("a|free|hug")["Os"], "0"),
        row("thin shell only", "baryons", T("a|0.05|hug"), "shell %.0f--%.0f" % (T("a|0.05|hug")["Rin_Gpc"], T("a|0.05|hug")["Rout_Gpc"]), 0.05, "0")]
def e(x):
    m, ex = ("%.2e" % x).split("e"); return r"$%s\times10^{%d}$" % (m, int(ex))
def rs(x): return ("{:,.0f}".format(x).replace(",", "{,}")) if x >= 100 else "%.2g" % x
def fmt(r):
    return (f'{r["lab"]} & {r["mat"]} & {r["DR"]:.2f} & {r["D"]:.0f} & {r["shape"]} & {e(r["M"])} & {r["Om"]:.2f} & {r["l0"]:.0f} & {r["R2"]:.3f} & {r["R2F"]:.3f} & '
            f'{r["delta"]} & {r["v2"]:.2f} & {r["age"]:.1f} & {rs(r["rsD"])}\\\\')
json.dump(rows, open("alt_rows.json", "w"), indent=1)
open("../alt_rows.tex", "w").write("\n".join(fmt(r) for r in rows) + "\n")
print(open("../alt_rows.tex").read())
# thin-shell table for uniform_bubble
th = []
for k, lab in [("a|free|hug", "shell only"), ("a|0.05|hug", "shell only"), ("b|free|hug", "shell + core"), ("b|0.05|hug", "shell + core"),
               ("a|free|marg", "shell only, margin"), ("b|free|marg", "shell + core, margin"), ("b|0.05|marg", "shell + core, margin"),
               ("a|free|hug|p2", "shell only, $p=2$"), ("b|0.05|hug|p2", "shell + core, $p=2$")]:
    d = T(k); h = d["Os"] * d["rin"]**3; mat = "baryons" if abs(d["Os"] - 0.05) < 1e-9 else "free"
    th.append(f'{lab} & {mat} & {d["DR"]:.2f} & {d["Os"]:.3f} & {d["mcc"]:.4f} & {h:.4f} & {d["mcc"] - h:+.4f} & {d["l0"]:.1f} & {d["R2_mig"]:.3f} & {d["R2_four"]:.3f} & {d["vesc2"]:.2f} & {d["age"]:.1f}\\\\'.replace("+-", "-").replace(" -", " $-$").replace(" +", " $+$"))
open("../thin_rows.tex", "w").write("\n".join(th) + "\n"); print("\n".join(th))
