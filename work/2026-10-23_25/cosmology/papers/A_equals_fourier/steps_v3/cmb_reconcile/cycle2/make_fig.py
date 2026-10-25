"""Cycle 2 figure: (left) exact-GR CMB quadrupole vs core mass fraction (far core, omb geometry + rel fit) against linear expectations;
(middle) quadrupole vs emission redshift (does the potential deepen into the past?); (right) R2 vs CMB quadrupole for B (far core + tilt)
and C (near-core variant). Reads A_ltb_cmb.json, A_full_*.json (if present), B_hybrid_fit*.json, C_*.json."""
import json, glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A = {r["name"]: r for r in json.load(open("A_ltb_cmb.json")) if "fail" not in r}
for fn in glob.glob("A_full_omb_core_x*.json"):
    r = json.load(open(fn)); A[f"omb|core-only|full x{r['mcc'] / 0.0017418224560860097:.2g}"] = dict(mcc=r["mcc"], cb_a=r["cb_a"])
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
om = sorted([(r["mcc"], abs(r["cb_a"][2])) for k, r in A.items() if k.startswith("omb")])
m = np.array([x[0] for x in om]); a = np.array([x[1] for x in om])
mm = np.geomspace(1e-7, 0.2, 50)
ax[0].loglog(m, a, "o-", color="C1", label="exact GR (LTB), far core, omb geometry (last point: hot spot dominated)")
r = A["rel|core-only"]; ax[0].loglog([r["mcc"]], [abs(r["cb_a"][2])], "s", color="C2", ms=9, label="exact GR, relativistic headline fit")
ax[0].loglog(mm, 5.6 * mm, "k--", label="linear, static potential (today's field)")
ax[0].loglog(mm, 3160 * mm, "k:", label="linear, constant core mass (potential grows into the past)")
ax[0].axhline(1e-5, color="r", ls="--", label="CMB bound 1e-5")
ax[0].axvspan(0.0016, 0.0027, color="C1", alpha=0.15); ax[0].axvspan(0.06, 0.1, color="C2", alpha=0.15)
ax[0].set_xlabel("core mass / critical mass within D"); ax[0].set_ylabel("|CMB quadrupole| (dipole removed)"); ax[0].legend(fontsize=7, loc="upper left"); ax[0].set_ylim(1e-7, 1e3)
ax[0].set_title("A. Far core: exact GR vs linear", fontsize=10)
for k, col in (("omb|core-only|Mc x0.01", "C1"), ("rel|core-only", "C2")):
    L = A[k]["ladder"]; z = np.array([v["zbar"] for v in L.values()]); a2 = np.array([abs(v["a2"]) for v in L.values()])
    ax[1].loglog(z, a2 / a2[-1], "o-", color=col, label=("far core, omb geometry (mass x0.01)" if col == "C1" else "far core, relativistic fit"))
zz = np.geomspace(2, 1100, 50); ax[1].loglog(zz, (1 + zz) / 1008, "k:", label="expected if potential grew as 1/a")
ax[1].set_xlabel("emission redshift"); ax[1].set_ylabel("quadrupole relative to z = 1000 emission"); ax[1].legend(fontsize=7.5); ax[1].set_ylim(1e-3, 3)
ax[1].set_title("Potential is constant before z ~ 20 (growing mode)", fontsize=10)
B = json.load(open("B_hybrid_fit.json"))["ladder_cavity"]
bx = [34 * b["mcc"] for b in B]; by = [b["R2"] for b in B]
ax[2].semilogx(bx, by, "o-", color="C0", label="B: far core + tilt, CMB from exact-GR slope (estimate)")
for c in json.load(open("C_near_core.json"))["cmb"] + json.load(open("C_small.json"))["cmb"]:
    pass
Cs = json.load(open("C_small.json")); Cn = json.load(open("C_near_core.json"))
pts = [(Cn["best"]["R2_s"], Cn["cmb"][0]["cb_a"][2], Cn["best"]["D"], Cn["best"]["z"])]
pts += [(f["R2_s"], c["cb_a"][2], f["D"], f["z"]) for f, c in zip(Cs["fits"], Cs["cmb"]) if f["mcc"] == 0.08]
dep = {d["D"]: d for d in json.load(open("C_depth.json")) if d["z"] == 0.075}
for R2, a2, D, z in pts:
    ax[2].semilogx([abs(a2)], [R2], "D", color="C3"); ax[2].annotate(f"D={D:g}", (abs(a2), R2), fontsize=8, xytext=(-8, 7), textcoords="offset points")
    if D in dep: ax[2].semilogx([abs(a2)], [dep[D]["R2_s"]], "x", color="C3", ms=8)
ax[2].plot([], [], "D", color="C3", label="C: near-core variant (D in Gpc), sampled at depth ~D/3")
ax[2].plot([], [], "x", color="C3", label="same variant sampled at z = 0.075")
ax[2].axvline(1e-5, color="r", ls="--"); ax[2].set_xlabel("|CMB quadrupole|"); ax[2].set_ylabel("R² vs Migkas (geometry + sine)")
ax[2].legend(fontsize=7, loc="lower left"); ax[2].set_ylim(0, 1); ax[2].set_title("B and C: fit quality vs CMB residual", fontsize=10)
plt.tight_layout(); plt.savefig("cycle2_summary.png", dpi=130)
