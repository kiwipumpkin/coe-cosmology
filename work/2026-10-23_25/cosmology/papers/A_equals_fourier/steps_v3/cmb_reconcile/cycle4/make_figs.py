import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 12})
A = [r for r in json.load(open("solveA.json")) if "err" not in r]; B = json.load(open("partB.json"))
LIM = 1e-5
PR = [("point", {}, "compact core (headline)", "k", "o"), ("pl", {"n": 2.0}, "(1) density rising toward CoE, n=2", "C0", "s"),
      ("pl", {"n": 1.0}, "(1b) rising, n=1", "C9", "v"), ("st", {}, "(2) thin Sedov-Taylor interior (no fit at D >= 49)", "C3", "^"),
      ("remnant", {"f": 0.5}, "(3) remnant core (50% of mass) + thin Sedov", "C2", "D"), ("remnant", {"f": 0.05}, "(3b) remnant core (5%) + thin Sedov", "C8", "d"),
      ("comp_st", {}, "(3c) remnant core + thin Sedov filling the cavity mass", "C1", "*")]
def rows(kind, a, s="all", ok=True):
    rr = sorted([r for r in A if r["slice"] == s and r["kind"] == kind and r["a"] == a], key=lambda r: r["D"])
    return [r for r in rr if (r["resid"] < 1e-6 and abs(r["CoE_axis_offset"]) < 5)] if ok else rr
fig, ax = plt.subplots(2, 2, figsize=(17, 12)); ax = ax.ravel()
for kind, a, lab, c, m in PR:
    rr = rows(kind, a); D = [r["D"] for r in rr]
    ax[0].semilogy(D, [abs(r["kin"][2]) for r in rr], marker=m, color=c, label=lab, ms=8)
    st = [abs(r["stat"][2]) for r in rr]
    if kind == "comp_st": ax[1].plot([], [], marker=m, color=c, label=lab + ": exactly zero for D >= 22.3")
    else: ax[1].semilogy(D, st, marker=m, color=c, label=lab, ms=8)
    ax[2].plot(D, [r["R2"] for r in rr], marker=m, color=c, label=lab, ms=8)
    rr2 = rows(kind, a, "z>0.07"); ax[2].plot([r["D"] for r in rr2], [r["R2"] for r in rr2], marker=m, color=c, ls="--", ms=6)
    ax[3].semilogy(D, [r["mcc"] for r in rr], marker=m, color=c, label=lab, ms=8)
for i in (0, 1):
    ax[i].axhline(LIM, color="0.3", lw=2, ls=":"); ax[i].set_ylim(1e-7, 5); ax[i].set_xscale("log")
    ax[i].text(14.8, 0.45 * LIM, "CMB quadrupole limit 1e-5", fontsize=11, color="0.3")
ax[0].set_title("Doppler from emitter motion on the 14.3 Gpc sphere: every profile fails by 3,000x or more")
ax[0].set_ylabel("|quadrupole of CMB dT/T|, emitter motion"); ax[0].legend(fontsize=9.5, loc="center right")
ax[1].set_ylim(1e-9, 5); ax[1].texts[0].set_position((14.8, 1.5 * LIM))
ax[1].set_title("Static gravitational redshift: set only by (interior mass - cavity mass)")
ax[1].set_ylabel("|quadrupole of CMB dT/T|, static potential"); ax[1].legend(fontsize=9.5, loc="lower right")
ax[2].set_title("Cluster fit after re-solving mass and axis (solid: all clusters, dashed: z>0.07)")
ax[2].set_ylabel("R2 against the Migkas slice profile"); ax[2].set_ylim(0.6, 0.95); ax[2].set_xscale("log"); ax[2].legend(fontsize=9.5, loc="lower left")
ax[3].set_title("Interior mass the cluster pattern needs (fit to all clusters)"); ax[3].set_ylabel("interior mass / critical mass within D")
ax[3].set_xscale("log"); ax[3].legend(fontsize=9.5, loc="upper left")
for a_ in ax:
    a_.set_xlabel("distance to the CoE, D (Gpc)"); a_.set_xticks([15, 20, 30, 50, 94]); a_.set_xticklabels(["15", "20", "30", "50", "94"]); a_.minorticks_off(); a_.grid(alpha=0.3)
fig.suptitle("Way out A: change the mass profile inside our radius (Chris's integrator, Omega_s = 0.049 shell, R = 14.3 Gpc)", fontsize=15)
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig("cycle4_A_profiles.png", dpi=110)
# ---------------- B
R = B["rows"]; D = np.array([r["D"] for r in R]); g = lambda k, l: np.array([r[k][l] for r in R])
fig, ax = plt.subplots(2, 2, figsize=(17, 12)); ax = ax.ravel()
r49 = [r for r in R if abs(r["D"] - 49.1746) < 1e-3][0]; M = B["multipliers"]
items = [("static redshift, core only", r49["stat_core"][2]), ("static redshift, core + shell cavity", r49["stat"][2]),
         ("adiabatic 1/3 of the static term", r49["stat"][2] / 3), ("growing-mode SW + ISW, linear (3.6x static)", M["growing_net_linear"] * r49["stat"][2]),
         ("growing-mode SW + ISW, exact-GR slope (6.1x)", M["growing_net_exactGR_cycle2"] * r49["stat"][2]),
         ("ISW / Rees-Sciama in the static-table picture", 0.0), ("Doppler from emitter motion", r49["kin"][2]),
         ("ionised-shell scattering (kSZ), change", -B["tau_shell"] / 2 * r49["kin"][2]), ("lensing of a uniform source", 0.0)]
y = np.arange(len(items))[::-1]
for yi, (lab, v) in zip(y, items):
    col = "C3" if v > 0 else ("C0" if v < 0 else "0.6"); ax[0].barh(yi, max(abs(v), 2e-7), color=col, height=0.4)
    ax[0].text(2.5e-7, yi + 0.27, lab + ("   (%+.1e)" % v if v != 0 else "   (zero)"), fontsize=10)
ax[0].set_xscale("log"); ax[0].set_xlim(2e-7, 3); ax[0].axvline(LIM, color="0.3", ls=":", lw=2); ax[0].set_yticks([])
ax[0].set_ylim(-0.6, len(items) - 0.2)
from matplotlib.patches import Patch; from matplotlib.lines import Line2D
ax[0].legend(handles=[Patch(color="C3", label="positive quadrupole"), Patch(color="C0", label="negative quadrupole"), Line2D([], [], ls=":", color="0.3", label="CMB limit 1e-5")], fontsize=10, loc="lower right"); ax[0].set_xlabel("|quadrupole of CMB dT/T| at D = 49 Gpc")
ax[0].set_title("Each CMB term at the headline D = 49 Gpc: only the static terms can change sign")
lin = lambda v: np.sign(v) * np.log10(1 + abs(v) / 1e-5)
for k, l, lab, c in (("stat_core", 2, "static, core only", "0.4"), ("stat_cavity", 2, "static, shell cavity", "C1"), ("stat", 2, "static, net", "C3"), ("kin", 2, "Doppler from emitter motion", "C0")):
    ax[1].plot(D, g(k, l), marker="o", color=c, label=lab)
ax[1].set_yscale("symlog", linthresh=1e-4); ax[1].axhline(0, color="k", lw=0.8); ax[1].axvline(B["static_null"]["D"], color="C3", ls="--", lw=1)
ax[1].text(27, 2e-4, "static term vanishes at D = %.2f Gpc\n(core mass = cavity mass; tolerance +/- %.0f Mpc)" % (B["static_null"]["D"], 1e3 * B["static_null"]["tol_D_Gpc"]), fontsize=10, color="C3")
ax[1].set_ylabel("CMB quadrupole (signed; linear inside +/- 1e-4)"); ax[1].set_title("Static redshift can be nulled by choosing D; the Doppler term cannot")
ax[1].legend(fontsize=10, loc="lower right")
cz = B["cancel"]; Dc = np.array([c["D"] for c in cz]); ms = np.array([c["m_star"] for c in cz])
ax[2].plot(Dc, ms, "ko-", label="factor on the static term needed to cancel the Doppler quadrupole")
for v, lab, c in ((1 / 3, "adiabatic 1/3", "C2"), (1.0, "static potential (Chris's integrator)", "C0")):
    ax[2].axhline(v, color=c, ls="--", label=lab)
ax[2].axhspan(M["growing_net_linear"], M["growing_net_exactGR_cycle2"], color="C3", alpha=0.25, label="growing-mode SW + ISW history (3.6 to 6.1)")
ax[2].set_yscale("symlog", linthresh=1); ax[2].set_ylabel("required factor (signed)"); ax[2].legend(fontsize=10, loc="lower right")
ax[2].set_title("SW + ISW + Doppler: l=2 cancels only for a growing-mode history, D ~ 40-60 Gpc")
ax[3].semilogy(Dc, [abs(c["a3"]) for c in cz], "o-", color="C4", label="octupole left after cancelling the quadrupole")
ax[3].semilogy(Dc, [abs(c["a4"]) for c in cz], "s-", color="C5", label="l = 4 left after cancelling the quadrupole")
ax[3].axhline(LIM, color="0.3", ls=":", lw=2, label="CMB limit 1e-5"); ax[3].set_ylim(3e-6, 1); ax[3].legend(fontsize=10, loc="lower left")
ax[3].set_ylabel("|residual multipole of CMB dT/T|"); ax[3].set_title("Left after cancelling l=2: 3e-3 to 0.16, i.e. 300x to 16,000x over")
for a_ in ax[1:]:
    a_.set_xscale("log"); a_.set_xlabel("distance to the CoE, D (Gpc)"); a_.set_xticks([15, 20, 30, 50, 100, 200]); a_.set_xticklabels(["15", "20", "30", "50", "100", "200"]); a_.minorticks_off(); a_.grid(alpha=0.3)
fig.suptitle("Way out B: CMB-specific terms on the headline compact core + Omega_s = 0.049 shell (cluster fit held fixed)", fontsize=15)
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig("cycle4_B_terms.png", dpi=110); print("ok")
