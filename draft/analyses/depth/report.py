import json, glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":15}); RV=14.3; BASE=237.27
W=json.load(open("web_shells.json")); CU=json.load(open("cumulative_curves.json"))
S={}
for f in glob.glob("fits_*.json"):
    d=json.load(open(f)); S[f[5:-5]]=d
NAMES={"k4_mv_coe":"4 objects, masses or voids","k4_mo_coe":"4 objects, masses only","k3_mv_coe":"3 objects, masses or voids",
       "k4_mv":"4 objects, masses or voids (CoE interior allowed)","k4_mo":"4 objects, masses only (CoE interior allowed)","k3_mv":"3 objects, masses or voids (CoE interior allowed)"}
def interp_cross(x,y,target):   # first x where y <= target (y decreasing), linear interpolation
    for i in range(len(x)):
        if y[i]<=target+1e-9:
            if i==0: return x[0]
            return float(x[i-1]+(x[i]-x[i-1])*(y[i-1]-target)/(y[i-1]-y[i]))
    return None
res={"A_fits":{}}
for tag,d in S.items():
    x=np.array([r["R_cut_gpc"] for r in d["scan"]]); y=np.array([r["rms_over_floor"] for r in d["scan"]]); y=np.minimum.accumulate(y)
    asym=y[-1]; gain=BASE-y; g_as=BASE-asym
    res["A_fits"][tag]=dict(label=NAMES.get(tag,tag),R_cut_gpc=x.tolist(),best_leftover_x=y.tolist(),asymptote_x=float(asym),
        flatten_within5pct_gpc=interp_cross(x,y,1.05*asym),
        R90_of_improvement_gpc=interp_cross(x,BASE-gain,BASE-0.9*g_as) if True else None,
        R95_of_improvement_gpc=interp_cross(x,BASE-gain,BASE-0.95*g_as),
        still_falling_last_step_pct=float(100*(y[-2]-y[-1])/y[-2]),
        positions_at_largest_cut=d["scan"][-1]["dist_gpc"],masses_at_largest_cut=d["scan"][-1]["M_kg"],
        coe_dist_at_largest_cut=d["scan"][-1]["dist_from_CoE_gpc"],max_coe_violation_gpc=max(r["max_coe_violation_gpc"] for r in d["scan"]) if d["coe_constraint"] else None)
res["B_and_C"]=W
json.dump(res,open("results.json","w"),indent=1)

# ---- figure 1
fig,ax=plt.subplots(figsize=(12.5,7.8))
for tag,col,mk,lab in (("k4_mv_coe","#1f77b4","o","4 objects placed (best fit)"),("k3_mv_coe","#ff7f0e","s","3 objects placed (best fit)")):
    a=res["A_fits"][tag]; ax.plot(a["R_cut_gpc"],a["best_leftover_x"],mk+"-",color=col,lw=2.8,ms=9,label=lab)
    xf=a["flatten_within5pct_gpc"]; yf=np.interp(xf,a["R_cut_gpc"],a["best_leftover_x"])
    ax.plot([xf,xf],[0,yf],ls="--",color=col,lw=1.8)
    ax.annotate(f"flattens at {xf:.0f} Gpc\n({xf/RV:.1f}× visible radius)",xy=(xf,yf),xytext=(xf+2.5,yf+(55 if tag=="k4_mv_coe" else 45)),color=col,fontsize=13,weight="bold",arrowprops=dict(arrowstyle="->",lw=1.5,color=col))
ax.axhline(BASE,color="0.5",ls=":",lw=1.2); ax.text(15.2,BASE+4,"start, nothing placed (237×)",fontsize=12,color="0.4")
ax.set_xlabel("objects allowed only within this distance of us (Gpc)"); ax.set_ylabel("best leftover CMB pattern (× the 10⁻⁵ target)")
ax.set_ylim(0,260); ax.set_xlim(14,59); ax.grid(alpha=.3); ax.legend(loc="upper right",fontsize=13)
sec=ax.secondary_xaxis("top",functions=(lambda x:x/RV,lambda x:x*RV)); sec.set_xlabel("in units of the visible radius (14.3 Gpc)")
ax.set_title("How far out do placed objects need to go?  ('flattens' = within 5% of the final value)",fontsize=15,pad=12)
fig.tight_layout(); fig.savefig("depth_cutoff.png",dpi=110); plt.close(fig)

# ---- figure 2 (extended web, objects start 15.0 Gpc = 14.3 + typical object radius)
D=W["B2_extended_web"]["inner_edge_15.0"]
fig,ax=plt.subplots(figsize=(12.5,7.5))
for g,col in (("m=1-2","#ff7f0e"),("m=3-5","#1f77b4")):
    c=CU[g]; ax.plot(c["R"],100*np.array(c["share"]),color=col,lw=3,label={"m=1-2":"lowest modes (m = 1–2)","m=3-5":"modes m = 3–5"}[g])
    for lvl,ls in ((90,"--"),(95,":")):
        r=D[g][f"R{lvl}_gpc"]; ax.plot([r],[lvl],"o",color=col,ms=10)
        if g=="m=3-5": ax.text(r-0.6,lvl+(1.5 if lvl==95 else -5),f"{r:.0f} Gpc",color=col,fontsize=13,weight="bold",ha="right")
        else: ax.text(r+0.6,lvl+(1.5 if lvl==95 else -5),f"{r:.0f} Gpc",color=col,fontsize=13,weight="bold",ha="left")
for lvl,ls in ((90,"--"),(95,":")): ax.axhline(lvl,color="0.4",ls=ls,lw=1.3); ax.text(61,lvl+0.8,f"{lvl}%",fontsize=12,color="0.3")
ax.set_xlim(14,65); ax.set_ylim(0,102); ax.grid(alpha=.3)
ax.set_xlabel("distance from us, R (Gpc)"); ax.set_ylabel("share of the random web's effect\ncoming from inside R (%)")
sec=ax.secondary_xaxis("top",functions=(lambda x:x/RV,lambda x:x*RV)); sec.set_xlabel("in units of the visible radius (14.3 Gpc)")
ax.set_title("Where does a random web of clusters and voids act from?",fontsize=16,pad=12)
ax.text(40,40,"random clusters and voids start 15 Gpc from us\n(just outside our 14.3 Gpc visible sphere,\nallowing for the objects' own size)",fontsize=12,color="0.3",ha="center")
ax.legend(loc="lower right",fontsize=13); fig.tight_layout(); fig.savefig("effect_vs_distance.png",dpi=110); plt.close(fig)
print(json.dumps(res["A_fits"],indent=1)[:4000])

# ---- summary table
A=res["A_fits"]; B=W
L=["# Tidal influence depth ('multipole horizon')\n","Visible radius d = 14.3 Gpc. Leftover = RMS of the CMB residual after SR + central mass, in units of the 1e-5 floor (start 237.3×).\n",
   "## A) Best leftover vs R_cut (objects allowed only within R_cut of us)\n",
   "|series|"+"|".join(f"{x:g}" for x in A["k4_mv_coe"]["R_cut_gpc"])+"|final|flattens (≤5% above final)|90% of improvement by|95% by|","|-|"+"-|"*(len(A["k4_mv_coe"]["R_cut_gpc"])+4)]
for tag in ("k4_mv_coe","k4_mo_coe","k3_mv_coe","k4_mv","k4_mo","k3_mv"):
    if tag not in A: continue
    a=A[tag]; f=lambda v:"–" if v is None else f"{v:.1f} Gpc ({v/RV:.2f} d)"
    L.append(f"|{a['label']}|"+"|".join(f"{y:.1f}" for y in a["best_leftover_x"])+"|"*(1+len(A['k4_mv_coe']['R_cut_gpc'])-len(a['best_leftover_x']))+f"{a['asymptote_x']:.1f}×|{f(a['flatten_within5pct_gpc'])}|{f(a['R90_of_improvement_gpc'])}|{f(a['R95_of_improvement_gpc'])}|")
L.append("\nSeries without '(CoE interior allowed)' keep every object ≥ D−d = 43.9 Gpc from the CoE (same region as the random web). Curves are made non-increasing (a larger R_cut can always reuse the smaller-R_cut solution).\n")
L.append("## B1) Random web, k=0, N=200, medium, 50% voids (1,000 layouts, centres > 14.3 Gpc), split into distance shells\n")
L.append("|shell (Gpc)|objects per layout|share of web effect (power)|share m=1–2|share m=3–5|pattern RMS per 1e51 kg (× floor)|m=1 / m=2 / m=3 / m=4 / m=5 amplitude per 1e51 kg (× floor)|median change in leftover from this shell alone (× floor)|\n|-|-|-|-|-|-|-|-|")
for sh in B["B1_ensemble_k0_N200_b_fv50"]["shells"]:
    L.append(f"|{sh['range_gpc'][0]}–{sh['range_gpc'][1]}|{sh['mean_objects']:.1f}|{100*sh['share_of_power']:.2f}%|{100*sh['share_m12']:.1f}%|{100*sh['share_m35']:.1f}%|{sh['rms_per_1e51kg_x_floor']:.1f}|"+" / ".join(f"{x:.2f}" for x in sh["mode_rms_per_1e51kg_x_floor"][:5])+f"|{sh['median_abs_change_in_leftover_x_floor']:.2f}|")
L.append("\n## B2) Depth = radius containing 90% / 95% of the random web's effect (web extended to 400 Gpc, uniform in volume, CoE interior excluded)\n")
L.append("|objects start at|all modes (leftover) R90 / R95|m=1–2 R90 / R95|m=3–5 R90 / R95|\n|-|-|-|-|")
for k,lab in (("inner_edge_14.3","14.3 Gpc (centres anywhere outside; near-ring objects dominate)"),("inner_edge_15.0","15.0 Gpc (14.3 + typical object radius) ← headline"),("inner_edge_15.5","15.5 Gpc (larger voids)")):
    d=B["B2_extended_web"][k]; g=lambda x:f"{d[x]['R90_gpc']:.1f} / {d[x]['R95_gpc']:.1f} Gpc ({d[x]['R90_gpc']/RV:.2f} / {d[x]['R95_gpc']/RV:.2f} d)"
    L.append(f"|{lab}|{g('all (leftover RMS)')}|{g('m=1-2')}|{g('m=3-5')}|")
d=B["B1_cumulative_depth_within_42.9"]; L.append(f"|B1 ensemble as run (14.3–42.9 Gpc)|{d['all (leftover RMS)']['R90_gpc']:.1f} / {d['all (leftover RMS)']['R95_gpc']:.1f}|{d['m=1-2']['R90_gpc']:.1f} / {d['m=1-2']['R95_gpc']:.1f}|{d['m=3-5']['R90_gpc']:.1f} / {d['m=3-5']['R95_gpc']:.1f}|")
L.append("\n## C) Analytic falloff vs measured\n\n|mode|expected amplitude falloff|measured slope 60–180 Gpc|measured slope 15–30 Gpc|far-field depth R90 / R95 (× inner edge r0)|\n|-|-|-|-|-|")
for m in range(1,6):
    an=B["C_analytic"][f"m={m}"]
    L.append(f"|m={m}|(R_vis/r)^{an['amp_falloff_power']}|{B['C_measured_power_law_slopes_60_180gpc'][f'm={m}']:.2f}|{B['C_measured_power_law_slopes_15_30gpc'][f'm={m}']:.2f}|{an['R90_over_r0']:.2f} / {an['R95_over_r0']:.2f}|")
open("summary_table.md","w").write("\n".join(L)+"\n")
