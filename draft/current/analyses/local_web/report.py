"""Collect raw_*.npz -> results.json, summary_table.md, helped_vs_hurt_local_web.png, example_web.png"""
import glob, json, re, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from web_core import *
H0h=0.7; MPC_PER_GPC=1000.0
VREG=region_volume()
# ---------- prior results (read-only) for comparison
prior=dict(
 surround_N200_medium=dict(helped=363,hurt=637,median=240.3,min=210.2,best_scaled_min=168.4,src="../shell_ensemble/surround/summary_table.md (bare, surround+gradient, N=200, level b 3e49-3e50 kg)"),
 surround_N500_medium=dict(helped=243,hurt=757,median=249.1,min=204.9,best_scaled_min=175.4,src="same, N=500 level b"),
 surround_N200_tiny=dict(helped=464,hurt=536,median=237.3,min=234.2,best_scaled_min=168.5,src="same, N=200 level a 1e48-1e49"),
 surround_N500_tiny=dict(helped=516,hurt=484,median=237.2,min=233.5,best_scaled_min=183.3,src="same, N=500 level a"),
 surround_N200_large=dict(helped=0,hurt=1000,median=990.5,min=262.7,best_scaled_min=169.8,src="same, N=200 level c 1e51-1e52"),
 with4_N200_medium=dict(helped=180,hurt=820,src="../shell_ensemble/with_clumps (k=4, baseline 88.8x, N=200 level b)"),
 with4_N500_medium=dict(helped=57,hurt=943,src="same, N=500 level b"),
 four_clumps=88.8, start=BASE0)
groups={}
for f in sorted(glob.glob("raw_*.npz")):
    b=f.split("/")[-1]
    m=re.match(r"raw_s([0-9.]+)_(full|tube|[abc]\.npz)",b); mg=re.match(r"raw_gap([0-9.]+)_s([0-9.]+)",b)
    if m: key=(float(m.group(1)),"check" if "check" in b else ("tube" if m.group(2).endswith(".npz") else m.group(2)))
    elif mg: key=(float(mg.group(2)),"gap"+mg.group(1).rstrip("."))
    else: continue
    z=np.load(f); groups.setdefault(key,[]).append({k:z[k] for k in z.files})
def cat(lst,k): return np.concatenate([d[k] for d in lst])
res={"start_x_floor":BASE0,"four_clumps_x_floor":BASE4,"region_volume_gpc3":VREG,"pack":PACK,"contrasts_fraction_of_cell_mass":contrasts(),"prior":prior,"runs":{}}
rows=[]
for (s,mode),lst in sorted(groups.items(),key=lambda t:(-t[0][0],t[0][1])):
    n=len(cat(lst,"f_rmsw")); r=dict(s_gpc=s,mode=mode,n=int(n))
    for p,name in (("f","foam"),("u","unstructured")):
        for j,base in (("0",BASE0),("4",BASE4)):
            L=cat(lst,f"{p}_left{j}"); B=cat(lst,f"{p}_best{j}"); Mu=np.abs(cat(lst,f"{p}_mult{j}"))
            r[f"{name}_vs{'237' if j=='0' else '88.8'}"]=dict(helped=int((L<base).sum()),hurt=int((L>=base).sum()),median=float(np.median(L)),p5=float(np.percentile(L,5)),
                best_as_drawn=float(L.min()),best_scaled_median=float(np.median(B)),best_scaled_min=float(B.min()),
                median_abs_multiplier=float(np.median(Mu)),multiplier_at_best_scaled=float(cat(lst,f"{p}_mult{j}")[np.argmin(B)]),
                median_change=float(np.median(L-base)),mean_abs_change=float(np.mean(np.abs(L-base))))
        r[f"{name}_web_rms_median"]=float(np.median(cat(lst,f"{p}_rmsw")))
        A=np.concatenate([np.array(d[f"{p}_amps0"]) for d in lst]); r[f"{name}_median_amps_m1_6"]=np.median(A,0)[:6].tolist()
    for k in ("n_cells","n_clusters","n_fil","n_wall","med_void_def","med_cluster","med_Rv","min_dist"):
        r[k]=float(np.median(cat(lst,k)))
    r["max_abs_net_rel"]=float(np.max(np.abs(cat(lst,"net_rel"))))
    r["full_region_equiv_voids"]=VREG/s**3; r["full_region_equiv_clusters"]=6.768*VREG/s**3
    r["void_R_eff_hinv_mpc"]=r["med_Rv"]*MPC_PER_GPC*H0h
    res["runs"][f"s{s:g}_{mode}"]=r; rows.append(r)
res["base_amps_m1_6"]=(amps(R)/FLOOR)[:6].tolist()
# power-law of web size vs cell size (tube/full production runs)
pr=[r for r in rows if r["mode"] in ("full","tube")]
ss=np.array([r["s_gpc"] for r in pr]); fw=np.array([r["foam_web_rms_median"] for r in pr]); uw=np.array([r["unstructured_web_rms_median"] for r in pr])
sel=ss<=1.0
if sel.sum()>=2:
    pf=np.polyfit(np.log(ss[sel]),np.log(fw[sel]),1); pu=np.polyfit(np.log(ss[sel]),np.log(uw[sel]),1)
    res["scaling_s_le_1"]=dict(foam_slope=float(pf[0]),unstr_slope=float(pu[0]),
        foam_web_rms_at_s0p0375=float(np.exp(np.polyval(pf,np.log(0.0375)))),unstr_web_rms_at_s0p0375=float(np.exp(np.polyval(pu,np.log(0.0375)))))
try: res["tube_check"]=json.load(open("check_tube.json"))
except Exception: pass
json.dump(res,open("results.json","w"),indent=1)

# ---------- tables
def lab(r):
    s=r["s_gpc"]; t=f"{s:g} Gpc"
    if abs(s-0.05)<1e-9: t+=" (real web: voids R≈20 h⁻¹Mpc)"
    if r["mode"].startswith("gap"): return t+f" [ring tube, cells kept ≥{r['mode'][3:]} Gpc outside the sphere]"
    return t+(" [full region]" if r["mode"]=="full" else " [ring tube]" if r["mode"]=="tube" else " [full-region check]")
md=[]
md.append(f"# Random local-web ('foam') layouts around the visible sphere\n\nStart: {BASE0:.1f}× the 1e-5 floor. Four fitted clumps (compensated): {BASE4:.1f}×. Region: 14.3–30 Gpc from us, CoE interior excluded (volume {VREG:,.0f} Gpc³).\n")
md.append("Packing used (shares of each cell's own mass and volume): voids 78% vol / 23% mass (δ = −0.70); walls 6% / 8%; filaments 15.5% / 40%; clusters (web junctions) 0.5% / 28.6%. Every cell is exactly mass-compensated, so net mass is zero.\n")
md.append("## Helped / hurt vs 237× (as drawn)\n\n|cell size|layouts|foam: helped / hurt|foam median|foam best|foam best-scaled (multiplier)|scattered: helped / hurt|scattered median|scattered best|scattered best-scaled|")
md.append("|-|-|-|-|-|-|-|-|-|-|")
for r in rows:
    f=r["foam_vs237"]; u=r["unstructured_vs237"]
    md.append(f"|{lab(r)}|{r['n']}|{f['helped']} / {f['hurt']}|{f['median']:.1f}×|{f['best_as_drawn']:.1f}×|{f['best_scaled_min']:.1f}× ({f['multiplier_at_best_scaled']:+.2g})|{u['helped']} / {u['hurt']}|{u['median']:.1f}×|{u['best_as_drawn']:.1f}×|{u['best_scaled_min']:.1f}× ({u['multiplier_at_best_scaled']:+.2g})|")
md.append("\nbest-scaled = one free strength multiplier on the whole web (least squares). A multiplier of +1 is the web as drawn; |multiplier| ≫ 1 means contrasts far beyond the accepted packing (voids would need δ < −1).\n")
md.append("## On top of the four clumps (start 88.8×)\n\n|cell size|layouts|foam: helped / hurt|foam median|foam best|foam best-scaled|scattered: helped / hurt|scattered median|scattered best|")
md.append("|-|-|-|-|-|-|-|-|-|")
for r in rows:
    f=r["foam_vs88.8"]; u=r["unstructured_vs88.8"]
    md.append(f"|{lab(r)}|{r['n']}|{f['helped']} / {f['hurt']}|{f['median']:.1f}×|{f['best_as_drawn']:.1f}×|{f['best_scaled_min']:.1f}×|{u['helped']} / {u['hurt']}|{u['median']:.1f}×|{u['best_as_drawn']:.1f}×|")
md.append("\n## Size of the web's own CMB pattern and what was placed\n\n|cell size|web pattern RMS, foam (× floor)|web pattern RMS, scattered|voids placed (median per layout)|cluster junctions placed|filament pieces|wall pieces|voids if the whole 14.3–30 Gpc region were filled|void R_eff (median)|void deficit (median)|cluster junction mass (median)|")
md.append("|-|-|-|-|-|-|-|-|-|-|-|")
for r in rows:
    md.append(f"|{lab(r)}|{r['foam_web_rms_median']:.1f}|{r['unstructured_web_rms_median']:.1f}|{r['n_cells']:,.0f}|{r['n_clusters']:,.0f}|{r['n_fil']:,.0f}|{r['n_wall']:,.0f}|{r['full_region_equiv_voids']:.2g}|{r['med_Rv']*1000:.0f} Mpc ({r['void_R_eff_hinv_mpc']:.0f} h⁻¹Mpc)|{r['med_void_def']:.2g} kg ({r['med_void_def']/MSUN:.1g} M☉)|{r['med_cluster']:.2g} kg ({r['med_cluster']/MSUN:.1g} M☉)|")
if "scaling_s_le_1" in res:
    sc=res["scaling_s_le_1"]
    md.append(f"\nWeb pattern RMS falls with cell size roughly as s^{sc['foam_slope']:.2f} (foam) and s^{sc['unstr_slope']:.2f} (scattered) for s ≤ 1 Gpc. Extrapolated to 15 h⁻¹Mpc voids (s = 0.0375 Gpc): foam {sc["foam_web_rms_at_s0p0375"]:.1f}×, scattered {sc["unstr_web_rms_at_s0p0375"]:.1f}× the floor.\n")
md.append("\n## Earlier random runs for comparison (bare, surround + gradient, 14.3–42.9 Gpc, 1,000 layouts)\n\n|run|helped / hurt vs 237×|median|best as drawn|best-scaled|")
md.append("|-|-|-|-|-|")
for k in ("surround_N200_tiny","surround_N500_tiny","surround_N200_medium","surround_N500_medium","surround_N200_large"):
    p=prior[k]; md.append(f"|{k.replace('surround_','').replace('_',' ')}|{p['helped']} / {p['hurt']}|{p['median']}×|{p['min']}×|{p['best_scaled_min']}×|")
md.append(f"|with four clumps, N200 medium|{prior['with4_N200_medium']['helped']} / {prior['with4_N200_medium']['hurt']} (vs 88.8×)||||")
md.append(f"|with four clumps, N500 medium|{prior['with4_N500_medium']['helped']} / {prior['with4_N500_medium']['hurt']} (vs 88.8×)||||")
open("summary_table.md","w").write("\n".join(md)+"\n")

# ---------- figure: helped vs hurt (foam vs scattered), per 1,000, vs 237x
prod=[r for r in rows if r["mode"] in ("full","tube")]
lbl=[("real web\n0.05 Gpc" if abs(r["s_gpc"]-0.05)<1e-9 else f"{r['s_gpc']:g} Gpc") for r in prod]
fig,axs=plt.subplots(1,2,figsize=(12,4.8),sharey=True)
for ax,(name,title) in zip(axs,(("foam","Cosmic-web foam (cells)"),("unstructured","Same pieces, scattered at random"))):
    x=np.arange(len(prod)); h=np.array([1000*r[f"{name}_vs237"]["helped"]/r["n"] for r in prod]); u=1000-h
    ax.bar(x-0.2,h,0.4,color="#2e9e44",label="helped"); ax.bar(x+0.2,u,0.4,color="#d43c3c",label="hurt")
    for xi,hv,uv in zip(x,h,u):
        ax.text(xi-0.2,hv+12,f"{hv:.0f}",ha="center",fontsize=9); ax.text(xi+0.2,uv+12,f"{uv:.0f}",ha="center",fontsize=9)
    ax.axhline(500,color="grey",ls="--",lw=1)
    ax.set_xticks(x); ax.set_xticklabels(lbl,fontsize=9); ax.set_title(title); ax.set_xlabel("cell size (bigger cells = coarser web)")
    ax.set_ylim(0,1100)
axs[0].set_ylabel("layouts out of 1,000"); axs[0].legend(loc="upper left",frameon=False)
fig.suptitle("Random local-web layouts: helped or hurt the 237× leftover",fontsize=13)
fig.tight_layout(); fig.savefig("helped_vs_hurt_local_web.png",dpi=130); plt.close(fig)

# ---------- example web slice (s = 2 Gpc, full region), |z| < 1 Gpc
X,M,rad,info,parts=foam(2.0,np.random.default_rng([2,0]),return_parts=True)
fig,ax=plt.subplots(figsize=(7.2,7.2)); zc=1.0
Pv,Mv,Rvv=parts["void"]; k=np.abs(Pv[:,2])<zc
for p,rv in zip(Pv[k],Rvv[k]): ax.add_patch(plt.Circle(p[:2],rv,color="#9cc3e6",alpha=0.55,lw=0))
Pf,Mf=parts["fil"]; k=np.abs(Pf[:,2])<zc; ax.scatter(Pf[k,0],Pf[k,1],s=2,color="#e8a33d",label="filament pieces")
Pc,Mc=parts["clu"]; k=np.abs(Pc[:,2])<zc; ax.scatter(Pc[k,0],Pc[k,1],s=4+30*Mc[k]/Mc.max(),color="#c0392b",label="cluster junctions")
ax.add_patch(plt.Circle((0,0),RMIN,fill=False,color="k",lw=1.5)); ax.add_patch(plt.Circle((0,0),30,fill=False,color="k",ls="--",lw=1))
ax.add_patch(plt.Circle(C_G[:2],A_IN_G,fill=False,color="grey",ls=":",lw=1.2))
ax.scatter([0],[0],marker="*",s=80,color="k"); ax.text(0,-2.2,"us",ha="center")
ax.scatter([],[],s=60,color="#9cc3e6",label="voids"); ax.legend(loc="upper right",frameon=False,fontsize=9)
ax.text(0,RMIN+0.6,"visible sphere 14.3 Gpc",ha="center",fontsize=9); ax.text(0,30.6,"30 Gpc",ha="center",fontsize=9)
u=C_G[:2]/np.linalg.norm(C_G[:2]); ax.annotate("toward CoE (interior left empty)",xy=tuple(22*u),xytext=tuple(5*u+np.array([0,-6])),fontsize=9,arrowprops=dict(arrowstyle="->",color="grey"))
ax.set_xlim(-33,33); ax.set_ylim(-33,33); ax.set_aspect("equal"); ax.set_xlabel("Gpc"); ax.set_ylabel("Gpc")
ax.set_title("Example web, 2 Gpc cells (slice through the CMB ring plane)")
fig.tight_layout(); fig.savefig("example_web.png",dpi=130); plt.close(fig)
print(open("summary_table.md").read())
