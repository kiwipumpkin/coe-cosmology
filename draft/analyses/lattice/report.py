"""Builds results.json, summary_table.md and the four figures from search_*.json and jitter.json."""
import json, os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("run_search.py").read().split("\nif __name__")[0])
buf=json.load(open("search_alone_with4_buf.json")); bare=json.load(open("search_alone_bare.json")) if os.path.exists("search_alone_bare.json") else []
jit=json.load(open("jitter.json"))
NICE={"SC_rocksalt":"Simple\ncubic","BCC_CsCl":"Body-\ncentred","FCC_L10":"Face-\ncentred","HCP_AB":"Hexagonal","shell_rings":"Rings\n(CoE axis)"}
ONE={"SC_rocksalt":"Simple cubic (rock-salt +/-)","BCC_CsCl":"Body-centred cubic","FCC_L10":"Face-centred cubic","HCP_AB":"Hexagonal close-packed","shell_rings":"Rings around the CoE axis"}
NPAR={k:(5 if k=="shell_rings" else 8) for k in KINDS}
NPAR_TXT={k:("spacing + 3 offsets (radial, polar, azimuthal) + amplitude; axis fixed to the CoE axis" if k=="shell_rings" else "spacing + 3 orientation angles + 3 phase offsets + amplitude") for k in KINDS}
REF=dict(start=BASE0,four_clumps=BASE4,four_clumps_uncompensated=83.4,
         best_random_alone_tuned=153.1,best_random_alone_as_drawn=175.2,best_random_with4_tuned=61.2,best_random_with4_as_drawn=62.8)
def best_of(rows,kind,t): r=[x for x in rows if x["kind"]==kind and x["target"]==t]; return min(r,key=lambda x:x["best"]) if r else None

res=dict(references=REF,param_counts=NPAR,param_notes=NPAR_TXT,
         conventions=dict(region="14.3-42.9 Gpc from us, every direction, CoE distance >= D-d = 43.9 Gpc (shell D-d..D+d and beyond D+d to D+3d); CoE interior empty",
            mass_rule="object mass = amp x rho_exp-lock(r_CoE) x V_site (amp = fraction of its cell's mass gathered into a cluster or removed from a void); voids rescaled so total void deficit = total cluster excess (net mass 0); |amp| <= 1/max(1,void_rescale) so no void below delta=-1",
            buffer_rule="main runs: every object (delta=+1 cluster / delta=-1 void sphere holding |M|) lies entirely outside the 14.3 Gpc vision sphere; 'bare' comparison: only the centre must be outside",
            tidal="field3 / residual / FLOOR from ../shell_ensemble (unchanged); point masses",
            search="per lattice type and spacing: Sobol samples over uniform rotations x phases (2048; 1024 for a<10 Gpc; 512 for rings), Nelder-Mead polish of the best starts, amplitude in closed form"),
         best_alone={},best_with4={},best_alone_bare={},all_buffer=buf,all_bare=bare,jitter=jit)
for k in KINDS:
    for t,key,rows in (("alone","best_alone",buf),("with4","best_with4",buf),("alone","best_alone_bare",bare)):
        b=best_of(rows,k,t)
        if b: res[key][k]=dict(leftover=b["best"],spacing_gpc=b["spacing_gpc"],euler_deg=b["euler_deg"],phase=b["phase"],amp=b["amp"],n_params=NPAR[k],
                              n_sites=b["n_sites"],n_clusters=b["n_clusters"],n_voids=b["n_voids"],M_obj_median_kg=b["M_obj_median_kg"],M_obj_max_kg=b["M_obj_max_kg"],
                              min_dist_from_us_gpc=b["min_dist_from_us_gpc"],frac_beyond_Dpd=b["frac_beyond_Dpd"],sobol_median=b["sobol_median"],
                              amps_over_floor=b["amps_over_floor"][:6])
# diagnostic: how much of each best lattice's effect comes from its few nearest objects (same masses, nothing refit)
def near_share(b,target):
    kind=b["kind"]; a=b["spacing_gpc"]/RV; r,info=evaluate(kind,a,b["x"],target)
    P,M=info["P"],info["M"]; d=np.linalg.norm(P,axis=1); o=np.argsort(d); out={}
    for n in (4,8):
        sel=o[:n]; Dd=d[sel]; w=(M[sel]/KG)@field3(np.arctan2(P[sel,1],P[sel,0]),np.arcsin(P[sel,2]/Dd),Dd)
        out[f"nearest_{n}_only"]=rms(target-w)/FLOOR
    out["nearest_dist_gpc"]=(d[o[:8]]*RV).round(1).tolist(); return out
for t,key,target in (("alone","best_alone",R),("with4","best_with4",R4)):
    for k in KINDS:
        b=best_of(buf,k,t); res[key][k].update(near_share(b,target))
json.dump(res,open("results.json","w"),indent=1)

def lev(M):
    if M<3e49: return "a (1e48–1e49)" if M>=1e48 else "below a"
    if M<1e51: return "b (3e49–3e50)" if M<=3.2e50 else "between b and c"
    return "c (1e51–1e52)" if M<=1e52 else "above c"
LINES=["# Crystal-lattice clusters + voids around the visible sphere: best leftover per lattice",""
   ,f"Start 237.3× floor; four fitted clumps (compensated, net mass zero) 88.8×. Leftover = RMS of the CMB residual after the model, in units of the 1e-5 floor.",""
   ,"## Lattice alone (main rule: every object fully outside our 14.3 Gpc vision sphere)",""
   ,"|lattice|best leftover|spacing a (Gpc)|orientation (Euler ZYZ°, rel. CoE axis)|fitted params|objects in region (+/−)|mass per object (median, kg)|mass level|closest object (Gpc from us)|nearest 4 objects only|typical orientation (median of search)|",
   "|-|-|-|-|-|-|-|-|-|-|-|"]
for k in KINDS:
    b=res["best_alone"].get(k)
    if not b: continue
    ori="axis fixed (CoE)" if k=="shell_rings" else "(%.0f, %.0f, %.0f)"%tuple(b["euler_deg"])
    LINES.append(f"|{ONE[k]}|**{b['leftover']:.1f}×**|{b['spacing_gpc']}|{ori}|{b['n_params']}|{b['n_sites']} ({b['n_clusters']}/{b['n_voids']})|{b['M_obj_median_kg']:.1e}|{lev(b['M_obj_median_kg'])}|{b['min_dist_from_us_gpc']:.1f}|{b['nearest_4_only']:.1f}×|{b['sobol_median']:.1f}×|")
LINES+=["","For comparison: four clumps 88.8× (5 fitted params: shared distance, 2 masses, 2 angles); best single random layout 153× with its strength tuned (1 fitted param, best of thousands), 175× as drawn.",""
   ,"## Lattice on top of the four clumps (start 88.8×)",""
   ,"|lattice|best leftover|spacing a (Gpc)|orientation (Euler ZYZ°)|fitted params (lattice only; clumps fixed)|mass per object (median, kg)|closest object (Gpc)|nearest 4 lattice objects only|",
   "|-|-|-|-|-|-|-|-|"]
for k in KINDS:
    b=res["best_with4"].get(k)
    if not b: continue
    ori="axis fixed (CoE)" if k=="shell_rings" else "(%.0f, %.0f, %.0f)"%tuple(b["euler_deg"])
    LINES.append(f"|{ONE[k]}|**{b['leftover']:.1f}×**|{b['spacing_gpc']}|{ori}|{b['n_params']} (+5 for the clumps)|{b['M_obj_median_kg']:.1e}|{b['min_dist_from_us_gpc']:.1f}|{b['nearest_4_only']:.1f}×|")
LINES+=["","For comparison: best single random layout on top of the four clumps 61× tuned / 63× as drawn.",""]
if bare:
    LINES+=["## Bare rule (only object centres outside 14.3 Gpc, as in the random 'bare' runs), lattice alone","","|lattice|best leftover|spacing|closest object (Gpc)|","|-|-|-|-|"]
    for k in KINDS:
        b=res["best_alone_bare"].get(k)
        if b: LINES.append(f"|{ONE[k]}|{b['leftover']:.1f}×|{b['spacing_gpc']}|{b['min_dist_from_us_gpc']:.2f}|")
    LINES.append("")
LINES+=["## Leftover vs spacing (lattice alone, best over orientation and phase; 7 fitted params at fixed spacing, 4 for rings)","","|spacing (Gpc)|"+"|".join(ONE[k] for k in KINDS)+"|","|-|"+"-|"*len(KINDS)]
for aG in SPACINGS:
    row=[]
    for k in KINDS:
        r=[x for x in buf if x["kind"]==k and x["spacing_gpc"]==aG and x["target"]=="alone"]; row.append(f"{r[0]['best']:.1f}" if r else "")
    LINES.append(f"|{aG}|"+"|".join(row)+"|")
LINES+=["","### same, on top of the four clumps (start 88.8×)","","|spacing (Gpc)|"+"|".join(ONE[k] for k in KINDS)+"|","|-|"+"-|"*len(KINDS)]
for aG in SPACINGS:
    row=[]
    for k in KINDS:
        r=[x for x in buf if x["kind"]==k and x["spacing_gpc"]==aG and x["target"]=="with4"]; row.append(f"{r[0]['best']:.1f}" if r else "")
    LINES.append(f"|{aG}|"+"|".join(row)+"|")
LINES+=["","## Jittered lattice (best lattice, every object shifted by a random 3D offset of RMS size f·a; 1,000 layouts per level; strength kept as fitted)",""]
for t,lab in (("alone","lattice alone"),("with4","lattice + four clumps")):
    j=jit[t]; LINES+=[f"**{lab}**: {ONE[j['kind']]}, a = {j['spacing_gpc']} Gpc; perfect lattice {j['perfect']:.1f}×, start {j['start']:.1f}×","",
        "|jitter (RMS shift)|helped / hurt vs start|median leftover|best|worst 5%|beat the perfect lattice|median with strength re-tuned|","|-|-|-|-|-|-|-|"]
    for f,o in j["levels"].items():
        if f=="0%": continue
        LINES.append(f"|{f} of a ({o['jitter_rms_gpc']:.1f} Gpc)|{o['helped_vs_start']} / {o['hurt_vs_start']}|{o['median']:.1f}×|{o['min']:.1f}×|{o['p95']:.1f}×|{o['better_than_perfect']}|{o['tuned_median']:.1f}×|")
    LINES.append("")
open("summary_table.md","w").write("\n".join(LINES))

# ---------- figures (minimal)
plt.rcParams.update({"font.size":12})
def refline(ax,y,txt,ls,c="k"):
    ax.axhline(y,ls=ls,c=c,lw=1.3); ax.text(4.55,y+2,txt,va="bottom",ha="left",fontsize=10,color=c)
fig,axs=plt.subplots(1,2,figsize=(15,6),sharey=True)
for ax,key,title,refs in ((axs[0],"best_alone","Lattice alone",[(BASE0,"start 237×","-"),(153.1,"best random layout 153×",":")]),
                          (axs[1],"best_with4","Lattice on top of the four clumps",[(BASE0,"start 237×","-"),(BASE4,"four clumps 88.8×","--"),(61.2,"best random + clumps 61×",":")])):
    ks=[k for k in KINDS if k in res[key]]; v=[res[key][k]["leftover"] for k in ks]
    bars=ax.bar(range(len(ks)),v,color="#4a7fb5",width=0.6)
    for i,(b,k) in enumerate(zip(v,ks)): ax.text(i,b+3,f"{b:.0f}×\n({res[key][k]['spacing_gpc']:g} Gpc)",ha="center",va="bottom",fontsize=10)
    for y,t,ls in refs: refline(ax,y,t,ls)
    if key=="best_alone": refline(ax,BASE4,"four clumps 88.8×","--")
    ax.set_xticks(range(len(ks))); ax.set_xticklabels([NICE[k] for k in ks],fontsize=10); ax.set_xlim(-0.6,len(ks)+2.4); ax.set_title(title)
axs[0].set_ylabel("Leftover CMB pattern (× the 1e-5 floor)\nlower is better"); axs[0].set_ylim(0,290)
fig.suptitle("Best leftover for each lattice type"); fig.tight_layout(); fig.savefig("best_by_lattice.png",dpi=130); plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(14,5.4),sharey=True)
for ax,t,title in ((axs[0],"alone","Lattice alone"),(axs[1],"with4","Lattice on top of the four clumps")):
    for k in KINDS:
        r=sorted([x for x in buf if x["kind"]==k and x["target"]==t],key=lambda x:x["spacing_gpc"])
        ax.plot([x["spacing_gpc"] for x in r],[x["best"] for x in r],marker="o",ms=4,label=ONE[k])
    ax.axhline(BASE0,c="k",lw=1.2); ax.axhline(BASE4,c="k",ls="--",lw=1.2)
    ax.text(40.5,BASE0," start",va="center",fontsize=10); ax.text(40.5,BASE4," four clumps",va="center",fontsize=10)
    ax.set_xlabel("Lattice spacing (Gpc)"); ax.set_title(title); ax.set_xlim(4,46)
axs[0].set_ylabel("Leftover (× floor), lower is better"); axs[0].set_ylim(0,260); axs[0].legend(fontsize=9,loc="lower left")
fig.suptitle("Leftover vs lattice spacing (best orientation at each spacing)"); fig.tight_layout(); fig.savefig("leftover_vs_spacing.png",dpi=130); plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(12,5),sharey=True)
for ax,t,title in ((axs[0],"alone","Lattice alone (vs start 237×)"),(axs[1],"with4","Lattice + four clumps (vs 88.8×)")):
    lv=[f for f in jit[t]["levels"] if f!="0%"]; h=[jit[t]["levels"][f]["helped_vs_start"] for f in lv]; u=[jit[t]["levels"][f]["hurt_vs_start"] for f in lv]
    x=np.arange(len(lv)); ax.bar(x-0.18,h,0.36,color="#3a9a4a",label="helped"); ax.bar(x+0.18,u,0.36,color="#c0392b",label="hurt")
    for i in x: ax.text(i-0.18,h[i]+10,str(h[i]),ha="center",fontsize=10); ax.text(i+0.18,u[i]+10,str(u[i]),ha="center",fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels([f"{f} of spacing" for f in lv]); ax.set_xlabel("How far each object is shifted (random)"); ax.set_title(title,fontsize=12)
axs[0].set_ylabel("Layouts out of 1,000"); axs[0].set_ylim(0,1100); axs[0].legend(loc="upper left")
fig.suptitle("Perfect lattice → shaken lattice → random: how often it still helps"); fig.tight_layout(); fig.savefig("jitter.png",dpi=130); plt.close(fig)

# example: best lattice alone, projected onto the ring plane (the plane of the CMB ring through us and the CoE)
b=min([x for x in buf if x["target"]=="alone"],key=lambda x:x["best"]); kind=b["kind"]; a=b["spacing_gpc"]/RV
r,info=evaluate(kind,a,b["x"],R); P=info["P"]*RV; M=info["M"]; Cg=C*RV
ex=np.cross([0,0,1.0],AX); ex/=np.linalg.norm(ex)   # in-plane axis perpendicular to CoE axis
X=P@AX; Y=P@ex
fig,ax=plt.subplots(figsize=(8.5,8))
sz=40*np.sqrt(np.abs(M)/np.abs(M).max())+4
cl=M>0; ax.scatter(X[cl],Y[cl],s=sz[cl],c="#c0392b",label="cluster (+)",alpha=0.85); ax.scatter(X[~cl],Y[~cl],s=sz[~cl],c="#2e6fb7",label="void (−)",alpha=0.85)
t=np.linspace(0,2*np.pi,400)
ax.plot(14.3*np.cos(t),14.3*np.sin(t),"k",lw=1.5,label="our visible sphere (14.3 Gpc)")
ax.plot(42.9*np.cos(t),42.9*np.sin(t),"k:",lw=1,label="edge of filled region (42.9 Gpc)")
for rr,ls in ((43.9,"--"),(72.5,"--")): ax.plot(58.2+rr*np.cos(t),rr*np.sin(t),c="gray",ls=ls,lw=1)
ax.plot([],[],c="gray",ls="--",label="shell D−d and D+d around the CoE")
ax.plot(0,0,"k*",ms=14); ax.text(1,-3,"us",fontsize=11); ax.plot(58.2,0,"kx",ms=10,mew=2); ax.text(54,3,"CoE",fontsize=11)
ax.set_xlim(-50,62); ax.set_ylim(-50,50); ax.set_aspect("equal"); ax.set_xlabel("Gpc along the line to the CoE"); ax.set_ylabel("Gpc")
ax.legend(loc="lower left",fontsize=9)
ax.set_title(f"Best lattice: {ONE[kind]}, spacing {b['spacing_gpc']:g} Gpc → {b['best']:.0f}× (start 237×)\n(seen from above: dots that look inside the black circle are above or below it, not inside it)",fontsize=11)
fig.tight_layout(); fig.savefig("example_lattice.png",dpi=130); plt.close(fig)
print("\n".join(LINES))
