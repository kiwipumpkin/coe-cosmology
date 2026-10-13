"""Masses-on-top-of-a-fixed-background report: web_add/*_omp0.json -> web_*.png + web_add/summary_table{SUFX}.json"""
import os,json,glob,re,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from common import *
plt.rcParams.update({"font.size":8,"savefig.dpi":200})
CAT=[("b_random","random growth (web)","#7f8c8d"),("b_foam_rim","foam-rule growth (web)","#5dade2"),("b_cmb_aware","CMB-aware growth","#7d3c98"),
     ("shallow_cmb_aware","CMB-aware, shallow voids","#b48fc9"),("wave","radial wave, 150 Mpc","#d35400"),("lattice","quasi-lattice, 120 Mpc","#16a085"),("none","no background (SR + central mass only)","black")]
VARN={"":"—","_core":"core filled","_inside":"our sphere filled","_tidal":"tidal-depth fill"}; VMK={"":"o","_core":"s","_inside":"^","_tidal":"o"}
SUFX=os.environ.get("WSUF","_ex1.5")
def parse(t):
    for c,_,_ in CAT:
        if t.startswith(c+"_s"): r=t[len(c)+2:]; m=re.match(r"(\d+)(.*)",r); return c,int(m.group(1)),m.group(2)
D=[]
for f in sorted(glob.glob(f"{OUT}/web_add/*_omp0{SUFX}.json")):
    d=json.load(open(f)); c,s,v=parse(d["tag"]); d.update(cat=c,seed=s,var=v); D.append(d)
cats=[c for c in CAT if any(d["cat"]==c[0] for d in D)]
# ---- 1: background leftover and bodies needed (stacked)
fig,axs=plt.subplots(2,1,figsize=(7,7.5),sharex=True)
for i,(c,nm,col) in enumerate(cats):
    for d in [d for d in D if d["cat"]==c]:
        j=i+0.08*(list(VARN).index(d["var"])-1)
        axs[0].scatter(j,d["background"]["total"],c=col,marker=VMK[d["var"]],s=26,edgecolor="k",lw=0.3)
        axs[1].scatter(j,d["N_live"],c=col,marker=VMK[d["var"]],s=26,edgecolor="k" if d["reached"] else "red",lw=0.3 if d["reached"] else 1.2)
axs[0].axhline(1,color="green",lw=1.5); axs[0].axhline(237.3,color="0.5",ls=":",lw=0.8); axs[0].set_yscale("log"); axs[0].set_ylabel("leftover before adding (× floor)\nSR + central mass + background")
axs[0].set_title("Background alone",fontsize=9); axs[1].set_title("Bodies added on top to reach under 1×",fontsize=9)
axs[1].set_yscale("log"); axs[1].set_ylabel("bodies needed (masses + voids)"); axs[1].set_xticks(range(len(cats))); axs[1].set_xticklabels([n for _,n,_ in cats],rotation=25,ha="right",fontsize=7)
hd=[Line2D([],[],marker=VMK[v],ls="",color="0.5",mec="k",label=n) for v,n in VARN.items()]+[Line2D([],[],color="green",lw=1.5,label="1× (CMB level)"),Line2D([],[],color="0.5",ls=":",label="237× (SR + central mass alone)"),Line2D([],[],marker="o",ls="",mfc="white",mec="red",mew=1.2,label="did not reach 1× (stopped)")]
fig.legend(handles=hd,loc="lower center",ncol=3,frameon=False,fontsize=7); fig.tight_layout(rect=(0,0.07,1,1)); fig.savefig(f"{OUT}/web_bodies_needed{SUFX}.png",bbox_inches="tight"); plt.close(fig)
# ---- 2: spectra before / after (medians per background type), stacked
mm=np.arange(1,13)
fig,axs=plt.subplots(2,1,figsize=(7,8),sharex=True)
for ax,key,ttl in ((axs[0],"background","Before adding: SR + central mass + background"),(axs[1],"final","After adding bodies (until under 1×, or until no allowed body helps): what is left")):
    for c,nm,col in cats:
        A=np.array([d[key]["amp"] for d in D if d["cat"]==c]); h=np.array([d[key]["rms_13up"] for d in D if d["cat"]==c])
        ax.plot(mm,np.median(A,0),"-o",color=col,ms=3,lw=1.2); ax.plot([14],[np.median(h)],"D",color=col,ms=5)
    ax.axhline(1,color="green",lw=1.5); ax.set_yscale("log"); ax.set_ylabel("amplitude (× floor)"); ax.set_title(ttl,fontsize=9)
axs[1].set_xticks(list(mm)+[14]); axs[1].set_xticklabels([str(m) for m in mm]+["13+\n(rms)"]); axs[1].set_xlabel("mode m around the ring (median over runs)")
hd=[Line2D([],[],color=col,marker="o",ms=3,label=nm) for c,nm,col in cats]+[Line2D([],[],color="green",lw=1.5,label="1× (CMB level)")]
fig.legend(handles=hd,loc="lower center",ncol=2,frameon=False,fontsize=7); fig.tight_layout(rect=(0,0.1,1,1)); fig.savefig(f"{OUT}/web_spectrum{SUFX}.png",bbox_inches="tight"); plt.close(fig)
# ---- 3: symmetric share before/after
fig,ax=plt.subplots(figsize=(7,4.2))
for i,(c,nm,col) in enumerate(cats):
    for d in [d for d in D if d["cat"]==c]:
        ax.scatter(i-0.15,d["background"]["sym"],c=col,marker=VMK[d["var"]],s=22,edgecolor="k",lw=0.3,alpha=0.6)
        ax.scatter(i+0.15,d["final"]["sym"],c=col,marker=VMK[d["var"]],s=26,edgecolor="k",lw=0.8)
ax.axhline(0.5,color="0.5",ls=":",lw=0.8); ax.axhspan(0.9,1.0,color="0.9",zorder=0); ax.text(-0.4,0.92,"near-symmetric (≥90%)",fontsize=6.5,color="0.4")
ax.set_ylim(0,1); ax.set_ylabel("share of m=1–12 power that is mirror-symmetric\nabout the CoE axis (ℓ 15°/195°)"); ax.set_xticks(range(len(cats))); ax.set_xticklabels([n for _,n,_ in cats],rotation=25,ha="right",fontsize=7)
hd=[Line2D([],[],marker="o",ls="",color="0.6",alpha=0.6,mec="k",label="before adding (left of each pair)"),Line2D([],[],marker="o",ls="",color="0.6",mec="k",mew=0.8,label="after adding to 1× (right)"),Line2D([],[],color="0.5",ls=":",label="50% = no preference")]
fig.legend(handles=hd,loc="lower center",ncol=3,frameon=False,fontsize=7); fig.tight_layout(rect=(0,0.08,1,1)); fig.savefig(f"{OUT}/web_symmetry{SUFX}.png",bbox_inches="tight"); plt.close(fig)
# ---- 4: side views + ring of what is left, one per background type (both-filled, first seed)
pick=[]
for c,nm,col in cats:
    cand=[d for d in D if d["cat"]==c]; cand=sorted(cand,key=lambda d:(d["var"] not in ("","_tidal"),d["seed"])); pick.append((cand[0],nm,col))
fig=plt.figure(figsize=(10.5,5.2*int(np.ceil(len(pick)/2)))); GO=fig.add_gridspec(int(np.ceil(len(pick)/2)),2,hspace=0.42,wspace=0.45)
for k,(d,nm,col) in enumerate(pick):
    sub=GO[k//2,k%2].subgridspec(2,1,height_ratios=[1.7,0.8],hspace=0.45); ax=fig.add_subplot(sub[0]); tag=d["tag"]
    f=np.load(f"{OUT}/final_{tag}.npz"); core_open=d["var"] in ("","_core")
    if "Xp" in f:
        v=f["kp"]=="void"; Xv,rv,Xc=f["Xp"][v],f["rp"][v],f["Xp"][~v]
    elif "Rv" in f: Xv=f["X"]; rv=np.full(len(Xv),float(f["Rv"])); Xc=np.zeros((0,3))
    else: Xv=np.zeros((0,3)); rv=np.zeros(0); Xc=np.zeros((0,3))
    sv=(np.abs(Xv[:,2])<1.5)&(np.linalg.norm(Xv,axis=1)<=30); sc=(np.abs(Xc[:,2])<1.5)&(np.linalg.norm(Xc,axis=1)<=30)
    ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(Xv[sv],rv[sv])],fc="#d6eaf8",ec="#85c1e9",lw=0.15,zorder=2))
    if len(Xc): ax.scatter(Xc[sc,0],Xc[sc,1],s=0.5,c="#e6b0aa",lw=0,zorder=3)
    if "zoom" in f:   # wave / lattice: inset showing the shells / lattice in the ring plane
        ins=ax.inset_axes([-0.42,0.55,0.4,0.4]); h=float(f["zoom_half"]); c0=f["zoom_center"]
        ins.imshow(f["zoom"].T,origin="lower",extent=[c0[0]-h,c0[0]+h,c0[1]-h,c0[1]+h],cmap="RdBu_r",vmin=-0.15,vmax=0.15)
        th=np.linspace(0,2*np.pi,400); ins.plot(RMIN*np.cos(th),RMIN*np.sin(th),color="#1f4e79",lw=1); ins.set_xlim(c0[0]-h,c0[0]+h); ins.set_ylim(c0[1]-h,c0[1]+h)
        ins.set_xticks([]); ins.set_yticks([]); ins.set_title("3 Gpc zoom at ℓ=195°: red over, blue under",fontsize=5)
        ax.text(-33,-31,"whole filled region (core + shell + our sphere)\ncarries the smooth %s field, contrast ±%.2f;\nown ring pattern %.1f×"%("150 Mpc shell" if "wave" in tag else "120 Mpc lattice",float(f["A"]),rms(f["w"])/FLOOR),fontsize=5.5,color=col,ha="right",va="bottom",bbox=dict(fc="white",ec="none",alpha=0.8,pad=0.5))
    X=np.array(d["X"]); M=np.array(d["M"]); live=np.abs(M)>=1e49; X,M=X[live],M[live]; far=np.linalg.norm(X,axis=1)>30; nfar=int(far.sum()); X,M=X[~far],M[~far]
    vm=M<0; rho=seq_core.rho_gpc(X)*seq_core.GPC**3; Rvd=(np.abs(M)/(0.7*rho*4/3*np.pi))**(1/3)
    ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(X[vm],Rvd[vm])],fc="none",ec="#1a5276",lw=1.0,zorder=6))
    ax.scatter(X[~vm,0],X[~vm,1],s=12+30*np.abs(M[~vm])/2e52,c="#c0392b",edgecolor="k",lw=0.3,zorder=7)
    ax.add_patch(Circle((0,0),RMIN,fill=False,ec="#1f4e79",lw=1.3,zorder=5)); ax.plot(0,0,"k+",ms=5)
    frame30(ax); ax.tick_params(labelsize=5.5)
    for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(2)
    ax.set_title(f"{nm} ({VARN[d['var']]}, seed {d['seed']})\nbefore {d['background']['total']:,.1f}×; added {d['n_mass']} masses + {d['n_void']} voids\n"+"(all ≥%.1f Gpc from the ring in 3D, projected here; %d beyond 30 Gpc not shown)"%(d.get("ring_clearance_gpc",0.1),nfar),fontsize=6,color=col,fontweight="bold")
    r_=fig.add_subplot(sub[1]); res=np.array(d["resid"])/FLOOR
    r_.axhspan(-1,1,color="green",alpha=0.25,lw=0); r_.axhline(0,color="0.6",lw=0.5); r_.plot(np.arange(360),res,color=col,lw=0.9)
    lim=1.3*max(np.abs(res).max(),2); r_.set_yscale('symlog',linthresh=1,linscale=0.6); r_.set_ylim(-lim,lim); r_.set_xlim(0,360); r_.set_xticks([0,90,180,270,360]); r_.tick_params(labelsize=5.5)
    fs=d["final"]; r_.set_title(f"left {d['left']:.2f}×: m1–5 {fs['rms_m1_5']:.2f}, m6–12 {fs['rms_m6_12']:.2f}, 13+ {fs['rms_13up']:.2f}; symmetric share {fs['sym']*100:.0f}%",fontsize=6.3,color=col)
    r_.set_xlabel("ℓ (°) around the ring",fontsize=5.8,labelpad=1); r_.set_ylabel("× $10^{-5}$",fontsize=6)
hd=[Line2D([],[],marker="o",ls="",mfc="#c0392b",mec="k",ms=5,label="added mass (size ∝ mass)"),Line2D([],[],marker="o",ls="",mfc="none",mec="#1a5276",ms=7,label="added void (to scale)"),
    Line2D([],[],marker="o",ls="",mfc="#d6eaf8",mec="#85c1e9",ms=6,label="background void (slab |z|<1.5)"),Line2D([],[],marker="o",ls="",color="#e6b0aa",ms=3,label="background cluster"),
    Line2D([],[],color="#1f4e79",lw=1.3,label="CMB ring / our sphere"),Line2D([],[],color="0.4",ls=":",label="30 Gpc from us"),Line2D([],[],color="#a93226",marker=">",lw=1.2,label="direction to the CoE"),plt.Rectangle((0,0),1,1,fc="green",alpha=0.25,label="±1× (CMB level)")]
fig.legend(handles=hd,loc="lower center",ncol=4,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.0)); fig.savefig(f"{OUT}/web_side_views{SUFX}.png",bbox_inches="tight"); plt.close(fig)
# ---- table
T=[dict(tag=d["tag"],cat=d["cat"],var=VARN[d["var"]],before=round(d["background"]["total"],2),before_sym=round(d["background"]["sym"],3),before_13up=round(d["background"]["rms_13up"],2),
        N=d["N_live"],masses=d["n_mass"],voids=d["n_void"],inside_us=d["inside_us"],core=d["in_core"],near_ring=d["near_ring_2gpc"],left=round(d["left"],3),reached=d["reached"],
        m1_5=round(d["final"]["rms_m1_5"],2),m6_12=round(d["final"]["rms_m6_12"],2),m13=round(d["final"]["rms_13up"],2),sym=round(d["final"]["sym"],3)) for d in D]
json.dump(T,open(f"{OUT}/web_add/summary_table{SUFX}.json","w"),indent=1)
for t in T: print(t)
