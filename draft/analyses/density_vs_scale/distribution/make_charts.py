import sys, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
sys.path.insert(0,"/workspace/clump_scan/inside_visible/smooth_field")
from lcdm import sigma_TH
S=json.load(open("distribution_stats.json")); P=json.load(open("../points_v2.json")); F=json.load(open("../fits_v2.json"))
rows={r["R"]:r for r in S["rows"]}; W=S["width_fit"]
RED="#d62728"; BLUE="#1f6fb2"
# ---------------- 1. histograms by size
Rs=[10,20,40,80,150]
fig,axs=plt.subplots(1,5,figsize=(17,4.6),sharey=True,gridspec_kw=dict(wspace=0.08))
for ax,R in zip(axs,Rs):
    a=np.load(f"cic_R{R}.npy").astype(float); r=rows[float(R)]
    lo,hi=min(a.min(),-0.8)*1.05,max(a.max(),0.3)*1.05
    if R==10: hi=10.5
    bins=np.linspace(-1,hi,70) if R<=20 else np.linspace(*( (-0.6,1.6) if R==40 else (-0.35,0.45) if R==80 else (-0.1,0.15)),60)
    c,e=np.histogram(a,bins); dens=c/a.size/np.diff(e)
    ax.stairs(np.where(dens>0,dens,np.nan),e,color="0.25",lw=1.6,fill=True,facecolor="0.82")
    s=r["sig_lin_presmoothed"]; sl=np.sqrt(np.log(1+s*s)); x=np.linspace(-0.999,e[-1],800)
    pdf=np.exp(-(np.log1p(x)+sl*sl/2)**2/(2*sl*sl))/((1+x)*sl*np.sqrt(2*np.pi))
    ax.plot(x,pdf,"--",color=BLUE,lw=2)
    ax.set_yscale("log"); ax.set_ylim(1e-3,2e2 if R==150 else 1e2)
    ax.axvline(0,color="0.6",lw=0.6)
    for p in P:   # published extremes of about this size
        if abs(np.log10(p["R"]/R))<0.17 and p["group"] not in ("tpp","tpp_min"):
            d=p["d"]*p["sign"]
            if e[0]<d<e[-1]: ax.plot(d,1.6e-3,"^",color="k",ms=9,zorder=5,clip_on=False)
    ax.plot(a.max(),1.6e-3,"^",color="k",mfc="white",ms=9,zorder=5); ax.plot(a.min(),1.6e-3,"^",color="k",mfc="white",ms=9,zorder=5)
    ax.set_title(f"spheres of radius {R} Mpc",fontsize=12)
    ax.text(0.97,0.95,f"spread (std) {r['std']:.2f}\n~{max(1,round(r['n_indep'])) if R<150 else '1–3'} independent",transform=ax.transAxes,ha="right",va="top",fontsize=9.5)
    ax.set_xlabel("density contrast δ",fontsize=11)
axs[0].set_ylabel("how common (fraction of volume per unit δ)",fontsize=11)
h=[Line2D([],[],color="0.25",lw=6,alpha=0.4),Line2D([],[],ls="--",color=BLUE,lw=2),
   Line2D([],[],ls="none",marker="^",mfc="white",mec="k",ms=9),Line2D([],[],ls="none",marker="^",color="k",ms=9)]
fig.legend(h,["2M++ local map (within 150 Mpc/h of us)","standard cosmology (linear σ, lognormal shape)","2M++ densest / emptiest sphere","published extreme objects of similar size"],
           loc="lower center",ncol=4,fontsize=10.5,frameon=False,bbox_to_anchor=(0.5,-0.08))
fig.suptitle("How common each density is, for regions of different size",fontsize=14,y=1.02)
fig.savefig("density_hist_by_size.png",dpi=120,bbox_inches="tight"); plt.close(fig)

# ---------------- 2. heat map
Hm=np.load("heat_2mpp.npz"); Rg=Hm["R"]; lb=Hm["lbins"]
le=np.log10(Rg); ed=np.r_[le[0]-(le[1]-le[0])/2,(le[1:]+le[:-1])/2,le[-1]+(le[-1]-le[-2])/2]; Re=10**ed
fig,(A,B)=plt.subplots(2,1,figsize=(11,10),sharex=True,gridspec_kw=dict(hspace=0.06))
norm=LogNorm(1e-4,0.3)
for ax,H in ((A,Hm["H_pos"]),(B,Hm["H_neg"])):
    Hc=np.where(H>0,H,np.nan)
    pc=ax.pcolormesh(Re,10**lb,Hc.T,cmap="YlOrBr",norm=norm,shading="flat",zorder=1)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1,20000); ax.grid(True,which="major",color="0.88",lw=0.6,zorder=0)
A.set_ylim(1e-4,3e3); B.set_ylim(1e-4,1.0); B.invert_yaxis()
A.set_ylabel("over-density  δ  (+)",fontsize=13); B.set_ylabel("under-density  δ  (−)",fontsize=13)
B.set_yticks([1e-4,1e-3,1e-2,0.1,1]); B.set_yticklabels(["−0.0001","−0.001","−0.01","−0.1","−1 (empty)"])
A.set_yticks([1e-4,1e-3,1e-2,0.1,1,10,100,1000]); A.set_yticklabels(["0.0001","0.001","0.01","0.1","1","10","100","1000"])
B.set_xlabel("region radius (Mpc)",fontsize=13); B.set_xticks([1,10,100,1000,10000]); B.set_xticklabels(["1","10","100","1,000","10,000"])
for p in P:
    ax=A if p["sign"]>0 else B; ax.plot(p["R"],p["d"],"o",ms=4.5,mfc="0.35",mec="none",alpha=0.8,zorder=4)
Rf=np.logspace(np.log10(2),np.log10(20000),200)
for ax,k in ((A,"over_noDef"),(B,"under_noDisp")):
    ax.plot(Rf,10**(F[k]["logA"]-F[k]["beta"]*np.log10(Rf)),"--",color="0.4",lw=1.2,zorder=3)
Rw=np.logspace(1,np.log10(20000),200); typ=10**(W["logA"]-W["beta"]*np.log10(Rw)); sg=np.array([sigma_TH(r/1000) for r in Rw])
for ax in (A,B):
    ax.plot(Rw,typ,"-",color=BLUE,lw=2,zorder=3)
    ax.plot(Rw,sg,":",color=BLUE,lw=2,zorder=3)
    ax.plot([1000,2000,4000],[0.068,0.055,0.059],"-*",color=RED,ms=16,lw=2,zorder=7)
A.annotate("required field\npeak ±5–7%",(2000,0.055),xytext=(2000,0.6),ha="center",fontsize=11,color=RED,weight="bold",
           arrowprops=dict(arrowstyle="-",color=RED,lw=0.8))
A.text(1.3,1.3e-4,"colour = share of the local 2M++ volume at that density\n(only measured out to 150 Mpc radius)",fontsize=9.5,color="#8a4b08")
h=[Line2D([],[],ls="none",marker="o",mfc="0.35",mec="none",ms=6),Line2D([],[],ls="--",color="0.4",lw=1.2),
   Line2D([],[],color=BLUE,lw=2),Line2D([],[],ls=":",color=BLUE,lw=2),Line2D([],[],color=RED,marker="*",ms=12,lw=2)]
L=A.legend(h,["published extremes (37)","best-fit trend of extremes","typical region (2M++ spread, extended)","typical region, standard cosmology","required field"],
           loc="upper right",fontsize=9.5,framealpha=0.95); L.set_zorder(20)
cb=fig.colorbar(pc,ax=[A,B],shrink=0.5,pad=0.02); cb.set_label("fraction of volume (per 0.1 dex of |δ|)",fontsize=10)
A.set_title("Where typical regions sit vs the extremes",fontsize=14,pad=12)
fig.savefig("density_vs_size_heatmap.png",dpi=120,bbox_inches="tight"); plt.close(fig)

# ---------------- 3. object counts
O=json.load(open("object_counts.json"))
fig,ax=plt.subplots(figsize=(9.5,6.5))
st={"clusters_redMaPPer":("o","#b2182b","galaxy clusters (redMaPPer, radius R200)"),
    "superclusters_Liivamagi2012":("s","#ef8a62","superclusters (SDSS, Liivamägi+2012)"),
    "voids_Mao2017":("D","#2166ac","voids (BOSS, Mao+2017)")}
for k,(m,c,l) in st.items():
    R=np.array(O[k]["R"]); n=np.array(O[k]["n_per_Gpc3_dex"],float); cnt=np.array(O[k]["counts"])
    g=cnt>=3; ax.plot(R[g],n[g],"-",marker=m,color=c,ms=6,lw=1.6,label=l)
for k,c,l in (("model_Tinker08_z0.2","#b2182b","model: halo mass function (Tinker+2008)"),("model_Jennings13_Vdn_z0.5","#2166ac","model: void size function (Jennings+2013)")):
    R=np.array(O[k]["R"]); n=np.array(O[k]["n_per_Gpc3_dex"],float); g=n>1e-3
    ax.plot(R[g],n[g],"--",color=c,lw=1.3,alpha=0.8,label=l)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1,1000); ax.set_ylim(1e-1,1e6)
ax.set_xlabel("object radius (Mpc)",fontsize=12); ax.set_ylabel("number per Gpc³ per factor-10 in radius",fontsize=12)
ax.set_xticks([1,10,100,1000]); ax.set_xticklabels(["1","10","100","1,000"])
ax.grid(True,which="major",color="0.88",lw=0.6); ax.legend(fontsize=9.5,loc="upper right",framealpha=0.95)
ax.text(32,0.13,"solid = observed catalogs (bins with ≥3 objects); dashed = theory.\nCluster drop below 2 Mpc = catalog cut (richness ≥20).\nTheory voids are matter-defined and smaller than\ngalaxy-defined BOSS voids, so the two do not overlap.",fontsize=8.5,color="0.35")
ax.set_title("How many objects of each size exist",fontsize=14)
fig.savefig("object_counts_vs_size.png",dpi=120,bbox_inches="tight"); plt.close(fig)
