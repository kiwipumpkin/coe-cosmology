import sys, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0,"/workspace/clump_scan/inside_visible/smooth_field")
from lcdm import sigma_TH
F=json.load(open("fits_v2.json")); P=json.load(open("points_v2.json"))
R=np.logspace(0,np.log10(20000),200); sig=np.array([sigma_TH(r/1000) for r in R])
def band(f,Rs):
    x=np.log10(Rs); y=f["logA"]-f["beta"]*x; c=np.array(f["cov"])
    sy=np.sqrt(c[0][0]+x*x*c[1][1]+2*x*c[0][1]); return 10**y,10**(y-sy),10**(y+sy)
fig,(a,b)=plt.subplots(2,1,figsize=(11,10.5),sharex=True,gridspec_kw=dict(hspace=0.06))
for ax in (a,b):
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1,20000)
    ax.fill_between(R,1e-7,3*sig,color="#cfe3f3",zorder=0)
    ax.axvspan(1,10,color="0.93",zorder=0.5)
    ax.grid(True,which="major",color="0.88",lw=0.6)
a.set_ylim(1e-4,3e3); b.set_ylim(1e-4,1.0); b.invert_yaxis()
a.set_ylabel("over-density  δ  (+)",fontsize=13); b.set_ylabel("under-density  δ  (−)",fontsize=13)
b.set_yticks([1e-4,1e-3,1e-2,0.1,1]); b.set_yticklabels(["−0.0001","−0.001","−0.01","−0.1","−1 (empty)"])
a.set_yticks([1e-4,1e-3,1e-2,0.1,1,10,100,1000]); a.set_yticklabels(["0.0001","0.001","0.01","0.1","1","10","100","1000"])
b.set_xlabel("region radius (Mpc)",fontsize=13)
b.set_xticks([1,10,100,1000,10000]); b.set_xticklabels(["1","10","100","1,000","10,000"])
Rf=np.logspace(np.log10(2),np.log10(20000),200)
def fl(ax,key,lab,ls="--",shade=True,col="k",lw=2):
    m,lo,hi=band(F[key],Rf)
    if shade: ax.fill_between(Rf,lo,hi,color="0.6",alpha=0.22,lw=0,zorder=1)
    ax.plot(Rf,m,ls,color=col,lw=lw,zorder=2,label=lab)
fo=F["over_noDef"]; fu=F["under_noDisp"]; f2=F["under_2mpp"]
fl(a,"over_noDef",f"Best fit, all measured peaks: δ ∝ R$^{{-{fo['beta']:.1f}}}$")
fl(b,"under_noDisp",f"Best fit, all voids (no disputed): |δ| ∝ R$^{{-{fu['beta']:.2f}}}$")
fl(b,"under_2mpp",f"2M++ map alone: |δ| ∝ R$^{{-{f2['beta']:.1f}}}$",ls=":",shade=False,col="#1f6fb2",lw=2)
sty={"def":dict(marker="^",mfc="white",mec="k",ms=8,label="Cluster edge R200 (by definition)"),
     "mass":dict(marker="o",mfc="k",mec="k",ms=7,label="Mass estimates (clusters, superclusters)"),
     "gal":dict(marker="o",mfc="0.6",mec="0.3",ms=7,label="Galaxy-count over-densities"),
     "tpp":dict(marker="s",mfc="#1f6fb2",mec="#1f6fb2",ms=6,label="2M++ public map: densest sphere"),
     "void_gal":dict(marker="o",mfc="k",mec="k",ms=7,label="Voids (galaxy counts)"),
     "void_mass":dict(marker="D",mfc="k",mec="k",ms=6,label="Voids (lensing / matter)"),
     "disputed":dict(marker="o",mfc="white",mec="k",ms=7,label="Local Hole claims (disputed)"),
     "tpp_min":dict(marker="s",mfc="#1f6fb2",mec="#1f6fb2",ms=6,label="2M++ public map: emptiest sphere")}
done=set()
for p in P:
    ax=a if p["sign"]>0 else b; s=dict(sty[p["group"]]); lab=s.pop("label")
    ax.plot(p["R"],p["d"],ls="none",zorder=5,**s); done.add(p["group"])
kw=dict(fontsize=10.5,zorder=6,va="center")
def lab(ax,R,d,t,dx=1.25,dy=1.0,ha="left"): ax.annotate(t,(R,d),xytext=(R*dx,d*dy),ha=ha,**kw)
lab(a,2.84,645,"Coma",dx=1.3); lab(a,7.2,11.5,"Virgo",dx=0.8,ha="right"); lab(a,11.4,20,"Shapley core",dx=0.8,dy=2.2,ha="right")
lab(a,36.6,1.62,"Saraswati",dx=0.8,dy=0.55,ha="right"); lab(a,79.4,1.0,"Sloan Great Wall",dx=1.25,dy=1.3)
lab(a,36.1,2.74,"BOSS Great Wall",dx=1.25,dy=1.9)
lab(b,23,0.8,"Local Void",dx=0.8,dy=0.85,ha="right"); lab(b,44,0.75,"Boötes",dx=1.25,dy=0.8)
lab(b,314,0.14,"Eridanus",dx=1.25,dy=0.9); lab(b,300,0.33,"Local Hole claims",dx=1.3,dy=1.15)
gx=np.array([1000,2000,4000]); pk=[0.068,0.055,0.059]; rm=[0.0105,0.0124,0.0141]
for ax in (a,b):
    ax.plot(gx,pk,"-*",color="#d62728",ms=16,lw=2,zorder=7)
    ax.plot(gx,rm,"--o",color="#d62728",mfc="white",ms=8,lw=1.5,zorder=7)
a.annotate("Required field\npeak ±5–7%",(1000,0.068),xytext=(700,0.3),ha="center",fontsize=11,color="#d62728",weight="bold")
a.annotate("rms ≈1%",(1000,0.0105),xytext=(800,0.006),ha="right",va="center",fontsize=10.5,color="#d62728")
a.plot(830,100,"s",color="#7b3fa0",ms=11,zorder=7)
a.annotate("Point-mass fit cores",(830,100),xytext=(1120,100),va="center",fontsize=10.5,color="#7b3fa0",weight="bold")
b.plot(3750,0.7,"s",color="#7b3fa0",ms=11,zorder=7)
b.annotate("Point-mass fit voids",(3750,0.7),xytext=(2800,0.7),ha="right",va="center",fontsize=10.5,color="#7b3fa0",weight="bold")
from matplotlib.lines import Line2D
def H(g): s=dict(sty[g]); l=s.pop("label"); return Line2D([],[],ls="none",**s),l
for ax,gs in ((a,["def","mass","gal","tpp"]),(b,["tpp_min","void_gal","void_mass","disputed"])):
    hh,ll=ax.get_legend_handles_labels(); hh=[h for h,l in zip(hh,ll) if l.startswith(("Best","2M++ map"))]; ll=[l for l in ll if l.startswith(("Best","2M++ map"))]
    for g in gs: h_,l_=H(g); hh.append(h_); ll.append(l_)
    L=ax.legend(hh,ll,loc="lower left" if ax is a else "upper right",fontsize=9,framealpha=0.95); L.set_zorder(20)
a.text(1.3,0.6,"linear theory\nunder-predicts here",fontsize=9,color="0.45")
b.text(1.3,0.0002,"pale blue = standard cosmology\n3σ range for that size",fontsize=9,color="#2c6ea3",va="top")
a.set_title("Observed extremes with best-fit trends (37 published / public-data points)",fontsize=14,pad=12)
fig.savefig("density_vs_scale_fits_radius_v2.png",dpi=130,bbox_inches="tight")
