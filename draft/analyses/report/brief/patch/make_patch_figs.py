"""Brief figures for the patch-ratio result. Data read-only from ../../../realistic_counts/patch_ratio (results.json,
thr_in10.npz/json); chart logic adapted from its plots.py / map.py at paper size, legends outside the plots."""
import os,sys,json,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
HERE=os.path.dirname(os.path.abspath(__file__)); FIG=os.path.join(HERE,"..","figures")
PR=os.path.abspath(os.path.join(HERE,"..","..","..","realistic_counts","patch_ratio"))
os.chdir(PR); sys.path.insert(0,PR); from ps_core import pinfo,CONFIGS
plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"savefig.dpi":300,"pdf.fonttype":42})
def save(fig,n): fig.savefig(f"{FIG}/{n}.pdf",bbox_inches="tight"); fig.savefig(f"{FIG}/{n}.png",bbox_inches="tight",dpi=200); plt.close(fig)
R=json.load(open("results.json")); T=R["thresholds"]
def fmapB(B,p):
    df=B/p["Delta"]; f0=p["f0"]; return 100*(((f0+df)/(1-f0-df))/(f0/(1-f0))-1)
fig,axs=plt.subplots(2,1,figsize=(7.1,9.2),gridspec_kw=dict(hspace=0.62))
fig.subplots_adjust(left=0.1,right=0.97,top=0.97,bottom=0.13)
ax=axs[0]
cv=[("in05","extreme 50 Mpc","inside, 0.5 Gpc domains,\n50 Mpc patches (+0.83 / −0.39)","#c0392b","-"),
    ("in05","typical ±0.3, 50 Mpc","inside, 0.5 Gpc domains,\n50 Mpc patches (±0.3)","#e67e22","--"),
    ("in20","extreme 50 Mpc","inside, 2 Gpc domains,\n50 Mpc patches (+0.83 / −0.39)","#2471a3","-"),
    ("both10","extreme 50 Mpc","inside + outside, no mirror, 1 Gpc,\n50 Mpc patches (+0.83 / −0.39)","#7d3c98","-")]
thr={}
for c,pn,lab,col,ls in cv:
    p=pinfo(pn); B=[x[0] for x in T[c]["curve"]]; L=[max(x[1],1e-3) for x in T[c]["curve"]]
    if T[c]["B_1x"]:
        k=[i for i,b in enumerate(B) if b<T[c]["B_lower"]+1e-9]; B=[B[i] for i in k]+[T[c]["B_1x"]]; L=[L[i] for i in k]+[1.0]
    x=[fmapB(b,p) for b in B]; ax.plot(x,L,ls,color=col,lw=1.6,marker="o",ms=2.5,label=lab)
    if T[c]["B_1x"]: xs=fmapB(T[c]["B_1x"],p); thr[(c,pn)]=xs; ax.plot(xs,1,"o",ms=7,mfc="none",mec=col,mew=1.4)
ax.axhline(1,color="green",lw=1.6); ax.text(1.1,1.2,"1× = CMB flat",color="green",fontsize=7.5)
ax.set_yscale("log"); ax.set_xscale("log"); ax.set_ylim(5e-1,400); ax.set_xlim(1,100)
ax.set_xlabel("allowed shift of the over : under patch ratio from the cosmic mean (%)"); ax.set_ylabel("best CMB leftover (× $10^{-5}$ floor)")
ax.set_title("(a) leftover vs allowed ratio shift"); ax.grid(alpha=0.25,which="both")
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.15),ncol=2,frameon=False,fontsize=7,title="open circle: first shift that reaches 1×",title_fontsize=7)
ax=axs[1]; tab=R["table"]
for pn,col,lab in [("extreme 30 Mpc","#922b21","30 Mpc patches (+2.1 / −0.55)"),("extreme 50 Mpc","#c0392b","50 Mpc patches (+0.83 / −0.39)"),
                   ("extreme 100 Mpc","#e59866","100 Mpc patches (+0.24 / −0.25)"),("typical ±0.3, 50 Mpc","#2471a3","50 Mpc patches, typical ±0.3")]:
    rr=sorted([t for t in tab if t["patch"]==pn and t["cfg"].startswith("in")],key=lambda t:t["H"])
    ax.plot([t["H"] for t in rr],[t["n_sigma_domain"] for t in rr],"-o",color=col,lw=1.6,ms=4,label=lab)
rr=sorted([t for t in tab if t["patch"]=="extreme 50 Mpc" and t["cfg"].startswith("in")],key=lambda t:t["H"])
ax.plot([t["H"] for t in rr],[t["n_sigma_LCDM"] for t in rr],"--s",color="0.3",lw=1.4,ms=4,label="ΛCDM scatter at domain size")
ax.axhline(1,color="green",lw=1.1); ax.text(2.3,1.06,"1σ",color="green",ha="right"); ax.axhline(3,color="green",lw=1.1,ls=":"); ax.text(2.3,3.15,"3σ",color="green",ha="right")
ax.set_xscale("log",base=2); ax.set_yscale("log"); ax.set_xticks([0.5,1,2]); ax.set_xticklabels(["0.5","1","2"]); ax.set_xlim(0.45,2.4)
ax.set_xlabel("domain size (Gpc)"); ax.set_ylabel("needed ratio shift per domain ÷ random scatter (σ)")
ax.set_title("(b) one domain at a time vs randomness"); ax.grid(alpha=0.25,which="both")
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.15),ncol=3,frameon=False,fontsize=7,title="inside, mirror allowed; one domain only",title_fontsize=7)
save(fig,"fig_patch_ratio")
print("thresholds % ratio shift up:",{k:round(v,1) for k,v in thr.items()})
for t in tab:
    if t["patch"]=="extreme 50 Mpc": print(t["cfg"],t["H"],round(t["ratio_up_pct"],1),round(t["ratio_dn_pct"],1),round(t["n_sigma_domain"],2),round(t.get("n_sigma_LCDM",np.nan),2))
# ---- map (side view = ring plane), inside, 1 Gpc domains, 50 Mpc extreme patches
g=CONFIGS["in10"]()[0]; z=np.load("thr_in10.npz"); d=z[z.files[0]]; P=g.P
p=pinfo("extreme 50 Mpc"); f0=p["f0"]; f=f0+d/p["Delta"]; rel=(f/(1-f))/(f0/(1-f0))
lob=json.load(open("pattern_in10.json"))["lobes"]
s=(P[:,2]>0)&(P[:,2]<0.6); v=100*(rel-1); lim=np.ceil(np.abs(v[s]).max())
fig=plt.figure(figsize=(3.4,3.9)); ax=fig.add_axes([0.13,0.25,0.8,0.7]); cax=fig.add_axes([0.2,0.1,0.66,0.025])
sc=ax.scatter(P[s,0],P[s,1],c=v[s],cmap="RdBu_r",vmin=-lim,vmax=lim,s=25,marker="s",linewidths=0,rasterized=True)
ax.add_patch(Circle((0,0),14.3,fill=False,ls="--",color="k",lw=0.9)); ax.plot(0,0,"k+",ms=7,mew=1.2)
for l in lob:
    a=np.radians(l["lon"]); ax.text(17.2*np.cos(a),17.2*np.sin(a),f"{l['rel_ratio']:.2f}:1",ha="center",va="center",fontsize=6.5,color="#922b21" if l["sign"]>0 else "#1a5276")
ax.set_xlim(-20,20); ax.set_ylim(-20,20); ax.set_aspect("equal"); ax.set_xlabel("Gpc",labelpad=1); ax.set_ylabel("Gpc",labelpad=1); ax.tick_params(labelsize=6.5)
cb=fig.colorbar(sc,cax=cax,orientation="horizontal"); cb.set_label("over : under patch ratio vs cosmic mean (%)",fontsize=6.8); cb.ax.tick_params(labelsize=6.5)
save(fig,"fig_patch_map")
print("lobes",[(l["lon"],round(l["rel_ratio"],3)) for l in lob])
