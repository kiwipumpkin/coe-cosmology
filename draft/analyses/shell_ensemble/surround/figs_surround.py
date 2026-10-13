import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D; from matplotlib.patches import Patch
plt.rcParams.update({"font.size":15})
c=json.load(open("counts.json")); C=c["counts"]; O=c["original_counts_recomputed"]
groups=[(N,lev) for lev in "abc" for N in (200,500)]
LV={"a":"tiny\n10⁴⁸–10⁴⁹ kg","b":"medium\n3×10⁴⁹–3×10⁵⁰ kg","c":"large\n10⁵¹–10⁵² kg"}
x=np.arange(len(groups)); w=0.36
help_=[C[f"grad_bare_N{N}_{l}"]["placed"]["rms_lower"] for N,l in groups]; hurt=[1000-h for h in help_]
oh=[O[f"bare_N{N}_{l}"]["placed"]["rms_lower"] for N,l in groups]; ohu=[1000-h for h in oh]
fig,ax=plt.subplots(figsize=(15,8))
b1=ax.bar(x-w/2,help_,w,color="#2a9d4a",label="helped: lowered the leftover CMB pattern")
b2=ax.bar(x+w/2,hurt,w,color="#d1495b",label="hurt: raised the leftover CMB pattern")
ax.bar(x-w/2,oh,w,fill=False,edgecolor="k",ls="--",lw=1.4); ax.bar(x+w/2,ohu,w,fill=False,edgecolor="k",ls="--",lw=1.4)
for xi,v in zip(x-w/2,help_): ax.text(xi,v+12,str(v),ha="center",fontsize=13,color="#1d6e34",weight="bold")
for xi,v in zip(x+w/2,hurt): ax.text(xi,v+12,str(v),ha="center",fontsize=13,color="#9b2335",weight="bold")
ax.set_xticks(x); ax.set_xticklabels([f"{N} objects\n{LV[l]}" for N,l in groups],fontsize=12.5)
for k in (1.5,3.5): ax.axvline(k,color="0.7",lw=1)
ax.axhline(500,color="0.5",ls=":",lw=1.2); ax.annotate("half\n(500)",xy=(1.005,500),xycoords=("axes fraction","data"),va="center",fontsize=11,color="0.4",annotation_clip=False)
ax.set_ylim(0,1100); ax.set_ylabel("number of random layouts (out of 1,000)")
ax.set_xlabel("number of clusters + voids per layout, and mass of each object")
ax.set_title("Random clusters and voids placed all around our visible sphere (incl. beyond D+d):\nhow many of 1,000 layouts helped vs hurt the leftover CMB pattern (237× the 10⁻⁵ target)",fontsize=16)
h=[b1,b2,Patch(fill=False,edgecolor="k",ls="--",label="dashed outline: earlier shell-only layout, for comparison")]
ax.legend(handles=h,loc="upper left",fontsize=12.5,framealpha=.95); ax.grid(axis="y",alpha=.3)
fig.tight_layout(); fig.savefig("helped_vs_hurt.png",dpi=110); plt.close(fig)

# example realization
o=json.load(open("ensemble.json")); key="grad_bare_N500_b"; ex=o["examples"][key]
P=np.array(ex["P"]); M=np.array(ex["M"]); B=np.array(ex["beyond"]); ell=np.radians(o["ell_coe"]); D=58.2
Cv=D*np.array([np.cos(ell),np.sin(ell)]); t=np.linspace(0,2*np.pi,400); cl=M>0
sz=12+60*np.abs(M)/np.abs(M).max()
CATS=((cl&~B,"o","red","cluster, in shell D−d..D+d"),(~cl&~B,"o","blue","void, in shell D−d..D+d"),
      (cl&B,"^","darkorange","cluster, beyond D+d"),(~cl&B,"^","deepskyblue","void, beyond D+d"))
fig=plt.figure(figsize=(26,9.5))
axs=[fig.add_subplot(1,3,1),fig.add_subplot(1,3,2)]
for a,lim in zip(axs,(None,48)):
    for r,ls,col,lab in ((43.9,"--","k","shell inner edge D−d (43.9 Gpc from CoE)"),(72.5,"-","k","shell outer edge D+d (72.5 Gpc)"),(101.1,":","purple","D+3d (101.1 Gpc), farthest point used")):
        a.plot(Cv[0]+r*np.cos(t),Cv[1]+r*np.sin(t),color=col,ls=ls,lw=1.6,label=lab)
    a.fill(14.3*np.cos(t),14.3*np.sin(t),color="C0",alpha=.25,label="our visible sphere (14.3 Gpc)"); a.plot(42.9*np.cos(t),42.9*np.sin(t),color="green",ls="-.",lw=1.4,label="placement region: up to 42.9 Gpc from us")
    keep=np.ones(len(M),bool) if lim is None else np.abs(P[:,2])<7
    for sel,mk,col,lab in CATS:
        sel=sel&keep; a.scatter(P[sel,0],P[sel,1],s=sz[sel]*(2.2 if lim else 1),marker=mk,c=col,alpha=.7,edgecolor="k",lw=.4,label=lab)
    a.scatter(*Cv,marker="x",s=250,c="k",lw=3,label="CoE (58.2 Gpc away)"); a.scatter(0,0,marker="*",s=450,c="gold",edgecolor="k",zorder=6,label="us (observer)")
    a.set_aspect("equal"); a.grid(alpha=.3); a.set_xlabel("x (Gpc)"); a.set_ylabel("y (Gpc)")
    if lim: a.set_xlim(-lim,lim); a.set_ylim(-lim,lim); a.set_title("Zoom on our visible sphere: thin slice\n(only objects within ±7 Gpc of this plane)",fontsize=15)
    else: a.set_title("Whole layout seen from above\n(all objects projected flat, so some appear to overlap the sphere)",fontsize=15); a.legend(loc="lower left",fontsize=10,framealpha=.9,ncol=2)
b=fig.add_subplot(1,3,3,projection="mollweide")
lon=np.arctan2(P[:,1],P[:,0]); lat=np.arcsin(P[:,2]/np.linalg.norm(P,axis=1))
lon=-((lon-ell+np.pi)%(2*np.pi)-np.pi)      # CoE direction at the centre of the map
for sel,mk,col,lab in CATS: b.scatter(lon[sel],lat[sel],s=sz[sel]*.8,marker=mk,c=col,alpha=.7,edgecolor="k",lw=.3)
b.scatter(0,0,marker="x",s=250,c="k",lw=3); b.text(0,np.radians(8),"toward CoE",ha="center",fontsize=13)
b.text(np.radians(150),np.radians(-12),"away from\nCoE",ha="center",fontsize=12)
b.set_xticklabels([]); b.grid(alpha=.4)
b.set_title("The sky as seen from us: the direction of every object\n(empty patch toward the CoE = CoE interior, inside D−d)",fontsize=15)
fig.suptitle("Example random layout: 500 clusters + voids (medium mass, 3×10⁴⁹–3×10⁵⁰ kg each) placed all around our visible sphere; circles = in the shell, triangles = beyond D+d",fontsize=17)
fig.tight_layout(rect=(0,0,1,.94)); fig.savefig("example_realization.png",dpi=95); plt.close(fig)
print("ok")
