import json,glob,os,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from ncw_core import *
plt.rcParams.update({"font.size":14,"axes.titlesize":15,"legend.fontsize":12})
KEY=[4,6,8,10,11,12]
OLD={4:35.62,6:15.87,8:4.46,10:1.26,11:0.58,12:0.72}           # hard wall, best of masses-only / masses+voids (solution_space)
OLD_ASYM={11:1.25,12:0.73}
def load(d):
    runs=[]
    for f in glob.glob(f"{d}/*.json"): runs+=json.load(open(f))
    return [r for r in runs if r["viol"]<=1e-3 and r.get("mirror_matches",0)==0 and r.get("on_plane",0)==0]
V={v:load(f"raw_{v}") for v in ("smooth","open")}; VA={v:load(f"raw_{v}_asym")+load(f"raw_{v}_asym_part1") for v in ("smooth","open")}
def best(runs):
    b={}
    for r in runs:
        N=len(r["kinds"])
        if N not in b or r["left"]<b[N]["left"]: b[N]=r
    return b
B={v:best(V[v]) for v in V}; BA={v:best(VA[v]) for v in VA if VA[v]}
fig,ax=plt.subplots(figsize=(10.5,6.8))
ax.plot(KEY,[OLD[n] for n in KEY],"o-",color="0.5",lw=2.5,ms=7,label="old: hard wall at the CoE interior")
cols={"smooth":"#1f77b4","open":"#e67e22"}; labs={"smooth":"no wall, smooth-density limit","open":"no wall, fully open"}
for v in ("smooth","open"):
    N=[n for n in KEY if n in B[v]]; ax.plot(N,[B[v][n]["left"] for n in N],"s-",color=cols[v],lw=2.5,ms=7,label=labs[v])
if "smooth" in BA:
    N=sorted(BA["smooth"]); ax.plot(N,[BA["smooth"][n]["left"] for n in N],"D:",color=cols["smooth"],lw=2,ms=7,label="smooth limit + no-mirror rule")
ax.axhline(1,color="green",lw=2); ax.text(12.25,1.08,"CMB flat\ntarget 1×",fontsize=12,color="green",va="bottom")
ax.set_yscale("log"); ax.set_ylim(0.01,60); ax.set_xlim(3.5,13.3); ax.set_xticks(KEY)
ax.set_yticks([0.01,0.1,1,10]); ax.set_yticklabels(["0.01","0.1","1","10"])
ax.set_xlabel("number of objects N (masses + deep voids)"); ax.set_ylabel("best leftover (× 10⁻⁵ floor, log scale)")
ax.set_title("Removing the hard wall at the CoE interior"); ax.legend(loc="upper right",frameon=False)
fig.tight_layout(); fig.savefig("a_best_vs_N_no_wall_logy.png",dpi=110); plt.close(fig)
# map: smooth, first N below 1x
below=[n for n in sorted(B["smooth"]) if B["smooth"][n]["left"]<1]; Nm=below[0] if below else min(B["smooth"],key=lambda n:B["smooth"][n]["left"])
m=B["smooth"][Nm]
k=np.array(list(m["kinds"])); X=np.array(m["X"]); M=np.array(m["M"]); Rv=np.array(m["Rv"]); d=np.linalg.norm(X,axis=1)
l=np.degrees(np.arctan2(X[:,1],X[:,0]))%360; b=np.degrees(np.arcsin(X[:,2]/d))
fig,(a1,a2)=plt.subplots(1,2,figsize=(15,7)); t=np.linspace(0,2*np.pi,300)
for rr,ls in ((14.3,"-"),(43,":")): a2.plot(rr*np.cos(t),rr*np.sin(t),"k",ls=ls,lw=1)
a2.plot(C_G[0]+A_IN_G*np.cos(t),C_G[1]+A_IN_G*np.sin(t),color="0.6",lw=1.2,ls="--")
ax_=C_G[0]+A_IN_G*np.cos(t); ay_=C_G[1]+A_IN_G*np.sin(t); ok_=(np.abs(ax_)<44)&(np.abs(ay_)<44); j_=np.argmax(np.where(ok_,ay_,-1e9))
a2.text(ax_[j_]+1.5,ay_[j_]-4,"old wall",color="0.4",fontsize=12)
ax_l=np.radians(15.16); a2.plot([-45*np.cos(ax_l),45*np.cos(ax_l)],[-45*np.sin(ax_l),45*np.sin(ax_l)],color="0.8",ls="--",lw=1)
for i in range(len(k)):
    if k[i]=="C":
        s=60+240*M[i]/2e52; a1.scatter(l[i],b[i],s=s,c="#c0392b"); a2.scatter(X[i,0],X[i,1],s=s,c="#c0392b",zorder=3)
    else:
        a1.scatter(l[i],b[i],s=60+240*Rv[i]/4,facecolors="none",edgecolors="#2471a3",lw=2.5)
        a2.add_patch(plt.Circle(X[i,:2],Rv[i],color="#2471a3",alpha=0.35))
a2.plot(0,0,"k+",ms=14); a2.set_aspect("equal"); a2.set_xlim(-45,45); a2.set_ylim(-45,45); a2.set_xlabel("Gpc")
a2.set_title("Top view (us = +, grey dashed arc = old wall,\nlight dashed line = axis ℓ 15°/195°)",fontsize=13)
a1.set_xlim(360,0); a1.set_ylim(-60,60); a1.set_xticks(range(0,361,45)); a1.axhline(0,color="0.8",lw=1)
a1.axvspan(270,300,color="0.9",zorder=-1); a1.text(285,55,"Migkas",ha="center",fontsize=11)
a1.set_xlabel("galactic longitude ℓ (deg)"); a1.set_ylabel("latitude b (deg)"); a1.set_title("Sky positions (masses red dots, voids blue rings)",fontsize=13)
fig.suptitle(f"No wall, smooth-density limit: {int((k=='C').sum())} masses + {int((k=='V').sum())} deep voids, leftover {m['left']:.2f}×"); fig.tight_layout()
fig.savefig("d_no_wall_smooth_map.png",dpi=105); plt.close(fig)
# printouts
def describe(tag,r):
    k=np.array(list(r["kinds"])); X=np.array(r["X"]); M=np.array(r["M"]); Rv=np.array(r["Rv"]); d=np.linalg.norm(X,axis=1); dep=core_depth(X)
    rho=rho_gpc(X); isC=k=="C"
    print(f"\n{tag} N={len(k)} {r['kinds']} left={r['left']:.3f}  clusters {M[isC].sum():.2e} voids {M[~isC].sum():.2e} net {M.sum():.1e}")
    print("  masses inside old wall: %d/%d (max depth %.1f Gpc); voids inside: %d/%d"%((dep[isC]>0).sum(),isC.sum(),dep[isC].max(),(dep[~isC]>0).sum(),(~isC).sum()))
    capC=(isC&(M>=0.98*MHI*rho/RHO)).sum(); print("  bounds: mass at cap",capC,"void Rv=4",int(((~isC)&(Rv>=3.98)).sum()),"d>=42.9",int((d>=42.9).sum()))
    for i in np.argsort(np.degrees(np.arctan2(X[:,1],X[:,0]))%360):
        print("   %s d=%.1f l=%.0f b=%+.0f M=%.2e Rv=%.2f depth_inside_wall=%+.1f"%(k[i],d[i],np.degrees(np.arctan2(X[i,1],X[i,0]))%360,np.degrees(np.arcsin(X[i,2]/d[i])),M[i],Rv[i],dep[i]))
for v in ("smooth","open"):
    for n in KEY:
        if n in B[v]:
            L=sorted(r["left"] for r in V[v] if len(r["kinds"])==n); print(v,n,B[v][n]["kinds"],"best %.3f"%L[0],"n=%d"%len(L),"below1x=%d"%sum(x<1 for x in L))
for v in BA:
    for n in sorted(BA[v]): print("asym",v,n,BA[v][n]["kinds"],round(BA[v][n]["left"],3),"odd",round(BA[v][n].get("odd_frac",0),4))
for v in ("smooth","open"):
    for n in (4,8,Nm): 
        if n in B[v]: describe(v,B[v][n])
# depth stats over all good runs (<2x)
for v in ("smooth","open"):
    deps=[np.array(r["core_depth"])[np.array(list(r["kinds"]))=="C"] for r in V[v] if r["left"]<2]
    if deps:
        a=np.concatenate(deps); print(v,"masses in fits <2x:",len(a),"inside wall %.0f%%"%(100*(a>0).mean()),"median depth of those %.1f, max %.1f"%(np.median(a[a>0]) if (a>0).any() else 0,a.max()))
print("MAP N",Nm)
json.dump({v:{n:{kk:vv for kk,vv in B[v][n].items() if kk!="p"} for n in B[v]} for v in B},open("ncw_best.json","w"))
