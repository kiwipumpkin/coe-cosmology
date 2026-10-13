import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":14,"axes.titlesize":15,"legend.fontsize":12})
res=json.load(open("results.json")); em=json.load(open("extend_mixed.json")); ec=json.load(open("extend_masses.json"))
MIX=["masses only","masses + deep voids","masses + shallow voids","deep voids only","shallow voids only"]
COL={"masses only":"#c0392b","masses + deep voids":"#8e44ad","masses + shallow voids":"#e67e22","deep voids only":"#2471a3","shallow voids only":"#5dade2"}
fig,ax=plt.subplots(figsize=(11,7.5))
for mix in MIX:
    g={}
    for r in res.values():
        if r["mix"]==mix: g.setdefault(r["N"],[]).append(r["best"])
    Ns=sorted(g); y=list(np.minimum.accumulate([min(g[n]) for n in Ns]))
    if mix=="masses + deep voids": Ns+= [s["N"] for s in em[1:]]; y+=[s["left"] for s in em[1:]]
    if mix=="masses only": Ns+= [s["N"] for s in ec[1:]]; y+=[s["left"] for s in ec[1:]]
    ax.plot(Ns,y,"o-",color=COL[mix],lw=2.5,ms=6,label=mix)
for v,t in [(237.3,"start 237×"),(88.8,"four clumps 89×"),(27.4,"27×")]:
    ax.axhline(v,color="0.5",ls="--",lw=1); ax.text(20.3,v,t,va="center",fontsize=11,color="0.35")
ax.axhline(1,color="green",lw=2); ax.text(20.3,1,"CMB flat\ntarget 1×",va="center",fontsize=11,color="green")
m=em[-1]; ax.annotate(f"{m['left']:.2f}× at N={m['N']}",(m["N"],m["left"]),xytext=(m["N"]+1.2,0.4),fontsize=12,color=COL["masses + deep voids"],arrowprops=dict(arrowstyle="->",color=COL["masses + deep voids"]))
c=ec[-1]; ax.annotate(f"{c['left']:.1f}× at N=20",(20,c["left"]),xytext=(16.3,12),fontsize=12,color=COL["masses only"],arrowprops=dict(arrowstyle="->",color=COL["masses only"]))
ax.set_yscale("log"); ax.set_ylim(0.3,350); ax.set_xlim(0.5,23.5); ax.set_xticks(range(1,21))
ax.set_yticks([0.5,1,2,5,10,20,50,100,200]); ax.set_yticklabels(["0.5","1","2","5","10","20","50","100","200"])
ax.set_xlabel("number of objects N"); ax.set_ylabel("best leftover (× floor, log scale)")
ax.set_title("Best leftover vs number of objects (N>8: add one object at a time, then refit all)")
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.11),ncol=3,frameon=False)
fig.tight_layout(); fig.savefig("a_best_vs_N_logy_extended.png",dpi=110); plt.close(fig)
# floor map
import sys; sys.path.insert(0,"."); from ss_core import C_G,A_IN_G
k=np.array(list(m["kinds"])); X=np.array(m["X"]); M=np.array(m["M"]); Rv=np.array(m["Rv"]); d=np.linalg.norm(X,axis=1)
l=np.degrees(np.arctan2(X[:,1],X[:,0]))%360; b=np.degrees(np.arcsin(X[:,2]/d))
fig,(a1,a2)=plt.subplots(1,2,figsize=(15,7)); t=np.linspace(0,2*np.pi,300)
for rr,ls in ((14.3,"-"),(43,":")): a2.plot(rr*np.cos(t),rr*np.sin(t),"k",ls=ls,lw=1)
a2.plot(C_G[0]+A_IN_G*np.cos(t),C_G[1]+A_IN_G*np.sin(t),color="0.5",lw=1)
ax_l=np.radians(15.16); a2.plot([-45*np.cos(ax_l),45*np.cos(ax_l)],[-45*np.sin(ax_l),45*np.sin(ax_l)],color="0.7",ls="--",lw=1)
for i in range(len(k)):
    if k[i]=="C":
        s=60+240*M[i]/2e52; a1.scatter(l[i],b[i],s=s,c="#c0392b"); a2.scatter(X[i,0],X[i,1],s=s,c="#c0392b",zorder=3)
    else:
        a1.scatter(l[i],b[i],s=60+240*Rv[i]/4,facecolors="none",edgecolors="#2471a3",lw=2.5)
        a2.add_patch(plt.Circle(X[i,:2],Rv[i],color="#2471a3",alpha=0.35))
a2.plot(0,0,"k+",ms=14); a2.set_aspect("equal"); a2.set_xlim(-45,45); a2.set_ylim(-45,45); a2.set_xlabel("Gpc")
a2.set_title("Top view (us = +, dashed = mirror axis ℓ 15°/195°,\ngrey arc = CoE interior edge)",fontsize=13)
a1.set_xlim(360,0); a1.set_ylim(-60,60); a1.set_xticks(range(0,361,45)); a1.axhline(0,color="0.8",lw=1)
a1.axvline(15.16,color="0.7",ls="--",lw=1); a1.axvline(195.16,color="0.7",ls="--",lw=1)
a1.axvspan(270,300,color="0.9",zorder=-1); a1.text(285,55,"Migkas",ha="center",fontsize=11)
a1.set_xlabel("galactic longitude ℓ (deg)"); a1.set_ylabel("latitude b (deg)"); a1.set_title("Sky positions (masses red dots, voids blue rings)",fontsize=13)
fig.suptitle(f"Floor-reaching layout: {int((k=='C').sum())} masses + {int((k=='V').sum())} deep voids, leftover {m['left']:.2f}×"); fig.tight_layout()
fig.savefig("d_floor_solution_map.png",dpi=105); plt.close(fig)
