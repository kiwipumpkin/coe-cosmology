import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":14,"axes.titlesize":15,"legend.fontsize":12})
res=json.load(open("results.json")); MIX=["masses only","masses + deep voids","masses + shallow voids","deep voids only","shallow voids only"]
COL={"masses only":"#c0392b","masses + deep voids":"#8e44ad","masses + shallow voids":"#e67e22","deep voids only":"#2471a3","shallow voids only":"#5dade2"}
def best_by(mix,key=lambda r:r["best"],agg=min):
    out={}
    for r in res.values():
        if r["mix"]!=mix: continue
        out.setdefault(r["N"],[]).append(r)
    return out
# (a) best leftover vs N
fig,ax=plt.subplots(figsize=(10,7.5))
for mix in MIX:
    g=best_by(mix); Ns=sorted(g); y=np.minimum.accumulate([min(r["best"] for r in g[n]) for n in Ns])
    ax.plot(Ns,y,"o-",color=COL[mix],lw=2.5,ms=7,label=mix)
for v,t in [(237.3,"start 237×"),(88.8,"four clumps 89×"),(27.4,"24-pair greedy fit 27×")]:
    ax.axhline(v,color="0.5",ls="--",lw=1); ax.text(8.1,v,t,va="center",fontsize=11,color="0.35")
ax.set_xlabel("number of objects N"); ax.set_ylabel("best leftover (× floor), best split"); ax.set_xlim(0.8,9.6); ax.set_ylim(0,250)
ax.set_xticks(range(1,9)); ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.12),ncol=3,frameon=False); ax.set_title("Best leftover vs number of fitted objects")
fig.tight_layout(); fig.savefig("a_best_vs_N.png",dpi=110); plt.close(fig)
# (b) distinct good solutions vs N (below 50x), best split per N for mixed
for thr in ("50","27"):
    fig,ax=plt.subplots(figsize=(10,6.5))
    for mix in MIX:
        g=best_by(mix); Ns=sorted(g); y=[max(r["distinct_below"][thr] for r in g[n]) for n in Ns]
        ax.plot(Ns,y,"o-",color=COL[mix],lw=2.5,ms=7,label=mix)
    ax.set_xlabel("number of objects N"); ax.set_ylabel(f"distinct solutions below {thr}× (best split)"); ax.set_xticks(range(1,9))
    ax.set_title(f"How many different arrangements reach below {thr}× (first 40 starts per class)"); ax.legend(loc="upper left")
    fig.tight_layout(); fig.savefig(f"b_distinct_below{thr}.png",dpi=110); plt.close(fig)
# (c) top-down + sky map of all good (<50x) basins: masses red, voids blue
fig,(a1,a2)=plt.subplots(1,2,figsize=(15,7))
for r in res.values():
    for b in r["basins"]:
        if b["left"]>=50: continue
        X=np.array(b["X"]).reshape(-1,3); k=np.array(list(b["kinds"]))
        if not len(X): continue
        l=np.degrees(np.arctan2(X[:,1],X[:,0]))%360; d=np.linalg.norm(X,axis=1); bb=np.degrees(np.arcsin(X[:,2]/d))
        for t,c in (("C","#c0392b"),("V","#2471a3"),("S","#5dade2")):
            s=k==t
            a1.scatter(l[s],bb[s],c=c,s=18,alpha=0.35,lw=0)
            a2.scatter(X[s,0],X[s,1],c=c,s=18,alpha=0.35,lw=0)
a1.set_xlim(360,0); a1.set_ylim(0,90); a1.set_xticks(range(0,361,45)); a1.set_xlabel("galactic longitude ℓ (deg)"); a1.set_ylabel("|latitude| (deg)")
a1.axvspan(270,300,color="0.85",zorder=-1); a1.text(285,85,"Migkas\nℓ 270–300",ha="center",va="top",fontsize=11)
a1.set_title("Sky position (masses red, voids blue)")
t=np.linspace(0,2*np.pi,300)
for rr,ls in ((14.3,"-"),(43,":")): a2.plot(rr*np.cos(t),rr*np.sin(t),"k",ls=ls,lw=1)
from ss_core import C_G,A_IN_G
a2.plot(C_G[0]+A_IN_G*np.cos(t),C_G[1]+A_IN_G*np.sin(t),color="0.5",lw=1); a2.plot(*C_G[:2],"kx")
a2.plot(0,0,"k+",ms=12); a2.set_aspect("equal"); a2.set_xlim(-45,45); a2.set_ylim(-45,45); a2.set_xlabel("Gpc"); a2.set_title("Top view (we are at +; CoE interior = grey arc)")
fig.suptitle("All distinct solutions below 50×, every class overlaid"); fig.tight_layout(); fig.savefig("c_where_objects_go.png",dpi=105); plt.close(fig)
