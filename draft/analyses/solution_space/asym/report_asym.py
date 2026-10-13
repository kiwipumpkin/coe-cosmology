import json,glob,os,sys,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from asym_core import *
plt.rcParams.update({"font.size":14,"axes.titlesize":15,"legend.fontsize":12})
P=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OLD={8:(4.46,"4+4"),9:(2.96,"5+4"),10:(1.26,"6+4"),11:(0.58,"6+5"),12:(0.72,"6+6")}
runs=[]
for f in glob.glob("raw/*.json"): runs+=json.load(open(f))["runs"]
ok=[r for r in runs if r["viol"]<=1e-3 and r["mirror_matches"]==0 and r["on_plane"]==0]
byk={}
for r in ok: byk.setdefault(r["kinds"],[]).append(r)
tab=[]
for k,v in sorted(byk.items(),key=lambda t:(len(t[0]),t[0])):
    L=sorted(x["left"] for x in v); n_all=sum(1 for r in runs if r["kinds"]==k)
    tab.append((len(k),k.count("C"),k.count("V"),len(v),n_all,L[0],np.median(L),sum(x<1 for x in L),sum(x<2 for x in L)))
best={}
for N in range(8,13):
    c=[r for r in ok if len(r["kinds"])==N]
    if c: best[N]=min(c,key=lambda r:r["left"])
json.dump(dict(table=tab,best={N:{k:v for k,v in r.items() if k!="p"} for N,r in best.items()}),open("asym_results.json","w"),indent=1)
# chart
fig,ax=plt.subplots(figsize=(10,6.5)); Ns=sorted(best)
ax.plot(list(OLD),[v[0] for v in OLD.values()],"o-",color="0.55",lw=2.5,ms=7,label="before (mirror layouts allowed)")
ax.plot(Ns,[best[N]["left"] for N in Ns],"s-",color="#8e44ad",lw=2.5,ms=8,label="now (no mirror layouts)")
for N in Ns:
    r=best[N]; ax.text(N,r["left"]*1.22,f"{r['left']:.2f}×",ha="center",fontsize=12,color="#8e44ad")
ax.axhline(1,color="green",lw=2); ax.text(12.35,1,"CMB flat\ntarget 1×",va="center",fontsize=12,color="green")
ax.set_yscale("log"); ax.set_ylim(0.3,10); ax.set_xlim(7.5,13.2); ax.set_xticks(range(8,13))
ax.set_yticks([0.5,1,2,5,10]); ax.set_yticklabels(["0.5","1","2","5","10"])
ax.set_xlabel("number of objects N (masses + deep voids)"); ax.set_ylabel("best leftover (× floor, log scale)")
ax.set_title("Best leftover with and without mirror-image layouts"); ax.legend(loc="upper right",frameon=False)
fig.tight_layout(); fig.savefig("a_best_vs_N_asym_logy.png",dpi=110); plt.close(fig)
# map of overall best
b=min(best.values(),key=lambda r:r["left"])
k=np.array(list(b["kinds"])); X=np.array(b["X"]); M=np.array(b["M"]); Rv=np.array(b["Rv"]); d=np.linalg.norm(X,axis=1)
l=np.degrees(np.arctan2(X[:,1],X[:,0]))%360; bb=np.degrees(np.arcsin(X[:,2]/d))
fig,(a1,a2)=plt.subplots(1,2,figsize=(15,7)); t=np.linspace(0,2*np.pi,300)
for rr,ls in ((14.3,"-"),(43,":")): a2.plot(rr*np.cos(t),rr*np.sin(t),"k",ls=ls,lw=1)
a2.plot(C_G[0]+A_IN_G*np.cos(t),C_G[1]+A_IN_G*np.sin(t),color="0.5",lw=1)
ax_l=np.radians(AXIS); a2.plot([-45*np.cos(ax_l),45*np.cos(ax_l)],[-45*np.sin(ax_l),45*np.sin(ax_l)],color="0.7",ls="--",lw=1)
for i in range(len(k)):
    if k[i]=="C":
        s=60+240*M[i]/2e52; a1.scatter(l[i],bb[i],s=s,c="#c0392b"); a2.scatter(X[i,0],X[i,1],s=s,c="#c0392b",zorder=3)
    else:
        a1.scatter(l[i],bb[i],s=60+240*Rv[i]/4,facecolors="none",edgecolors="#2471a3",lw=2.5)
        a2.add_patch(plt.Circle(X[i,:2],Rv[i],color="#2471a3",alpha=0.35))
a2.plot(0,0,"k+",ms=14); a2.set_aspect("equal"); a2.set_xlim(-45,45); a2.set_ylim(-45,45); a2.set_xlabel("Gpc")
a2.set_title("Top view (us = +, dashed = mirror axis ℓ 15°/195°,\ngrey arc = CoE interior edge)",fontsize=13)
a1.set_xlim(360,0); a1.set_ylim(-60,60); a1.set_xticks(range(0,361,45)); a1.axhline(0,color="0.8",lw=1)
a1.axvline(AXIS,color="0.7",ls="--",lw=1); a1.axvline(AXIS+180,color="0.7",ls="--",lw=1)
a1.axvspan(270,300,color="0.9",zorder=-1); a1.text(285,55,"Migkas",ha="center",fontsize=11)
a1.set_xlabel("galactic longitude ℓ (deg)"); a1.set_ylabel("latitude b (deg)"); a1.set_title("Sky positions (masses red dots, voids blue rings)",fontsize=13)
fig.suptitle(f"Best layout with no mirror images: {int((k=='C').sum())} masses + {int((k=='V').sum())} deep voids, leftover {b['left']:.2f}×"); fig.tight_layout()
fig.savefig("d_asym_solution_map.png",dpi=105); plt.close(fig)
# printout
for row in tab: print("N=%d %dC+%dV ok %d/%d best %.2f med %.2f <1x:%d <2x:%d"%row)
for N in sorted(best):
    r=best[N]; k=np.array(list(r["kinds"])); M=np.array(r["M"]); Rv=np.array(r["Rv"]); X=np.array(r["X"]); d=np.linalg.norm(X,axis=1)
    print(f"\nN={N} {r['kinds']} left={r['left']:.3f} odd={r['odd_frac']:.3f} closest pair (ang,dlog,slog)={np.round(r['closest_pair'],3)} minplane={r['min_plane_deg']:.1f}")
    print(" clusters %.3e voids %.3e net %.2e (%.2f%%)"%(M[k=='C'].sum(),M[k=='V'].sum(),M.sum(),100*M.sum()/M[k=='C'].sum()))
    dc=np.linalg.norm(X-C_G,axis=1)
    print(" bounds: M cap",int(((k=='C')&(M>=0.98*MHI)).sum()),"M floor",int(((k=='C')&(M<=1.05*MLO)).sum()),"Rv cap",int(((k=='V')&(Rv>=3.98)).sum()),"Rv floor",int(((k=='V')&(Rv<=0.52)).sum()),
          "d at 43",int((d>=42.9).sum()),"touch sphere",int((d-np.where(k=='V',Rv,0)<=14.35).sum()),"touch CoE",int((dc-np.where(k=='V',Rv,0)<=A_IN_G+0.05).sum()))
    for i in np.argsort(np.degrees(np.arctan2(X[:,1],X[:,0]))%360):
        print("  %s d=%.1f l=%.0f b=%+.0f M=%.2e Rv=%.2f"%(k[i],d[i],np.degrees(np.arctan2(X[i,1],X[i,0]))%360,np.degrees(np.arcsin(X[i,2]/d[i])),M[i],Rv[i]))
    print(" top modes",np.round(r["amps"],2))
