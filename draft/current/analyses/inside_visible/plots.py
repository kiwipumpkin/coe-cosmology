import json, glob, os, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
plt.rcParams.update({"font.size":14,"axes.titlesize":15,"legend.fontsize":12})
SS="../solution_space"
def best_out():
    m={}; x={}
    for f in glob.glob(f"{SS}/raw/*.json"):
        d=json.load(open(f)); k=d["kinds"]
        if "S" in k: continue
        L=min(r["left"] for r in d["runs"] if r["viol"]<=1e-3); N=len(k)
        if set(k)=={"C"}: m[N]=min(m.get(N,1e9),L)
        elif set(k)=={"C","V"}: x[N]=min(x.get(N,1e9),L)
    for f,dd in (("extend_masses.json",m),("extend_mixed.json",x)):
        for s in json.load(open(f"{SS}/{f}")):
            if s["N"]>8 and s["N"]<=12: dd[s["N"]]=s["left"]
    return m,x
def best_in(tag):
    r=json.load(open(f"results_{tag}.json")); m={}; x={}; v={}
    for c in r.values():
        dd={"masses only":m,"masses + deep voids":x,"deep voids only":v}[c["mix"]]
        dd[c["N"]]=min(dd.get(c["N"],1e9),c["best"])
    return m,x,v,r
om,ox=best_out(); im,ix,iv,res=best_in("cap2e52")
fig,ax=plt.subplots(figsize=(10.5,7.5))
def line(d,col,ls,lab,mk="o"):
    N=sorted(d); ax.plot(N,[d[n] for n in N],ls,color=col,lw=2.6,ms=7,marker=mk,label=lab)
line(im,"#c0392b","-","inside: masses only")
line(ix,"#8e44ad","-","inside: masses + deep voids")
im53,_,_,_=best_in("cap1e53"); im53={n:v for n,v in im53.items()}
N=sorted(im53); ax.plot(N,[im53[n] for n in N],":",color="#e59866",lw=2,marker="^",ms=6,label="inside: masses only, cap 1e53 kg")
line(om,"#c0392b","--","outside: masses only",mk="s")
line(ox,"#8e44ad","--","outside: masses + deep voids",mk="s")
for v,t,c,ls in [(237.3,"start 237×","0.4","--"),(88.8,"four clumps 89×","0.4","--"),(1.0,"CMB flat target 1×","green","-")]:
    ax.axhline(v,color=c,ls=ls,lw=1.8 if c=="green" else 1); ax.text(12.3,v*1.12,t,va="bottom",fontsize=11,color=c if c=="green" else "0.35")
ax.set_yscale("log"); ax.set_xlim(0.7,14.6); ax.set_ylim(0.01,400); ax.set_xticks(range(1,13))
ax.set_xlabel("number of objects N"); ax.set_ylabel("best leftover (× 10⁻⁵ floor)")
ax.set_title("Objects inside our visible sphere vs outside it")
ax.legend(loc="lower left",frameon=False); fig.tight_layout(); fig.savefig("best_vs_N_inside_vs_outside_logy.png",dpi=110); plt.close(fig)
# top view of best inside solution: chosen = fewest objects reaching < 1x (else overall best)
pick=res[sys.argv[1] if len(sys.argv)>1 else "CCCCVVVV"]
b=pick["best_run"]; X=np.array(b["X"]); M=np.array(b["M"]); Rv=np.array(b["Rv"]); Rc=np.array(b["Rc"]); k=np.array(list(pick["kinds"]))
fig,axs=plt.subplots(1,2,figsize=(15,7.6))
for ax,(i,j,lab) in zip(axs,[(0,1,"top view (looking down on the CMB ring plane)"),(0,2,"side view (ring plane edge-on)")]):
    ax.add_patch(Circle((0,0),14.3,fill=False,ls="--",color="k",lw=1.5))
    for n in range(len(k)):
        if k[n]=="C":
            ax.add_patch(Circle((X[n,i],X[n,j]),Rc[n],color="#c0392b",alpha=0.9,zorder=3))
            if j==1: ax.text(X[n,i]+Rc[n]+0.2,X[n,j]+0.2,f"{M[n]:.1e}".replace("e+","e"),fontsize=9,color="#922b21")
        else: ax.add_patch(Circle((X[n,i],X[n,j]),Rv[n],color="#2471a3",alpha=0.45))
    ax.plot(0,0,"k+",ms=16,mew=2.5); ax.set_xlim(-16,16); ax.set_ylim(-16,16); ax.set_aspect("equal")
    if j==2: ax.axhline(0,color="0.6",lw=1)
    ax.set_xlabel("Gpc"); ax.set_ylabel("Gpc"); ax.set_title(lab,fontsize=13)
# CoE direction arrow
from inside_core import C_G
u=C_G/np.linalg.norm(C_G); axs[0].annotate("",xy=(u[0]*15.5,u[1]*15.5),xytext=(u[0]*12.5,u[1]*12.5),arrowprops=dict(arrowstyle="->",color="0.3"))
axs[0].text(u[0]*13.0,u[1]*13.0+1.3,"to CoE",fontsize=10,color="0.3",ha="center")
fig.suptitle(f"Inside-sphere solution: {pick['n_C']} masses, {pick['n_V']} voids, leftover {pick['best']:.2f}× (masses kg labelled)",fontsize=14)
fig.tight_layout(); fig.savefig(f"inside_universe_{pick['kinds']}.png",dpi=110); plt.close(fig)
print("picked",pick["kinds"],pick["best"]); print("inside m",im); print("inside x",ix); print("inside v",iv); print("out m",om); print("out x",ox)
