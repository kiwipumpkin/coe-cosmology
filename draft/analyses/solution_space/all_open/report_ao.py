import json,glob,os,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from ao_core import *
plt.rcParams.update({"font.size":14,"axes.titlesize":15,"legend.fontsize":12})
KEY=[4,6,8,10,11,12]
def load(d):
    runs=[]
    for f in glob.glob(f"{d}/*.json"): runs+=json.load(open(f))
    return [r for r in runs if r["viol"]<=1e-3 and r.get("mirror_matches",0)==0 and r.get("on_plane",0)==0]
def bestN(runs):
    b={}
    for r in runs:
        N=len(r["kinds"])
        if N not in b or r["left"]<b[N]["left"]: b[N]=r
    return b
A={v:load(f"raw_{v}") for v in ("smooth","open")}; BA={v:bestN(A[v]) for v in A}
AS=load("raw_smooth_asym"); BAS=bestN(AS)
ins=json.load(open(f"{CS}/inside_visible/results_cap2e52.json")); IN={}
for c in ins.values(): IN[c["N"]]=min(IN.get(c["N"],1e9),c["best"])
OUT={n:r["left"] for n,r in bestN([r for f in glob.glob(f"{SS}/no_core_wall/raw_smooth/*.json") for r in json.load(open(f)) if r["viol"]<=1e-3]).items()}
fig,ax=plt.subplots(figsize=(10.5,6.8))
def ln(d,fmt,col,lab,lw=2.5):
    N=[n for n in KEY if n in d]; ax.plot(N,[d[n] if not isinstance(d[n],dict) else d[n]["left"] for n in N],fmt,color=col,lw=lw,ms=7,label=lab)
ln(OUT,"s--","#7f8c8d","outside only (no wall)")
ln(IN,"o-","#c0392b","inside our sphere only")
ln(BA["smooth"],"D-","#1f77b4","all open (smooth-density limit)")
ln(BA["open"],"^:","#5dade2","all open (no density limit)",lw=2)
if BAS:
    N_=sorted(BAS); ax.plot(N_,[BAS[n]["left"] for n in N_],"P",color="#17a589",ms=13,label="all open + no-mirror rule")
ax.axhline(1,color="green",lw=2); ax.text(12.25,1.1,"CMB flat\ntarget 1×",fontsize=12,color="green",va="bottom")
ax.set_yscale("log"); ax.set_ylim(0.005,60); ax.set_xlim(3.5,13.3); ax.set_xticks(KEY)
ax.set_yticks([0.01,0.1,1,10]); ax.set_yticklabels(["0.01","0.1","1","10"])
ax.set_xlabel("number of objects N (masses and/or deep voids)"); ax.set_ylabel("best leftover (× 10⁻⁵ floor, log scale)")
ax.set_title("Everything open vs inside-only vs outside-only"); ax.legend(loc="upper right",frameon=False)
fig.tight_layout(); fig.savefig("a_best_vs_N_all_open_logy.png",dpi=110); plt.close(fig)
below=[n for n in sorted(BA["smooth"]) if BA["smooth"][n]["left"]<1]; Nm=below[0] if below else min(BA["smooth"],key=lambda n:BA["smooth"][n]["left"])
m=BA["smooth"][Nm]; k=np.array(list(m["kinds"])); X=np.array(m["X"]); M=np.array(m["M"]); Rv=np.array(m["Rv"]); d=np.linalg.norm(X,axis=1)
l=np.degrees(np.arctan2(X[:,1],X[:,0]))%360; b=np.degrees(np.arcsin(X[:,2]/d))
fig,(a1,a2)=plt.subplots(1,2,figsize=(15,7)); t=np.linspace(0,2*np.pi,300)
a2.plot(RV*np.cos(t),RV*np.sin(t),"k",lw=1.3); a2.text(0,RV+1,"our 14.3 Gpc sphere",ha="center",fontsize=11)
a2.plot(43*np.cos(t),43*np.sin(t),"k:",lw=1)
ax_=C_G[0]+A_IN_G*np.cos(t); ay_=C_G[1]+A_IN_G*np.sin(t); a2.plot(ax_,ay_,color="0.6",ls="--",lw=1.2)
ok_=(np.abs(ax_)<44)&(np.abs(ay_)<44); j_=np.argmax(np.where(ok_,ay_,-1e9)); a2.text(ax_[j_]+1.5,ay_[j_]-4,"old wall",color="0.4",fontsize=12)
for i in range(len(k)):
    if k[i]=="C":
        s=50+240*M[i]/2e52; a1.scatter(l[i],b[i],s=s,c="#c0392b",edgecolors="k" if d[i]<RV else "none",lw=1.5); a2.scatter(X[i,0],X[i,1],s=s,c="#c0392b",zorder=3)
    else:
        a1.scatter(l[i],b[i],s=60+240*Rv[i]/4,facecolors="none",edgecolors="#2471a3",lw=2.5); a2.add_patch(plt.Circle(X[i,:2],Rv[i],color="#2471a3",alpha=0.35))
a2.plot(0,0,"k+",ms=14); a2.set_aspect("equal"); a2.set_xlim(-45,45); a2.set_ylim(-45,45); a2.set_xlabel("Gpc")
a2.set_title("Top view (us = +; dotted = 43 Gpc reach)",fontsize=13)
a1.set_xlim(360,0); a1.set_ylim(-90,90); a1.set_xticks(range(0,361,45)); a1.axhline(0,color="0.8",lw=1)
a1.axvspan(270,300,color="0.9",zorder=-1); a1.text(285,82,"Migkas",ha="center",fontsize=11)
a1.set_xlabel("galactic longitude ℓ (deg)"); a1.set_ylabel("latitude b (deg)")
a1.set_title("Sky positions (masses red, black edge = inside our sphere;\nvoids blue rings)",fontsize=13)
fig.suptitle(f"All open, smooth-density limit: {int((k=='C').sum())} masses + {int((k=='V').sum())} deep voids, leftover {m['left']:.2f}×"); fig.tight_layout()
fig.savefig("d_all_open_map.png",dpi=105); plt.close(fig)
# text
print("inside-only",{n:round(IN.get(n,np.nan),3) for n in KEY}); print("outside-only",{n:round(OUT.get(n,np.nan),3) for n in KEY})
for v in ("smooth","open"):
    for n in KEY:
        L=sorted([r for r in A[v] if len(r["kinds"])==n],key=lambda r:r["left"])
        if not L: continue
        r=L[0]; d=np.linalg.norm(r["X"],axis=1); kk=np.array(list(r["kinds"]))
        print(f"{v} N={n} best {r['left']:.3f} {r['kinds']} ({r['start']}) runs={len(L)} below1x={sum(x['left']<1 for x in L)} inside {int((d<RV).sum())}/{n} (masses {int(((d<RV)&(kk=='C')).sum())}/{int((kk=='C').sum())})")
for n,r in sorted(BAS.items()): print(f"asym smooth N={n} best {r['left']:.3f} {r['kinds']} odd {r.get('odd_frac',0):.4f} runs={sum(len(x['kinds'])==n for x in AS)}")
for v in ("smooth",):
    G=[r for r in A[v] if r["left"]<2]
    for t_ in ("C","V"):
        dd=np.concatenate([np.linalg.norm(np.array(r["X"])[np.array(list(r["kinds"]))==t_],axis=1) for r in G if t_ in r["kinds"]])
        print(v,t_,"in fits<2x:",len(dd),"inside sphere %.0f%%"%(100*(dd<RV).mean()),"median d inside %.1f"%np.median(dd[dd<RV]) if (dd<RV).any() else "","median d outside %.1f"%np.median(dd[dd>=RV]) if (dd>=RV).any() else "")
    dep=np.concatenate([A_IN_G-np.linalg.norm(np.array(r["X"])-C_G,axis=1) for r in G]); print("objects inside old wall: %.0f%%"%(100*(dep>0).mean()))
d=np.linalg.norm(X,axis=1); isC=k=="C"; rho=rho_gpc(X); dep=A_IN_G-np.linalg.norm(X-C_G,axis=1)
print(f"\nMAP N={Nm} {m['kinds']} left={m['left']:.3f} masses {M[isC].sum():.2e} voids {M[~isC].sum():.2e} net {M.sum():.1e}")
print(" bounds: mass at density cap",int((isC&(M>=0.98*MHI*rho/RHO)).sum()),"mass floor",int((isC&(M<=1.05*MLO)).sum()),"void Rv=4",int(((~isC)&(Rv>=3.98)).sum()),"d>=42.9",int((d>=42.9).sum()),
      "ring margin",int((np.abs(d-RV)-np.where(isC,np.array(m['Rc']),Rv)<=0.15).sum()),"near us",int((d-np.where(isC,np.array(m['Rc']),Rv)<=0.55).sum()))
for i in np.argsort(d): print("  %s d=%.1f l=%.0f b=%+.0f M=%.2e Rv=%.2f %s inside_wall=%+.1f"%(k[i],d[i],l[i],b[i],M[i],Rv[i],"INSIDE" if d[i]<RV else "outside",dep[i]))
print(" top modes",np.round(m["amps"],2))

# by start type: best per N for warm-inside / warm-outside / mixed / cold
for v in ("smooth","open"):
    for n in KEY:
        row=[]
        for tg in ("warm_inside","warm_outside","mixed","cold"):
            L=[r for r in A[v] if len(r["kinds"])==n and r["start"].startswith(tg)]
            if L:
                r=min(L,key=lambda r:r["left"]); dd=np.linalg.norm(r["X"],axis=1); row.append(f"{tg}:{r['left']:.3f}({int((dd<RV).sum())}/{n} in)")
        print(v,n," ".join(row))
