"""Charts + numbers for the core-open rerun. -> best_vs_N_core_open_logy.png, side_view_<scen>.png, side_views_all.png, report_out.txt, best.json"""
import os,sys,json,glob,numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import run_co as rc
from co_core import *
plt.rcParams.update({"font.size":13})
OLD=SS
def ld(f):
    if not os.path.exists(f): return []
    r=json.load(open(f)); return r["runs"] if isinstance(r,dict) else r
runs={}
for f in glob.glob("raw/*.json"):
    for r in json.load(open(f)): runs.setdefault(r["scen"],[]).append(r)
def nb(r):  # number of bodies (partners + ordinary)
    if r["scen"]=="ep": return len(r["Xp"])+len(r["Xo"])
    return len(r["name"])
def best_by(scen,key):
    b={}
    for r in runs.get(scen,[]):
        if not r["ok"]: continue
        k=key(r)
        if k not in b or r["left"]<b[k]["left"]: b[k]=r
    return b
BN={s:best_by(s,nb) for s in ("out","asym","ao","ep")}
# ---- walled / earlier comparisons (from the earlier runs' tables and raw files)
WALL={"out":{4:35.6,6:15.9,8:4.46,9:2.96,10:1.26,11:0.58,12:0.72},
      "asym":{8:6.19,9:3.90,10:2.30,11:1.25,12:0.73}}
NOWALL43={"out":{4:33.7,6:14.2,8:2.97,10:0.46,11:0.25,12:0.074}}
ao_old={}
for f in glob.glob(f"{OLD}/all_open/raw_open/*.json"):
    for r in ld(f):
        if r["viol"]<=1e-3: n=len(r["kinds"]); ao_old[n]=min(ao_old.get(n,1e9),r["left"])
WALL["ao"]=ao_old    # all-open earlier: no wall, but reach 43 Gpc (core only partly reachable)
ep_old={}
for f in glob.glob(f"{OLD}/extreme_partners/raw*/*.json"):
    for r in ld(f):
        if r["viol"]<=1e-3 and abs(r["net"])<=1e51: n=r["npart"]+len(r["kinds"]); ep_old[n]=min(ep_old.get(n,1e9),r["left"])
WALL["ep"]=ep_old
LAB={"out":"outside, masses + deep voids","asym":"outside, no mirror layouts","ep":"giant partners (+ ordinary bodies)","ao":"all open (inside, shell, core)"}
COL={"out":"#c0392b","asym":"#8e44ad","ep":"#d68910","ao":"#1f77b4"}
out=[]
def P(*a): s=" ".join(str(x) for x in a); print(s); out.append(s)
# ---------------- chart 1: best vs N
fig,ax=plt.subplots(figsize=(10,7.5))
SKIP={"out":{9},"ao":{7}}   # undersampled mixes with no warm start (5+4, 7 masses); listed in report_out.txt, left off the chart
CO_=json.load(open("coreonly_best.json")) if os.path.exists("coreonly_best.json") else {}
for s in ("out","asym","ep","ao"):
    N=sorted(n for n in BN[s] if n not in SKIP.get(s,())); ax.plot(N,[BN[s][n]["left"] for n in N],"-o",color=COL[s],lw=2.5,ms=8,label=f"{LAB[s]}: core open")
    W=WALL[s]; Nw=sorted(n for n in W if n>=min(N)-1 and n<=max(N)+0)
    lab={"out":"walled","asym":"walled","ep":"walled","ao":"earlier (43 Gpc reach)"}[s]
    if Nw: ax.plot(Nw,[W[n] for n in Nw],"--",marker="s",mfc="white",color=COL[s],lw=1.6,ms=7,label=f"{LAB[s]}: {lab}")
if CO_:
    Nc=sorted(int(k) for k in CO_); ax.plot(Nc,[CO_[str(n)] for n in Nc],"-^",color="#6e2c00",lw=2,ms=8,label="all bodies forced inside the core")
ax.axhline(1,color="green",lw=2.5,label="target 1× (CMB flat, 10⁻⁵ floor)")
ax.set_yscale("log"); ax.set_ylim(0.03,400); ax.set_yticks([0.1,1,10,100]); ax.set_yticklabels(["0.1","1","10","100"])
ax.set_xlabel("number of bodies N"); ax.set_ylabel("best leftover (× floor, log scale)"); ax.set_title("Best leftover vs number of bodies, CoE core open")
ax.grid(True,which="major",color="0.9")
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.13),ncol=2,frameon=False,fontsize=11)
fig.tight_layout(); fig.savefig("best_vs_N_core_open_logy.png",dpi=110); plt.close(fig)
# ---------------- tables
P("# best leftover by N  (core open | walled/earlier)  and bodies in the core in that best layout")
for s in ("out","asym","ep","ao"):
    P(f"\n## {s}: {LAB[s]}")
    for n in sorted(BN[s]):
        r=BN[s][n]; cd=np.array(r["core_depth"]); ic=cd>0
        P(f"N={n:2d} {r['name']:22s} open {r['left']:8.3f}  walled/earlier {WALL[s].get(n,float('nan')):8.3f}"+(f"  nowall43 {NOWALL43['out'].get(n,float('nan')):.3f}" if s=="out" else "")+
          f"  in-core {ic.sum()}/{len(cd)}  depths {np.round(cd[ic],1).tolist()}  start={r['start']}")
# ---------------- core usage pooled
P("\n# where bodies go: pooled over feasible fits with leftover < 2x  (and < 1.5 x best of their class)")
pool_all=[]
for s in ("out","asym","ep","ao"):
    Xs=[];Ms=[];K=[]
    for r in runs.get(s,[]):
        if not r["ok"] or r["left"]>2: continue
        if s=="ep":
            X=np.r_[np.array(r["Xp"]).reshape(-1,3),np.array(r["Xo"]).reshape(-1,3)]; M=np.r_[r["Mp"],r["Mo"]]; k=["P"]*len(r["Xp"])+list(r["name"].split("+")[1] if "+" in r["name"] else "")
        else: X=np.array(r["X"]); M=np.array(r["M"]); k=list(r["name"])
        Xs.append(X); Ms.append(M); K+=k
    if not Xs: P(s,"no fits below 2x"); continue
    X=np.vstack(Xs); M=np.concatenate(Ms); K=np.array(K); cd=core_depth(X); ic=cd>0; d=np.linalg.norm(X,axis=1)
    lon=(np.degrees(np.arctan2(X[:,1],X[:,0])))%360
    P(f"{s}: {len(Xs)} fits, {len(K)} bodies; in core {ic.sum()} ({100*ic.mean():.0f}%)  masses in core {((K=='C')&ic).sum()}/{(K=='C').sum()}  voids {((K=='V')&ic).sum()}/{(K=='V').sum()}"+
      (f"  partners {((K=='P')&ic).sum()}/{(K=='P').sum()}" if s=="ep" else ""))
    if ic.any():
        P(f"   core bodies: depth median {np.median(cd[ic]):.1f} max {cd[ic].max():.1f} Gpc; distance from us median {np.median(d[ic]):.1f} max {d[ic].max():.1f} Gpc; lon median {np.median(lon[ic]):.0f}")
        P(f"   core bodies |M| median {np.median(np.abs(M[ic])):.2e} kg; deeper than 10 Gpc: {(cd>10).sum()}; farther than 50 Gpc from us: {(ic&(d>50)).sum()}")
# ---------------- do the core bodies matter? remove them from each scenario's best layout
P("\n# effect of the core bodies in the overall best layout per scenario (remove them, no refit)")
def model_of(r): return rc.model(r["scen"],r["name"])
def parts(r):
    """per-body ring patterns (n,360) and total; uses the same pattern code"""
    m=model_of(r); p=np.array(r["p"])
    if r["scen"]=="ep":
        du,Xp,Mp,Rp,dp,Xo,Mo,Rvo,Rco=m.decode(p); X=np.r_[Xp,Xo]; M=np.r_[Mp,Mo]
        Wb=M[:,None]*unit_field(X); base=du*ep_core.S1; return X,M,Wb,base
    X,M,Rv,Rc=m.decode(p); Wb=M[:,None]*unit_field(X)
    if r["scen"]=="ao": base=np.zeros(360); comp=lambda Mm:-Mm.sum()*ao_core.FLOC
    else: base=np.zeros(360); comp=lambda Mm:-(Mm.sum()/KG)*ss_core.FC
    return X,M,Wb,(base,comp)
BEST={}
for s in ("out","asym","ep","ao"):
    r=min([x for x in runs.get(s,[]) if x["ok"]],key=lambda x:x["left"]); BEST[s]=r
    X,M,Wb,b=parts(r); cd=core_depth(X); ic=cd>0
    if s=="ep": w=Wb.sum(0)+b; wn=Wb[~ic].sum(0)+b   # net-mass rule broken by removal; reported as is
    else:
        base,comp=b; w=Wb.sum(0)+comp(M); wn=Wb[~ic].sum(0)+comp(M[~ic])
    m=model_of(r); chk=m.leftover(np.array(r["p"]))
    a_core=amps(Wb[ic].sum(0))/FLOOR if ic.any() else np.zeros(12)
    P(f"{s} best {r['name']} {r['left']:.3f}x (check {chk:.3f}, recon {rms(R-w)/FLOOR:.3f}); core bodies {ic.sum()}: without them {rms(R-wn)/FLOOR:.2f}x;"
      f" their own pattern rms {rms(Wb[ic].sum(0))/FLOOR if ic.any() else 0:.2f}x, m1-6 {np.round(a_core[:6],2).tolist()}")
    BEST[s]["_ic"]=ic.tolist()
# per-body detail of each best
P("\n# best layouts body by body (lon deg, lat deg, dist Gpc, core depth Gpc, mass kg / radius Gpc, own ring rms x)")
for s,r in BEST.items():
    X,M,Wb,b=parts(r); cd=core_depth(X); d=np.linalg.norm(X,axis=1)
    if s=="ep": k=["P"]*len(r["Xp"])+list(r["name"].split("+")[1] if "+" in r["name"] else ""); sz=list(r["Rp"])+list(r["Rvo"])
    else: k=list(r["name"]); sz=r["Rv"]
    P(f"## {s} {r['name']} {r['left']:.3f}x"+(f" delta_us={r['delta_us']:+.3f}" if s=="ep" else ""))
    for i in range(len(k)):
        P(f"  {k[i]} lon {np.degrees(np.arctan2(X[i,1],X[i,0]))%360:6.1f} lat {np.degrees(np.arcsin(X[i,2]/d[i])):+6.1f} d {d[i]:6.1f} core {cd[i]:+6.1f}  M {M[i]:+.2e}"+
          (f" R {sz[i]:.2f}" if k[i]!="C" else "")+f"  own {rms(Wb[i])/FLOOR:8.1f}")
# ---------------- seeded starts: did seeded core bodies stay?
P("\n# core-seeded starts: bodies in the core at the end, and how those fits scored vs warm/cold (median of class-best ratio)")
for s in ("out","asym","ep","ao"):
    R_=[x for x in runs.get(s,[]) if x["ok"]]; byc={}
    for x in R_: byc.setdefault(x["name"],[]).append(x)
    rat={"warm":[],"cold":[],"core":[]}; endcore={"warm":[],"cold":[],"core":[]}
    for nm,L in byc.items():
        b=min(y["left"] for y in L)
        for y in L:
            t="core" if y["start"].startswith("core") else ("warm" if y["start"].startswith("warm") else "cold")
            rat[t].append(y["left"]/b); endcore[t].append(np.mean(np.array(y["core_depth"])>0))
    P(f"{s}: "+"; ".join(f"{t}: n={len(rat[t])} median left/best {np.median(rat[t]) if rat[t] else float('nan'):.2f}, wins {sum(1 for nm,L in byc.items() if min(L,key=lambda y:y['left'])['start'].startswith(t))}, mean share of bodies ending in core {100*np.mean(endcore[t]) if endcore[t] else 0:.0f}%" for t in rat))
# ---------------- kernel: one 2e52 kg mass toward the CoE
P("\n# one 2e52 kg mass on the us-CoE line: its ring pattern by mode (x floor); target m3..m5 = 200,185,173")
for dd in (20,30,45,60,80,100):
    a=amps(2e52*unit_field((CU*dd)[None])[0])/FLOOR; P(f"  d={dd:3d} Gpc (core depth {A_IN_G-abs(CN-dd):+.1f}): m1..m6 {np.round(a[:6],2).tolist()}")
P("start amps m1..12",np.round(amps(R)/FLOOR,1).tolist())
open("report_out.txt","w").write("\n".join(out)+"\n")
json.dump({s:{k:v for k,v in r.items() if k!="resid"} for s,r in BEST.items()},open("best.json","w"))
# ---------------- side views
t=np.linspace(0,2*np.pi,400); deg=np.degrees(TH)
def side(ax,r,s):
    ax.set_aspect("equal")
    ax.add_patch(Circle(C_G[:2],A_IN_G,color="#f5cba7",alpha=0.55,zorder=0,lw=0))
    ax.plot(C_G[0]+A_IN_G*np.cos(t),C_G[1]+A_IN_G*np.sin(t),color="#ca6f1e",lw=1.2)
    ax.plot(C_G[0],C_G[1],"x",color="#ca6f1e",ms=10,mew=2); ax.text(C_G[0],C_G[1]-6,"CoE",color="#ca6f1e",ha="center",va="top",fontsize=11)
    ax.text(C_G[0],C_G[1]+A_IN_G-6,"core (CoE interior)",color="#a04000",ha="center",va="top",fontsize=11)
    ax.plot(RV*np.cos(t),RV*np.sin(t),"k",lw=1.3); ax.plot(0,0,"k+",ms=14,mew=2)
    if s=="ep":
        Xp=np.array(r["Xp"]).reshape(-1,3); Rp=r["Rp"]; dp=r["dp"]
        for i in range(len(Xp)):
            ax.add_patch(Circle(Xp[i,:2],Rp[i],fill=False,ec="#2471a3" if dp[i]<0 else "#c0392b",lw=2,ls="-")); ax.text(Xp[i,0],Xp[i,1],f"P{i+1}",ha="center",va="center",fontsize=11)
        X=np.array(r["Xo"]).reshape(-1,3); M=np.array(r["Mo"]); Rv=np.array(r["Rvo"]); k=list(r["name"].split("+")[1]) if "+" in r["name"] else []
    else: X=np.array(r["X"]); M=np.array(r["M"]); Rv=np.array(r["Rv"]); k=list(r["name"])
    for i in range(len(k)):
        if k[i]=="C": ax.scatter(X[i,0],X[i,1],s=40+260*M[i]/2e52,c="#c0392b",zorder=3,edgecolors="k",lw=0.5)
        else: ax.add_patch(Circle(X[i,:2],Rv[i],color="#2471a3",alpha=0.45,zorder=2))
    ax.set_xlim(-110,50); ax.set_ylim(-70,70); ax.set_xlabel("Gpc"); ax.set_ylabel("Gpc")
def ring(ax,r,s):
    y=np.array(r["resid"])/FLOOR; ax.plot(deg,y,color=COL[s],lw=2); ax.axhline(0,color="0.8",lw=1)
    ax.axhspan(-1,1,color="green",alpha=0.25,lw=0,zorder=0)
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle around the ring (galactic longitude ℓ, degrees)"); ax.set_ylabel("leftover (× floor)")
def title(r,s):
    ic=np.array(r["_ic"]); n=nb(r)
    return f"{LAB[s]}: best layout, N={n}, leftover {r['left']:.2f}×, bodies in core: {ic.sum()}"
for s,r in BEST.items():
    fig,axs=plt.subplots(2,1,figsize=(10,12),gridspec_kw=dict(height_ratios=[2.6,1]))
    side(axs[0],r,s); axs[0].set_title(title(r,s)+"\nside view (CMB-ring plane seen face-on; off-plane bodies projected)",fontsize=12)
    ring(axs[1],r,s); axs[1].set_title("leftover around the ring (own y-scale; green band = ±1× CMB level)",fontsize=12)
    fig.text(0.5,0.005,"+ = us; black circle = our 14.3 Gpc sphere; shaded = CoE core, open to bodies;\nred dots = masses (sized by mass); blue = voids at true radius; P = giant partner",ha="center",fontsize=10)
    fig.tight_layout(rect=(0,0.04,1,1)); fig.savefig(f"side_view_{s}.png",dpi=100); plt.close(fig)
print("charts done")
