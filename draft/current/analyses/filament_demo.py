import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import least_squares
exec(open("scan2.py").read().split("out={")[0])
S=json.load(open("scan.json")); p=S["sym"]["4"]["p"]; D0=p[0]; KG=7.7e51/0.0039
nodes0=[(phi+p[2],p[1]),(phi-p[2],p[1]),(phi+p[4],p[3]),(phi-p[4],p[3])]
def ptfield(a_deg,rho): return field(np.radians(a_deg),rho)
def arc(a1,a2,rho,n=None):
    n=n or max(8,int(abs(a2-a1)/0.25))  # unit total mass spread along arc a1->a2 (deg) at radius rho
    a=np.linspace(a1,a2,n); return sum(ptfield(x,rho) for x in a)/n
def radial(a,rho,half,n=40):
    r=np.linspace(rho-half,rho+half,n); return sum(ptfield(a,x) for x in r)/n
amps=lambda x:np.array([np.hypot(2*np.mean(x*np.cos(k*th)),2*np.mean(x*np.sin(k*th))) for k in range(1,16)])
def ring_order(nodes): return sorted(nodes,key=lambda t:t[0]%360)
def links_model(nodes,rho,Mlink):
    ns=ring_order(nodes); out=np.zeros_like(th); segs=[]
    angs=[n[0]%360 for n in ns]; L=[(angs[(i+1)%4]-angs[i])%360 for i in range(4)]; Ltot=sum(L)
    ML=Mlink if np.ndim(Mlink) else [Mlink*L[i]/Ltot for i in range(4)]
    for i in range(4):
        a1=angs[i]; a2=a1+L[i]; out+=ML[i]*arc(a1,a2,rho); segs.append((a1,a2,ML[i]))
    return out,segs
def links_from_nodes(nodes,rho,f):  # each node hands fraction f of its mass, half to each neighbouring link
    ns=ring_order(nodes); k=[n[1] for n in ns]; ML=[f*(k[i]+k[(i+1)%4])/2 for i in range(4)]
    return links_model(nodes,rho,ML)
Mtot=sum(k for _,k in nodes0)
res={}
pts=sum(k*ptfield(a,D0) for a,k in nodes0); res["4 points (fit)"]=pts
for f in (0.5,1.0):
    lm,segs=links_from_nodes(nodes0,D0,f); res[f"{int(f*100)}% of mass moved to links"]=(1-f)*pts+lm
for h in (10,20):
    res[f"each point → ±{h}° arc"]=sum(k*arc(a-h,a+h,D0) for a,k in nodes0)
res["each point → radial ±0.3R"]=sum(k*radial(a,D0,0.3) for a,k in nodes0)
# refit: nodes+links web, symmetric pairs, D, node mass per pair, link mass, all >=0
def web(q):
    D,k1,o1,k2,o2=q[:5]; ML=q[5:9]; ns=[(phi+o1,k1),(phi-o1,k1),(phi+o2,k2),(phi-o2,k2)]
    m=sum(k*ptfield(a,D) for a,k in ns); lm,_=links_model(ns,D,list(ML)); return m+lm
best=None
for ML0 in (0,0.004,0.008,0.016):
  for D_ in (1.4,1.7,2.2,3):
    f_=least_squares(lambda q:web(q)-R,[D_,p[1],p[2],p[3],p[4]]+[ML0*np.random.default_rng(int(D_*10)).uniform(0.2,1.8) for _ in range(4)],bounds=([1.05,0,1,0,1,0,0,0,0],[20,1,179,1,179,1,1,1,1]),max_nfev=1500)
    if best is None or f_.cost<best.cost: best=f_
q=best.x; res["refit: nodes + links (all positive)"]=web(q)
# links-only refit (no point mass at nodes)
def webL(q):
    D,o1,o2=q[:3]; ns=[(phi+o1,0),(phi-o1,0),(phi+o2,0),(phi-o2,0)]; return links_model(ns,D,list(q[3:7]))[0]
bestL=None
for D_ in (1.3,1.7,2.2,3,5):
  for o in ((147,60),(120,40),(160,90),(100,20)):
    f_=least_squares(lambda r:webL(r)-R,[D_,*o,0.004,0.004,0.004,0.004],bounds=([1.05,1,1,0,0,0,0],[20,179,179,5,5,5,5]),max_nfev=1500)
    if bestL is None or f_.cost<bestL.cost: bestL=f_
qL=bestL.x; res["refit: links only (positive)"]=webL(qL)
floor=1e-5
rows=[("post-SR residual (nothing added)",rms(R),amps(R))]
for k,v in res.items(): left=R-v; rows.append((k,rms(left),amps(left)))
for k,r,a in rows: print(f"{k:45s} RMS={r:.3e} ({r/floor:.0f}x)  m1..8 x floor:",np.round(a[:8]/floor).astype(int))
print("web refit D Gpc",q[0]*14.3,"node kg",q[1]*KG,q[3]*KG,"offs",q[2],q[4],"links kg",np.array(q[5:9])*KG)
print("links-only D Gpc",qL[0]*14.3,"offs",qL[1],qL[2],"links kg",np.array(qL[3:7])*KG)
print("4pt total",Mtot*KG)
json.dump({"rows":[(k,r,list(a)) for k,r,a in rows],"web":list(q),"linksonly":list(qL)},open("filament_demo.json","w"),indent=1)
# ---- figure
plt.rcParams.update({"font.size":12})
fig=plt.figure(figsize=(15,15))
def geo(ax,nodes,rho,segs,title,showpts=True,ML=None):
    ax.plot(np.linspace(0,2*np.pi,400),np.ones(400),"C0",lw=2)
    ax.plot([np.radians(phi)]*2,[0,rho*1.2],"k",lw=1); ax.plot([np.radians(phi+180)]*2,[0,rho*1.2],"k--",lw=1)
    for a1,a2,mm in segs:
        a=np.radians(np.linspace(a1,a2,100)); ax.plot(a,np.full(100,rho),color="C3",lw=1+8*min(1,mm/(max(x[2] for x in segs)+1e-30)),alpha=.8)
    if showpts:
        for a,k in nodes:
            if k>0: ax.scatter(np.radians(a),rho,s=60+400*k/0.004,color="C3",edgecolor="k",zorder=3)
    ax.set_ylim(0,rho*1.3); ax.set_yticklabels([]); ax.set_title(title,pad=18,fontsize=12)
a1=fig.add_subplot(3,2,1,projection="polar"); geo(a1,nodes0,D0,[],f"A. Four point clumps at {D0*14.3:.1f} Gpc\n{Mtot*KG:.1e} kg total")
_,segs=links_from_nodes(nodes0,D0,1.0)
a2=fig.add_subplot(3,2,2,projection="polar"); geo(a2,nodes0,D0,segs,"B. Same mass, same nodes, moved onto\nfilaments linking them along the shell",showpts=False)
ns=[(phi+q[2],q[1]),(phi-q[2],q[1]),(phi+q[4],q[3]),(phi-q[4],q[3])]; _,segs2=links_model(ns,q[0],list(q[5:9]))
a3=fig.add_subplot(3,2,3,projection="polar"); geo(a3,ns,q[0],segs2,f"C. Refit nodes + links (all positive)\nD = {q[0]*14.3:.1f} Gpc, links {np.sum(q[5:9])*KG:.1e} kg")
ax=fig.add_subplot(3,2,4); ax.plot(ELL,R*1e3,"C3",lw=2,label="post-SR residual")
for k,c in (("4 points (fit)","k"),("100% of mass moved to links","C1"),("refit: nodes + links (all positive)","C2")):
    ax.plot(ELL,res[k]*1e3,color=c,lw=1.6,ls="--",label=k)
ax.set_xlim(0,360); ax.set_xlabel("Galactic longitude ℓ (deg)"); ax.set_ylabel("ΔT/T ×10³"); ax.grid(alpha=.3); ax.legend(fontsize=9,loc="lower left")
ax=fig.add_subplot(3,1,3); ms=np.arange(1,16); sel=[0,1,2,3,7,6]; w=0.14
labs=["post-SR residual","after 4 points","after 50% to links","after 100% to links","links-only refit","nodes+links refit"]
for j,i in enumerate(sel):
    ax.bar(ms+(j-2.5)*w,rows[i][2]/floor,w,label=f"{labs[j]} (RMS {rows[i][1]/floor:.0f}×)")
ax.axhline(1,ls="--",color="purple",label="10⁻⁵ floor"); ax.set_yscale("log"); ax.set_xticks(ms); ax.set_xlabel("multipole m"); ax.set_ylabel("amplitude ÷ 10⁻⁵"); ax.legend(fontsize=9,ncol=2); ax.grid(alpha=.3,axis="y")
fig.suptitle("From points to links: the four-clump CMB sketch with filaments",fontsize=15)
fig.tight_layout(); fig.savefig("filament_demo.png",dpi=90); print("saved")
