"""One quadrant of the red CoE universe, 0..D+R from the CoE, seed-2000 web style.
Template: seed-2000 vmn_gen web (Mmin 3e15, lcdm), inner cube side L tiled periodically; bodies kept with
prob shp(r_CoE/42)/pmax (densweb/gen.py rule: number density ∝ red rho, masses unchanged)."""
import sys,numpy as np,matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection,PatchCollection
sys.path.insert(0,'/workspace/coe-cosmology-gh/draft/current/analyses/void_mass_network');import vmn_gen as G
g=G.generate(Mmin=3e15,seed=2000,comp_mode="lcdm");S=g["S"]/1000;XC=g["XC"]/1000;Rv=g["Rv"]/1000;L=g["links"];MC=g["MC"]
CC=[-0.041,-0.035,0.1353,0.0];Rb=1.5*(1+1/2.937);D=42.0;R=14.3;RM=D+R
shp=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)))
xs=np.linspace(0,RM/D,4000);pmax=shp(xs).max()
Lc=23.0;h=Lc/2  # cube inside 20 Gpc ball (half-diag 19.9)
iv=np.all(np.abs(S)<h,1);ic=np.all(np.abs(XC)<h,1)
HZ=4.0 # slab kept for profile
rng=np.random.default_rng(7)
# tile so that z=0 is a tile boundary (template's empty 0.5 Gpc hole never in z=0 slab)
Vs,Vr,Cs,Cm,Lk=[],[],[],[],[]
for i in range(-1,int(RM//Lc)+2):
 for j in range(-int(RM//Lc)-2,int(RM//Lc)+2):
  for k in (-1,0):
   off=np.array([i*Lc+h,j*Lc+h,k*Lc+h])
   sv=S[iv]+off;sc=XC[ic]+off
   av=rng.random(len(sv))<shp(np.linalg.norm(sv,axis=1)/D)/pmax
   ac=rng.random(len(sc))<shp(np.linalg.norm(sc,axis=1)/D)/pmax
   inq=lambda P:(P[:,0]>-1)&(np.abs(P[:,1])<P[:,0]+2)&(np.linalg.norm(P[:,:2],axis=1)<RM+1)&(np.abs(P[:,2])<HZ)
   mv=av&inq(sv);mc=ac&inq(sc)
   gv=-np.ones(len(S),int);gv[np.where(iv)[0][mv]]=np.arange(mv.sum())+sum(len(a) for a in Vs)
   gc=-np.ones(len(XC),int);gc[np.where(ic)[0][mc]]=np.arange(mc.sum())+sum(len(a) for a in Cs)
   Vs.append(sv[mv]);Vr.append(Rv[iv][mv]);Cs.append(sc[mc]);Cm.append(MC[ic][mc])
   l=np.stack([gc[L[:,0]],gv[L[:,1]]],1);Lk.append(l[(l>=0).all(1)])
Vs,Vr,Cs,Cm,Lk=map(np.vstack if 0 else np.concatenate,(Vs,Vr,Cs,Cm,Lk))
US=np.array([D,0,0])
kv=np.linalg.norm(Vs-US,axis=1)>0.5;kc=np.linalg.norm(Cs-US,axis=1)>0.5  # 0.5 Gpc empty zone around us
# ---- profile (3D slab |z|<HZ, quadrant x,y>0)
rb=np.linspace(0,RM,29);rc=0.5*(rb[1:]+rb[:-1])
def prof(P,w=None):
  q=(P[:,0]>0)&(np.abs(P[:,1])<=P[:,0]);r=np.linalg.norm(P[q],axis=1);n,_=np.histogram(r,rb,weights=None if w is None else w[q])
  # volume of quarter-annulus x slab (MC)
  return n
rr=rng.random((800000,3))*[RM,2*RM,2*HZ]-[0,RM,HZ];rr=rr[np.abs(rr[:,1])<=rr[:,0]];vol,_=np.histogram(np.linalg.norm(rr,axis=1),rb);vol=vol/800000*RM*2*RM*2*HZ
nC=prof(Cs)/vol;nV=prof(Vs)/vol;mC=prof(Cs,Cm.astype(float))/vol
ok=vol>0.5
# ---- figure
fig,ax=plt.subplots(2,1,figsize=(12,26),gridspec_kw=dict(height_ratios=[1.75,1]))
a=ax[0];hz=0.4
wd=lambda P:(np.abs(P[:,1])<=P[:,0])&(np.linalg.norm(P[:,:2],axis=1)<=RM)
sv=(np.abs(Vs[:,2])<hz)&kv&wd(Vs);sc=(np.abs(Cs[:,2])<hz)&kc&wd(Cs)
rcut=np.sqrt(np.maximum(Vr[sv]**2-Vs[sv,2]**2,0))
a.add_collection(PatchCollection([plt.Circle(p[:2],r) for p,r in zip(Vs[sv],rcut) if r>0],facecolor="#2471a3",alpha=0.35,edgecolor="none",zorder=1))
lk=Lk[sc[Lk[:,0]]&sv[Lk[:,1]]];nl=len(lk);dec=1
if nl>15000: dec=int(np.ceil(nl/15000));lk=lk[::dec]
a.add_collection(LineCollection(np.stack([Cs[lk[:,0],:2],Vs[lk[:,1],:2]],1),colors="0.3",linewidths=0.25,alpha=0.5,zorder=2))
m=0.01*Cm[sc];a.scatter(Cs[sc,0],Cs[sc,1],s=1.5+25*(m/m.max())**1.5,c="#c0392b",edgecolors="none",zorder=3)
t=np.linspace(0,2*np.pi,800);a.plot(D+R*np.cos(t),R*np.sin(t),"k",lw=1.6,zorder=4)
tw=np.linspace(-np.pi/4,np.pi/4,300)
for rr_ in (D-R,D,D+R): a.plot(rr_*np.cos(tw),rr_*np.sin(tw),"k:",lw=0.9,zorder=4)
for sg in (1,-1): a.plot([0,RM*np.cos(np.pi/4)],[0,sg*RM*np.sin(np.pi/4)],"k-",lw=1.2,zorder=4)
a.plot(0,0,"k*",ms=18,zorder=5);a.text(0.6,0.6,"CoE",fontsize=13,weight="bold",zorder=6,bbox=dict(fc="w",alpha=.8))
a.plot(D,0,"k+",ms=16,mew=2.5,zorder=5);a.text(D+0.5,0.6,"us",fontsize=13,weight="bold",zorder=6,bbox=dict(fc="w",alpha=.8))
a.set_xlim(0,RM);a.set_ylim(-39.8,39.8);a.set_aspect("equal");a.set_xlabel("Gpc from CoE along CoE→us axis (l 0°, b 15°)");a.set_ylabel("Gpc from CoE (perpendicular, in the CoE–us plane)")
a.set_title(f"Red CoE universe, wedge ±45° about CoE→us axis, r = 0 to D+R = {RM:.1f} Gpc: slab |z| < {hz} Gpc, seed-2000 web at 1% mass\n{sc.sum()} clusters, {sv.sum()} void cores, {nl} void–cluster links"+(f" (every {dec}th link drawn)" if dec>1 else "")+f"; black arc = our R = {R} Gpc sphere, dotted = r = D−R, D, D+R",fontsize=12)
a.text(0.01,0.99,"blue = void cores, red = clusters (size ~ mass), grey = bipartite void↔cluster links\nnumber density ∝ red ρ(r_CoE) (poly ln ρ, C = (−0.041, −0.035, 0.1353, 0)), web tiled from seed-2000 template",transform=a.transAxes,va="top",fontsize=10.5,bbox=dict(fc="w"),zorder=7)
b=ax[1];rho=shp(rc/D)/shp(1.0)
nC0=np.nanmean(nC[ok][np.abs(rc[ok]-D)<4]);nV0=np.nanmean(nV[ok][np.abs(rc[ok]-D)<4]);mC0=np.nanmean(mC[ok][np.abs(rc[ok]-D)<4])
b.plot(rc[ok],nC[ok]/nC0,"o-",c="#c0392b",label="cluster number density (÷ value at D)")
b.plot(rc[ok],nV[ok]/nV0,"s-",c="#2471a3",label="void number density (÷ value at D)")
b.plot(rc[ok],mC[ok]/mC0,"^-",c="k",alpha=.6,label="cluster mass density (÷ value at D)")
rf=np.linspace(0,RM,500);b.plot(rf,shp(rf/D)/shp(1),"m",lw=2.5,label="red ρ(r)/ρ(D)")
for v_,s_ in ((D-R,"D−R"),(D,"D"),(D+R,"D+R")): b.axvline(v_,c="0.4",ls="--");b.text(v_+0.3,1.17,s_,fontsize=12)
b.axvspan(D-20,D+20,color="y",alpha=.12,label="r covered by the fitted 20 Gpc web")
b.set_xlim(0,RM);b.set_ylim(0.6,1.25);b.set_xlabel("r from CoE (Gpc)");b.set_ylabel("relative density");b.grid(alpha=.3)
b.set_title(f"Density vs r from CoE (counted in the 3D wedge |angle| ≤ 45°, slab |z| < {HZ} Gpc, 2 Gpc bins)",fontsize=12)
b.legend(loc="lower center",ncol=2,fontsize=10.5)
fig.tight_layout();fig.savefig("wedge_web.png",dpi=100)
np.savez("wedge_profile.npz",r=rc[ok],nC=nC[ok],nV=nV[ok],mC=mC[ok],rho=rho[ok])
print("pmax",pmax,"argmax r",xs[shp(xs).argmax()]*D,"rho(0)",shp(0)/shp(1),"rho(D+R)",shp(RM/D),"links",nl,dec,sc.sum(),sv.sum())
for x,y,z,w in zip(rc[ok],nC[ok]/nC0,nV[ok]/nV0,rho[ok]):print(f"{x:5.1f} {y:.3f} {z:.3f} {w:.3f}")
