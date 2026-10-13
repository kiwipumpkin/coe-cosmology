"""FCC lattice analogue vs foam s2 edited web (Chris, 10 Oct).
FCC clusters, cubic cell edge a = 171 Mpc (Einasto period), voids in the tetrahedral holes (4 clusters around each void) = fluorite-type
arrangement. 3D Delaunay tetrahedralization of FCC box and of the whole foam s2 web (ball 30 Gpc) BEFORE/AFTER the entropy-rule edits.
Masses: FCC uses the web's median cluster mass and median void deficit (negative mass). Outputs in fcc_lattice/."""
import sys,json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import edit_before_after as EB
from common import *
import seq_core
OUTD="fcc_lattice"; G=6.674e-11; GPC=seq_core.GPC
d=json.load(open("minedit_ent/b_foam_rim_s2_tidalcore_entropy.json")); base=d["tag"].replace("_entropy","")
LAY={"web BEFORE":EB.layout(base),"web AFTER":EB.layout(base,d)}
Xb,rb,vb,Mb,_=LAY["web BEFORE"]; MC=float(np.median(Mb[~vb])); MV=float(np.median(Mb[vb]))
A=0.171  # Gpc
# ---------------- FCC + tetrahedral-hole voids
def fcc(ncell,voids=True,octa=False):
    B=np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]]); T=np.array([[i,j,k] for i in (.25,.75) for j in (.25,.75) for k in (.25,.75)])
    O=np.array([[.5,.5,.5],[.5,0,0],[0,.5,0],[0,0,.5]])
    g=np.array([[i,j,k] for i in range(ncell) for j in range(ncell) for k in range(ncell)],float)
    P=[(g[:,None,:]+B[None]).reshape(-1,3)]; kind=[np.zeros(len(P[0]),int)]
    if voids: P.append((g[:,None,:]+T[None]).reshape(-1,3)); kind.append(np.ones(len(P[-1]),int))
    if octa:  P.append((g[:,None,:]+O[None]).reshape(-1,3)); kind.append(np.full(len(P[-1]),2))
    return np.vstack(P)*A,np.concatenate(kind)
def analyse(X,M,inner,name,vol_region):
    """inner(points)->bool: interior mask (stats only from tets / links / nodes fully interior, avoids hull slivers)."""
    rng=np.random.default_rng(1); Xj=X+rng.normal(0,1e-7*np.ptp(X),X.shape)   # break exact co-sphericity of the lattice
    D=Delaunay(Xj); S=D.simplices; ins=inner(X)
    ok=ins[S].all(1); S=S[ok]; q=X[S]
    vol=np.abs(np.einsum("ij,ij->i",q[:,1]-q[:,0],np.cross(q[:,2]-q[:,0],q[:,3]-q[:,0])))/6
    nz=vol>1e-6*np.median(vol); n_flat=int((~nz).sum()); S=S[nz]; vol=vol[nz]   # flat slivers from co-spherical lattice points (zero volume)
    E=np.unique(np.sort(np.vstack([S[:,[i,j]] for i in range(4) for j in range(i+1,4)]),axis=1),axis=0)
    L=np.linalg.norm(X[E[:,0]]-X[E[:,1]],axis=1)
    deg=np.bincount(E.ravel(),minlength=len(X)); core=inner(X*1.0) & ins
    # degree only for nodes well inside (all their tets kept): nodes inside a smaller interior
    nin=int(ins.sum()); dsp=(vol_region/nin)**(1/3)
    U=-G*M[E[:,0]]*M[E[:,1]]/(L*GPC)
    p=vol/vol.sum(); Svol=float(-(p*np.log(p)).sum())
    lv=np.log10(vol/vol.mean()); bins=np.arange(-6,3.01,0.1); h,_=np.histogram(lv,bins); ph=h[h>0]/h.sum(); Shist=float(-(ph*np.log(ph)).sum())
    kinds=np.where(M[E[:,0]]>0,0,1)+np.where(M[E[:,1]]>0,0,1)
    return dict(name=name,n_nodes=nin,n_tets=int(len(S)),n_links=int(len(E)),mean_spacing_Gpc=dsp,mean_tet_vol=float(vol.mean()),
                S_vol=Svol,lnN=float(np.log(len(vol))),evenness=Svol/np.log(len(vol)),S_hist_logvol=Shist,cv_vol=float(vol.std()/vol.mean()),
                U_J=float(U.sum()),U_per_node_J=float(U.sum()/nin),u_dimless=float(U.sum()/nin/(G*MC**2/(dsp*GPC))),
                U_cc=float(U[kinds==0].sum()/nin),U_cv=float(U[kinds==1].sum()/nin),U_vv=float(U[kinds==2].sum()/nin),
                frac_links_cc=float((kinds==0).mean()),frac_links_cv=float((kinds==1).mean()),frac_links_vv=float((kinds==2).mean()),
                mean_link_over_spacing=float(L.mean()/dsp),n_flat_dropped=n_flat),dict(vol=vol/vol.mean(),L=L/dsp,deg=None,E=E,S=S)
RES={}; RAW={}
# FCC box: 10 cells per side (1.71 Gpc), stats from the inner 6x6x6 cells
NC=10; lo,hi=2*A,8*A
for nm,vo in (("FCC + tetra-hole voids",True),("FCC clusters only",False)):
    X,k=fcc(NC,voids=vo); M=np.where(k==0,MC,MV)
    inner=lambda P: np.all((P>=lo-1e-9)&(P<hi-1e-9),axis=1)
    RES[nm],RAW[nm]=analyse(X,M,inner,nm,(hi-lo)**3)
    E=RAW[nm]["E"]; deg=np.bincount(E.ravel(),minlength=len(X)); ii=np.all((X>=lo+A)&(X<hi-A),axis=1); RAW[nm]["deg"]=deg[ii]
    RES[nm]["mean_degree"]=float(deg[ii].mean()); RES[nm]["n_voids_per_cluster"]=float((k==1).sum()/(k==0).sum())
# web: whole modeled ball, interior = |x| < 28 Gpc (drops hull slivers at the 30 Gpc edge)
RIN=28.0; Vin=4/3*np.pi*RIN**3
for nm,(X,r,v,M,_) in LAY.items():
    inner=lambda P: np.linalg.norm(P,axis=1)<RIN
    RES[nm],RAW[nm]=analyse(X,M,inner,nm,Vin)
    E=RAW[nm]["E"]; deg=np.bincount(E.ravel(),minlength=len(X)); ii=np.linalg.norm(X,axis=1)<RIN-3; RAW[nm]["deg"]=deg[ii]
    RES[nm]["mean_degree"]=float(deg[ii].mean()); RES[nm]["n_voids_per_cluster"]=float(v.sum()/(~v).sum())
RES["_conventions"]=dict(a_Gpc=A,M_cluster_kg=MC,M_void_kg=MV,web_region="ball |x|<28 Gpc of the 30 Gpc model (whole modeled region, hull slivers dropped)",
    fcc_region="inner 6x6x6 cells of a 10^3-cell box (1.03 Gpc cube)",u_dimless="U per node / (G M_cluster^2 / mean spacing)",
    S_vol="-sum p ln p, p_i = V_i/sum V (size evenness: equal cells -> max)",S_hist="Shannon entropy of the histogram of log10(V/mean V), 0.1-dex bins (spread of cell sizes: delta spikes -> low)")
json.dump(RES,open(f"{OUTD}/fcc_vs_web.json","w"),indent=1)
for k,o in RES.items():
    if not k.startswith("_"): print(k,{kk:(round(vv,4) if isinstance(vv,float) else vv) for kk,vv in o.items()})
np.savez(f"{OUTD}/fcc_vs_web_raw.npz",**{f"{k}|{q}":RAW[k][q] for k in RAW for q in ("vol","L","deg")})
