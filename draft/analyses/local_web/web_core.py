"""Local cosmic-web ("foam") model around the visible sphere, scored with the SAME tidal pattern, residual and floor as
../shell_ensemble and ../lattice (imported read-only via ../lattice/lattice_core.py; nothing upstream is modified).

Geometry (Gpc unless noted): region = 14.3 < |x| <= RMAX (default 30) from us, every direction, with CoE distance >= D-d = 43.9
(CoE interior excluded, as in all earlier surround/lattice/depth runs). Positions are converted to R_vis units for the tidal code."""
import os, sys, numpy as np, numba as nb
from scipy.spatial import Delaunay
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(os.path.dirname(HERE),"lattice"))
_cwd=os.getcwd(); import lattice_core as lc; os.chdir(_cwd)
R,R4,BASE0,BASE4,FLOOR,KG,RV,rms,amps=lc.R,lc.R4,lc.BASE0,lc.BASE4,lc.FLOOR,lc.KG,lc.RV,lc.rms,lc.amps
C_G=lc.C*RV                          # CoE position in Gpc (58.2 Gpc from us)
A_IN_G=lc.A_IN*RV                    # D-d = 43.9 Gpc
L_G=lc.L*RV; D_G=lc.D*RV             # exp-lock scale length 922.65 Gpc, D = 58.2 Gpc
GPC=lc.GPC; RHO=lc.RHO               # 3.0857e25 m, 2.9e-27 kg/m^3 (mean matter density, same constant as earlier runs)
TH=lc.field3.__globals__["th"]       # 360 in-plane sightlines (same ring as every earlier run)
MSUN=1.989e30
RMIN=RV                              # 14.3 Gpc visible radius

# ---------------- packing "quality" (fractions of each cell's own mass M_cell = rho*V_cell) ----------------
# Volume / mass shares chosen inside the accepted ranges (Cautun+2014 NEXUS+, SDSS void catalogues):
#   voids 78% vol / 23.4% mass (delta = -0.70), walls 6% / 8%, filaments 15.5% / 40%, clusters(knots) 0.5% / 28.6%.
PACK=dict(f_void=0.78, delta_void=-0.70, f_wall=0.06, m_wall=0.08, f_fil=0.155, m_fil=0.40, f_clu=0.005, m_clu=0.286)
def contrasts(p=PACK):
    """mass moved per cell (fraction of M_cell): void deficit (<0) and excess in walls / filaments / clusters (>0); sums to 0."""
    dv=p["f_void"]*p["delta_void"]                        # -0.546
    ew=p["m_wall"]-p["f_wall"]; ef=p["m_fil"]-p["f_fil"]; ec=p["m_clu"]-p["f_clu"]
    tot=dv+ew+ef+ec
    assert abs(tot)<1e-9, tot
    assert p["delta_void"]>=-1
    return dict(void=dv,wall=ew,fil=ef,clu=ec)

def rho_gpc(xg):                     # exp-lock density at CoE distance r, normalised to the mean density at D (kg/m^3)
    r=np.linalg.norm(xg-C_G,axis=-1); return RHO*np.exp(-(r-D_G)/L_G)

def in_region(xg,rmax):
    s=np.linalg.norm(xg,axis=-1); return (s>RMIN)&(s<=rmax)&(np.linalg.norm(xg-C_G,axis=-1)>=A_IN_G)

# ---------------- fast tidal pattern: identical formula to field3 (checked in check_kernel) ----------------
@nb.njit(parallel=True,fastmath=False,cache=True)
def _pattern(X,Y,Z,K,cth,sth):
    n=X.shape[0]; nt=cth.shape[0]; out=np.zeros(nt)
    for t in nb.prange(nt):
        x=cth[t]; y=sth[t]; acc=0.0
        for i in range(n):
            dx=X[i]-x; dy=Y[i]-y; dz=Z[i]
            r2=dx*dx+dy*dy+dz*dz; D2=X[i]*X[i]+Y[i]*Y[i]+Z[i]*Z[i]
            acc+=K[i]*((dx*x+dy*y)/(r2*np.sqrt(r2))-(X[i]*x+Y[i]*y)/(D2*np.sqrt(D2)))
        out[t]=acc
    return out-out.mean()

def pattern(xg,M):
    """xg: (n,3) positions in Gpc; M: kg (signed). Returns Delta T/T pattern on the 360-point ring (same units as R)."""
    P=np.ascontiguousarray(xg/RV)
    return _pattern(P[:,0].copy(),P[:,1].copy(),P[:,2].copy(),np.ascontiguousarray(M/KG),np.cos(TH),np.sin(TH))

def check_kernel(n=300,seed=0):
    rng=np.random.default_rng(seed); u=rng.normal(size=(n,3)); u/=np.linalg.norm(u,axis=1)[:,None]
    P=u*rng.uniform(1.02,3,n)[:,None]; M=rng.normal(size=n)*1e50
    Dd=np.linalg.norm(P,axis=1); ref=(M/KG)@lc.field3(np.arctan2(P[:,1],P[:,0]),np.arcsin(P[:,2]/Dd),Dd)
    w=pattern(P*RV,M); return float(np.max(np.abs(w-ref))/np.max(np.abs(ref)))

# ---------------- foam (Poisson-Voronoi web built from the Delaunay dual) ----------------
def _ball_shell(n,rng,r0,r1):
    u=rng.normal(size=(n,3)); u/=np.linalg.norm(u,axis=1)[:,None]
    return u*(rng.uniform(r0**3,r1**3,n)**(1/3))[:,None]

def _circumcenters(P,T):
    a=P[T[:,0]]; b=P[T[:,1]]-a; c=P[T[:,2]]-a; d=P[T[:,3]]-a
    A=np.stack([b,c,d],1); rhs=0.5*np.stack([(b*b).sum(1),(c*c).sum(1),(d*d).sum(1)],1)
    det=np.linalg.det(A); ok=np.abs(det)>1e-12*np.max(np.abs(det))
    x=np.zeros_like(a); x[ok]=np.linalg.solve(A[ok],rhs[ok][...,None])[...,0]
    return a+x, np.abs(det)/6.0, ok

def foam(s,rng,rmax=30.0,pack=PACK,return_parts=False,sel=None,gen=None,gap=0.0):
    """One Poisson-Voronoi foam realisation with mean cell spacing s (Gpc) = n_seed^(-1/3).
    Each kept cell (seed in region, whole cell outside the 14.3 Gpc sphere) is exactly mass-compensated:
      void  : -0.546 M_cell at the seed (cell centre); radius R_v = (3 f_void V_cell / 4pi)^(1/3)
      walls : +0.02 M_cell shared equally over the cell's Voronoi faces (Delaunay edges), placed at the face (seed-pair midpoint)
      fils  : +0.245 M_cell shared over the cell's Voronoi edges (Delaunay triangles), placed at the edge midpoint (between 2 vertices)
      clust : +0.281 M_cell shared over the cell's Voronoi vertices (Delaunay tetra circumcentres)
    Shared faces/edges/vertices collect from every neighbouring cell, so junction masses vary naturally.
    M_cell = rho_exp-lock(r_CoE) * V_cell, V_cell = barycentric dual volume (1/4 of each incident tetrahedron; partitions space)."""
    cc=contrasts(pack); pad=3.0*s
    r0=max(0.0,RMIN-pad); r1=rmax+pad
    Vsh=4/3*np.pi*(r1**3-r0**3); n=rng.poisson(Vsh/s**3)
    if gen is None:
        S=_ball_shell(n,rng,r0,r1)
        S=S[np.linalg.norm(S-C_G,axis=1)>=A_IN_G-pad]
    else: S=gen(s,rng,pad)          # custom seed generator (e.g. only a tube around the CMB ring), returns Gpc positions
    tri=Delaunay(S); T=tri.simplices; NB=tri.neighbors
    cen,vol,okt=_circumcenters(S,T)
    nS=len(S)
    Vc=np.bincount(T.ravel(),weights=np.repeat(vol/4,4),minlength=nS)
    # which cells are kept: seed in region, cell entirely outside the vision sphere (all its vertices and its void sphere)
    dcen=np.linalg.norm(cen,axis=1); dcen[~okt]=0.0
    minv=np.full(nS,np.inf); np.minimum.at(minv,T.ravel(),np.repeat(dcen,4))
    Rv=(3*pack["f_void"]*Vc/(4*np.pi))**(1/3)
    keep=in_region(S,rmax)&(minv>RMIN+gap)&(np.linalg.norm(S,axis=1)>=RMIN+gap+Rv)
    if sel is not None: keep&=sel(S)
    # cells whose tets are all well-formed only (hull cells never satisfy region test because of padding)
    Mcell=np.where(keep,rho_gpc(S)*Vc*GPC**3,0.0)
    # clusters on tets (Voronoi vertices)
    ntet=np.bincount(T.ravel(),minlength=nS)
    share=np.where(keep,cc["clu"]*Mcell/np.maximum(ntet,1),0.0)
    Mclu=share[T].sum(1)
    # filaments on triangles (Voronoi edges): faces of tets shared with a neighbour (count once: t<nb) -- interior only
    t_idx,j_idx=np.nonzero(NB>T.shape[0]*0-1)        # all (t,j) with neighbour >= 0
    nbj=NB[t_idx,j_idx]; once=t_idx<nbj; t_idx,j_idx,nbj=t_idx[once],j_idx[once],nbj[once]
    mask=np.ones((len(t_idx),4),bool); mask[np.arange(len(t_idx)),j_idx]=False
    TRI=T[t_idx][mask].reshape(-1,3)
    ntri=np.bincount(TRI.ravel(),minlength=nS)
    share=np.where(keep,cc["fil"]*Mcell/np.maximum(ntri,1),0.0)
    Mfil=share[TRI].sum(1); Pfil=0.5*(cen[t_idx]+cen[nbj])
    # walls on Delaunay edges (Voronoi faces)
    E=np.concatenate([T[:,[a,b]] for a in range(4) for b in range(a+1,4)]); E.sort(1)
    E=np.unique(E,axis=0); ne=np.bincount(E.ravel(),minlength=nS)
    share=np.where(keep,cc["wall"]*Mcell/np.maximum(ne,1),0.0)
    Mwal=share[E].sum(1); Pwal=0.5*(S[E[:,0]]+S[E[:,1]])
    Mvoid=cc["void"]*Mcell
    parts=dict(void=(S[keep],Mvoid[keep],Rv[keep]),clu=(cen[Mclu>0],Mclu[Mclu>0]),fil=(Pfil[Mfil>0],Mfil[Mfil>0]),wall=(Pwal[Mwal>0],Mwal[Mwal>0]))
    info=dict(s=s,n_cells=int(keep.sum()),n_clusters=int((Mclu>0).sum()),n_fil=int((Mfil>0).sum()),n_wall=int((Mwal>0).sum()),
              mean_Vcell=float(Vc[keep].mean()),net_mass_rel=float((Mvoid[keep].sum()+Mclu.sum()+Mfil.sum()+Mwal.sum())/np.abs(Mvoid[keep]).sum()),
              median_void_deficit_kg=float(np.median(-Mvoid[keep])) if keep.any() else 0.0,median_cluster_kg=float(np.median(Mclu[Mclu>0])) if (Mclu>0).any() else 0.0,
              median_void_R_eff_gpc=float(np.median(Rv[keep])) if keep.any() else 0.0,min_dist_gpc=float(min(np.linalg.norm(cen[Mclu>0],axis=1).min(initial=np.inf),np.linalg.norm(S[keep],axis=1).min(initial=np.inf))))
    X=np.concatenate([parts[k][0] for k in ("void","clu","fil","wall")]); M=np.concatenate([parts[k][1] for k in ("void","clu","fil","wall")])
    rad=np.concatenate([parts["void"][2],np.zeros(len(M)-len(parts["void"][2]))])
    if return_parts: return X,M,rad,info,parts
    return X,M,rad,info

def unstructured(M,rad,rng,rmax=30.0):
    """Control: the SAME masses (same voids, clusters, filament and wall pieces, net zero), but every piece placed independently
    at random in the region (no cells, no local compensation). Number density follows the gradient as in the surround runs:
    positive pieces ∝ rho(r_CoE), voids ∝ 1/rho. Voids keep their whole sphere outside 14.3 Gpc (centre >= 14.3 + R_v)."""
    n=len(M); out=np.empty((n,3)); need=np.arange(n)
    rmin=RMIN+rad; A_OUT=D_G+rmax
    while need.size:
        k=need.size; P=_ball_shell(int(k*1.6)+16,rng,RMIN,rmax)[:k] if False else None
        P=_ball_shell(k,rng,RMIN,rmax)
        sP=np.linalg.norm(P,axis=1); r=np.linalg.norm(P-C_G,axis=1)
        pos=M[need]>0
        acc=np.where(pos,np.exp(-(r-A_IN_G)/L_G),np.exp((r-A_OUT)/L_G))
        ok=(r>=A_IN_G)&(sP>=rmin[need])&(rng.uniform(size=k)<acc)
        out[need[ok]]=P[ok]; need=need[~ok]
    return out

def score(w,target):
    left=target-w; ww=float(np.dot(w,w)); sb=float(np.dot(target,w)/ww) if ww>0 else 0.0
    return rms(left)/FLOOR, rms(target-sb*w)/FLOOR, sb

# ---------------- tube around the CMB ring (where essentially all of the effect comes from at small cell sizes) ----------------
def ring_dist(X):
    """distance (Gpc) from the circle of in-plane sightline points at 14.3 Gpc (the ring the residual is measured on)"""
    return np.sqrt((np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2)

def tube_gen(T):
    """seed generator: Poisson seeds (density 1/s^3) in the torus of radius T+pad around the ring (CoE-interior part removed)."""
    def gen(s,rng,pad):
        a=T+pad; Vslab=np.pi*((RMIN+a)**2-max(RMIN-a,0)**2)*2*a; n=rng.poisson(Vslab/s**3)
        rc=np.sqrt(rng.uniform(max(RMIN-a,0)**2,(RMIN+a)**2,n)); ph=rng.uniform(0,2*np.pi,n); z=rng.uniform(-a,a,n)
        S=np.c_[rc*np.cos(ph),rc*np.sin(ph),z]
        return S[(ring_dist(S)<a)&(np.linalg.norm(S-C_G,axis=1)>=A_IN_G-pad)]
    return gen

def foam_tube(s,rng,T,**kw):
    return foam(s,rng,gen=tube_gen(T),sel=lambda S:ring_dist(S)<T,**kw)

def unstructured_tube(M,rad,rng,T,rmax=30.0,gap=0.0):
    """unstructured control inside the same tube volume: same pieces, independent random positions in region ∩ tube."""
    n=len(M); out=np.empty((n,3)); need=np.arange(n); rmin=RMIN+gap+rad; A_OUT=D_G+rmax
    while need.size:
        k=need.size
        rc=np.sqrt(rng.uniform(max(RMIN-T,0)**2,(RMIN+T)**2,k)); ph=rng.uniform(0,2*np.pi,k); z=rng.uniform(-T,T,k)
        P=np.c_[rc*np.cos(ph),rc*np.sin(ph),z]
        sP=np.linalg.norm(P,axis=1); r=np.linalg.norm(P-C_G,axis=1); pos=M[need]>0
        acc=np.where(pos,np.exp(-(r-A_IN_G)/L_G),np.exp((r-A_OUT)/L_G))
        ok=(ring_dist(P)<T)&(r>=A_IN_G)&(sP>=rmin[need])&(sP<=rmax)&(rng.uniform(size=k)<acc)
        out[need[ok]]=P[ok]; need=need[~ok]
    return out

def region_volume(rmax=30.0,n=4_000_000,seed=0):
    rng=np.random.default_rng(seed); P=_ball_shell(n,rng,RMIN,rmax); V=4/3*np.pi*(rmax**3-RMIN**3)
    return float(V*np.mean(np.linalg.norm(P-C_G,axis=1)>=A_IN_G))
