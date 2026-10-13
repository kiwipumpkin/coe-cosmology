"""Imperfect-crystal ladder (Chris, 10 Oct): perfect FCC (clusters) + voids in tetrahedral holes, thermal jitter, point defects
(vacancies, interstitials, voids in octahedral holes), polycrystal (Voronoi grains, random orientation, spacing +-5%, disordered
boundaries), fully random (Poisson). Two compositions: 'fluorite' (voids in ALL tetra holes, 2 voids/cluster) and 'web ratio'
(voids in 1/12 of tetra holes, 1 void/6 clusters like foam s2). Same 3D Delaunay statistics as fcc_compare.py. Outputs fcc_lattice/."""
import json,numpy as np
from scipy.spatial import Delaunay
from scipy.spatial.transform import Rotation
exec(open("fcc_compare.py").read().split("RES={}; RAW={}")[0])   # analyse(), fcc(), masses MC MV, A, LAY
rng=np.random.default_rng(7)
B=np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]]); T=np.array([[i,j,k] for i in (.25,.75) for j in (.25,.75) for k in (.25,.75)])
O=np.array([[.5,.5,.5],[.5,0,0],[0,.5,0],[0,0,.5]])
def sites(n0,n1,a):
    g=np.array([[i,j,k] for i in range(n0,n1) for j in range(n0,n1) for k in range(n0,n1)],float)
    f=lambda U:(g[:,None,:]+U[None]).reshape(-1,3)*a
    return f(B),f(T),f(O)
def crystal(L,a,fill,rot=None,shift=None,vac=0.,inter=0.,octa=0.):
    """lattice points (kind 0 cluster, 1 void) covering the cube [0,L]^3; fill = fraction of tetra holes holding a void"""
    n=int(np.ceil(L*1.8/a))+1; c,t,o=sites(-n,n,a)
    rot=Rotation.identity() if rot is None else rot; sh=np.zeros(3) if shift is None else shift
    tv=t[rng.random(len(t))<fill]
    if octa>0:   # move a fraction of the voids from tetra to (random) octahedral holes
        k=rng.random(len(tv))<octa; tv=np.vstack([tv[~k],o[rng.choice(len(o),k.sum(),replace=False)]])
    c=c[rng.random(len(c))>=vac]
    P=np.vstack([c,tv]); K=np.r_[np.zeros(len(c),int),np.ones(len(tv),int)]
    P=rot.apply(P)+sh; m=np.all((P>=0)&(P<L),axis=1); P,K=P[m],K[m]
    if inter>0:   # extra clusters at random positions
        ni=int(inter*(K==0).sum()); P=np.vstack([P,rng.uniform(0,L,(ni,3))]); K=np.r_[K,np.zeros(ni,int)]
    return P,K
def spacing(P,L): return (L**3/len(P))**(1/3)
def jitter(P,f,d): return P+rng.normal(0,f*d/np.sqrt(3),P.shape)
def poisson(L,n,fv): P=rng.uniform(0,L,(n,3)); return P,(rng.random(n)<fv).astype(int)
def poly(L,a,fill,ng=10,w=0.6,dsp=0.05,jit=0.05):
    seeds=rng.uniform(0,L,(ng,3)); P_=[];K_=[];G_=[]
    for g in range(ng):
        ag=a*(1+rng.uniform(-dsp,dsp)); P,K=crystal(L,ag,fill,rot=Rotation.random(random_state=int(rng.integers(1e9))),shift=rng.uniform(0,ag,3)+seeds[g]-seeds[g])
        P=jitter(P,jit,ag/12**(1/3)); m=np.all((P>=0)&(P<L),axis=1); P,K=P[m],K[m]
        dd=np.linalg.norm(P[:,None,:]-seeds[None],axis=2); o=np.argsort(dd,1); d1=dd[np.arange(len(P)),o[:,0]]; d2=dd[np.arange(len(P)),o[:,1]]
        keep=(o[:,0]==g)&((d2-d1)/2>w*a/2); P_.append(P[keep]); K_.append(K[keep]); G_.append(np.full(keep.sum(),g))
    P=np.vstack(P_); K=np.concatenate(K_); Gr=np.concatenate(G_)
    # disordered grain boundaries: Poisson points at the mean crystal density in the boundary zone
    Q=rng.uniform(0,L,(400000,3)); dd=np.linalg.norm(Q[:,None,:]-seeds[None],axis=2); dd.sort(1); bz=(dd[:,1]-dd[:,0])/2<=w*a/2
    vb=L**3*bz.mean(); rho=len(P)/(L**3-vb); nb=int(rho*vb); Qb=Q[bz][rng.choice(bz.sum(),nb,replace=False)]
    fv=(K==1).mean(); P=np.vstack([P,Qb]); K=np.r_[K,(rng.random(nb)<fv).astype(int)]; Gr=np.r_[Gr,np.full(nb,-1)]
    return P,K,Gr,seeds
def run(P,K,L,name,marg):
    M=np.where(K==0,MC,MV); inner=lambda X: np.all((X>=marg)&(X<L-marg),axis=1)
    r,raw=analyse(P,M,inner,name,(L-2*marg)**3)
    deg=np.bincount(raw["E"].ravel(),minlength=len(P)); ii=np.all((P>=marg+0.2)&(P<L-marg-0.2),axis=1); raw["deg"]=deg[ii]
    r["mean_degree"]=float(deg[ii].mean()); r["n_voids_per_cluster"]=float((K==1).sum()/max((K==0).sum(),1)); return r,raw
RES={}; RAW={}; L0=14*A; MG=2*A
for comp,fill in (("fluorite",1.0),("web ratio",1/12)):
    P0,K0=crystal(L0,A,fill); d0=spacing(P0,L0); C={}
    C["perfect"]=(P0,K0)
    for f in (0.05,0.15,0.30): C[f"jitter {int(f*100)}%"]=(jitter(P0,f,d0),K0)
    for f in (0.03,0.10): C[f"defects {int(f*100)}%"]=crystal(L0,A,fill,vac=f,inter=f,octa=f)
    Lp=18*A; Pp,Kp,Gp,seeds=poly(Lp,A,fill); C["polycrystal"]=(Pp,Kp)
    C["random"]=poisson(L0,len(P0),(K0==1).mean())
    if comp=="fluorite": np.savez("fcc_lattice/polycrystal_fluorite.npz",P=Pp,K=Kp,G=Gp,seeds=seeds,L=Lp)
    for nm,(P,K) in C.items():
        L=Lp if nm=="polycrystal" else L0; key=f"{comp} | {nm}"
        RES[key],RAW[key]=run(P,K,L,key,MG); o=RES[key]
        print(f"{key:28s} N={o['n_nodes']:6d} tets={o['n_tets']:6d} Shist={o['S_hist_logvol']:.3f} even={o['evenness']:.4f} cv={o['cv_vol']:.3f} u={o['u_dimless']:+8.2f} deg={o['mean_degree']:.2f} L/d={o['mean_link_over_spacing']:.3f}")
W=json.load(open("fcc_lattice/fcc_vs_web.json"))
for k in ("web BEFORE","web AFTER"): RES[k]=W[k]
Wr=np.load("fcc_lattice/fcc_vs_web_raw.npz")
for k in ("web BEFORE","web AFTER"): RAW[k]={q:Wr[f"{k}|{q}"] for q in ("vol","L","deg")}
json.dump(RES,open("fcc_lattice/ladder.json","w"),indent=1,default=float)
np.savez("fcc_lattice/ladder_raw.npz",**{f"{k}|{q}":RAW[k][q] for k in RAW for q in ("vol","L","deg")})
