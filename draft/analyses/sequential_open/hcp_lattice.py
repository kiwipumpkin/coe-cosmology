"""HCP lattice analogue (Chris, 11 Oct). Same conventions as fcc_compare.py / fcc_ladder.py:
clusters on HCP sites (ABAB stacking), nearest-neighbour spacing d = 171/sqrt(2) = 120.9 Mpc (same as the FCC version, cubic cell 171 Mpc),
c = sqrt(8/3) d. Voids in tetrahedral holes (4 clusters each) or, as a variant, octahedral holes (6 clusters each).
Masses: cluster = foam s2 median 1.40e50 kg, void = median deficit -8.39e50 kg (negative mass); 'web ratio' = voids in 1/12 of tetra holes
(1 void per 6 clusters; 6*1.40-8.39 ~ 0, so this mix is already almost exactly compensated).
Outputs hcp_lattice/: hcp_ladder.json, hcp_ladder_raw.npz, hcp_shells.json, hcp_ewald.json. Figures: hcp_figs.py."""
import json,numpy as np
from scipy.spatial import Delaunay,cKDTree
from scipy.spatial.transform import Rotation
from scipy.special import erfc
exec(open("fcc_ladder.py").read().split("RES={}; RAW={}; L0")[0])   # analyse(), run(), poly(), jitter(), poisson(), crystal() (FCC), MC, MV, A, LAY
crystal_fcc=crystal; OD="hcp_lattice/"
# ---------------- HCP generator (fractional hex coords)
HB=np.array([[0,0,0],[1/3,2/3,.5]]); HT=np.array([[0,0,3/8],[0,0,5/8],[1/3,2/3,1/8],[1/3,2/3,7/8]]); HO=np.array([[2/3,1/3,.25],[2/3,1/3,.75]])
def hex_cart(F,d): c=np.sqrt(8/3)*d; return np.c_[d*(F[:,0]-F[:,1]/2),d*np.sqrt(3)/2*F[:,1],c*F[:,2]]
def hcp_sites(lo,hi,d):
    """all HCP sites / tetra holes / octa holes inside the box lo..hi (cartesian, Gpc)"""
    c=np.sqrt(8/3)*d; h=np.sqrt(3)/2*d
    v=np.arange(np.floor(lo[1]/h)-1,np.ceil(hi[1]/h)+2); w=np.arange(np.floor(lo[2]/c)-1,np.ceil(hi[2]/c)+2)
    u=np.arange(np.floor(lo[0]/d+v.min()/2)-1,np.ceil(hi[0]/d+v.max()/2)+2)
    g=np.stack(np.meshgrid(u,v,w,indexing="ij"),-1).reshape(-1,3)
    out=[]
    for Bs in (HB,HT,HO):
        X=hex_cart((g[:,None,:]+Bs[None]).reshape(-1,3),d); m=np.all((X>=lo)&(X<hi),axis=1); out.append(X[m])
    return out
def crystal_hcp(L,a,fill,rot=None,shift=None,vac=0.,inter=0.,octa=0.,hole="tetra"):
    """same interface as fcc crystal(): a = equivalent FCC cube edge (d = a/sqrt2); kind 0 cluster, 1 void.
    hole='tetra': voids in a fraction `fill` of tetra holes; hole='octa': fraction `fill` of octa holes."""
    d=a/np.sqrt(2)
    if rot is None: lo,hi=np.full(3,-d),np.full(3,L+d)
    else: R=np.sqrt(3)*L+2*a; lo,hi=np.full(3,-R),np.full(3,R)
    c,t,o=hcp_sites(lo,hi,d)
    if rot is not None:
        k=lambda X: X[np.linalg.norm(X,axis=1)<R]; c,t,o=k(c),k(t),k(o)
    rot=Rotation.identity() if rot is None else rot; sh=np.zeros(3) if shift is None else shift
    H=t if hole=="tetra" else o
    tv=H[rng.random(len(H))<fill]
    if octa>0:
        kk=rng.random(len(tv))<octa; tv=np.vstack([tv[~kk],o[rng.choice(len(o),kk.sum(),replace=False)]])
    c=c[rng.random(len(c))>=vac]
    P=np.vstack([c,tv]); K=np.r_[np.zeros(len(c),int),np.ones(len(tv),int)]
    P=rot.apply(P)+sh; m=np.all((P>=0)&(P<L),axis=1); P,K=P[m],K[m]
    if inter>0:
        ni=int(inter*(K==0).sum()); P=np.vstack([P,rng.uniform(0,L,(ni,3))]); K=np.r_[K,np.zeros(ni,int)]
    return P,K
if __name__=="__main__":
    d0=A/np.sqrt(2); print(f"d = {d0*1e3:.1f} Mpc, c = {np.sqrt(8/3)*d0*1e3:.1f} Mpc")
    # ---- sanity: hole coordination
    c,t,o=hcp_sites(np.full(3,-.6),np.full(3,.6),d0); TK=cKDTree(c)
    for nm,H,k in (("tetra",t,4),("octa",o,6)):
        dd,_=TK.query(H[np.linalg.norm(H,axis=1)<.3],k+1); print(nm,"hole: nearest clusters /d =",np.round(dd[0]/d0,4))
    # ---- coordination shells (exact, all pairs) for clusters, cluster->void and void->void, HCP vs FCC
    SH={}
    def shells(X,Y,ctr,rmax=2.3,self_=True):
        TY=cKDTree(Y); i=np.argsort(np.linalg.norm(X-ctr,axis=1))[:40]; out={}
        for j in i:
            dd,_=TY.query(X[j],200,distance_upper_bound=rmax*d0+1e-9); dd=dd[np.isfinite(dd)]; dd=dd[dd>1e-9]/d0
            r=np.round(dd,3); u_,n_=np.unique(r,return_counts=True)
            for a_,b_ in zip(u_,n_): out.setdefault(float(a_),[]).append(int(b_))
        return {k:float(np.median(v)) for k,v in sorted(out.items())}
    Lb=1.2
    for nm in ("HCP","FCC"):
        if nm=="HCP": c,t,o=hcp_sites(np.zeros(3),np.full(3,Lb),d0)
        else: c,t,o=sites(-1,int(Lb/A)+2,A); kk=lambda X:X[np.all((X>=0)&(X<Lb),axis=1)]; c,t,o=kk(c),kk(t),kk(o)
        ctr=np.full(3,Lb/2)
        SH[nm]={"cluster->cluster":shells(c,c,ctr),"tetra void->cluster":shells(t,c,ctr,1.2),"octa void->cluster":shells(o,c,ctr,1.2),
                "tetra void->tetra void":shells(t,t,ctr,1.5),"octa void->octa void":shells(o,o,ctr,1.5),
                "per cluster: tetra holes, octa holes":[len(t)/len(c),len(o)/len(c)]}
        print(nm,json.dumps(SH[nm]))
    json.dump(SH,open(OD+"hcp_shells.json","w"),indent=1)
    # ---- ladder (same recipe and box sizes as fcc_ladder.py)
    RES={}; RAW={}; L0=14*A; MG=2*A
    for comp,fill,hole in (("web ratio",1/12,"tetra"),("web ratio octa",1/6,"octa"),("fluorite",1.0,"tetra")):
        P0,K0=crystal_hcp(L0,A,fill,hole=hole); dsp=spacing(P0,L0); C={}
        C["perfect"]=(P0,K0)
        if comp=="web ratio":
            for f in (0.05,0.15,0.30): C[f"jitter {int(f*100)}%"]=(jitter(P0,f,dsp),K0)
            for f in (0.03,0.10): C[f"defects {int(f*100)}%"]=crystal_hcp(L0,A,fill,vac=f,inter=f,octa=f)
            crystal=crystal_hcp; Lp=18*A; Pp,Kp,Gp,seeds=poly(Lp,A,fill); crystal=crystal_fcc; C["polycrystal"]=(Pp,Kp)
            np.savez(OD+"hcp_polycrystal_webratio.npz",P=Pp,K=Kp,G=Gp,seeds=seeds,L=Lp)
        for nm,(P,K) in C.items():
            L=Lp if nm=="polycrystal" else L0; key=f"HCP {comp} | {nm}"
            RES[key],RAW[key]=run(P,K,L,key,MG); o_=RES[key]
            print(f"{key:34s} N={o_['n_nodes']:6d} Shist={o_['S_hist_logvol']:.3f} cv={o_['cv_vol']:.3f} u={o_['u_dimless']:+8.2f} deg={o_['mean_degree']:.2f} L/d={o_['mean_link_over_spacing']:.3f}",flush=True)
    # HCP and FCC clusters only (Delaunay)
    for nm,fn in (("HCP clusters only",crystal_hcp),("FCC clusters only",crystal_fcc)):
        P,K=fn(L0,A,0.0); RES[nm],RAW[nm]=run(P,K,L0,nm,MG); print(nm,RES[nm]["S_hist_logvol"],RES[nm]["u_dimless"],RES[nm]["mean_degree"])


    # compensated web seed 3 (BEFORE = random balanced start, AFTER = entropy-rule edits reaching 0.998x)
    dC=json.load(open("minedit_comp/comp_s3_entropy.json")); bC=dC["tag"].replace("_entropy","")
    RIN=28.0; Vin=4/3*np.pi*RIN**3
    for nm,lay in (("comp web s3 BEFORE",EB.layout(bC)),("comp web s3 AFTER",EB.layout(bC,dC))):
        X,r,v,M,_=lay; inner=lambda Q: np.linalg.norm(Q,axis=1)<RIN
        RES[nm],RAW[nm]=analyse(X,M,inner,nm,Vin); E=RAW[nm]["E"]; deg=np.bincount(E.ravel(),minlength=len(X)); ii=np.linalg.norm(X,axis=1)<RIN-3
        RAW[nm]["deg"]=deg[ii]; RES[nm]["mean_degree"]=float(deg[ii].mean()); RES[nm]["n_voids_per_cluster"]=float(v.sum()/(~v).sum())
        print(nm,RES[nm]["S_hist_logvol"],RES[nm]["u_dimless"],RES[nm]["mean_degree"])
    json.dump(RES,open(OD+"hcp_ladder.json","w"),indent=1,default=float)
    np.savez(OD+"hcp_ladder_raw.npz",**{f"{k}|{q}":RAW[k][q] for k in RAW for q in ("vol","L","deg")})
