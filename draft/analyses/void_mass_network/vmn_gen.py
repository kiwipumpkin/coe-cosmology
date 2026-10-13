"""Random void+cluster network generator (bipartite void<->cluster links) following rules.md.
Construction (van de Weygaert & Icke 1989 Voronoi picture; Poisson-Voronoi geometry):
  1. void seeds: Poisson (or hard-core) points in a ball of radius R_out (+margin); number set so that the mean number of
     clusters bordering a void equals kv_mean (ours) given the Tinker08 cluster abundance above Mmin.
  2. void cells = Voronoi cells of the seeds; volume and centroid by Monte Carlo (nearest-seed assignment).
     Void core = fraction fV of the cell volume at density (1+delta_v) rho  (Cautun+14 NEXUS+: 77% vol, 1+delta=0.2).
     Drawn/void radius Rv = min( Req=(3 fV Vcell/4pi)^(1/3), d_ij Req_i/(Req_i+Req_j) over neighbours j )  -> void spheres never overlap.
  3. web nodes = Voronoi vertices (Delaunay circumcentres; each touches exactly 4 cells). Vertices closer than
     merge_lam x (mean seed spacing) are merged (union-find) -> nodes touching 4,5,6.. voids (target mean from Codis+18).
  4. clusters = the n_cl highest-scoring nodes, n_cl = n_Tinker(>Mmin) x V; score = (node empty-sphere radius)^1 x
     lognormal scatter (ours); masses by abundance matching to Tinker08 (rank k -> n(>M)V = k-1/2).
     Links: cluster <-> every void cell it touches (bipartite; no void-void, no cluster-cluster links).
  5. mass bookkeeping for the CMB kernel (point masses, kg, positions Gpc):
       void deficit  D_i = -(-delta_v) rho fV Vcell_i  at the cell centroid
       clusters      C_j at the node
       ridge (walls+filaments+smaller haloes) W_i = f_comp(Rv_i) |D_i| - sum_j C_j/k_c(j)  at the cell centroid
     f_comp modes: 'lcdm' (default, ours): net cell contrast drawn N(0, sigma_lin(Reff)) so each cell is compensated up to
     a LCDM-sized fluctuation; 'hamaus': f = 1 - s ln(Reff/R_comp) (Hamaus+14: compensated at R_comp=17.6 Mpc/h, bigger
     voids under-, smaller over-compensated; slope s and its extrapolation to these scales are ours); 'exact': f=1.
Units: Mpc inside the generator; exported X in Gpc, M in kg."""
import numpy as np, lss
from scipy.spatial import Delaunay, cKDTree
MPC_M=3.0857e22; MSUN=1.989e30
RHO_MSUN=lss.RHO_M            # Msun / Mpc^3
DEF=dict(R_out=20000.0, margin=1500.0, Mmin=3e15, kv_mean=8.0, kc_target=4.5, merge_lam=None,
         fV=0.77, dv=-0.8, comp_mode="lcdm", s_comp=0.05, R_comp=17.6/0.70, score_scatter=0.5, hardcore=0.7,
         gap_us=500.0, mc_per_cell=60, seed=0)
LAM_DEFAULT=0.15      # calibrated in check_stats.py so that <voids per cluster> ~ 4.5 (see rules.md)
def _seeds(rng,n,R,hard):
    pts=[]; need=n
    while need>0:
        u=rng.normal(size=(int(need*1.3)+10,3)); u/=np.linalg.norm(u,axis=1)[:,None]
        p=u*(R*rng.uniform(size=len(u))**(1/3))[:,None]; pts.append(p); need-=len(p)
    P=np.vstack(pts)[:n]
    if hard>0:   # crude hard-core thinning + refill (ours, optional)
        d0=hard*(4/3*np.pi*R**3/n)**(1/3)
        for _ in range(30):
            t=cKDTree(P); pr=t.query_pairs(d0,output_type="ndarray")
            if len(pr)==0: break
            bad=np.unique(pr[:,1]); u=rng.normal(size=(len(bad),3)); u/=np.linalg.norm(u,axis=1)[:,None]
            P[bad]=u*(R*rng.uniform(size=len(bad))**(1/3))[:,None]
    return P
def _circum(P,simp):
    A=P[simp[:,0]]; B=P[simp[:,1]]-A; C=P[simp[:,2]]-A; D=P[simp[:,3]]-A
    M=np.stack([B,C,D],1); rhs=0.5*np.stack([(B*B).sum(1),(C*C).sum(1),(D*D).sum(1)],1)
    ok=np.abs(np.linalg.det(M))>1e-9*np.linalg.norm(B,axis=1)**3
    x=np.zeros_like(B); x[ok]=np.linalg.solve(M[ok],rhs[ok][...,None])[...,0]
    c=A+x; return c,np.linalg.norm(x,axis=1),ok
def _uf_merge(n,pairs):
    par=np.arange(n)
    def f(i):
        while par[i]!=i: par[i]=par[par[i]]; i=par[i]
        return i
    for a,b in pairs:
        ra,rb=f(a),f(b)
        if ra!=rb: par[max(ra,rb)]=min(ra,rb)
    return np.array([f(i) for i in range(n)])
def generate(**kw):
    p=dict(DEF); p.update(kw); rng=np.random.default_rng(p["seed"]); lam=p["merge_lam"] if p["merge_lam"] is not None else LAM_DEFAULT
    R=p["R_out"]; Rg=R+p["margin"]; Vin=4/3*np.pi*R**3
    n_cl_dens=lss.cluster_n_above(p["Mmin"]); nv_dens=n_cl_dens*p["kc_target"]/p["kv_mean"]
    nseed=rng.poisson(nv_dens*4/3*np.pi*Rg**3); S=_seeds(rng,nseed,Rg,p["hardcore"]); spacing=nv_dens**(-1/3)
    # cells: MC volumes and centroids
    tree=cKDTree(S); nmc=int(p["mc_per_cell"]*nseed); u=rng.normal(size=(nmc,3)); u/=np.linalg.norm(u,axis=1)[:,None]
    Q=u*(Rg*rng.uniform(size=nmc)**(1/3))[:,None]; _,own=tree.query(Q,workers=-1)
    cnt=np.bincount(own,minlength=nseed); Vcell=cnt*(4/3*np.pi*Rg**3/nmc)
    cen=np.zeros((nseed,3))
    for k in range(3): cen[:,k]=np.bincount(own,weights=Q[:,k],minlength=nseed)
    cen=np.where(cnt[:,None]>0,cen/np.maximum(cnt,1)[:,None],S)
    Req=(3*p["fV"]*Vcell/(4*np.pi))**(1/3)
    dk,ik=tree.query(S,k=16,workers=-1); dk=dk[:,1:]; ik=ik[:,1:]
    # non-overlap: each neighbour pair splits its separation in proportion to the two equivalent radii
    Rv=np.minimum(Req,(dk*Req[:,None]/(Req[:,None]+Req[ik])).min(1))
    # nodes
    tri=Delaunay(S); simp=tri.simplices; cc,rad,ok=_circum(S,simp)
    keep=ok&(np.linalg.norm(cc,axis=1)<Rg); newid=np.full(len(keep),-1); newid[keep]=np.arange(keep.sum())
    nb=tri.neighbors; ii=np.repeat(np.arange(len(nb)),4); jj=nb.ravel(); m=(jj>ii)&(jj>=0)
    ii,jj=newid[ii[m]],newid[jj[m]]; m=(ii>=0)&(jj>=0); tet_edges=np.c_[ii[m],jj[m]]      # foam edges (Voronoi edges)
    e6=np.vstack([simp[:,[a,b]] for a in range(4) for b in range(a+1,4)]); e6.sort(1); vv_edges=np.unique(e6,axis=0)  # shared faces
    simp=simp[keep]; cc=cc[keep]; rad=rad[keep]
    if lam>0:
        t2=cKDTree(cc); pr=t2.query_pairs(lam*spacing,output_type="ndarray"); lab=_uf_merge(len(cc),pr)
    else: lab=np.arange(len(cc))
    ul,inv=np.unique(lab,return_inverse=True); nn=len(ul)
    fe=inv[tet_edges]; fe=fe[fe[:,0]!=fe[:,1]]; fe.sort(1); foam_edges=np.unique(fe,axis=0)          # node-node edges after merging
    npos=np.zeros((nn,3))
    for k in range(3): npos[:,k]=np.bincount(inv,weights=cc[:,k],minlength=nn)
    npos/=np.bincount(inv,minlength=nn)[:,None]
    nrad=np.zeros(nn); np.maximum.at(nrad,inv,rad)
    # node -> void sets
    pairs=np.unique(np.c_[np.repeat(inv,4),simp.ravel()],axis=0)      # (node, void)
    # inside region only (node inside R, voids may be in margin)
    inR=np.linalg.norm(npos,axis=1)<R
    n_cl=rng.poisson(n_cl_dens*Vin)
    score=np.log(nrad)+p["score_scatter"]*rng.normal(size=nn); score[~inR]=-np.inf
    order=np.argsort(-score)[:n_cl]
    Mcl=lss.cluster_masses_ranked(n_cl,p["Mmin"],Vin)                  # descending
    cl_node=order; XC=npos[cl_node]; MC=Mcl
    isC=np.full(nn,-1); isC[cl_node]=np.arange(n_cl)
    lk=pairs[isC[pairs[:,0]]>=0]; links=np.c_[isC[lk[:,0]],lk[:,1]]     # (cluster index, void index)
    kc=np.bincount(links[:,0],minlength=n_cl)
    vin=np.linalg.norm(S,axis=1)<R
    kv=np.bincount(links[:,1],minlength=nseed)
    # mass bookkeeping (Msun)
    D=(-p["dv"])*RHO_MSUN*p["fV"]*Vcell
    Reff=(3*Vcell/(4*np.pi))**(1/3)
    if p["comp_mode"]=="hamaus":      # literature trend (Hamaus+14) extrapolated with slope s_comp (ours)
        fcomp=np.clip(1-p["s_comp"]*np.log(np.maximum(Reff,1e-3)/p["R_comp"]),0,1.2)
    elif p["comp_mode"]=="exact": fcomp=np.ones(nseed)
    else:   # 'lcdm' (ours): net cell contrast ~ N(0, sigma_lin(Reff)) -> f = 1 + delta_cell*Vcell/|D|
        dcell=lss.sigma(Reff)*rng.normal(size=nseed); fcomp=1+dcell*RHO_MSUN*Vcell/D
    share=np.bincount(links[:,1],weights=MC[links[:,0]]/kc[links[:,0]],minlength=nseed)
    W=fcomp*D-share
    out=dict(params=p,lam=lam,spacing=spacing,S=S,cen=cen,Vcell=Vcell,Rv=Rv,Reff=(3*Vcell/(4*np.pi))**(1/3),void_in=vin,kv=kv,
             XC=XC,MC=MC,kc=kc,links=links,node_rad=nrad[cl_node],D=D,W=W,fcomp=fcomp,share=share,
             n_nodes_in=int(inR.sum()),npos=npos,node_in=inR,node_cl=isC,foam_edges=foam_edges,vv_edges=vv_edges,node_void_pairs=pairs)
    return out
def kernel_bodies(g,R_lim=None):
    """point masses for the CMB kernel: positions Gpc, masses kg. Cells whose seed is outside R_out are dropped,
    and every body within gap_us of us is dropped (same 0.5 Gpc exclusion as the joint fits)."""
    p=g["params"]; R_lim=R_lim or p["R_out"]; v=np.linalg.norm(g["S"],axis=1)<R_lim
    X=np.vstack([g["cen"][v],g["cen"][v],g["XC"]]); M=np.concatenate([-g["D"][v],g["W"][v],g["MC"]])
    far=np.linalg.norm(X,axis=1)>p["gap_us"]; return X[far]/1000.0,M[far]*MSUN
