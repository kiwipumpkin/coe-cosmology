"""Non-crossing ("planar-in-3D") link set, per Chris's rule:
  * links are tubes of radius r_t (ours: rt_frac x median void radius); two links CROSS if their axes come closer than
    2 r_t, unless they share an end node;
  * a link must not pass through another node's body: void core sphere (radius = drawn non-overlapping core Rv) or
    cluster body (radius R200m); its own two end nodes are exempt;
  * cluster-cluster links must not pass through any void core (covered by the rule above, since cc ends are clusters);
  * greedy: candidates sorted shortest first; a candidate that crosses an already-kept link is dropped
    (= of two crossing links the longer one goes); a candidate through a body is dropped outright.
Candidates: cc (direct foam edges), cc2 (cluster-cluster via one empty web node, straight), cv (bipartite), vv (shared walls).
Works in a sub-ball (R_sub) to keep the pair search tractable; statistics are read for nodes inside R_count."""
import numpy as np, numba
from scipy.spatial import cKDTree
import lss
@numba.njit(cache=True)
def segseg(p1,q1,p2,q2):
    d1=q1-p1; d2=q2-p2; r=p1-p2; a=d1@d1; e=d2@d2; f=d2@r; c=d1@r; b=d1@d2; den=a*e-b*b
    s=0.0
    if den>1e-12: s=min(max((b*f-c*e)/den,0.0),1.0)
    t=(b*s+f)/e
    if t<0.0: t=0.0; s=min(max(-c/a,0.0),1.0)
    elif t>1.0: t=1.0; s=min(max((b-c)/a,0.0),1.0)
    dv=p1+d1*s-(p2+d2*t); return np.sqrt(dv@dv)
@numba.njit(cache=True)
def greedy(P0,P1,EA,EB,order,nb_ptr,nb_idx,ok_body,thr):
    n=len(P0); kept=np.zeros(n,np.bool_); why=np.zeros(n,np.int8)   # 0 kept,1 body,2 crossing
    for ii in range(n):
        i=order[ii]
        if not ok_body[i]: why[i]=1; continue
        bad=False
        for k in range(nb_ptr[i],nb_ptr[i+1]):
            j=nb_idx[k]
            if not kept[j]: continue
            if EA[i]==EA[j] or EA[i]==EB[j] or EB[i]==EA[j] or EB[i]==EB[j]: continue
            if segseg(P0[i],P1[i],P0[j],P1[j])<thr: bad=True; break
        if bad: why[i]=2
        else: kept[i]=True
    return kept,why
@numba.njit(cache=True)
def body_check(P0,P1,EA,EB,b_ptr,b_idx,BX,BR):
    n=len(P0); ok=np.ones(n,np.bool_)
    for i in range(n):
        d=P1[i]-P0[i]; L2=d@d
        for k in range(b_ptr[i],b_ptr[i+1]):
            j=b_idx[k]
            if j==EA[i] or j==EB[i]: continue
            t=0.0
            if L2>0: t=min(max((BX[j]-P0[i])@d/L2,0.0),1.0)
            v=P0[i]+d*t-BX[j]
            if np.sqrt(v@v)<BR[j]: ok[i]=False; break
    return ok
def _csr(lists):
    ptr=np.zeros(len(lists)+1,np.int64); ptr[1:]=np.cumsum([len(l) for l in lists]); idx=np.fromiter((x for l in lists for x in l),np.int64,count=ptr[-1]); return ptr,idx
def candidates(g,cc,cc2,R_sub):
    S=g["S"]; XC=g["XC"]; ncl=len(XC); cv=g["links"]; vv=g["vv_edges"]
    pos=np.vstack([XC,S]); body_R=np.r_[lss.R_of_M(g["MC"])*(1/200)**(1/3),g["Rv"]]     # R200m for clusters, core sphere for voids
    E=[]; T=[]
    for t,(a,b) in ((0,(cc[:,0],cc[:,1])),(1,(cc2[:,0],cc2[:,1])),(2,(cv[:,0],ncl+cv[:,1])),(3,(ncl+vv[:,0],ncl+vv[:,1]))):
        E.append(np.c_[a,b]); T.append(np.full(len(a),t))
    E=np.vstack(E); T=np.concatenate(T)
    inside=(np.linalg.norm(pos[E[:,0]],axis=1)<R_sub)&(np.linalg.norm(pos[E[:,1]],axis=1)<R_sub)
    return E[inside],T[inside],pos,body_R
def planarize(g,cc,cc2,R_sub=6500.0,rt_frac=0.05):
    E,T,pos,BRall=candidates(g,cc,cc2,R_sub); P0=pos[E[:,0]].copy(); P1=pos[E[:,1]].copy(); L=np.linalg.norm(P1-P0,axis=1); mid=0.5*(P0+P1)
    rt=rt_frac*np.median(g["Reff"][g["void_in"]]); thr=2*rt
    # bodies near each segment
    bsel=np.where(np.linalg.norm(pos,axis=1)<R_sub+500)[0]; bt=cKDTree(pos[bsel]); Rmax=BRall[bsel].max()
    bl=bt.query_ball_point(mid,L/2+Rmax,workers=-1); bl=[bsel[np.array(x,int)] for x in bl]; bptr,bidx=_csr(bl)
    okb=body_check(P0,P1,E[:,0],E[:,1],bptr,bidx,pos,BRall)
    # neighbouring segments
    mt=cKDTree(mid); nl=mt.query_ball_point(mid,L+thr,workers=-1); nptr,nidx=_csr(nl)   # kept partners are never longer (shortest-first)
    order=np.argsort(L,kind="stable")
    kept,why=greedy(P0,P1,E[:,0],E[:,1],order,nptr,nidx,okb,thr)
    return dict(E=E,T=T,L=L,kept=kept,why=why,rt=rt,pos=pos)
