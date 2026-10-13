"""'Minimal edit' test: start from a randomly grown web (our sphere + shell to 30 Gpc incl. core cap), no new bodies.
Find the fewest EXISTING voids/clusters to move and/or resize so the ring leftover (with SR + central mass) drops to <= 1x.
Per edited body: new position (x,y,z) and size. Voids: radius 0.5-4 Gpc (delta -0.7 of the local exp-lock density).
  Clusters: mass within a factor 10 of its grown mass (radius scales as M^1/3, same density contrast).
Rules for an edited body (penalties driven to zero and checked): inside 30 Gpc, >= 0.5 Gpc from us, surface >= 3 Gpc from the
  CMB ring, no overlap (0.02 Gpc gap) with any other body (40 nearest checked) or other edited bodies. Unedited bodies stay as grown.
Kernel: the growth runs' own seq_core.unit_pattern_fast (pattern = sum M K(X)); net mass is allowed to change and is reported.
Method: (1) leverage of every body = best linearised reduction of chi^2 from changing that body alone (Jacobian of its pattern
  w.r.t. position and size, ridge-damped, step capped at 5 Gpc); (2) take the top few, (3) fit them jointly (bounded
  nonlinear least squares, with the bodies already edited re-fitted together every few rounds), (4) repeat until <= 1x.
usage: python3 minedit.py TAG   -> minedit/<TAG>.json (+ ckpt)"""
import os,sys,json,time,pickle,numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from common import *
from fit_on_web import spec
K=seq_core.unit_pattern_fast; GPCm=seq_core.GPC; RT=30.0; CLR=3.0; GAP=0.02; MOVE=5.0
tag=sys.argv[1]; NEAR=os.environ.get("NEAR","0")=="1"; SF="_nr" if NEAR else ""; OUTF=f"{OUT}/minedit/{tag}{SF}.json"; CK=f"{OUT}/ckpt/minedit_{tag}{SF}.pkl"; KMAX=int(os.environ.get("KMAX",600))
f=np.load(f"{OUT}/final_{tag}.npz"); X0=f["Xp"].astype(float); r0=f["rp"].astype(float); isV=(f["kp"]=="void"); M0=f["Mp"].astype(float); n=len(X0)
w0=f["w"]; T=R-w0                       # leftover before any edit
P0=K(X0)*1.0                            # (n,360) unit patterns
def dring(X): X=np.atleast_2d(X); return np.sqrt((np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2)
def mass_of(i,X,s):     # s: void radius (Gpc) or log10(M/M0) for clusters
    rho=seq_core.rho_gpc(np.atleast_2d(X))*GPCm**3
    Mv=-0.70*rho*4/3*np.pi*s**3; Mc=M0[i]*10**s
    return np.where(isV[i],Mv,Mc), np.where(isV[i],s,r0[i]*10**(s/3))
def s0(i): return np.where(isV[i],r0[i],0.0)
def slo(i): return np.where(isV[i],0.5,-1.0)
def shi(i): return np.where(isV[i],4.0,1.0)
def contrib(i,q):
    q=q.reshape(-1,4); X=q[:,:3]; M,rr=mass_of(i,X,q[:,3]); return M@K(X)-M0[i]@P0[i],X,M,rr
def penalties(i,q,tree,others):
    q=q.reshape(-1,4); X=q[:,:3]; M,rr=mass_of(i,X,q[:,3]); d=np.linalg.norm(X,axis=1); dr=dring(X)
    v=[np.maximum(0,d+rr-RT),np.maximum(0,rr+0.5-d),np.maximum(0,CLR-(dr-rr))]
    dd,jj=tree.query(X,k=40); rj=others[jj]; v.append(np.maximum(0,rr[:,None]+rj+GAP-dd).ravel())
    if len(i)>1:
        a,b=np.triu_indices(len(i),1); v.append(np.maximum(0,rr[a]+rr[b]+GAP-np.linalg.norm(X[a]-X[b],axis=1)))
    return np.concatenate(v)
def build_tree(edited,Xe,re_):
    mask=np.ones(n,bool); mask[edited]=False
    return cKDTree(X0[mask]),r0[mask]
def refit(sel,q,fixedw):
    """jointly fit bodies sel (params q) against T - fixedw (fixedw = contribution of other edited bodies)"""
    sel=np.array(sel); tree,others=build_tree(sel,None,None); lo=np.zeros((len(sel),4)); hi=np.zeros((len(sel),4))
    lo[:,:3]=X0[sel]-MOVE*3; hi[:,:3]=X0[sel]+MOVE*3; lo[:,3]=slo(sel); hi[:,3]=shi(sel); q=np.clip(q.reshape(-1,4),lo+1e-9,hi-1e-9).ravel(); lam=30.0
    tgt=T-fixedw
    h=1e-4
    def body_c(x):              # per-body pattern contributions (K,360)
        x=x.reshape(-1,4); M,rr=mass_of(sel,x[:,:3],x[:,3]); return M[:,None]*K(x[:,:3])
    def body_p(x):              # per-body penalty rows (region + neighbours), (K, 3+40)
        x=x.reshape(-1,4); X=x[:,:3]; M,rr=mass_of(sel,X,x[:,3]); d=np.linalg.norm(X,axis=1); dr=dring(X)
        dd,jj=tree.query(X,k=40); return np.c_[np.maximum(0,d+rr-RT),np.maximum(0,rr+0.5-d),np.maximum(0,CLR-(dr-rr)),np.maximum(0,rr[:,None]+others[jj]+GAP-dd)]
    a_,b_=np.triu_indices(len(sel),1)
    def jac_all(x,lam):
        Kn=len(sel); x=x.reshape(-1,4); c0=body_c(x); p0=body_p(x); npair=len(a_)
        J=np.zeros((360+Kn*43+npair,4*Kn))
        for a in range(4):
            xp=x.copy(); xp[:,a]+=h; dc=(body_c(xp)-c0)/h; dp=(body_p(xp)-p0)/h
            J[:360,a::4]=-dc.T/FLOOR/np.sqrt(360)
            for k in range(Kn): J[360+k*43:360+(k+1)*43,4*k+a]=lam*dp[k]
        if npair:
            X=x[:,:3]; M,rr=mass_of(sel,X,x[:,3]); dv=X[a_]-X[b_]; dn=np.maximum(np.linalg.norm(dv,axis=1),1e-9); act=(rr[a_]+rr[b_]+GAP-dn)>0
            drs=np.where(isV[sel],1.0,rr*np.log(10)/3); rows=360+Kn*43+np.arange(npair)
            for c in range(3):
                J[rows,4*a_+c]=lam*act*(-dv[:,c]/dn); J[rows,4*b_+c]=lam*act*(dv[:,c]/dn)
            J[rows,4*a_+3]=lam*act*drs[a_]; J[rows,4*b_+3]=lam*act*drs[b_]
        return J
    for it in range(4):
        fn=lambda x: np.concatenate([(tgt-contrib(sel,x)[0])/FLOOR/np.sqrt(360), lam*body_p(x).ravel(), lam*(penalties(sel,x,tree,others)[len(sel)*43:] if len(sel)>1 else np.zeros(0))])
        r=least_squares(fn,q,jac=lambda x,l=lam: jac_all(x,l),bounds=(lo.ravel(),hi.ravel()),method="trf",x_scale="jac",max_nfev=int(os.environ.get("NFEV",8))*20,xtol=1e-8,ftol=1e-8); q=r.x
        if penalties(sel,q,tree,others).max(initial=0)<=1e-3: break
        lam*=10
    return q,penalties(sel,q,tree,others).max(initial=0)
def leverage(res,Xc,sc,edited):
    """linearised one-body chi^2 reduction for every unedited body"""
    h=0.05; base=P0; J=np.zeros((n,4,360))
    for a in range(3):
        Xp=X0.copy(); Xp[:,a]+=h; J[:,a]=M0[:,None]*(K(Xp)-base)/h
    J[:,3]=np.where(isV[:,None],(3*M0/r0)[:,None]*base,(np.log(10)*M0)[:,None]*base)      # d/dR (void) or d/dlog10M (cluster)
    G=np.einsum("nai,nbi->nab",J,J); g=np.einsum("nai,i->na",J,res)
    mu=1e-3*np.trace(G,axis1=1,axis2=2)[:,None,None]/4*np.eye(4)[None]
    dl=np.linalg.solve(G+mu,g[...,None])[...,0]
    sc_=np.minimum(1,MOVE/np.maximum(np.linalg.norm(dl[:,:3],axis=1),1e-12)); dl*=sc_[:,None]
    gain=2*np.einsum("na,na->n",dl,g)-np.einsum("na,nab,nb->n",dl,G,dl); gain[edited]=-np.inf; gain[BAN]=-np.inf; gain[((dring(X0)-r0)<CLR)&(not NEAR)]=-np.inf      # bodies grown within 3 Gpc of the ring stay as grown
    return gain/FLOOR**2
st=dict(sel=[],q=np.zeros(0),curve=[(0,float(rms(T)/FLOOR))],rounds=0,t=0.0)
if os.path.exists(CK): st=pickle.load(open(CK,"rb")); print("resumed",tag,len(st["sel"]),flush=True)
t0=time.time()-st["t"]; BAN=np.zeros(n,bool); BAN[st.get("ban",[])]=True
def total_edit(sel,q):
    return contrib(np.array(sel),q)[0] if len(sel) else np.zeros(360)
while True:
    sel=st["sel"]; q=st["q"]; we=total_edit(sel,q); res=T-we; L=rms(res)/FLOOR
    if L<=1 or len(sel)>=KMAX: break
    cv=st["curve"]
    if len(sel)>=400 and [c for c in cv if c[0]<=len(sel)-300] and L>0.99*[c for c in cv if c[0]<=len(sel)-300][-1][1]: print("plateau",flush=True); break
    gain=leverage(res,None,None,np.array(sel,int))
    nb=int(np.clip(len(sel)//6,2,24)); new=list(np.argsort(-gain)[:nb]); new=[j for j in new if np.isfinite(gain[j]) and gain[j]>0]
    if not new: break
    qn=np.c_[X0[new],s0(np.array(new))].ravel()
    st["rounds"]+=1
    if st["rounds"]%6==0 and len(sel)+len(new)<=150:      # joint refit of everything edited so far + the new ones
        allsel=sel+new; qa,v=refit(allsel,np.r_[q,qn],np.zeros(360)); 
        Lq=rms(T-contrib(np.array(allsel),qa)[0])/FLOOR
        if v<=1e-3 and Lq<rms(T-we-contrib(np.array(new),qn)[0])/FLOOR: st["sel"],st["q"]=allsel,qa
        else: qn2,v2=refit(new,qn,we); st["sel"],st["q"]=(sel+new,np.r_[q,qn2]) if v2<=1e-3 else (sel,q)
    else:
        qn2,v2=refit(new,qn,we)
        if v2<=1e-3: st["sel"],st["q"]=sel+new,np.r_[q,qn2]
        else:   # drop the new ones that could not be placed legally; ban them by marking as edited-but-unchanged
            st["sel"],st["q"]=sel+new,np.r_[q,qn]          # unchanged params = original body (it is legal where it was grown)
    L=rms(T-total_edit(st["sel"],st["q"]))/FLOOR
    if L>=st["curve"][-1][1]-1e-6:          # no improvement: undo this round, ban the new candidates
        st["sel"],st["q"]=sel,q; BAN[new]=True; st["ban"]=list(np.where(BAN)[0]); L=st["curve"][-1][1]
        if BAN.sum()>5000: break
    else: st["curve"].append((len(st["sel"]),float(L)))
    st["t"]=time.time()-t0
    pickle.dump(st,open(CK+".tmp","wb")); os.replace(CK+".tmp",CK)
    print(f"  {tag} edited={len(st['sel'])} left={L:.2f} t={st['t']:.0f}s",flush=True)
sel=np.array(st["sel"],int); q=st["q"].reshape(-1,4) if len(sel) else np.zeros((0,4)); res=T-total_edit(list(sel),st["q"])
Xn=q[:,:3]; Mn,rn=mass_of(sel,Xn,q[:,3]) if len(sel) else (np.zeros(0),np.zeros(0))
changed=(np.linalg.norm(Xn-X0[sel],axis=1)>1e-3)|(np.abs(rn-r0[sel])/r0[sel]>1e-3)
out=dict(tag=tag+SF,near_ring_editable=NEAR,n_frozen_near_ring=int(((dring(X0)-r0)<CLR).sum()) if not NEAR else 0,n_bodies=n,n_voids=int(isV.sum()),before=spec(T),after=spec(res),left=float(rms(res)/FLOOR),reached=bool(rms(res)/FLOOR<=1),
  n_edited=int(changed.sum()),frac_edited=float(changed.sum()/n),curve=st["curve"],idx=sel.tolist(),changed=changed.tolist(),is_void=isV[sel].tolist(),
  X_old=X0[sel].tolist(),X_new=Xn.tolist(),r_old=r0[sel].tolist(),r_new=rn.tolist(),M_old=M0[sel].tolist(),M_new=Mn.tolist(),
  move=np.linalg.norm(Xn-X0[sel],axis=1).tolist(),dist_us_old=np.linalg.norm(X0[sel],axis=1).tolist(),dist_us_new=np.linalg.norm(Xn,axis=1).tolist(),
  dist_ring_old=dring(X0[sel]).tolist(),dist_ring_new=dring(Xn).tolist(),net_mass_change=float(Mn.sum()-M0[sel].sum()),resid=res.tolist(),secs=time.time()-t0)
json.dump(out,open(OUTF,"w")); print("DONE",tag,"edited",out["n_edited"],"of",n,"left %.3f"%out["left"],flush=True)
