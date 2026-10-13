"""Minimal edit v3: sparse joint refit over EXISTING bodies (group OMP + bounded joint least squares), Gauss-Newton outer loop.
Edits per body: rescale/delete (mass scale a in [0,10]; a below the realistic minimum -> deleted) and small moves.
Every body, incl. near-ring ones, may be resized or deleted in place; moves are only allowed so a body keeps >=3 Gpc surface clearance
  from the ring, stays inside 30 Gpc and >=0.5 Gpc from us, and cannot reach a neighbour (move bound per axis = min(1 Gpc, half the gap to the
  nearest neighbour /sqrt3, (ring clearance-3)/sqrt3, ...)); growth is capped so the new radius cannot reach half the nearest gap -> no overlaps by construction.
Selection: at each step the bodies with the largest group gain (their 4 columns: size + x,y,z) against the current residual are added (batch),
  then ALL selected bodies are re-fitted jointly (bounded least squares, BVLS/TRF). The exact nonlinear pattern is recomputed every few steps
  (re-linearised at the current positions). Stop at exact leftover <= 1x or stall.
usage: python3 minedit3.py TAG -> minedit3/<TAG>.json, ckpt/minedit3_<TAG>.pkl"""
import os,sys,json,time,pickle,numpy as np
from scipy.optimize import lsq_linear
from scipy.spatial import cKDTree
from common import *
from fit_on_web import spec
K=seq_core.unit_pattern_fast; RT=30.0; CLR=3.0; GAP=0.02; AMAX=10.0; MB=1.0
tag=sys.argv[1]; OUTF=f"{OUT}/minedit3/{tag}.json"; CK=f"{OUT}/ckpt/minedit3_{tag}.pkl"; KMAX=int(os.environ.get("KMAX",4000))
f=np.load(f"{OUT}/final_{tag}.npz"); X0=f["Xp"].astype(float); r0=f["rp"].astype(float); isV=(f["kp"]=="void"); M0=f["Mp"].astype(float); n=len(X0)
T=R-f["w"]; rho0=seq_core.rho_gpc(X0)
XS=X0.copy(); AS=np.ones(n)            # start state: as grown, or the end state of the move/resize/delete pass (minedit2) with START=v2
if os.environ.get("START")=="v2":
    _d=json.load(open(f"{OUT}/minedit2/{tag}.json")); _i=np.array(_d["idx"],int); XS[_i]=np.array(_d["X_new"]); AS[_i]=np.array(_d["a_new"])
    OUTF=f"{OUT}/minedit3/{tag}_v2.json"; CK=f"{OUT}/ckpt/minedit3_{tag}_v2.pkl"
def dring(X): X=np.atleast_2d(X); return np.sqrt((np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2)
AMIN=np.where(isV,(0.5/r0)**3,0.1); NEAR=(dring(X0)-r0)<CLR
rS=r0*np.maximum(AS,1e-9)**(1/3); alS=AS>0; trS=cKDTree(XS[alS]); rA=rS[alS]; dd,jj=trS.query(XS,k=41); _self=dd<1e-9; dd=np.where(_self,np.inf,dd); gapn=np.min(dd-rA[jj]-rS[:,None]-GAP,axis=1)     # free gap to nearest surface
gapn=np.maximum(gapn,0)
amax=np.maximum(AS,np.minimum(AMAX,((rS+gapn/2)/r0)**3))                                          # growth cannot use more than half the gap
s=np.linalg.norm(XS,axis=1); clr=dring(XS)-r0*amax**(1/3)
DXS=XS-X0
mb=np.minimum.reduce([np.full(n,MB),gapn/2/np.sqrt(3),np.maximum(0,clr-CLR)/np.sqrt(3),np.maximum(0,RT-s-r0*amax**(1/3))/np.sqrt(3),np.maximum(0,s-0.5-r0*amax**(1/3))/np.sqrt(3)])
mb[(dring(XS)-rS)<CLR]=0.0; NEAR=NEAR          # near-ring bodies: resize or delete in place (a move out would be > 1 Gpc)
def mass(i,X,a):
    return np.where(isV[i],M0[i]*a*seq_core.rho_gpc(np.atleast_2d(X))/rho0[i],M0[i]*a)
def exact(Xc,ac):     # exact pattern of the web with current positions/scales
    W=np.zeros(360)
    for i in range(0,n,3000): sl=slice(i,i+3000); W+=mass(np.arange(n)[sl],Xc[sl],ac[sl])@K(Xc[sl])
    return W
def columns(Xc,ac):   # linearisation at the current state: d pattern / d a, d pattern / d x,y,z   -> (n,4,360)
    h=0.02; P=K(Xc); J=np.zeros((n,4,360)); m1=mass(np.arange(n),Xc,np.ones(n)); J[:,0]=m1[:,None]*P
    for a in range(3):
        Xp=Xc.copy(); Xp[:,a]+=h; J[:,1+a]=(mass(np.arange(n),Xp,ac)[:,None]*K(Xp)-(m1*ac)[:,None]*P)/h
    return J
st=dict(sel=np.zeros(0,int),X=XS.copy(),a=AS.copy(),curve=[(0,float(rms(T)/FLOOR))],it=0,t=0.0)
if os.path.exists(CK): st=pickle.load(open(CK,"rb")); print("resumed",tag,len(st["sel"]),flush=True)
t0=time.time()-st["t"]
def save(): st["t"]=time.time()-t0; pickle.dump(st,open(CK+".tmp","wb")); os.replace(CK+".tmp",CK)
def n_edit(X,a): return int(((np.linalg.norm(X-X0,axis=1)>1e-3)|(np.abs(a-1)>1e-3)).sum())
stall=0
while True:
    Xc,ac=st["X"],st["a"]; W=exact(Xc,ac); res=R-W; L=rms(res)/FLOOR
    if L<=1 or len(st["sel"])>=KMAX or stall>=3: break
    J=columns(Xc,ac); J[:,1:]*=(mb>0)[:,None,None]                  # bodies without move room keep only the size column
    sel=list(st["sel"])
    for inner in range(8):                                             # add a batch, joint refit (linear), repeat; then re-linearise
        # group gain of each unselected body against the linear residual
        G=np.einsum("nai,nbi->nab",J,J)+1e-12*np.eye(4)[None]*np.einsum("nai,nai->n",J,J)[:,None,None]; g=np.einsum("nai,i->na",J,res)
        gain=np.einsum("na,na->n",g,np.linalg.solve(G,g[...,None])[...,0]); gain[sel]=-1
        nb=int(np.clip(len(sel)//5,4,60)); new=list(np.argsort(-gain)[:nb]); sel=sel+new
        si=np.array(sel); A=J[si].transpose(2,0,1).reshape(360,-1)/FLOOR   # columns ordered (body, param)
        lo=np.c_[-ac[si],np.repeat(-mb[si,None],3,1)-(Xc[si]-XS[si])].ravel(); hi=np.c_[amax[si]-ac[si],np.repeat(mb[si,None],3,1)-(Xc[si]-XS[si])].ravel()
        lo=np.minimum(lo,0); hi=np.maximum(hi,0)
        b0=(R-W)/FLOOR
        act=(hi-lo)>1e-9; xs=np.zeros(A.shape[1])
        sol=lsq_linear(A[:,act],b0,bounds=(lo[act],hi[act]),method="bvls" if act.sum()<=1600 else "trf",lsmr_tol="auto",max_iter=400); xs[act]=sol.x
        d=xs.reshape(-1,4); res=(b0-A@xs)*FLOOR; Llin=rms(res)/FLOOR
        if Llin<=0.7: break
    # apply, snap sizes below the realistic minimum to deletion, recompute exactly
    Xn=Xc.copy(); an=ac.copy(); an[si]=np.clip(ac[si]+d[:,0],0,amax[si]); Xn[si]=Xc[si]+d[:,1:]
    an[(an<AMIN)&(an>0)]=0.0
    Ln=rms(R-exact(Xn,an))/FLOOR
    if Ln<L: st["X"],st["a"]=Xn,an; stall=0
    else: stall+=1
    st["sel"]=si; st["it"]+=1; st["curve"].append((n_edit(st["X"],st["a"]),float(min(Ln,L)))); save()
    print(f"  {tag} it={st['it']} selected={len(si)} edited={st['curve'][-1][0]} lin={Llin:.2f} exact={Ln:.2f} t={st['t']:.0f}s",flush=True)
# ---- prune: try to undo the smallest edits while staying <= max(1, L) (makes the edited set as small as possible)
Xc,ac=st["X"],st["a"]; L=rms(R-exact(Xc,ac))/FLOOR; target=max(1.0,L)*1.0
ed=np.where((np.linalg.norm(Xc-X0,axis=1)>1e-3)|(np.abs(ac-1)>1e-3))[0]
J=columns(Xc,ac); eff=np.array([np.linalg.norm((mass(np.array([i]),Xc[i:i+1],ac[i:i+1])[0]*K(Xc[i:i+1])[0]-M0[i]*K(X0[i:i+1])[0])) for i in ed])
for i in ed[np.argsort(eff)]:
    Xt,at=Xc.copy(),ac.copy(); Xt[i]=X0[i]; at[i]=1.0
    dW=M0[i]*K(X0[i:i+1])[0]-mass(np.array([i]),Xc[i:i+1],ac[i:i+1])[0]*K(Xc[i:i+1])[0]
    Lt=rms(R-(st.get("_W",None) if False else exact(Xc,ac))-dW)/FLOOR if False else None
    break
# (exact incremental prune)
W=exact(Xc,ac); pr=0
for i in ed[np.argsort(eff)]:
    if os.environ.get('START')=='v2': break
    dW=M0[i]*K(X0[i:i+1])[0]-mass(np.array([i]),Xc[i:i+1],ac[i:i+1])[0]*K(Xc[i:i+1])[0]
    if rms(R-W-dW)/FLOOR<=target:
        W=W+dW; Xc[i]=X0[i]; ac[i]=1.0; pr+=1
st["X"],st["a"]=Xc,ac; save(); print("pruned",pr,flush=True)
res=R-exact(Xc,ac); L=rms(res)/FLOOR
ed=np.where((np.linalg.norm(Xc-X0,axis=1)>1e-3)|(np.abs(ac-1)>1e-3))[0]
dele=ac[ed]<=0; mv=np.linalg.norm(Xc[ed]-X0[ed],axis=1); moved=(~dele)&(mv>1e-3); resized=(~dele)&(np.abs(ac[ed]-1)>1e-3)
# legality check of the final layout
rr=r0*np.maximum(ac,0)**(1/3); alive=ac>0; Xa=Xc[alive]; ra=rr[alive]; t2=cKDTree(Xa); pairs=np.array(list(t2.query_pairs(2*ra.max()+GAP)))
ov=0 if not len(pairs) else int(((ra[pairs[:,0]]+ra[pairs[:,1]]+GAP-np.linalg.norm(Xa[pairs[:,0]]-Xa[pairs[:,1]],axis=1))>1e-6).sum())
ovo=0 if not len(pairs) else int(((ra[pairs[:,0]]+ra[pairs[:,1]]-np.linalg.norm(Xa[pairs[:,0]]-Xa[pairs[:,1]],axis=1))>1e-6).sum())
mvd=ed[moved]; ringbad=int(((dring(Xc[mvd])-rr[mvd])<CLR-1e-6).sum())
Mn=mass(ed,Xc[ed],ac[ed])
out=dict(tag=tag,method="group-OMP + joint bounded LSQ + prune",n_bodies=n,n_voids=int(isV.sum()),n_near_ring=int(NEAR.sum()),before=spec(T),after=spec(res),left=float(L),reached=bool(L<=1),
  n_edited=int(len(ed)),n_deleted=int(dele.sum()),n_moved=int(moved.sum()),n_resized=int(resized.sum()),n_near_ring_edited=int(NEAR[ed].sum()),frac_edited=float(len(ed)/n),
  overlaps_original_rule=ovo,overlaps_with_gap=ov,moved_within_3gpc_of_ring=ringbad,pruned=pr,
  curve=st["curve"],idx=ed.tolist(),is_void=isV[ed].tolist(),near_ring=NEAR[ed].tolist(),deleted=dele.tolist(),moved=moved.tolist(),resized=resized.tolist(),a_new=ac[ed].tolist(),
  X_old=X0[ed].tolist(),X_new=Xc[ed].tolist(),r_old=r0[ed].tolist(),r_new=np.where(dele,0,rr[ed]).tolist(),M_old=M0[ed].tolist(),M_new=np.where(dele,0,Mn).tolist(),move=mv.tolist(),
  dist_us_old=np.linalg.norm(X0[ed],axis=1).tolist(),dist_us_new=np.linalg.norm(Xc[ed],axis=1).tolist(),dist_ring_old=dring(X0[ed]).tolist(),dist_ring_new=dring(Xc[ed]).tolist(),
  net_mass_change=float(np.where(dele,0,Mn).sum()-M0[ed].sum()),resid=res.tolist(),secs=time.time()-t0,
  delete_all_near_ring=dict(left=float(rms(T+(M0[NEAR][:,None]*K(X0[NEAR])).sum(0))/FLOOR)))
json.dump(out,open(OUTF,"w")); print("DONE",tag,"edited",len(ed),"del",int(dele.sum()),"moved",int(moved.sum()),"left %.3f"%L,"overlaps",ovo,flush=True)
