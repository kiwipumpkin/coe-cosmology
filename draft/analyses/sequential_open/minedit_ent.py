"""ENTROPY VARIANT (minedit_ent.py, Chris 10 Oct): move + resize only, no deletion, no additions, NO 3 Gpc ring backstop
(only the growth rule that a body may not straddle the ring: surface >= 0.1 Gpc off it). Instead: every entropy component of
entropy_core.Score must stay >= its starting value (hard check on every accepted round; one-sided soft penalty inside the fit).
--- original v2 docstring follows ---
Minimal edit v2 (Chris, 9 Oct late): edit EXISTING bodies only -- move, rescale, or DELETE (no additions).
Per edited body: (x,y,z,a), a = mass scale relative to the grown body (void: deficit scale, R = R0 a^1/3; cluster: M = a M0, r = r0 a^1/3).
  a in [0,10] (void R up to ~3.2 Gpc for R0 = 1.5; cluster mass up to x10). Realistic minimum: void R >= 0.5 Gpc, cluster >= 0.1 M0;
  anything fitted below that is snapped to a = 0 = DELETED (counted as an edit).
Rules for a body that remains: inside 30 Gpc, >= 0.5 Gpc from us, surface >= 3 Gpc from the ring, no overlap (0.02 Gpc gap). Near-ring bodies
  (grown within 3 Gpc of the ring) are editable: they must move out, shrink, or be deleted. Rules fade out as a body is deleted (a -> 0).
Selection by leverage: for every unedited body, the better of (a) linearised move/rescale chi^2 gain and (b) exact deletion gain.
Fit: the selected batch is fitted jointly with bounded nonlinear least squares from several starts (as found / deleted / relocated at random
  uniform-in-volume within 5 Gpc), best legal start kept; every 5 rounds the most recent <=60 edited bodies are re-fitted jointly together.
  A round is kept only if the leftover goes down. Stop at <= 1x, plateau (<1% gain over 300 edits) or KMAX edits.
usage: python3 minedit2.py TAG -> minedit2/<TAG>.json, ckpt/minedit2_<TAG>.pkl"""
import os,sys,json,time,pickle,numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from common import *
from fit_on_web import spec
K=seq_core.unit_pattern_fast; RT=30.0; CLR=0.1; CLR3=3.0; GAP=0.02; MOVE=5.0; AMAX=10.0
tag=sys.argv[1]; CLEAR=float(os.environ.get("CLEAR","-1")); NODEL=True; NMODE="entropy"; SF="_entropy"; OUTF=f"{OUT}/minedit_ent/{tag}{SF}.json"; CK=f"{OUT}/ckpt/minedit_ent_{tag}{SF}.pkl"; KMAX=int(os.environ.get("KMAX",5000))
f=np.load(f"{OUT}/final_{tag}.npz"); X0=f["Xp"].astype(float); r0=f["rp"].astype(float); isV=(f["kp"]=="void"); M0=f["Mp"].astype(float); n=len(X0)
T=R-f["w"]; P0=K(X0)*1.0; C0=M0[:,None]*P0
rho0=seq_core.rho_gpc(X0)
AMIN=np.where(isV,(0.5/r0)**3,0.1)
def dring(X): X=np.atleast_2d(X); return np.sqrt((np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2)
NEARRING=(dring(X0)-r0)<CLR3        # reporting only: no special treatment
import entropy_core as EC
ESC=EC.Score(X0,r0,isV); S0=ESC.S0; EMARG=float(os.environ.get("EMARG","2e-5")); LAME=float(os.environ.get("LAME","3000"))
SOFTK=("radial","ring","angular")
def Hsoft(base,Xs):
    cr,cg,ca=EC.counts_soft(Xs); return np.array([EC.Hvol(base[0]+cr,EC.W_R),EC.Hvol(base[1]+cg,EC.W_G),EC.Hvol(base[2]+ca,EC.W_A)])
def ent_layout(): return ESC.components(CUR["X"],CUR["r"],isV)
def mass_of(i,X,a):
    a=np.maximum(a,0); Mv=M0[i]*a*seq_core.rho_gpc(np.atleast_2d(X))/rho0[i]       # void deficit follows local density
    return np.where(isV[i],Mv,M0[i]*a), r0[i]*np.maximum(a,1e-9)**(1/3)
def body_c(i,x):
    x=x.reshape(-1,4); M,_=mass_of(i,x[:,:3],x[:,3]); return M[:,None]*K(x[:,:3])
def body_p(i,x,tree,others):
    x=x.reshape(-1,4); X=x[:,:3]; M,rr=mass_of(i,X,x[:,3]); pf=np.minimum(1,np.maximum(x[:,3],0)/AMIN[i])
    d=np.linalg.norm(X,axis=1); dr=dring(X); dd,jj=tree.query(X,k=40)
    return pf[:,None]*np.c_[np.maximum(0,d+rr-RT),np.maximum(0,rr+0.5-d),np.maximum(0,CLR-(dr-rr)),np.maximum(0,rr[:,None]+others[jj]+GAP-dd)]
def pair_p(i,x):
    x=x.reshape(-1,4); X=x[:,:3]; _,rr=mass_of(i,X,x[:,3]); a,b=np.triu_indices(len(i),1)
    return np.maximum(0,rr[a]+rr[b]+GAP-np.linalg.norm(X[a]-X[b],axis=1))
CUR={"X":X0.copy(),"r":r0.copy()}          # current layout (edited bodies at their new positions/sizes; deleted -> r = 0)
def set_cur(sel,q):
    X=X0.copy(); r=r0.copy()
    if len(sel):
        qq=q.reshape(-1,4); si=np.array(sel); X[si]=qq[:,:3]; r[si]=np.where(qq[:,3]>0,r0[si]*np.maximum(qq[:,3],1e-9)**(1/3),0.0)
    CUR["X"],CUR["r"]=X,r
def tree_excl(sel):
    m=CUR["r"]>0; m[sel]=False; return cKDTree(CUR["X"][m]),CUR["r"][m]
def viol(sel,q,tree,others):
    v=body_p(sel,q,tree,others).max(initial=0)
    return max(v,pair_p(sel,q).max(initial=0)) if len(sel)>1 else v
def fitset(sel,q,tgt):
    """joint bounded LSQ of bodies sel (params q) to target tgt (= T minus all other edits, original patterns of sel removed already)"""
    sel=np.array(sel); tree,others=tree_excl(sel); Kn=len(sel); a_,b_=np.triu_indices(Kn,1); h=1e-4
    m_=CUR["r"]>0; m_[sel]=False; m_&=np.linalg.norm(CUR["X"],axis=1)<=RT; base=EC.counts_hard(CUR["X"][m_]); Xstart=CUR["X"][sel].copy()
    Shard=ent_layout(); Hc=np.array([Shard[k] for k in SOFTK]); Hs0=Hsoft(base,Xstart); Smin=np.array([S0[k] for k in SOFTK])+EMARG
    def erow(x): return LAME*np.maximum(0,Smin-(Hc+Hsoft(base,x.reshape(-1,4)[:,:3])-Hs0))
    lo=np.c_[X0[sel]-15,AMIN[sel] if NODEL else np.zeros(Kn)]; hi=np.c_[X0[sel]+15,np.full(Kn,AMAX)]; q=np.clip(q.reshape(-1,4),lo+1e-9,hi-1e-9).ravel()
    def jac(x,lam):
        x=x.reshape(-1,4); c0=body_c(sel,x); p0=body_p(sel,x,tree,others); J=np.zeros((360+Kn*43+len(a_)+3,4*Kn)); e0=erow(x)
        for k in range(Kn):
            for a in range(3):
                xq=x.copy(); xq[k,a]+=h; J[-3:,4*k+a]=(erow(xq)-e0)/h
        for a in range(4):
            xp=x.copy(); xp[:,a]+=h; J[:360,a::4]=-((body_c(sel,xp)-c0)/h).T/FLOOR/np.sqrt(360); dp=(body_p(sel,xp,tree,others)-p0)/h
            for k in range(Kn): J[360+k*43:360+(k+1)*43,4*k+a]=lam*dp[k]
            if len(a_): J[360+Kn*43:,a::4]=0
        if len(a_):   # pair rows by finite difference per coordinate pair (only active pairs matter)
            pp=pair_p(sel,x); act=np.where(pp>0)[0]
            for a in range(4):
                xp=x.copy(); xp[:,a]+=h; dp=(pair_p(sel,xp)-pp)/h          # all bodies shifted: wrong for positions -> do per body for active pairs
                for t_ in act:
                    for bb in (a_[t_],b_[t_]):
                        xq=x.copy(); xq[bb,a]+=h; J[360+Kn*43+t_,4*bb+a]=lam*(pair_p(sel,xq)[t_]-pp[t_])/h
        return J
    lam=30.0
    for it in range(4):
        fn=lambda x: np.concatenate([(tgt-body_c(sel,x).sum(0))/FLOOR/np.sqrt(360),lam*body_p(sel,x,tree,others).ravel(),lam*(pair_p(sel,x) if Kn>1 else np.zeros(0)),erow(x)])
        q=least_squares(fn,q,jac=lambda x,l=lam: jac(x,l),bounds=(lo.ravel(),hi.ravel()),method="trf",x_scale="jac",max_nfev=int(os.environ.get("NFEV",60)),xtol=1e-8,ftol=1e-8).x
        if viol(sel,q,tree,others)<=1e-3: break
        lam*=10
    # snap sub-minimum sizes to deletion, then a short polish of the survivors
    q=q.reshape(-1,4)
    if not NODEL: q[q[:,3]<AMIN[sel],3]=0.0
    return q.ravel(),viol(sel,q.ravel(),tree,others)
def leverage(res,edited):
    h=0.05; J=np.zeros((n,4,360))
    for a in range(3):
        Xp=X0.copy(); Xp[:,a]+=h; J[:,a]=M0[:,None]*(K(Xp)-P0)/h
    J[:,3]=C0
    G=np.einsum("nai,nbi->nab",J,J); g=np.einsum("nai,i->na",J,res); mu=1e-3*np.trace(G,axis1=1,axis2=2)[:,None,None]/4*np.eye(4)[None]
    dl=np.linalg.solve(G+mu,g[...,None])[...,0]; dl*=np.minimum(1,MOVE/np.maximum(np.linalg.norm(dl[:,:3],axis=1),1e-12))[:,None]
    gl=2*np.einsum("na,na->n",dl,g)-np.einsum("na,nab,nb->n",dl,G,dl); # (entropy variant: near-ring bodies are free like any other)
    gd=np.sum(res**2)-np.sum((res[None]+C0)**2,axis=1)                                      # exact deletion gain
    gain=gl if NODEL else np.maximum(gl,gd); gain[edited]=-np.inf; gain[BAN]=-np.inf; return gain
def total(sel,q):
    if not len(sel): return np.zeros(360)
    return body_c(np.array(sel),q).sum(0)-C0[np.array(sel)].sum(0)
st=dict(sel=[],q=np.zeros(0),curve=[(0,float(rms(T)/FLOOR))],rounds=0,t=0.0,ban=[])
NCLR=0
if CLEAR>=0:          # start by deleting every body whose surface is within CLEAR Gpc of the ring (removes the fine-scale texture)
    cl=np.where((dring(X0)-r0)<CLEAR)[0]; NCLR=len(cl)
    st["sel"]=cl.tolist(); st["q"]=np.c_[X0[cl],np.zeros(NCLR)].ravel(); st["curve"].append((NCLR,float(rms(T+C0[cl].sum(0))/FLOOR)))
def relocate(ids,Xc,rc):
    """move each body straight away from the nearest ring point to surface clearance 3.05 Gpc; if blocked, step further / rotate"""
    Xc=Xc.copy(); out=np.zeros((len(ids),3)); rr=np.random.default_rng(7)
    for k,i in enumerate(ids):
        x=X0[i]; rho=np.hypot(x[0],x[1]); p=np.array([x[0]/rho*RMIN,x[1]/rho*RMIN,0.0]); u=x-p; u/=max(np.linalg.norm(u),1e-9)
        al=rc>0; al[i]=False; tr=cKDTree(Xc[al]); ra=rc[al]; done=False
        for R_ in (4.0,6.0,9.0,13.0):     # random candidates around the body; keep the legal one closest to where it was
            C=x+rr.normal(size=(6000,3))*R_/2; C=np.r_[C,p+u*(CLR+r0[i]+0.05+np.linspace(0,R_,40))[:,None]]
            s_=np.linalg.norm(C,axis=1); ok=(dring(C)-r0[i]>=CLR+0.01)&(s_+r0[i]<=RT-0.01)&(s_-r0[i]>=0.51)
            C=C[ok]
            if not len(C): continue
            dd,jj=tr.query(C,k=30); ok=np.all(dd>=ra[jj]+r0[i]+GAP+0.005,axis=1); C=C[ok]
            if len(C): c=C[np.argmin(np.linalg.norm(C-x,axis=1))]; done=True; break
        if not done: raise RuntimeError(f"no legal spot for {i}")
        Xc[i]=c; out[k]=c
    return out
if False:      # stage 1: all near-ring bodies leave the 3 Gpc band (relocate-first: fixed nearest legal spots; forced: fitted, starting from their own spots)
    nrb=np.where(NEARRING)[0]; nrb=nrb[np.argsort(dring(X0[nrb])-r0[nrb])]
    Xrel=relocate(nrb,X0,r0.copy()); NCLR=len(nrb)
    sel_=[]; q_=np.zeros(0)
    if NMODE=="reloc":
        sel_=list(nrb); q_=np.c_[Xrel,np.ones(NCLR)].ravel(); set_cur(sel_,q_)
        for c0 in range(0,NCLR,40):         # joint refit (move/resize) of the relocated bodies, chunk by chunk
            ch=list(nrb[c0:c0+40]); oth=[b for b in sel_ if b not in ch]; qq=q_.reshape(-1,4); pos={b:k for k,b in enumerate(sel_)}
            qo=np.array([qq[pos[b]] for b in oth]).ravel(); qc=np.array([qq[pos[b]] for b in ch]).ravel()
            tg=T-(total(oth,qo) if oth else 0)+C0[np.array(ch)].sum(0); set_cur(sel_,q_)
            qn,v=fitset(ch,qc,tg)
            if v<=1e-3 and rms(tg-body_c(np.array(ch),qn).sum(0))<rms(tg-body_c(np.array(ch),qc).sum(0)):
                qq=qq.copy(); qn=qn.reshape(-1,4)
                for k,b in enumerate(ch): qq[pos[b]]=qn[k]
                q_=qq.ravel()
            print(f"  stage1 reloc chunk {c0} left={rms(T-total(sel_,q_))/FLOOR:.2f}",flush=True)
    else:
        for c0 in range(0,NCLR,40):         # forced out: fit the chunk from its own (illegal) spots + from the relocation spots; best legal start
            ch=np.array(nrb[c0:c0+40]); set_cur(sel_,q_)
            tg=T-(total(sel_,q_) if sel_ else 0)+C0[ch].sum(0); best=None
            for s0_ in (np.c_[X0[ch],np.ones(len(ch))],np.c_[Xrel[c0:c0+40],np.ones(len(ch))]):
                qn,v=fitset(list(ch),s0_.ravel(),tg)
                if v<=1e-3:
                    Lq=rms(tg-body_c(ch,qn).sum(0))
                    if best is None or Lq<best[0]: best=(Lq,qn)
            if best is None: raise RuntimeError("forced chunk could not be made legal")
            sel_=sel_+list(ch); q_=np.r_[q_,best[1]]
            print(f"  stage1 forced chunk {c0} left={rms(T-total(sel_,q_))/FLOOR:.2f}",flush=True)
    st["sel"]=sel_; st["q"]=q_; st["curve"].append((NCLR,float(rms(T-total(sel_,q_))/FLOOR))); st["nclr"]=NCLR
    pickle.dump(st,open(CK+".tmp","wb")); os.replace(CK+".tmp",CK)
if NODEL: NCLR=st.get("nclr",NCLR)
if os.path.exists(CK): st=pickle.load(open(CK,"rb")); print("resumed",tag,len(st["sel"]),flush=True)
BAN=np.zeros(n,bool); BAN[st["ban"]]=True; t0=time.time()-st["t"]; rng=np.random.default_rng(int(time.time()))
while True:
    sel=st["sel"]; q=st["q"]; set_cur(sel,q); we=total(sel,q); res=T-we; L=rms(res)/FLOOR; cv=st["curve"]
    if L<=1 or len(sel)>=KMAX: break
    old=[c for c in cv if c[0]<=len(sel)-300]
    if len(sel)>=300 and old and L>0.99*old[-1][1]: print("plateau",flush=True); break
    gain=leverage(res,np.array(sel,int)); nb=int(np.clip(len(sel)//6,2,24)); new=[j for j in np.argsort(-gain)[:nb] if np.isfinite(gain[j]) and gain[j]>0]
    if not new: break
    st["rounds"]+=1; new=np.array(new)
    tgt=T-we+C0[new].sum(0)               # what the new bodies (fitted) must supply
    starts=[np.c_[X0[new],np.ones(len(new))]]+([] if NODEL else [np.c_[X0[new],np.zeros(len(new))]])
    for _ in range(2):
        u=rng.normal(size=(len(new),3)); u/=np.linalg.norm(u,axis=1)[:,None]; starts.append(np.c_[X0[new]+u*MOVE*rng.uniform(size=(len(new),1))**(1/3),np.ones(len(new))])
    U=X0[new]+rng.uniform(-15,15,(len(new),3)); bad=np.linalg.norm(U,axis=1)>RT-1
    while bad.any(): U[bad]=X0[new][bad]+rng.uniform(-15,15,(bad.sum(),3)); bad=np.linalg.norm(U,axis=1)>RT-1
    starts.append(np.c_[U,np.ones(len(new))])
    best=None
    for s0_ in starts:
        set_cur(sel,q); qq,v=fitset(new,s0_.ravel(),tgt)
        if v<=1e-3:
            set_cur(sel+list(new),np.r_[q,qq]); Sx=ent_layout(); set_cur(sel,q)
            if not ESC.ok(Sx): continue                       # hard: no entropy component may drop below the start
            Lq=rms(tgt-body_c(new,qq).sum(0))
            if best is None or Lq<best[0]: best=(Lq,qq)
    if best is None: BAN[new]=True; st["ban"]=list(np.where(BAN)[0]); continue
    nsel=sel+list(new); nq=np.r_[q,best[1]]; set_cur(nsel,nq)
    if st["rounds"]%5==0:                  # joint refit of the most recent <=60 edited bodies
        blk=nsel[max(NCLR,len(nsel)-60):]; qb=nq.reshape(-1,4)[-len(blk):].ravel(); oth=nsel[:-len(blk)]; othq=nq.reshape(-1,4)[:len(nsel)-len(blk)].ravel()
        tg=T-total(oth,othq)+C0[np.array(blk)].sum(0); qb2,v=fitset(blk,qb,tg)
        if v<=1e-3 and rms(tg-body_c(np.array(blk),qb2).sum(0))<rms(tg-body_c(np.array(blk),qb).sum(0)):
            set_cur(nsel,np.r_[othq,qb2])
            if ESC.ok(ent_layout()): nq=np.r_[othq,qb2]
            set_cur(nsel,nq)
    Ln=rms(T-total(nsel,nq))/FLOOR
    if Ln>=L-1e-6: BAN[new]=True; st["ban"]=list(np.where(BAN)[0])
    else: st["sel"],st["q"]=nsel,nq; st["curve"].append((len(nsel),float(Ln)))
    st["t"]=time.time()-t0; pickle.dump(st,open(CK+".tmp","wb")); os.replace(CK+".tmp",CK)
    print(f"  {tag} edited={len(st['sel'])} left={min(Ln,L):.2f} t={st['t']:.0f}s",flush=True)
    if BAN.sum()>8000: break
# ---- stage 2 (FINE=1): fine-scale (m13+) pass. Linear leverage cannot see that a near-ring body should move several Gpc,
# so rank the unedited bodies within 6 Gpc of the ring by their own m13+ power and fit them jointly in batches of 12 from spread starts
# (in place, pushed 2-8 Gpc away from the ring, uniform in volume within 15 Gpc). Same acceptance: leftover must drop, entropy rule hard.
def hp(c): F=np.fft.rfft(c,axis=-1); F[...,:13]=0; return np.fft.irfft(F,n=360,axis=-1)
if os.environ.get("FINE")=="1":
    ban2=set(st.get("ban2",[])); fails=st.get("fails2",0); st["stage2_start"]=st.get("stage2_start",len(st["sel"]))
    pw=np.linalg.norm(hp(C0),axis=1); pool=[i for i in np.argsort(-pw) if dring(X0[i])-r0[i]<6]
    while fails<40:
        sel=st["sel"]; q=st["q"]; set_cur(sel,q); we=total(sel,q); L=rms(T-we)/FLOOR
        if L<=1 or len(sel)>=KMAX: break
        es=set(sel); new=np.array([i for i in pool if i not in es and i not in ban2][:12])
        if not len(new): break
        tgt=T-we+C0[new].sum(0); rr_=np.hypot(X0[new,0],X0[new,1]); P_=np.c_[X0[new,0]/rr_*RMIN,X0[new,1]/rr_*RMIN,np.zeros(len(new))]
        u_=X0[new]-P_; u_/=np.maximum(np.linalg.norm(u_,axis=1),1e-9)[:,None]
        starts=[np.c_[X0[new],np.ones(len(new))]]+[np.c_[X0[new]+u_*rng.uniform(2,8,(len(new),1)),np.ones(len(new))] for _ in range(2)]
        U=X0[new]+rng.uniform(-15,15,(len(new),3)); bad=np.linalg.norm(U,axis=1)>RT-1
        while bad.any(): U[bad]=X0[new][bad]+rng.uniform(-15,15,(bad.sum(),3)); bad=np.linalg.norm(U,axis=1)>RT-1
        starts.append(np.c_[U,np.ones(len(new))]); best=None; diag=[]
        for s0_ in starts:
            set_cur(sel,q); qq,v=fitset(new,s0_.ravel(),tgt)
            set_cur(sel+list(new),np.r_[q,qq]); Sx=ent_layout(); okE=ESC.ok(Sx); set_cur(sel,q)
            diag.append((round(float(v),4),okE,round(float(rms(tgt-body_c(new,qq).sum(0))/FLOOR),3),[k for k in S0 if Sx[k]<S0[k]]))
            if v<=1e-3 and okE:
                # absorb the large-scale change: 24 more bodies picked by leverage on the new residual, fitted jointly together with the batch
                s1=sel+list(new); q1=np.r_[q,qq]; set_cur(s1,q1); r1=T-total(s1,q1)
                g_=leverage(r1,np.array(s1,int)); ab=[j for j in np.argsort(-g_)[:24] if np.isfinite(g_[j]) and g_[j]>0]
                if ab:
                    ab=np.array(ab); blk=list(new)+list(ab); qb=np.r_[qq,np.c_[X0[ab],np.ones(len(ab))].ravel()]
                    set_cur(sel,q); tg=T-we+C0[np.array(blk)].sum(0); qb2,v2=fitset(blk,qb,tg)
                    if v2<=1e-3:
                        set_cur(sel+blk,np.r_[q,qb2]); ok2=ESC.ok(ent_layout()); set_cur(sel,q)
                        if ok2: qq_,bl_=qb2,blk
                        else: qq_,bl_=qq,list(new)
                    else: qq_,bl_=qq,list(new)
                else: qq_,bl_=qq,list(new)
                Lq=rms(T-total(sel+bl_,np.r_[q,qq_]))/FLOOR; diag.append(("absorbed",len(bl_),round(float(Lq),3)))
                if best is None or Lq<best[0]: best=(Lq,qq_,bl_)
        Ln=best[0] if best else np.inf
        if Ln<L-1e-4: st["sel"],st["q"]=sel+best[2],np.r_[q,best[1]]; st["curve"].append((len(st["sel"]),float(Ln))); fails=0
        else: ban2.update(new.tolist()); fails+=1
        st["ban2"]=list(ban2); st["fails2"]=fails; st["t"]=time.time()-t0; pickle.dump(st,open(CK+".tmp","wb")); os.replace(CK+".tmp",CK)
        print(f"  stage2 {tag} edited={len(st['sel'])} left={min(Ln,L):.3f} fails={fails} base={rms(tgt-C0[new].sum(0))/FLOOR:.3f} starts={diag}",flush=True)
sel=np.array(st["sel"],int); q=st["q"].reshape(-1,4) if len(sel) else np.zeros((0,4)); res=T-total(list(sel),st["q"])
Xn=q[:,:3]; an=q[:,3]; Mn,rn=mass_of(sel,Xn,an) if len(sel) else (np.zeros(0),np.zeros(0))
dele=an<=0; mv=np.linalg.norm(Xn-X0[sel],axis=1); moved=(~dele)&(mv>0.05); resized=(~dele)&(np.abs(an-1)>0.01)
# reference: delete every near-ring body and nothing else
Tdel=T+C0[NEARRING].sum(0)
out=dict(no_delete=NODEL,nodel_mode=NMODE if NODEL else None,n_stage1=NCLR if NODEL else 0,tag=tag+SF,clear_D=CLEAR,n_cleared=NCLR,n_bodies=n,n_voids=int(isV.sum()),n_near_ring=int(NEARRING.sum()),before=spec(T),after=spec(res),left=float(rms(res)/FLOOR),reached=bool(rms(res)/FLOOR<=1),
  delete_all_near_ring=dict(left=float(rms(Tdel)/FLOOR),spec=spec(Tdel)),
  n_edited=int((dele|moved|resized).sum()),n_deleted=int(dele.sum()),n_moved=int(moved.sum()),n_resized=int(resized.sum()),n_near_ring_edited=int(NEARRING[sel].sum()),
  frac_edited=float((dele|moved|resized).sum()/n),curve=st["curve"],idx=sel.tolist(),is_void=isV[sel].tolist(),near_ring=NEARRING[sel].tolist(),
  deleted=dele.tolist(),moved=moved.tolist(),resized=resized.tolist(),a_new=an.tolist(),
  X_old=X0[sel].tolist(),X_new=Xn.tolist(),r_old=r0[sel].tolist(),r_new=np.where(dele,0,rn).tolist(),M_old=M0[sel].tolist(),M_new=np.where(dele,0,Mn).tolist(),move=mv.tolist(),
  dist_us_old=np.linalg.norm(X0[sel],axis=1).tolist(),dist_us_new=np.linalg.norm(Xn,axis=1).tolist(),dist_ring_old=dring(X0[sel]).tolist(),dist_ring_new=dring(Xn).tolist(),
  net_mass_change=float(np.where(dele,0,Mn).sum()-M0[sel].sum()),resid=res.tolist(),secs=time.time()-t0)
set_cur(list(sel),st["q"]); S1=ent_layout(); Xa=CUR["X"]; ra=CUR["r"]
L0=EC.Score.links(X0,r0>0,isV); L1=EC.Score.links(Xa,ra>0,isV)
from scipy.stats import ks_2samp
VOLR=4/3*np.pi*RT**3; edges=np.arange(0,15.01,1.0); vb=np.histogram(EC.dring(EC._P),edges)[0]/len(EC._P)*VOLR
out.update(stage2_from=st.get("stage2_start"),entropy_before=S0,entropy_after=S1,entropy_ok=bool(ESC.ok(S1)),ks_before_after={k:float(ks_2samp(L0[k],L1[k]).statistic) for k in L0},
  ring_profile=dict(edges=edges.tolist(),vol=vb.tolist(),before=np.histogram(dring(X0),edges)[0].tolist(),after=np.histogram(dring(Xa),edges)[0].tolist()),
  n_end_within3_ring=int(((dring(Xa)-ra)<3).sum()),n_start_within3_ring=int(NEARRING.sum()))
json.dump(out,open(OUTF,"w")); print("DONE",tag,"edited",out["n_edited"],"del",out["n_deleted"],"left %.3f"%out["left"],flush=True)
