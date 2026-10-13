"""COMPENSATED-WEB minimal edit (Chris, 10 Oct), same frame / kernel / residual / entropy rule as minedit_ent.py.
Web: final_<tag>.npz from compweb_gen.py (each void + its wall clusters = one compensated group).
Edits = move + resize only (no deletion, no additions), no 3 Gpc wall. Two kinds of edit unit, both keep every void's compensation exact:
  GROUP unit g   (4 params d_x,d_y,d_z,a): the void and all its clusters move rigidly by d and scale by a
                  (all masses x a x rho(new centre)/rho(old centre), all radii and wall offsets x a^(1/3)) -> net mass of the group stays 0.
  CLUSTER unit k (3 params u): one cluster slides to a new spot inside its own void's wall (|u| from R+r_c+gap to 1.3R, in the group's
                  frame), mass unchanged -> compensation unchanged.
Body positions: x_k = c0_g + d_g + a_g^(1/3) u_k.  Rules for every body: inside 30 Gpc, >= 0.5 Gpc from us, surface >= 0.1 Gpc from the ring,
no overlaps (gap 0.02). ENTROPY RULE: every component of entropy_core.Score must stay >= its start (hard check every accepted round,
one-sided soft penalty inside the fit). Leverage selection, joint bounded least squares from several starts, refit of the recent units
every 5 rounds, then a fine-scale (m13+) pass on near-ring units. Checkpoint ckpt/minedit_comp_<tag>.pkl (resumes).
usage: python3 minedit_comp.py comp_s1 -> minedit_comp/comp_s1_entropy.json"""
import os,sys,json,time,pickle,numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from common import *
from fit_on_web import spec
import entropy_core as EC
K=seq_core.unit_pattern_fast; RT=30.0; CLR=0.1; GAP=0.02; MOVE=5.0; MOVEC=2.0; AMAX=10.0; WALL=1.3
tag=sys.argv[1]; OUTF=f"{OUT}/minedit_comp/{tag}_entropy.json"; CK=f"{OUT}/ckpt/minedit_comp_{tag}.pkl"
KMAX=int(os.environ.get("KMAX",5000)); TMAX=float(os.environ.get("TMAX",7200)); TMAX2=float(os.environ.get("TMAX2",3600))
f=np.load(f"{OUT}/final_{tag}.npz"); X0=f["Xp"].astype(float); r0=f["rp"].astype(float); isV=(f["kp"]=="void"); M0=f["Mp"].astype(float); grp=f["grp"].astype(int); n=len(X0)
T=R-f["w"]; NG=grp.max()+1; vb=np.zeros(NG,int); vb[grp[isV]]=np.where(isV)[0]; c0=X0[vb]; RV0=r0[vb]; rhoc0=seq_core.rho_gpc(c0)
U0=X0-c0[grp]; MEM=[np.where(grp==g)[0] for g in range(NG)]; CL=np.where(~isV)[0]
NU=NG+len(CL)                                   # units: 0..NG-1 groups, NG.. cluster units
def unit_bodies(u): return MEM[u] if u<NG else np.array([CL[u-NG]])
ESC=EC.Score(X0,r0,isV); S0=ESC.S0; EMARG=float(os.environ.get("EMARG","2e-5")); LAME=float(os.environ.get("LAME","3000")); SOFTK=("radial","ring","angular")
def dring(X): X=np.atleast_2d(X); return np.sqrt((np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2)
AMIN=np.maximum((0.5/RV0)**3,0.1)               # void radius >= 0.5 Gpc, masses >= 0.1x
# ---------------- state
def fresh(): return dict(D=np.zeros((NG,3)),A=np.ones(NG),U=U0.copy(),edited=[],curve=[(0,float(rms(T)/FLOOR))],rounds=0,t=0.0,ban=[],stage=1)
def layout_of(D,A,U,idx=None):
    idx=np.arange(n) if idx is None else idx; g=grp[idx]; s=A[g]**(1/3); c=c0[g]+D[g]
    X=c+s[:,None]*U[idx]; M=M0[idx]*A[g]*seq_core.rho_gpc(c)/rhoc0[g]; r=r0[idx]*s; return X,M,r
st=pickle.load(open(CK,"rb")) if os.path.exists(CK) else fresh()
if os.path.exists(CK): print("resumed",tag,len(st["edited"]),flush=True)
def cur(): return layout_of(st["D"],st["A"],st["U"])
def ent_of(X,r): return ESC.components(X,r,isV)
def Hsoft(base,Xs):
    cr,cg,ca=EC.counts_soft(Xs); return np.array([EC.Hvol(base[0]+cr,EC.W_R),EC.Hvol(base[1]+cg,EC.W_G),EC.Hvol(base[2]+ca,EC.W_A)])
# ---------------- batch fit
class Batch:
    def __init__(s,units,D,A,U):
        s.units=list(units); s.D,s.A,s.U=D.copy(),A.copy(),U.copy()
        s.np_=[4 if u<NG else 3 for u in s.units]; s.off=np.r_[0,np.cumsum(s.np_)]
        s.bod=np.unique(np.concatenate([unit_bodies(u) for u in s.units])); s.pos={b:i for i,b in enumerate(s.bod)}
        s.ub=[np.array([s.pos[b] for b in unit_bodies(u)]) for u in s.units]
        X,M,r=layout_of(D,A,U); s.Xcur,s.Mcur=X[s.bod],M[s.bod]; s.Ccur=(s.Mcur[:,None]*K(s.Xcur)).sum(0)
        al=np.ones(n,bool); al[s.bod]=False; s.tree=cKDTree(X[al]); s.oth=r[al]
        inR=al&(np.linalg.norm(X,axis=1)<=RT); s.base=EC.counts_hard(X[inR])
        s.Hc=np.array([ent_of(X,r)[k] for k in SOFTK]); s.Hs0=Hsoft(s.base,s.Xcur); s.Smin=np.array([S0[k] for k in SOFTK])+EMARG
        s.wall=[(i,u) for i,u in enumerate(s.units) if u>=NG]
    def q0(s):
        q=[]
        for u in s.units:
            if u<NG: q+= list(s.D[u])+[s.A[u]]
            else: q+= list(s.U[CL[u-NG]])
        return np.array(q,float)
    def bounds(s):
        lo=[];hi=[]
        for u in s.units:
            if u<NG: lo+=[-15,-15,-15,AMIN[u]]; hi+=[15,15,15,AMAX]
            else: g=grp[CL[u-NG]]; w=WALL*RV0[g]; lo+=[-w]*3; hi+=[w]*3
        return np.array(lo),np.array(hi)
    def apply(s,q):
        D,A,U=s.D.copy(),s.A.copy(),s.U.copy()
        for i,u in enumerate(s.units):
            p=q[s.off[i]:s.off[i+1]]
            if u<NG: D[u]=p[:3]; A[u]=p[3]
            else: U[CL[u-NG]]=p
        return D,A,U
    def rows(s,q):
        D,A,U=s.apply(q); X,M,r=layout_of(D,A,U,s.bod)
        cmb=(T_tgt_holder[0]-(M[:,None]*K(X)).sum(0))/FLOOR/np.sqrt(360)
        d=np.linalg.norm(X,axis=1); dr=dring(X); dd,jj=s.tree.query(X,k=40)
        bp=np.c_[np.maximum(0,d+r-RT),np.maximum(0,r+0.5-d),np.maximum(0,CLR-(dr-r)),np.maximum(0,r[:,None]+s.oth[jj]+GAP-dd)].ravel()
        wr=[]
        for i,u in s.wall:
            k=CL[u-NG]; g=grp[k]; uu=np.linalg.norm(U[k]); wr+=[max(0,RV0[g]+r0[k]+GAP-uu),max(0,uu-WALL*RV0[g])]
        a_,b_=np.triu_indices(len(s.bod),1); pr=np.maximum(0,r[a_]+r[b_]+GAP-np.linalg.norm(X[a_]-X[b_],axis=1))
        er=np.maximum(0,s.Smin-(s.Hc+Hsoft(s.base,X)-s.Hs0))
        return cmb,bp,np.array(wr),pr,er,X,M,r
    def viol(s,q):
        _,bp,wr,pr,_,_,_,_=s.rows(q); return max(bp.max(initial=0),wr.max(initial=0),pr.max(initial=0))
    def fit(s,q,tgt,nfev=None):
        T_tgt_holder[0]=tgt; lo,hi=s.bounds(); q=np.clip(q,lo+1e-9,hi-1e-9); lam=30.0; h=1e-4
        def fn(x,lam): c,bp,wr,pr,er,_,_,_=s.rows(x); return np.concatenate([c,lam*bp,lam*wr,lam*pr,LAME*er])
        nb_=len(s.bod); a_,b_=np.triu_indices(nb_,1)
        def parts(x):
            D,A,U=s.apply(x); X,M,r=layout_of(D,A,U,s.bod); return X,M,r,U
        def brow(X,r):
            d=np.linalg.norm(X,axis=1); dr=dring(X); dd,jj=s.tree.query(X,k=40)
            return np.c_[np.maximum(0,d+r-RT),np.maximum(0,r+0.5-d),np.maximum(0,CLR-(dr-r)),np.maximum(0,r[:,None]+s.oth[jj]+GAP-dd)]
        nw=len(s.wall); o_b=360; o_w=360+43*nb_; o_p=o_w+2*nw; o_e=o_p+len(a_)
        widx={i:t for t,(i,u) in enumerate(s.wall)}
        pidx=[np.where(np.isin(a_,bi)|np.isin(b_,bi))[0] for bi in s.ub]
        def jac(x,lam):
            X,M,r,U=parts(x); P=M[:,None]*K(X); B=brow(X,r); f0=fn(x,lam); J=np.zeros((len(f0),len(x)))
            for i,u in enumerate(s.units):
                bi=s.ub[i]; pi=pidx[i]
                for j in range(s.off[i],s.off[i+1]):
                    xp=x.copy(); xp[j]+=h; Xp,Mp,rp,Up=parts(xp)
                    J[:360,j]=-((Mp[bi,None]*K(Xp[bi])).sum(0)-P[bi].sum(0))/FLOOR/np.sqrt(360)/h
                    Bp=brow(Xp[bi],rp[bi])
                    for t_,b in enumerate(bi): J[o_b+43*b:o_b+43*(b+1),j]=lam*(Bp[t_]-B[b])/h
                    if i in widx:
                        k=CL[u-NG]; g=grp[k]; uu=np.linalg.norm(Up[k]); w2=np.array([max(0,RV0[g]+r0[k]+GAP-uu),max(0,uu-WALL*RV0[g])])
                        J[o_w+2*widx[i]:o_w+2*widx[i]+2,j]=(lam*w2-f0[o_w+2*widx[i]:o_w+2*widx[i]+2])/h
                    if len(pi):
                        pp=np.maximum(0,rp[a_[pi]]+rp[b_[pi]]+GAP-np.linalg.norm(Xp[a_[pi]]-Xp[b_[pi]],axis=1)); J[o_p+pi,j]=(lam*pp-f0[o_p+pi])/h
                    er=np.maximum(0,s.Smin-(s.Hc+Hsoft(s.base,Xp)-s.Hs0)); J[o_e:,j]=(LAME*er-f0[o_e:])/h
            return J
        for it in range(4):
            q=least_squares(lambda x: fn(x,lam),q,jac=lambda x: jac(x,lam),bounds=(lo,hi),method="trf",x_scale="jac",max_nfev=nfev or int(os.environ.get("NFEV",30)),xtol=1e-8,ftol=1e-8).x
            if s.viol(q)<=1e-3: break
            lam*=10
        return q,s.viol(q)
T_tgt_holder=[None]
def total_now(D,A,U):
    X,M,_=layout_of(D,A,U); return np.einsum("n,ni->i",M,K(X))
# fast incremental bookkeeping of the web pattern
X_,M_,r_=cur(); PAT=M_[:,None]*K(X_); W=PAT.sum(0)
def resid(): return T-(W-f["w"])
def commit(D,A,U,bod):
    global W
    X,M,_=layout_of(D,A,U,bod); new=M[:,None]*K(X); W=W-PAT[bod].sum(0)+new.sum(0); PAT[bod]=new
    st["D"],st["A"],st["U"]=D,A,U
def leverage(res,excl):
    X,M,r=cur(); P=K(X); h=0.05; Jx=np.zeros((n,3,360))
    for a in range(3):
        Xp=X.copy(); Xp[:,a]+=h; Jx[:,a]=M[:,None]*(K(Xp)-P)/h
    g_=grp; s=st["A"][g_]**(1/3); ha=0.05; c=c0[g_]+st["D"][g_]
    Xa=c+(s*(1+ha)**(1/3))[:,None]*st["U"]; Ja=(M[:,None]*(1+ha)*K(Xa)-M[:,None]*P)/ha*st["A"][g_][:,None]   # d/d a (a in absolute units)
    JG=np.zeros((NG,4,360)); np.add.at(JG,g_,np.concatenate([Jx,Ja[:,None]/st["A"][g_][:,None,None]],1))
    JC=Jx[CL]
    def gain(J,mv):
        G=np.einsum("nai,nbi->nab",J,J); gg=np.einsum("nai,i->na",J,res); mu=1e-3*np.trace(G,axis1=1,axis2=2)[:,None,None]/J.shape[1]*np.eye(J.shape[1])[None]
        dl=np.linalg.solve(G+mu,gg[...,None])[...,0]; dl*=np.minimum(1,mv/np.maximum(np.linalg.norm(dl[:,:3],axis=1),1e-12))[:,None]
        return 2*np.einsum("na,na->n",dl,gg)-np.einsum("na,nab,nb->n",dl,G,dl)
    gn=np.r_[gain(JG,MOVE),gain(JC,MOVEC)]; gn[list(excl)]=-np.inf
    gn[st["ban"]]=-np.inf
    for u in st["edited"]:      # a cluster unit is not picked while its group unit is being fitted in the same batch (handled in pick)
        pass
    return gn
def pick(gn,nb):
    out=[]; groups=set()
    for u in np.argsort(-gn):
        if len(out)>=nb or not np.isfinite(gn[u]) or gn[u]<=0: break
        g=u if u<NG else grp[CL[u-NG]]
        if (u<NG and any((v>=NG and grp[CL[v-NG]]==u) for v in out)) or (u>=NG and g in groups): continue
        out.append(int(u)); groups.add(g) if u<NG else None
    return out
def starts_for(units,rng,far=False):
    """start parameter vectors: as is; random nudges; (far) pushed away from the ring / random within 15 Gpc"""
    b=Batch(units,st["D"],st["A"],st["U"]); q=b.q0(); S=[q.copy()]
    for k in range(2 if far else 1):
        qq=q.copy()
        for i,u in enumerate(units):
            o=b.off[i]
            if u<NG:
                c=c0[u]+st["D"][u]
                if far and k==0:
                    rho=np.hypot(c[0],c[1]); p=np.array([c[0]/rho*RMIN,c[1]/rho*RMIN,0]); v=c-p; v/=max(np.linalg.norm(v),1e-9); qq[o:o+3]+=v*rng.uniform(2,8)
                elif far and k==1: qq[o:o+3]+=rng.uniform(-15,15,3)
                else: v=rng.normal(size=3); v/=np.linalg.norm(v); qq[o:o+3]+=v*MOVE*rng.uniform()**(1/3)
            else:
                kk=CL[u-NG]; g=grp[kk]; v=rng.normal(size=3); v/=np.linalg.norm(v); qq[o:o+3]=v*rng.uniform(RV0[g]+r0[kk]+GAP+0.01,WALL*RV0[g]-0.01)
        lo,hi=b.bounds(); S.append(np.clip(qq,lo+1e-6,hi-1e-6))
    return b,S
def try_batch(units,rng,far=False):
    """fit units jointly from several starts; returns (left, D,A,U, bodies) of the best legal + entropy-ok start, or None"""
    b,S=starts_for(units,rng,far); tgt=resid()+b.Ccur; best=None; DIAG.clear()
    for q0 in S:
        q,v=b.fit(q0,tgt)
        if v>1e-3: DIAG.append(("viol",round(float(v),4))); continue
        D,A,U=b.apply(q); X,M,r=layout_of(D,A,U); Sx=ent_of(X,r)
        if not ESC.ok(Sx): DIAG.append(("ent",[k for k in S0 if Sx[k]<S0[k]])); continue
        c,_,_,_,_,_,_,_=b.rows(q); Lq=np.sqrt(np.sum(c**2))*np.sqrt(360)/np.sqrt(360)
        Lq=rms(tgt-(layout_of(D,A,U,b.bod)[1][:,None]*K(layout_of(D,A,U,b.bod)[0])).sum(0))/FLOOR
        DIAG.append(("ok",round(float(Lq),2)))
        if best is None or Lq<best[0]: best=(Lq,D,A,U,b.bod)
    return best
DIAG=[]
rng=np.random.default_rng(int(time.time())); t0=time.time()-st["t"]
def save():
    st["t"]=time.time()-t0; pickle.dump(st,open(CK+".tmp","wb")); os.replace(CK+".tmp",CK)
# ---------------- interleaved search: each step tries a fine-scale batch (near-ring units ranked by their own m13+ power, starts
# include pushes 2-8 Gpc away from the ring and random spots within 15 Gpc) and, if that fails, a leverage batch (linearised gain).
# Every accepted fine batch is followed by a leverage 'absorb' batch; every 5 accepted steps the last 30 units are refitted jointly.
def hp(c): F_=np.fft.rfft(c,axis=-1); F_[...,:13]=0; return np.fft.irfft(F_,n=360,axis=-1)
st.setdefault("ban2",[]); st.setdefault("fails",0); st.setdefault("acc",0); FAILMAX=int(os.environ.get("FAILMAX",40))
def accept(best,units):
    commit(best[1],best[2],best[3],best[4]); st["edited"]+=[u for u in units if u not in st["edited"]]
while st["stage"]<3:
    L=rms(resid())/FLOOR; ne=len(st["edited"]); cv=st["curve"]
    if L<=1 or ne>=KMAX or time.time()-t0>TMAX or st["fails"]>=FAILMAX: break
    old=[c for c in cv if c[0]<=ne-300]
    if ne>=300 and old and L>0.99*old[-1][1]: print("plateau",flush=True); break
    st["rounds"]+=1; ok=False
    X,M,r=cur(); near=(dring(X)-r)<6
    gp=np.zeros(NG); np.add.at(gp,grp,np.linalg.norm(hp(PAT),axis=1)); cp=np.linalg.norm(hp(PAT[CL]),axis=1)
    pw=np.r_[gp,cp]; nu=np.r_[np.bincount(grp[near],minlength=NG)>0,near[CL]]; pw[~nu]=-1; pw[st["ban2"]]=-1
    new=pick(pw,8); kind=""
    if new:
        best=try_batch(new,rng,far=True)
        if best is not None and best[0]<L-1e-4:
            accept(best,new); ok=True; kind="fine"
            ab=pick(leverage(resid(),new),16)
            if ab:
                b2=try_batch(ab,rng)
                if b2 is not None and b2[0]<rms(resid())/FLOOR-1e-6: accept(b2,ab); kind+="+absorb"
        else: st["ban2"]=sorted(set(st["ban2"])|set(new))
    if not ok:
        nb=int(np.clip(ne//6,2,24)); new=pick(leverage(resid(),st["edited"]),nb)
        if new:
            best=try_batch(new,rng)
            if best is not None and best[0]<L-1e-6: accept(best,new); ok=True; kind="leverage"
            else: st["ban"]=sorted(set(st["ban"])|set(new))
    if ok:
        st["fails"]=0; st["acc"]+=1
        if st["acc"]%5==0:
            blk=st["edited"][-30:]; gs=set(u for u in blk if u<NG); blk=[u for u in blk if u<NG or grp[CL[u-NG]] not in gs]
            b2=try_batch(blk,rng)
            if b2 is not None and b2[0]<rms(resid())/FLOOR-1e-6: commit(b2[1],b2[2],b2[3],b2[4]); kind+="+refit"
        st["curve"].append((len(st["edited"]),float(rms(resid())/FLOOR)))
    else: st["fails"]+=1
    save(); print(f"  {tag} step {st['rounds']} {kind or 'rejected'} units={len(st['edited'])} left={rms(resid())/FLOOR:.3f} fails={st['fails']} t={time.time()-t0:.0f}s {'' if ok else DIAG}",flush=True)
st["stage"]=3; save()
# ---------------- output (same schema as minedit_ent, so edit_before_after / tri_view can draw it)
Xa,Ma,ra=cur(); res=resid()
mv=np.linalg.norm(Xa-X0,axis=1); moved=mv>0.05; resized=np.abs(Ma/M0-1)>0.01; ch=moved|resized; sel=np.where(ch)[0]
S1=ent_of(Xa,ra); L0=EC.Score.links(X0,r0>0,isV); L1=EC.Score.links(Xa,ra>0,isV)
from scipy.stats import ks_2samp
VOLR=4/3*np.pi*RT**3; edges=np.arange(0,15.01,1.0); vbv=np.histogram(EC.dring(EC._P),edges)[0]/len(EC._P)*VOLR
gnet=np.bincount(grp,weights=Ma,minlength=NG); gdef=-Ma[vb]
eu=np.array(st["edited"],int)
out=dict(tag=tag+"_entropy",n_bodies=n,n_voids=int(isV.sum()),n_groups=int(NG),before=spec(T),after=spec(res),left=float(rms(res)/FLOOR),reached=bool(rms(res)/FLOOR<=1),
  n_units_edited=int(len(eu)),n_group_units=int((eu<NG).sum()),n_cluster_units=int((eu>=NG).sum()),
  n_edited=int(ch.sum()),n_deleted=0,n_moved=int(moved.sum()),n_resized=int(resized.sum()),frac_edited=float(ch.sum()/n),
  n_voids_moved=int((moved&isV).sum()),n_clusters_moved=int((moved&~isV).sum()),curve=st["curve"],
  idx=sel.tolist(),is_void=isV[sel].tolist(),deleted=[False]*len(sel),moved=moved[sel].tolist(),resized=resized[sel].tolist(),
  X_old=X0[sel].tolist(),X_new=Xa[sel].tolist(),r_old=r0[sel].tolist(),r_new=ra[sel].tolist(),M_old=M0[sel].tolist(),M_new=Ma[sel].tolist(),move=mv[sel].tolist(),
  dist_ring_old=dring(X0[sel]).tolist(),dist_ring_new=dring(Xa[sel]).tolist(),move_median_void=float(np.median(mv[moved&isV])) if (moved&isV).any() else 0.0,
  move_median_cluster=float(np.median(mv[moved&~isV])) if (moved&~isV).any() else 0.0,mass_factor_median=float(np.median((Ma/M0)[resized])) if resized.any() else 1.0,
  compensation_max_rel=float(np.max(np.abs(gnet)/gdef)),resid=res.tolist(),secs=time.time()-t0,
  entropy_before=S0,entropy_after=S1,entropy_ok=bool(ESC.ok(S1)),ks_before_after={k:float(ks_2samp(L0[k],L1[k]).statistic) for k in L0},
  ring_profile=dict(edges=edges.tolist(),vol=vbv.tolist(),before=np.histogram(dring(X0),edges)[0].tolist(),after=np.histogram(dring(Xa),edges)[0].tolist()),
  n_end_within3_ring=int(((dring(Xa)-ra)<3).sum()),n_start_within3_ring=int(((dring(X0)-r0)<3).sum()))
json.dump(out,open(OUTF,"w")); print("DONE",tag,"edited",out["n_edited"],"units",len(eu),"left %.3f"%out["left"],flush=True)
