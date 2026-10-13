"""Bodies on top of a FIXED background with a uniform-in-volume placement prior.
Modeled region (Chris standard, 9 Oct): our sphere + shell to 30 Gpc from us, INCLUDING the CoE-core cap inside 30 Gpc.
Prior: body positions a priori uniform in volume over that region. The candidate bank is drawn that way (no oversampling near the ring).
Closeness penalty (stated): for a body whose centre is d from the CMB ring, F(d) = fraction of the region's volume within d of the ring
  (empirical, from 3 M uniform points; F ~ 0.025 d^2 for small d, so F=1 at d0 ~ 6.3 Gpc). -log F(d) is minus the log probability that a
  uniformly placed body lands that close; it equals 2 log(d0/d) near the ring, i.e. the (d0/d)^2 law in log form with d0 set by the
  geometry, not chosen. Under the prior each -log F is Exp(1) distributed, so sum over N bodies ~ Gamma(N,1) -> a p-value for the layout.
Objective: chi^2/2 + LAM * sum(-log F), chi^2 = sum over the 360 ring samples of (residual / 1x)^2 (1x per sample as the noise).
  LAM=1 is the plain Bayesian posterior with that noise; LAM=0 is unpenalized (same bank); LAM>1 is a sensitivity test.
Greedy: add the candidate with the best (chi^2 gain/2 - LAM*penalty); keep it only if the joint refit improves; jointly refit all masses (bounded LSQ, clusters <= 2e52 kg,
  voids delta -0.7 R 0.5-4 Gpc, no straddling the ring, >= 0.5 Gpc from us, nothing under 1e49 kg); stop at < 1x or when no candidate
  improves the penalized objective. Checkpointed every 25 picks.
usage: LAM=1 python3 prior_add.py TAG -> web_add/<TAG>_prior_lam<LAM>.json"""
import os,sys,json,time,pickle,numpy as np
from scipy.optimize import lsq_linear
from scipy.stats import gamma
WEBDIR=os.path.dirname(os.path.abspath(__file__)); CO="/workspace/clump_scan/core_open"; sys.path.insert(0,CO)
_c=os.getcwd(); os.chdir(CO); from co_core import *; import ao_core; os.chdir(_c)
sys.path.insert(0,WEBDIR); from fit_on_web import spec
tag=sys.argv[1]; LAM=float(os.environ.get("LAM",1)); NB=int(os.environ.get("NBANK",80000)); NMAX=int(os.environ.get("NMAX",400)); RT=30.0
out_fn=f"{WEBDIR}/web_add/{tag}_prior_lam{LAM:g}.json"; ck=f"{WEBDIR}/ckpt/prior_{tag}_lam{LAM:g}.pkl"
if os.path.exists(out_fn): sys.exit(0)
rng=np.random.default_rng(7)
T=R-np.load(f"{WEBDIR}/final_{tag}.npz")["w"]
def dring_of(P): return np.sqrt((np.hypot(P[:,0],P[:,1])-RMIN)**2+P[:,2]**2)
def uniform(n):
    out=[]
    while sum(len(o) for o in out)<n:
        P=rng.uniform(-RT,RT,(4*n,3)); s=np.linalg.norm(P,axis=1); out.append(P[(s<=RT)&(s>=0.5)])
    return np.concatenate(out)[:n]
Uq=np.sort(dring_of(uniform(3_000_000)))
def Fd(dd): return np.clip(np.searchsorted(Uq,dd)/len(Uq),1/len(Uq),1.0)
P=uniform(NB); dr=dring_of(P); d=np.linalg.norm(P,axis=1)
room=np.minimum(dr-0.1,d-0.5); rho=rho_gpc(P)*GPC**3
cR_max=np.minimum(MHI,4/3*np.pi*DELTA_C*rho*np.clip(room,0,None)**3)
vR=np.clip(np.minimum(4.0,room),0,None); v_max=np.where(vR>=0.5,0.7*rho*4/3*np.pi*vR**3,0.0)
ok=(cR_max>MLO)|(v_max>0); P,dr,cR_max,v_max,rho=P[ok],dr[ok],cR_max[ok],v_max[ok],rho[ok]
pen=-np.log(Fd(dr))
S=1e52; A_=(unit_field(P)-ao_core.FLOC[None]).T*S/FLOOR      # columns: pattern per 1e52 kg in units of 1x
Tn=T/FLOOR; lo=-v_max/S; hi=np.maximum(cR_max,0)/S; n2=np.sum(A_*A_,0); amin=MLO/S
st=dict(sel=[],x=np.zeros(0),ban=np.zeros(len(lo),bool),npick=0,traj=[],t=0.0)
if os.path.exists(ck): st=pickle.load(open(ck,"rb")); print("resumed",tag,len(st["sel"]),flush=True)
sel,x,ban=st["sel"],st["x"],st["ban"]; res=Tn-(A_[:,sel]@x if len(sel) else 0); t0=time.time()-st["t"]; stop="max bodies"
while len(sel)<NMAX and st["npick"]<4000:
    c=A_.T@res; a=np.clip(c/n2,lo,hi); a=np.where(np.abs(a)>=amin,a,0.0)
    gain=(2*a*c-a*a*n2)/2; score=gain-LAM*pen; score[sel]=-np.inf; score[ban]=-np.inf
    j=int(np.argmax(score))
    if score[j]<=0: stop="no body improves the penalized objective"; break
    sel.append(j); st["npick"]+=1; sg=st.setdefault("sg",[]); sg.append(1 if c[j]>0 else -1)
    # each body keeps the type (cluster/void) it was picked as, with a real floor (>= 1e49 kg); no dropping/banning
    sgv=np.array(sg); si=np.array(sel); lb=np.where(sgv>0,amin,lo[si]); ub=np.where(sgv>0,hi[si],-amin)
    xn=lsq_linear(A_[:,si],Tn,bounds=(lb,ub),method="bvls",max_iter=5000).x
    xg=np.r_[x,a[j]]; Lg=rms(Tn-A_[:,si]@xg); Lnew=rms(Tn-A_[:,si]@xn); Lold=rms(res)
    if Lg<Lnew: xn,Lnew=xg,Lg                       # bvls can stall; the greedy step (others unchanged, new body at its own best mass) is always feasible
    if Lnew>=Lold-1e-9: sel.pop(); sg.pop(); ban[j]=True; continue
    x=xn
    res=Tn-(A_[:,sel]@x if len(sel) else 0); L=rms(res)
    st["traj"].append((len(sel),float(L)))
    if L<1: stop="reached 1x"; break
    if st["npick"]%25==0:
        st.update(sel=sel,x=x,ban=ban,t=time.time()-t0); pickle.dump(st,open(ck+".tmp","wb")); os.replace(ck+".tmp",ck)
        print(f"  {tag} lam={LAM:g} N={len(sel)} left={L:.2f} t={time.time()-t0:.0f}s",flush=True)
X=P[sel]; M=x*S; dX=dr[sel]; lp=float(np.sum(np.log(Fd(dX)))); N=len(sel)
pval=float(gamma.sf(-lp,N)) if N else 1.0
res1=res*FLOOR
out=dict(tag=tag,lam=LAM,stop=stop,N=N,left=float(rms(res)),reached=bool(rms(res)<1),n_mass=int((M>0).sum()),n_void=int((M<0).sum()),
  inside_us=int((np.linalg.norm(X,axis=1)<RMIN).sum()),in_core=int((core_depth(X)>0).sum()) if N else 0,
  within_3gpc=int((dX<3).sum()),within_3gpc_expected=float(N*Fd(np.array([3.0]))[0]),within_1_5gpc=int((dX<1.5).sum()),
  logp_closeness=lp,logp_expected=-float(N),p_value=pval,dist_ring=dX.tolist(),F_of_d=Fd(dX).tolist(),
  background=spec(T),final=spec(res1),resid=res1.tolist(),X=X.tolist(),M=M.tolist(),traj=st["traj"][::5]+st["traj"][-1:],nbank=int(len(P)),secs=time.time()-t0,
  F_ref={str(q):float(Fd(np.array([q]))[0]) for q in (0.5,1,1.5,2,3,4,5,6,8,10,15)})
os.makedirs(f"{WEBDIR}/web_add",exist_ok=True); json.dump(out,open(out_fn,"w"))
if os.path.exists(ck): os.remove(ck)
print("DONE",tag,"lam",LAM,stop,"N",N,"left %.3f"%out["left"],"m/v",out["n_mass"],out["n_void"],"within3",out["within_3gpc"],"exp %.1f"%out["within_3gpc_expected"],"logp %.1f (exp %d) p=%.2g"%(lp,-N,pval),"t %.0f"%out["secs"],flush=True)
