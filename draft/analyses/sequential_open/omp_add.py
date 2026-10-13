"""Add bodies on top of a FIXED background until leftover < 1x, any type, anywhere (core_open/all-open kernel + FLOC compensation).
Fast linear method: a bank of candidate positions (half within 2 Gpc of the CMB ring, a quarter inside our sphere, a quarter anywhere in the
filled region incl. the core); each position is ONE body whose signed mass is fitted: + = cluster (<= 2e52 kg), - = deep void
(delta -0.7, radius 0.5-4 Gpc -> deficit bound). Bounds also keep every body from straddling the ring (radius + 0.1 Gpc clearance) and
0.5 Gpc from us. Greedy: add the candidate that best matches what's left, then re-fit ALL chosen masses together (bounded least squares);
repeat until < 1x. Bodies are not checked against each other or against the background's spheres (positions are spread; noted).
usage: python3 omp_add.py TAG [SEED] -> web_add/<TAG>.json"""
import os,sys,json,time,numpy as np
from scipy.optimize import lsq_linear
WEBDIR=os.path.dirname(os.path.abspath(__file__)); CO="/workspace/clump_scan/core_open"; sys.path.insert(0,CO)
_c=os.getcwd(); os.chdir(CO); from co_core import *; import ao_core; os.chdir(_c)
sys.path.insert(0,WEBDIR); from fit_on_web import spec
tag=sys.argv[1]; seed=int(sys.argv[2]) if len(sys.argv)>2 else 0
_o=f"{WEBDIR}/web_add/{tag}_omp{seed}{os.environ.get('OSUF','')}.json"; _l=_o+".running"
if os.path.exists(_o) or os.path.exists(_l): sys.exit(0)
open(_l,"w").close(); rng=np.random.default_rng(seed)
T=R-np.load(f"{WEBDIR}/final_{tag}.npz")["w"]; RX=float(os.environ.get("RING_EX",0.1));   # clearance between a body's surface and the ring (Gpc)
NB=int(os.environ.get("NBANK",40000)); NMAX=int(os.environ.get("NMAX",1500))
def dirs(n): u=rng.normal(size=(n,3)); return u/np.linalg.norm(u,axis=1,keepdims=True)
# bank
n1,n2=NB//2,NB//4; n3=NB-n1-n2
ph=rng.uniform(0,2*np.pi,n1); a=rng.uniform(RX,RX+3.0,n1)*np.sqrt(rng.uniform(0,1,n1))**0; b=rng.uniform(0,2*np.pi,n1)
t1=np.c_[(RMIN+a*np.cos(b))*np.cos(ph),(RMIN+a*np.cos(b))*np.sin(ph),a*np.sin(b)]
t2=dirs(n2)*(rng.uniform(0.5**3,RMIN**3,n2)**(1/3))[:,None]
t3=[]
while len(t3)<n3:
    if rng.uniform()<0.4: x=dirs(1)[0]*42.9*rng.uniform()**(1/3)
    else: x=C_G+dirs(1)[0]*A_IN_G*rng.uniform()**(1/3)
    if np.linalg.norm(x)>0.5: t3.append(x)
P=np.r_[t1,t2,np.array(t3)]; d=np.linalg.norm(P,axis=1); dring=np.sqrt((np.hypot(P[:,0],P[:,1])-RMIN)**2+P[:,2]**2)
room=np.minimum(dring-RX,d-0.5)                       # largest radius a body here may have
rho=rho_gpc(P)*GPC**3
cR_max=np.minimum(MHI,4/3*np.pi*DELTA_C*rho*np.clip(room,0,None)**3)          # cluster radius must fit in 'room'
vR=np.clip(np.minimum(4.0,room),0,None); v_max=np.where(vR>=0.5,0.7*rho*4/3*np.pi*vR**3,0.0)   # void needs R>=0.5
ok=(cR_max>MLO)|(v_max>0); P,cR_max,v_max,rho=P[ok],cR_max[ok],v_max[ok],rho[ok]
F=(unit_field(P)-ao_core.FLOC[None])                    # per kg; (n,360)
S=1e52; A_=F.T*S                                        # columns in units of 1e52 kg
lo=-v_max/S; hi=np.maximum(cR_max,0)/S
nrm=np.linalg.norm(A_,axis=0)
if os.environ.get("METHOD","drop")=="floor":
    vmin=0.7*rho*4/3*np.pi*0.5**3/S                          # smallest void deficit (R = 0.5 Gpc)
    sel=[]; sgn=[]; x=np.zeros(0); res=T.copy(); traj=[]; t0=time.time(); stall=False
    # each picked body gets its type (cluster / void) at pick time and keeps a real floor: cluster >= 1e49 kg, void >= 0.5 Gpc
    while True:
        c=(A_.T@res)/np.maximum(nrm,1e-300); c[sel]=0
        c=np.where((c>0)&(hi>MLO/S),c,np.where((c<0)&(lo<-vmin),-c,0)); j=int(np.argmax(c))
        if c[j]<=0: stall=True; break
        sel.append(j); sgn.append(1 if (A_[:,j]@res)>0 else -1); sg=np.array(sgn); si=np.array(sel)
        lb=np.where(sg>0,MLO/S,lo[si]); ub=np.where(sg>0,hi[si],-vmin[si])
        r=lsq_linear(A_[:,si],T,bounds=(lb,ub),method="bvls" if len(sel)<300 else "trf",lsmr_tol="auto",max_iter=3000); x=r.x
        Lold=L if len(sel)>1 else 1e9; res=T-A_[:,si]@x; L=rms(res)/FLOOR
        if (len(sel) in (1,2,4,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024) or L<1): traj.append((len(sel),float(L)))
        if len(sel)%25==0: print(f'  {tag} N={len(sel)} left={L:.2f} t={time.time()-t0:.0f}s',flush=True)
        if L<1 or len(sel)>=NMAX: break
    npick=len(sel)
else:
    sel=[]; x=np.zeros(0); res=T.copy(); traj=[]; t0=time.time(); ban=np.zeros(len(lo),bool); npick=0; stall=False; bestL=1e9; bsel=[]; bx=np.zeros(0)
    # greedy pick + bounded joint refit; a body whose fitted mass falls below 1e49 kg is removed and not offered again
    while True:
        c=(A_.T@res)/np.maximum(nrm,1e-300); c[sel]=0; c[ban]=0
        c=np.where((c>0)&(hi>0),c,np.where((c<0)&(lo<0),-c,0)); j=int(np.argmax(c))
        if c[j]<=0: stall=True; break
        sel.append(j); r=lsq_linear(A_[:,sel],T,bounds=(lo[sel],hi[sel]),method="bvls",max_iter=2000)
        x=r.x; keep=np.abs(x*S)>=MLO; npick+=1
        ban[np.array(sel)[~keep]]=True; sel=list(np.array(sel)[keep]); x=x[keep]
        for _ in range(6):       # refit after removals until every kept body is above the floor
            if not ((~keep).any() and len(sel)): break
            x=lsq_linear(A_[:,sel],T,bounds=(lo[sel],hi[sel]),method="bvls").x; keep=np.abs(x*S)>=MLO
            ban[np.array(sel)[~keep]]=True; sel=list(np.array(sel)[keep]); x=x[keep]
        res=T-(A_[:,sel]@x if len(sel) else 0); L=rms(res)/FLOOR
        if L<bestL: bestL,bsel,bx=L,list(sel),x.copy()
        if npick>6000: stall=True; break
        if (len(sel) in (1,2,4,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024) or L<1) and (not traj or traj[-1][0]!=len(sel)): traj.append((len(sel),float(L)))
        if L<1 or len(sel)>=NMAX: break
if os.environ.get('METHOD','drop')!='floor' and bestL<1e9: sel,x=bsel,bx; res=T-A_[:,sel]@x   # report the best state reached
X=P[sel]; M=x*S; live=np.abs(M)>=MLO
out=dict(tag=tag,seed=seed,nbank=int(len(P)),background=spec(T),traj=traj,N=int(len(sel)),N_live=int(live.sum()),left=float(rms(res)/FLOOR),reached=bool(rms(res)/FLOOR<1),
         n_mass=int((M>=MLO).sum()),n_void=int((M<=-MLO).sum()),inside_us=int((np.linalg.norm(X[live],axis=1)<RMIN).sum()),in_core=int((core_depth(X[live])>0).sum()),
         near_ring_2gpc=int((np.sqrt((np.hypot(X[live,0],X[live,1])-RMIN)**2+X[live,2]**2)<2).sum()),
         M_total_pos=float(M[M>0].sum()),M_total_neg=float(M[M<0].sum()),final=spec(res),resid=res.tolist(),X=X.tolist(),M=M.tolist(),secs=time.time()-t0)
os.makedirs(f"{WEBDIR}/web_add",exist_ok=True); out["ring_clearance_gpc"]=RX; json.dump(out,open(f"{WEBDIR}/web_add/{tag}_omp{seed}{os.environ.get('OSUF','')}.json","w"))   # seed 0 = first pass (drop rule), seed 1 = floor-respecting
os.remove(_l) if os.path.exists(_l) else None
print("DONE",tag,"stalled" if stall or len(sel)>=NMAX else "","N",out["N"],"live",out["N_live"],"left %.3f"%out["left"],"masses",out["n_mass"],"voids",out["n_void"],"inside",out["inside_us"],"core",out["in_core"],"nearring",out["near_ring_2gpc"],"traj",traj,"t %.0f"%out["secs"],flush=True)
