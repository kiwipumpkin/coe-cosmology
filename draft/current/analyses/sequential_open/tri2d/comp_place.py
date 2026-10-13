"""Compensated placement (Chris, 10 Oct): each void (radius R, contrast delta_v) has deficit M_def = |delta_v| rho (4/3 pi R^3);
clusters (radius r_c, contrast delta_c, excess mass M_c = delta_c rho (4/3 pi r_c^3)) are added around it until their mass equals
M_def (last one trimmed), at radius R..1.3R (compensation wall, Hamaus / Clampitt-Jain profile). 2D periodic box (a layer through
the voids). Then a constrained ground-state search: Metropolis anneal of U = -G sum M_i M_j / r_ij (all pairs, minimum image)
with moves = cluster anywhere in its own void's wall, or whole void + its clusters translated; no overlaps."""
import json,numpy as np,numba
from scipy.spatial import Delaunay
G=6.674e-11; GPC=3.0857e25; V=lambda r: 4/3*np.pi*r**3
RHO=8.39e50/(0.70*V(1.5)*GPC**3)          # kg/m^3, from the web's median void (R = 1.5 Gpc, delta_v = -0.70)
rng=np.random.default_rng(21); L=85.0
def draw_void(): return float(np.exp(rng.uniform(np.log(0.75),np.log(3.0)))),float(rng.uniform(-0.9,-0.5))
def draw_cluster(): return float(np.exp(rng.uniform(np.log(0.10),np.log(0.25)))),float(np.exp(rng.uniform(np.log(30),np.log(300))))   # delta_c 30-300: median M_c ~ web's 1.4e50 kg
def mind(a,b): d=a-b; return d-L*np.round(d/L)
VO=[]; CC=np.zeros((0,2)); RR=np.zeros(0)   # random sequential addition of void+wall discs (radius 1.3R, gap 0.5 Gpc)
for t in range(60000):
    R,dv=draw_void(); c=rng.uniform(0,L,2); d=c-CC; d-=L*np.round(d/L)
    if np.all(np.hypot(d[:,0],d[:,1])>1.3*R+1.3*RR+0.5): VO.append((c,R,dv)); CC=np.vstack([CC,c]); RR=np.r_[RR,R]
P=[];K=[];M=[];RAD=[];OWN=[];rows=[]
for vi,(c,R,dv) in enumerate(VO):
    Md=abs(dv)*RHO*V(R)*GPC**3; P.append(c); K.append(1); M.append(-Md); RAD.append(R); OWN.append(vi); acc=0.0; n=0; mine=[]; tries=0
    while acc<Md*(1-1e-9) and tries<3000:
        tries+=1; rc,dc=draw_cluster(); mc=dc*RHO*V(rc)*GPC**3
        if acc+mc>Md: mc=Md-acc; dc=mc/(RHO*V(rc)*GPC**3)
        th=rng.uniform(0,2*np.pi); rr=rng.uniform(R,1.3*R); x=(c+rr*np.array([np.cos(th),np.sin(th)]))%L
        if all(np.hypot(*mind(x,P[j]))>rc+RAD[j] for j in mine): P.append(x); K.append(0); M.append(mc); RAD.append(rc); OWN.append(vi); mine.append(len(P)-1); acc+=mc; n+=1
    rows.append(dict(R=R,delta_v=dv,M_def=Md,M_held=acc,n_clusters=n,mean_Mc=acc/max(n,1),wall_r_mean=float(np.mean([np.hypot(*mind(P[j],c)) for j in mine])/R)))
P=np.array(P);K=np.array(K);M=np.array(M);RAD=np.array(RAD);OWN=np.array(OWN); VC=np.array([v[0] for v in VO]); VR=np.array([v[1] for v in VO])
print("voids",len(VO),"clusters",int((K==0).sum()),"nodes",len(P),"compensation max resid",max(abs(sum(M[OWN==i])) for i in range(len(VO)))/RHO)
@numba.njit(cache=True)
def field_row(i,P,M,L):
    s=0.0
    for j in range(len(P)):
        if j==i: continue
        dx=P[j,0]-P[i,0]; dy=P[j,1]-P[i,1]; dx-=L*np.round(dx/L); dy-=L*np.round(dy/L); s+=M[j]/np.sqrt(dx*dx+dy*dy)
    return s
@numba.njit(cache=True)
def U_all(P,M,L):
    u=0.0
    for i in range(len(P)): u+=M[i]*field_row(i,P,M,L)
    return -0.5*u
@numba.njit(cache=True)
def ok_pos(i,x,P,RAD,K,OWN,L,skip_own):
    for j in range(len(P)):
        if j==i or (skip_own and OWN[j]==OWN[i]): continue
        dx=P[j,0]-x[0]; dy=P[j,1]-x[1]; dx-=L*np.round(dx/L); dy-=L*np.round(dy/L); d=np.sqrt(dx*dx+dy*dy)
        lim=RAD[i]+RAD[j]
        if K[j]==1 and K[i]==1: lim=1.3*RAD[i]+1.3*RAD[j]   # walls do not overlap
        if d<lim: return False
    return True
@numba.njit(cache=True)
def anneal(P,M,RAD,K,OWN,VIDX,L,T,nsteps,step,seed):
    np.random.seed(seed); N=len(P); nv=len(VIDX); acc=0
    for t in range(len(T)):
        beta=1.0/T[t]
        for s in range(nsteps):
            if np.random.random()<0.85:      # cluster move inside its own wall
                i=np.random.randint(N)
                if K[i]==1: continue
                v=VIDX[OWN[i]]; R=RAD[v]; th=2*np.pi*np.random.random(); rr=R+0.3*R*np.random.random()
                x=np.empty(2); x[0]=(P[v,0]+rr*np.cos(th))%L; x[1]=(P[v,1]+rr*np.sin(th))%L
                if not ok_pos(i,x,P,RAD,K,OWN,L,False): continue
                old=P[i].copy(); f0=field_row(i,P,M,L); P[i]=x; f1=field_row(i,P,M,L); dE=-M[i]*(f1-f0)
                if dE<=0 or np.random.random()<np.exp(-beta*dE): acc+=1
                else: P[i]=old
            else:                              # translate void + its clusters
                k=np.random.randint(nv); mem=np.where(OWN==k)[0]; d=step*(2*np.random.random(2)-1)
                old=P[mem].copy(); E0=0.0
                for i in mem: E0+=M[i]*field_row(i,P,M,L)
                for i in mem: P[i,0]=(P[i,0]+d[0])%L; P[i,1]=(P[i,1]+d[1])%L
                good=True
                for i in mem:
                    if not ok_pos(i,P[i],P,RAD,K,OWN,L,True): good=False; break
                if not good: P[mem]=old; continue
                E1=0.0
                for i in mem: E1+=M[i]*field_row(i,P,M,L)
                # internal pairs unchanged; external change = -(E1-E0) minus the internal double counting (cancels)
                dE=-(E1-E0)
                if dE<=0 or np.random.random()<np.exp(-beta*dE): acc+=1
                else: P[mem]=old
    return acc
VIDX=np.where(K==1)[0]; Mn=M/1.40e50       # mass in units of M_c(web median)
P0=P.copy(); U0=U_all(P,Mn,L); print("U start (G Mc^2/Gpc per node)",U0/len(P))
T=np.geomspace(20.0,0.05,30)
import time; t0=time.time(); acc=anneal(P,Mn,RAD,K,OWN,VIDX,L,T[:1],200,1.0,4); print('jit',time.time()-t0); t0=time.time(); acc=anneal(P,Mn,RAD,K,OWN,VIDX,L,T,15*len(P),1.0,5); print('anneal s',time.time()-t0); U1=U_all(P,Mn,L); print("accepted",acc,"U end per node",U1/len(P))
np.savez("comp_layout.npz",P0=P0,P=P,K=K,M=M,RAD=RAD,OWN=OWN,L=L)
json.dump(dict(rows=rows,RHO=RHO,U_start=U0,U_end=U1,n_nodes=len(P),n_voids=len(VO),n_clusters=int((K==0).sum())),open("comp_rows.json","w"),indent=1)
