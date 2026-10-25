# Direct-summation leapfrog (KDK) N-body, Plummer softening eps_ij^2=r_i^2+r_j^2, signed masses (voids negative, a_i independent of own sign).
import sys,os,json,time,numpy as np,numba as nbm
sys.path.insert(0,'/workspace/clump_scan/sequential');H=os.getcwd();os.chdir('/workspace/clump_scan/sequential');import seq_core as S;os.chdir(H)
th=np.radians(np.arange(360)) if not hasattr(S,'TH') else S.TH;F=S.FLOOR
G=6.674e-11;GPC=S.GPC;MCOE=float(sys.argv[1]);TAG=sys.argv[2];DT=float(os.environ.get('DT','0.005'));TEND=float(os.environ.get('TEND','3.0'))
def mm(v):m=np.arange(1,61)[:,None];a=np.hypot(2*np.mean(v*np.cos(m*th),-1),2*np.mean(v*np.sin(m*th),-1));return a/np.sqrt(2)/F
f=np.load('/workspace/clump_scan/sequential_open/final_comp_s3.npz');X=f['Xp'].astype(float).copy();r=f['rp'].astype(float);M=f['Mp'].astype(float);n=len(M)
K=G/GPC**2*M  # acceleration units Gpc-scaled: a [m/s^2] = K_j * d[Gpc]/q[Gpc^3]... (G M / GPC^2)
@nbm.njit(parallel=True,fastmath=True)
def acc(X,K,r2):
    n=X.shape[0];A=np.zeros((n,3))
    for i in nbm.prange(n):
        ax=ay=az=0.0
        for j in range(n):
            if j==i:continue
            dx=X[j,0]-X[i,0];dy=X[j,1]-X[i,1];dz=X[j,2]-X[i,2];q=dx*dx+dy*dy+dz*dz+r2[i]+r2[j];w=K[j]/(q*np.sqrt(q))
            ax+=w*dx;ay+=w*dy;az+=w*dz
        A[i,0]=ax;A[i,1]=ay;A[i,2]=az
    return A
@nbm.njit(parallel=True,fastmath=True)
def pot(X,K,r2,M):
    n=X.shape[0];U=np.zeros(n)
    for i in nbm.prange(n):
        s=0.0
        for j in range(i+1,n):
            dx=X[j,0]-X[i,0];dy=X[j,1]-X[i,1];dz=X[j,2]-X[i,2];s-=M[i]*K[j]/np.sqrt(dx*dx+dy*dy+dz*dz+r2[i]+r2[j])
        U[i]=s
    return U.sum()
C=np.array(S.C_G);KC=G*MCOE/GPC**2;g0=C/np.linalg.norm(C)**3
def acoe(X):d=C-X;return KC*(d/np.linalg.norm(d,axis=1)[:,None]**3-g0)
def Ucoe(X):return float(np.sum(M*(-KC/np.linalg.norm(C-X,axis=1)+KC*(X@g0))))
r2=r*r;A=lambda X:acc(X,K,r2)+(acoe(X) if MCOE>0 else 0)
a0=A(X);arms=np.sqrt(np.mean(np.sum(a0**2,1)))
TAU=np.sqrt(1.0/arms)   # time unit: rms body moves ~1 Gpc in ~1 tau under initial accel (Gpc, s-scaled)
ti=0.1;V=a0*ti*TAU      # Zel'dovich-consistent start: v = g * t_i (displacement ∝ g), t_i = 0.1 tau
dt=DT*TAU
def energy(X,V):return float(0.5*np.sum(M*np.sum(V*V,1))+pot(X,K,r2,M)+(Ucoe(X) if MCOE>0 else 0))
def mom(V):return np.sum(M[:,None]*V,0)
E0=energy(X,V);P0=mom(V);Pscale=np.sum(np.abs(M)[:,None]*np.abs(V),0).sum()
X0=X.copy();out=dict(MCOE=MCOE,tau_s=TAU,tau_Gyr=TAU/3.156e16,t_i=ti,dt_tau=DT,rows=[])
def rec(t,X,V):
    v=S.R-S.pattern(X,M);a=mm(v);d=np.linalg.norm(X-X0,axis=1)
    from scipy.spatial import cKDTree;p=np.array(list(cKDTree(X).query_pairs(2*r.max())));ov=int(np.sum(np.linalg.norm(X[p[:,0]]-X[p[:,1]],axis=1)<r[p[:,0]]+r[p[:,1]])) if len(p) else 0
    E=energy(X,V);P=mom(V);row=dict(t_tau=t,total=float(S.rms(v)/F),m2=float(a[1]),per_m=a[:12].tolist(),m13_60=float(np.sqrt(np.sum(a[12:]**2))),rms_disp=float(np.sqrt(np.mean(d*d))),max_disp=float(d.max()),overlaps=ov,dE_rel=float((E-E0)/abs(E0)),dP_rel=float(np.linalg.norm(P-P0)/Pscale))
    out['rows'].append(row);print(TAG,'t=%.2f disp rms %.3f max %.2f  total %.1f m2 %.2f m1 %.1f m3-12max %.1f ov %d dE %.2e dP %.2e'%(t,row['rms_disp'],row['max_disp'],row['total'],a[1],a[0],np.max(a[2:12]),ov,row['dE_rel'],row['dP_rel']),flush=True)
rec(0,X,V);a=a0;t=0.0;nst=int(round(TEND/DT));every=int(round(0.1/DT))
for k in range(1,nst+1):
    V=V+0.5*dt*a;X=X+dt*V;a=A(X);V=V+0.5*dt*a;t=k*DT
    if k%every==0:rec(t,X,V)
json.dump(out,open(f'nb_{TAG}.json','w'),indent=1);np.savez(f'nb_{TAG}_final.npz',X=X,V=V)
