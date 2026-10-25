import sys,os,json,numpy as np
MODE=sys.argv[1] if len(sys.argv)>1 else 'A'
from scipy.spatial import cKDTree
sys.path.insert(0,'/workspace/clump_scan/sequential');H=os.getcwd();os.chdir('/workspace/clump_scan/sequential');import seq_core as S;os.chdir(H)
th=np.radians(np.arange(360)) if not hasattr(S,'TH') else S.TH;F=S.FLOOR
def mm(v):m=np.arange(1,61)[:,None];a=np.hypot(2*np.mean(v*np.cos(m*th),-1),2*np.mean(v*np.sin(m*th),-1));return a/np.sqrt(2)/F
c2,s2=np.cos(2*th),np.sin(2*th)
f=np.load('/workspace/clump_scan/sequential_open/final_comp_s3.npz');X=f['Xp'].copy();r=f['rp'].copy();M=f['Mp'].copy();isV=f['kp']=='void'
res=S.R-S.pattern(X,M)
MCmax=np.abs(M[~isV]).max();MVmax=np.abs(M[isV]).max();DV=S.DELTA_V if hasattr(S,'DELTA_V') else None
print('cluster Mmax',MCmax,'void |M|max',MVmax,'DELTA_V',DV)
def rad(x,m,void):
    rho=S.rho_gpc(x);return (abs(m)/(abs(DV) if void else 1.0)/rho/S.GPC**3*3/(4*np.pi))**(1/3)   # void: |M|=|delta_v| rho V ; cluster: delta=+1 sphere (cluster_radius)
rng=np.random.default_rng(11);NC=300000
P=S.rand_dir(rng,NC)*(rng.uniform(0.5**3,30**3,NC)**(1/3))[:,None]
P=P[np.linalg.norm(P-S.C_G,axis=1)>=S.A_IN_G+0.5]
U=np.vstack([S.unit_pattern_fast(P[i:i+20000]) for i in range(0,len(P),20000)])
uc=2*U@c2/360;us=2*U@s2/360   # m=2 cos/sin coefficients per kg
def step(res,X,r):
    tc,ts=2*np.mean(res*c2),2*np.mean(res*s2);nn=uc**2+us**2;Mopt=(tc*uc+ts*us)/nn   # residual - M*u : want M*u = res part
    best=None
    for void,lo,hi in ((True,-MVmax,0.0),(False,0.0,MCmax)):
        Mb=np.clip(Mopt,lo,hi);gain=(tc**2+ts**2)-((tc-Mb*uc)**2+(ts-Mb*us)**2)
        order=np.argsort(-gain)[:4000];tree=cKDTree(X)
        for i in order:
            if gain[i]<=0:break
            x=P[i];rr=rad(x,Mb[i],void);d=np.linalg.norm(x);dr=np.sqrt((np.hypot(x[0],x[1])-14.3)**2+x[2]**2)
            if d-rr<0.5 or d+rr>30 or np.linalg.norm(x-S.C_G)-rr<S.A_IN_G or dr-rr<0.1:continue
            near=tree.query_ball_point(x,rr+r.max()+0.02)
            if near and np.any(np.linalg.norm(X[near]-x,axis=1)<rr+r[near]+0.02):continue
            if MODE=='B':
                cand=res-Mb[i]*U[i];ao=np.delete(mm(cand)[:12],1);a0=np.delete(mm(res)[:12],1)
                if np.max(ao)>max(np.max(a0),1.0)+1e-9 or S.rms(cand)>S.rms(res):continue
            if best is None or gain[i]>best[0]:best=(gain[i],i,Mb[i],rr,void)
            break
    return best
rows=[];a=mm(res);rows.append(dict(N=0,m2=float(a[1]),total=float(S.rms(res)/F),worst_other_1_12=float(np.max(np.delete(a[:12],1))),per_m=a[:12].tolist()))
print(0,rows[-1]['m2'],rows[-1]['total'],rows[-1]['worst_other_1_12'],flush=True);added=[];res0=res.copy()
for N in range(1,41):
    b=step(res,X,r)
    if b is None:print('no valid improving body');break
    g,i,m,rr,void=b;x=P[i];res=res-m*U[i];X=np.vstack([X,x]);r=np.r_[r,rr];added.append(dict(x=x.tolist(),M=float(m),r=float(rr),kind='void' if void else 'cluster'))
    a=mm(res);rows.append(dict(N=N,m2=float(a[1]),total=float(S.rms(res)/F),worst_other_1_12=float(np.max(np.delete(a[:12],1))),per_m=a[:12].tolist(),kind=added[-1]['kind'],M=float(m)))
    print(N,added[-1]['kind'],'%.2e'%m,'m2 %.3f total %.2f worst %.2f'%(rows[-1]['m2'],rows[-1]['total'],rows[-1]['worst_other_1_12']),flush=True)
    if a[1]<0.05:break
json.dump(dict(rows=rows,added=added),open(f'rta_{MODE}.json','w'),indent=1);np.savez(f'rta_{MODE}.npz',res0=res0,res=res,R=S.R)
