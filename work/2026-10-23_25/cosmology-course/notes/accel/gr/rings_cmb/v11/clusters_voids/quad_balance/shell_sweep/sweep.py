import sys,os,json,numpy as np
sys.path.insert(0,'/workspace/clump_scan/sequential');H=os.getcwd();os.chdir('/workspace/clump_scan/sequential');import seq_core as S;os.chdir(H)
th=np.radians(np.arange(360));F=S.FLOOR
def mm(v):m=np.arange(1,13)[:,None];a=np.hypot(2*np.mean(v*np.cos(m*th),-1),2*np.mean(v*np.sin(m*th),-1));return a/np.sqrt(2)/F
f=np.load('/workspace/clump_scan/sequential_open/final_comp_s3.npz');X=f['Xp'];M=f['Mp']
d=np.linalg.norm(X,axis=1);inn=d<S.RV;out=~inn
print('N',len(X),'inner',inn.sum(),'outer',out.sum(),'Mout net',M[out].sum(),'Minn net',M[inn].sum(),'outer d range',d[out].min(),d[out].max())
Pin=S.pattern(X[inn],M[inn]);Pout=S.pattern(X[out],M[out])
def vmf(mu,kappa,K,rng):   # Wood sampler, mu (3,)
    if kappa==0: v=rng.normal(size=(K,3));return v/np.linalg.norm(v,axis=1,keepdims=True)
    u=rng.uniform(size=K);w=1+np.log(u+(1-u)*np.exp(-2*kappa))/kappa
    v=rng.normal(size=(K,3));v-=np.outer(v@mu,mu);v/=np.linalg.norm(v,axis=1,keepdims=True)
    return w[:,None]*mu+np.sqrt(np.clip(1-w*w,0,None))[:,None]*v
def fib(K):
    i=np.arange(K)+.5;z=1-2*i/K;p=np.pi*(1+5**.5)*i;s=np.sqrt(1-z*z);return np.c_[s*np.cos(p),s*np.sin(p),z]
FB=fib(400)
def smoothed(center,sig,K=400,seed=0):
    rng=np.random.default_rng(seed);Y=X[out]-center;s=np.linalg.norm(Y,axis=1);U=Y/s[:,None];pts=[];ms=[]
    for j in range(len(Y)):
        if sig is None:
            R_=np.linalg.qr(rng.normal(size=(3,3)))[0];D=FB@R_.T     # randomly rotated Fibonacci sphere: exact uniform shell
        else: D=vmf(U[j],1/np.radians(sig)**2,K,rng)
        pts.append(center+s[j]*D);ms.append(np.full(len(D),M[out][j]/len(D)))
    return S.pattern(np.vstack(pts),np.concatenate(ms))
cen={'us':np.zeros(3),'CoE':np.asarray(S.C_G,float)}
SIG=[0,5,10,20,40,80,'uniform'];rows=[]
base=dict(full_web=mm(S.R-Pin-Pout),no_web=mm(S.R),inner_only=mm(S.R-Pin))
for k,v in base.items():print(k,'m2 %.2f'%v[1])
for cn,c in cen.items():
    for sg in SIG:
        if sg==0: Po=Pout;Po2=Pout
        else:
            sgv=None if sg=='uniform' else sg
            Po=smoothed(c,sgv,seed=1);Po2=smoothed(c,sgv,seed=2)
        res=S.R-Pin-Po;a=mm(res);ao=mm(Po);noise=mm(Po-Po2)[1]/np.sqrt(2)
        r=dict(center=cn,sigma_deg=sg,m2=a[1],total=S.rms(res)/F,per_m=a.tolist(),outer_only_m2=ao[1],outer_only_rms=S.rms(Po)/F,mc_noise_m2=noise)
        rows.append(r);print(cn,sg,'m2 %.3f total %.2f | outer-only m2 %.3f rms %.3f (MC noise m2 %.3f)'%(r['m2'],r['total'],ao[1],r['outer_only_rms'],noise),flush=True)
json.dump(dict(base={k:v.tolist() for k,v in base.items()},rows=rows,counts=dict(inner=int(inn.sum()),outer=int(out.sum()))),open('shell_sweep.json','w'),indent=1,default=float)
