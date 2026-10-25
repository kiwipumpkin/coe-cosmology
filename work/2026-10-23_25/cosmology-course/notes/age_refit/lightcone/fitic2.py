import json,numpy as np; from ic2 import *; from lc import F,CC,sig
from scipy.optimize import differential_evolution,minimize
ELL6=np.arange(0,360,6.); G6=CC.G_FULL[::6]
def cmbR(p,step=6):
    ELL6=np.arange(0,360,float(step)); G6=CC.G_FULL[::step]; n=len(ELL6)
    b,q,mu,s,U,phi=p; r=S_ic2(b,q,mu,s,U,(ELL6-phi)%360)
    if r is None: return None
    th=np.radians(ELL6); X=np.c_[np.ones(n),np.cos(th),np.sin(th)]
    def rm(v): c=np.linalg.lstsq(X,v,rcond=None)[0]; return v-X[:,1:]@c[1:]
    S=r['S']; R=rm(G6*S/np.mean(G6*S)-1); base=rm(G6/G6.mean()-1); m=rm(S/S.mean()-1)
    rms=lambda v: float(np.sqrt(np.mean(v**2)))
    return dict(resid=rms(R),cancel=100*(1-rms(R)/rms(base)),dTT=rms(m),corr=float(np.corrcoef(m,base)[0,1]),S=S,R=R,m=m,base=base,T=2.725*S/S.mean(),Tabs_ratio=float(S.mean()),
                zU=r['zU'],Dm=r['Dm'],Dg=r['Dg'])
def migkas(p,step=6):
    b,q,mu,s,U=p[:5]; afs=[]
    for zk in F.ZMIG:
        r=S_ic2(b,zk/U,mu,s,U,np.arange(0,360,float(step)))
        if r is None: return None
        z=1/r['S']-1; d=r['tf']; afs.append((z/d)/np.mean(z/d))
    af=np.mean(afs,0); af=np.interp(np.arange(360.),np.append(np.arange(0,360,float(step)),360),np.append(af,af[0])); return af/af.mean()
def cm(p): r=cmbR(p); return 1e3 if r is None else r['resid']
def mg(p): af=migkas(p); return 1e3 if af is None or not np.all(np.isfinite(af)) or np.any(af<=0) else F.chi2_l0(af)[0]
BOX=[(1e-4,0.999),(1e-4,3),(-1,1),(0.2,20),(1.0001,30),(0,360)]   # DE search box only; NM polish unbounded
out={}
for name,f,bx in (('cmb',cm,BOX),('migkas',mg,BOX[:5])):
    r=differential_evolution(f,bx,seed=5,popsize=20,maxiter=300,tol=1e-10,polish=False)
    r2=minimize(f,r.x,method='Nelder-Mead',options=dict(maxiter=20000,xatol=1e-10,fatol=1e-12)); x=r2.x if r2.fun<r.fun else r.x
    out[name]=dict(x=list(map(float,x)),fun=float(min(r.fun,r2.fun))); json.dump(out,open('fitic2.json','w'),indent=1); print(name,out[name],flush=True)
# cross-evaluate at 6 deg and 1 deg
for name in ('cmb','migkas'):
    x=np.array(out[name]['x']); p6=list(x[:5])+[x[5] if len(x)>5 else 0.]
    if name=='migkas': chi,l0=F.chi2_l0(migkas(x,6)); p6[5]=l0
    for st in (6,1):
        c=cmbR(p6,st); af=migkas(p6,st); chi,l0=F.chi2_l0(af) if af is not None else (np.nan,np.nan)
        out[name][f'{st}deg']=dict(chi2=chi,sigma=sig(chi,6),l0=l0,cmb=None if c is None else {k:c[k] for k in ('resid','cancel','dTT','corr','Tabs_ratio')},
          T_obs_3000_range=None if c is None else [float(3000*c['S'].min()),float(3000*c['S'].max())],T_obs_3000_mean=None if c is None else float(3000*c['S'].mean()),
          zU_range=None if c is None else [float(c['zU'].min()),float(c['zU'].max())],curves=None if c is None else dict(m=c['m'].tolist(),base=c['base'].tolist(),R=c['R'].tolist()),af=None if af is None else af.tolist())
out['note']='optimiser on 6-deg grid; best fits re-scored at 6 and 1 deg'; json.dump(out,open('fitic2.json','w'),indent=1,default=float); print('DONE',flush=True)
