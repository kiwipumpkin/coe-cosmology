"""rec_fit6: visible circle proper R=c t* (NOT stretched). Stretch a(now)/a(emit), a=D/R_rec, enters ONLY as wavelength factor.
c=1,R_rec=1. D''=-mu/D^2, D(0)=R_rec, D'(0)=U0>1 (t=0: we cross R_rec). Emitters on circle |x-D n0|=R, emitted when they crossed R_rec,
travel time t* every direction. Variants: A stretch on / B stretch off; TD on (grav clock+gamma) / off. LOS Doppler always on."""
import json,sys,numpy as np
sys.path.insert(0,'.'); from lc import F,sig
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution,minimize
from scipy.special import expit
from functools import lru_cache
@lru_cache(maxsize=64)
def traj(U0,mu,ts):
    s=solve_ivp(lambda t,y:[y[1],-mu/y[0]**2],(0,ts),[1.,U0],dense_output=True,rtol=1e-10,atol=1e-12); return s.sol
def terms(p,ell,d=None,var=('A',1)):
    lU,fm,lt=p; U0=1+np.exp(lU); ts=np.exp(lt);
    # mu as fraction of the max allowed: 2mu < min emitter radius (D-R)
    so=traj(U0,0.,ts); Dg=so(ts)[0]
    mu=expit(fm)*min(0.5*(U0**2-1),0.5*max(Dg-ts,1e-9))  # U>c throughout and 2mu<r_e are the stated rules; checked below
    sol=traj(U0,mu,ts); Do,Uo=sol(ts); R=ts if d is None else d
    if Do-R<=2*mu or Uo<=1 or U0**2-2*mu*(1-1/Do)<=0: return None
    th=np.radians(ell); n=np.c_[np.cos(th),np.sin(th)]; E=np.array([Do,0])+R*n; re=np.linalg.norm(E,axis=1)
    if np.any(2*mu/re>=1): return None
    u=np.sqrt(U0**2-2*mu*(1-1/re)); vp=((u-Uo)/re)[:,None]*E; vn=np.sum(vp*n,1); v2=np.sum(vp**2,1)
    if np.any(v2>=1): return None
    te=ts-R; st=Do/sol(te)[0] if var[0]=='A' else 1.0
    g=np.sqrt((1-2*mu/Do)/(1-2*mu/re)) if var[1] else np.ones_like(re)
    gam=1/np.sqrt(1-v2) if var[1] else np.ones_like(re); dop=1+vn
    return dict(stretch=st*np.ones_like(re),grav=g,gamma=gam,dop=dop,onez=st*g*gam*dop,mu=mu,U0=U0,ts=ts,D=Do,U=Uo,H=Uo/Do,R=R,r_e=re)
def harm(T,ell):
    th=np.radians(ell); X=np.c_[np.ones(len(th)),np.cos(th),np.sin(th)]; f=T/T.mean()-1
    c=np.linalg.lstsq(X,f,rcond=None)[0]; r=f-X[:,1:]@c[1:]
    return dict(dip=float(np.hypot(c[1],c[2])),quad=float(2*np.mean(f*np.cos(2*th))),rms_nodip=float(np.sqrt(np.mean(r**2))),f=f,r=r)
def migkas(p,var,step):
    ell=np.arange(0,360,float(step)); r0=terms(p,ell,None,var)
    if r0 is None: return None
    dg=np.linspace(1e-3,r0['ts']*0.999,120); Z=[]
    for d in dg:
        q=terms(p,ell,d,var)
        if q is None: return None
        Z.append(q['onez']-1)
    Z=np.array(Z); zm=Z.mean(1); afs=[]
    for zk in F.ZMIG:
        if zm.max()>=zk and np.all(np.diff(zm)>0): i=int(np.searchsorted(zm,zk)); i=min(max(i,1),len(dg)-1)
        else: i=max(1,int(round(zk/F.ZMIG[-1]*(len(dg)-1))))   # fallback when z never reaches zk (variant B)
        A=Z[i]/dg[i]; afs.append(A/np.mean(A))
    af=np.mean(afs,0); af=np.interp(np.arange(360.),np.append(ell,360),np.append(af,af[0])); return af/af.mean()
WT,WC=0.01,0.01
def obj(p,var):
    try:
        q=terms(p,np.arange(0,360,6.),None,var); af=migkas(p,var,6)
        if q is None or af is None or np.any(af<=0) or not np.all(np.isfinite(af)): return 1e6
        T=3000/q['onez']; h=harm(T,np.arange(0,360,6.))
        return F.chi2_l0(af)[0]+(np.log(T.mean()/2.725)/WT)**2+(h['rms_nodip']/WC)**2
    except Exception: return 1e6
BOX=[(-3,4),(-6,6),(np.log(10),np.log(3000))]
out={}
for var in (('A',1),('A',0),('B',1),('B',0)):
    k=f"{var[0]}_{'TD' if var[1] else 'noTD'}"
    r=differential_evolution(obj,BOX,args=(var,),seed=3,popsize=12,maxiter=60,tol=1e-8,polish=False)
    r2=minimize(obj,r.x,args=(var,),method='Nelder-Mead',options=dict(maxiter=3000)); x=r2.x if r2.fun<r.fun else r.x
    res=dict(x=x.tolist(),obj=float(min(r.fun,r2.fun)))
    for st in (6,1):
        ell=np.arange(0,360,float(st)); q=terms(x,ell,None,var); af=migkas(x,var,st); c2,l0=F.chi2_l0(af)
        T=3000/q['onez']; h=harm(T,ell)
        res[f'{st}deg']=dict(chi2=c2,sigma=sig(c2,6),l0=l0,Tmean=float(T.mean()),Tmin=float(T.min()),Tmax=float(T.max()),dipole=h['dip'],quad=h['quad'],rms_nodip=h['rms_nodip'],
          T=T.tolist(),f=h['f'].tolist(),r=h['r'].tolist(),af=af.tolist(),
          terms={t:[float(q[t].min()),float(q[t].mean()),float(q[t].max())] for t in ('stretch','grav','gamma','dop','onez')})
    res['params']={t:float(q[t]) for t in ('U0','mu','ts','D','U','H','R')}; res['params']['2mu_over_min_re']=float(2*q['mu']/q['r_e'].min())
    out[k]=res; print(k,res['obj'],res['params'],{s:(res[s]['chi2'],res[s]['sigma'],res[s]['Tmean'],res[s]['dipole'],res[s]['quad']) for s in ('6deg','1deg')},flush=True)
json.dump(out,open('rec_fit6.json','w'),indent=1,default=float)
