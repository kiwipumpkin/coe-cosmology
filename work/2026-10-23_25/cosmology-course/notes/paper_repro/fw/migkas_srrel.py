"""SR framework, relativistic: front layer launched at beta_f (fixed per run), our layer at beta_o (free, relativistic).
Layers: proper speed p = q*U, U = p(beta_o); q from 0.05 to q_front = p(beta_f)/U. Sky circle must lie inside the front."""
import numpy as np,json,sys
from scipy.optimize import minimize
from multiprocessing import Pool
import core,dyn
TH6=np.arange(0,360,6.0); TH1=np.arange(0,360,1.0)
BF=float(sys.argv[1]) if len(sys.argv)>1 else 0.99
RMIN=float(sys.argv[2]) if len(sys.argv)>2 else 0.0   # minimum launch beta_obs/beta_front (observer near front)
SUF='' if RMIN==0 else f'_near{RMIN}'
pb=lambda b: b/np.sqrt(1-b*b)
def model(p,BF=BF):
    try: MT,LD,Rb=core.profile(p['DR'],p['alpha'])
    except Exception: return None
    U=pb(p['bo']); qf=pb(BF)/U
    if qf<=1: return None
    to=dyn.sr_tobs(U,p['K'],MT)
    if to is None: return None
    QG=np.geomspace(0.05,qf,400); r,pp=dyn.sr_integrate(QG,U,p['K'],MT,to); v=pp/np.sqrt(1+pp*pp)
    rho=1/p['DR']
    if not (np.all(np.isfinite(v)) and np.all(np.diff(r)>0) and v.min()>0) or r[0]>1-rho or r[-1]<1+rho: return None
    return (lambda x: np.interp(x,r,v)),dict(t_obs=to,qf=qf,U=U,R_front=float(r[-1]),beta_front_now=float(v[-1]),beta_obs_now=float(np.interp(1,r,v)),LD=LD)
def chi2(p,th):
    m=model(p)
    if m is None: return 800.0,None
    af=core.af_from_profile(m[0],p['DR'],th,relativistic=True); c,mb,full=core.score(af,th,p['l0']); return c,(mb,full,m[1])
def unpack(u):
    p=dict(DR=10**u[0],K=10**u[1],bo=1/(1+np.exp(-u[2])),alpha=10**u[3],l0=u[4]%360)
    if not (1.05<p['DR']<200 and 0.3<p['alpha']<8 and p['bo']>=max(0.5,RMIN*BF) and p['bo']<BF): return None
    return p
def obj(u,th):
    p=unpack(u); return 800.0 if p is None else chi2(p,th)[0]
def _run(u0):
    r=minimize(lambda u: obj(u,TH6),u0,method='Nelder-Mead',options=dict(maxiter=600,xatol=1e-4,fatol=1e-4,adaptive=True)); return r.fun,r.x.tolist()
if __name__=='__main__':
    lg=lambda b: np.log(b/(1-b))
    S=[[np.log10(DR),np.log10(K),lg(bo),np.log10(al),20.] for DR in ((4,10,30) if RMIN==0 else (12,25,60)) for K in (0.01,0.3) for bo in ((0.6,0.85) if RMIN==0 else (RMIN*BF+0.3*(BF-RMIN*BF),RMIN*BF+0.7*(BF-RMIN*BF))) for al in (0.5,)]
    S=[s for s in S if unpack(s) is not None]
    with Pool(8) as P: R=P.map(_run,S)
    R.sort(key=lambda z:z[0]); b6=R[0]
    r1=minimize(lambda u: obj(u,TH1),b6[1],method='Nelder-Mead',options=dict(maxiter=800,xatol=1e-5,fatol=1e-5,adaptive=True))
    out=dict(beta_front_launch=BF,min_beta_ratio=RMIN)
    for tag,u,th in (('6deg',b6[1],TH6),('1deg',r1.x,TH1)):
        p=unpack(u)
        if p is None: out[tag]=dict(error='no valid fit'); continue
        c,res=chi2(p,th)
        if res is None: out[tag]=dict(error='no valid fit',chi2=c); continue
        mb,full,info=res; sig,_=core.sigma_from_delta(core.CHI2_NULL-c,5)
        out[tag]=dict(params=p,chi2=c,R2=core.r2u(core.Y,mb),sigma=float(sig),n_free=5,D_Gpc=p['DR']*core.RV,curve=full.tolist(),
                      age_Gyr=info['t_obs']*p['DR']*core.RV/core.C_GPC_GYR,**{k:v for k,v in info.items()})
    json.dump(out,open(f'migkas_srrel_{BF}{SUF}.json','w'),indent=1,default=float)
    print(BF,{t:{k:out[t].get(k) for k in ('R2','chi2','sigma','D_Gpc','params','beta_obs_now','beta_front_now','R_front','age_Gyr')} for t in ('6deg','1deg')},flush=True)
