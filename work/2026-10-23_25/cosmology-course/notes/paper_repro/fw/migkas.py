import numpy as np,json,sys
from scipy.optimize import minimize
from multiprocessing import Pool
import core,dyn
QG=np.geomspace(0.05,3.0,300)
TH6=np.arange(0,360,6.0); TH1=np.arange(0,360,1.0)
def model(fw,p):
    """returns (vprof, info) or None. p dict."""
    try: MT,LD,Rb=core.profile(p['DR'],p['alpha'])
    except Exception: return None
    if fw=='classical':
        V0=dyn.cl_V0(p['mu'],MT); r,v=dyn.cl_integrate(QG,V0,p['mu'],MT)
        info=dict(V0=V0,t_obs=1.0)
    elif fw=='sr':
        to=dyn.sr_tobs(p['U'],p['K'],MT)
        if to is None: return None
        r,pp=dyn.sr_integrate(QG,p['U'],p['K'],MT,to); v=pp/np.sqrt(1+pp*pp); info=dict(t_obs=to)
    else:
        r=np.linspace(0.3*(1-1/p['DR']),1.5*(1+1/p['DR']),400); M=(p['K']/2)*np.interp(r,core.RG,MT)
        E,ok=dyn.ltb_E(r,M,p['tau'])
        if not ok.all() or np.any(1+2*E<=0): return None
        v=np.sqrt(2*M/r+2*E); info=dict(t_obs=p['tau'],min_1p2E=float((1+2*E).min()))
    if not (np.all(np.isfinite(v)) and np.all(np.diff(r)>0) and v.min()>0): return None
    rho=1/p['DR']
    if r[0]>1-rho or r[-1]<1+rho: return None
    info.update(MT=MT,LD=LD,Rb=Rb); return (lambda x: np.interp(x,r,v)),info
def chi2(fw,p,th):
    m=model(fw,p)
    if m is None: return 800.0,None
    af=core.af_from_profile(m[0],p['DR'],th,relativistic=(fw=='sr'))
    c,mb,full=core.score(af,th,p['l0']); return c,(mb,full,m[1])
NAMES=dict(classical=['DR','mu','alpha','l0'],sr=['DR','K','U','alpha','l0'],gr=['DR','K','tau','alpha','l0'])
def unpack(fw,u):
    p={}
    for k,x in zip(NAMES[fw],u):
        p[k]= (x%360) if k=='l0' else 10**x
    if not (1.05<p['DR']<10 and 0.3<p['alpha']<8): return None
    return p
def obj(fw,u,th):
    p=unpack(fw,u)
    if p is None: return 800.0
    return chi2(fw,p,th)[0]
def _run(args):
    fw,u0=args
    r6=minimize(lambda u: obj(fw,u,TH6),u0,method='Nelder-Mead',options=dict(maxiter=500,xatol=1e-4,fatol=1e-4,adaptive=True))
    return r6.fun,r6.x.tolist()
def starts(fw):
    S=[]
    for DR in (1.5,2.0,2.7,3.5,4.5):
        for l0 in (15.,25.):
            for al in (0.5,1.5):
                if fw=='classical':
                    for mu in (0.3,2.0): S.append([np.log10(DR),np.log10(mu),np.log10(al),l0])
                elif fw=='sr':
                    for K,U in ((0.3,0.4),(1.,0.8)): S.append([np.log10(DR),np.log10(K),np.log10(U),np.log10(al),l0])
                else:
                    for K,tau in ((0.3,1.5),(2.,0.5)): S.append([np.log10(DR),np.log10(K),np.log10(tau),np.log10(al),l0])
    return S
if __name__=='__main__':
    fw=sys.argv[1]; S=starts(fw)
    with Pool(8) as P: R=P.map(_run,[(fw,u) for u in S])
    R.sort(key=lambda z:z[0]); best6=R[0]
    r1=minimize(lambda u: obj(fw,u,TH1),best6[1],method='Nelder-Mead',options=dict(maxiter=800,xatol=1e-5,fatol=1e-5,adaptive=True))
    out={}
    nf=len(NAMES[fw])
    for tag,u,th in (('6deg',best6[1],TH6),('1deg',r1.x,TH1)):
        p=unpack(fw,u); c,(mb,full,info)=chi2(fw,p,th); sig,pv=core.sigma_from_delta(core.CHI2_NULL-c,nf)
        out[tag]=dict(params=p,chi2=c,R2=core.r2u(core.Y,mb),sigma=float(sig),n_free=nf,L_over_D=info['LD'],t_obs=info['t_obs'],
                      D_Gpc=p['DR']*core.RV,curve=full.tolist(),extra={k:v for k,v in info.items() if k in ('V0','min_1p2E')})
    out['top_starts_6deg']=[z[0] for z in R[:5]]
    json.dump(out,open(f'migkas_{fw}.json','w'),indent=1,default=float)
    print(fw,{t:(out[t]['R2'],out[t]['chi2'],out[t]['sigma'],out[t]['params']) for t in ('6deg','1deg')})
