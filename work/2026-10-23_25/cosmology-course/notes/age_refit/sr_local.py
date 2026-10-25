import json,numpy as np,os,sys
os.environ['VMAX']=os.environ.get('VMAX','1'); law=sys.argv[1]; sys.argv=['x',law]
exec(open('fit_forward.py').read().split('x=fit()')[0])
from scipy.optimize import minimize_scalar
S0=CC.rms(CC.resid(CC.G_FULL,np.ones(360)))
def score(S,phi):
    Sl=np.interp((CC.ELL-phi)%360,np.arange(361.),np.append(S,S[0])); r=CC.resid(CC.G_FULL,Sl); return CC.rms(r),r
def sr_eval(e,q):
    c=F.cmb_sr_local(e['o'],q); S=c['S']
    if not np.all(np.isfinite(S)): return dict(ok=False,nfound=int(np.isfinite(S).sum()))
    rl,res=score(S,e['l0']); 
    L=np.arange(0,360,1.); v=[score(S,p)[0] for p in L]; pb=L[int(np.argmin(v))]
    m=minimize_scalar(lambda p:score(S,p)[0],bounds=(pb-1,pb+1),method='bounded'); rf,resf=score(S,m.x)
    fo=CC.fourier(res,6); lam=e['L']['Hbar']/H70
    return dict(ok=True,cancel_l0=100*(1-rl/S0),res_l0=rl/1e-5,cancel_phifree=100*(1-rf/S0),res_phifree=rf/1e-5,phi_best=float(m.x%360),
        modes_l0={m_['m']:m_['A'] for m_ in fo['modes']},modes_phifree={m_['m']:m_['A'] for m_ in CC.fourier(resf,6)['modes']},
        T_mean=float(2.725*np.mean(S)/np.mean(S)),S_mean=float(np.mean(S)),S_dipole=float(CC.fourier(S/S.mean()-1,2)['modes'][0]['A']),
        tLS_0=float(c['te'][0]*lam),tLS_180=float(c['te'][180]*lam),dtLS_0m180=float((c['te'][0]-c['te'][180])*lam))
if __name__=='__main__':
    out={}
    for s in ('_vc',''):
        f=f'fwd/{law}{s}_ckpt_free.json'
        if not os.path.exists(f): continue
        if s=='': os.environ['VMAX']='0'; VMAX=0.0
        z=np.array(json.load(open(f))['x']); e=ev(z,True)
        out[s or 'nolimit']={str(q):sr_eval(e,q) for q in (0.1,0.23,0.5,1.0,2.0)}
        json.dump(out,open(f'fwd/sr_{law}.json','w'),indent=1,default=float); print(law,s,{q:(round(d.get('cancel_l0',np.nan),1),round(d.get('cancel_phifree',np.nan),1)) for q,d in out[s or 'nolimit'].items()},flush=True)
