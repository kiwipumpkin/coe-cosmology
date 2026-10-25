import json,numpy as np,os,sys
os.environ['VMAX']='1'; sys.argv=['x','eta']
exec(open('fit_forward.py').read().split('x=fit()')[0])
from scipy.optimize import differential_evolution,minimize
S0=CC.rms(CC.resid(CC.G_FULL,np.ones(360)))
def full(z):
    e=ev(np.array(z),True)
    if e is None: return None
    c=F.cmb(e['o']); S=c['S']
    if not np.all(np.isfinite(S)): return dict(e=e,rms=np.inf)
    Sl=np.interp((CC.ELL-e['l0'])%360,np.arange(361.),np.append(S,S[0])); rms=CC.rms(CC.resid(CC.G_FULL,Sl)); return dict(e=e,rms=rms,cancel=100*(1-rms/S0),Tm=float(np.mean(c['Tobs'])))
def cc(z):
    r=full(z); return 1e3 if r is None or not np.isfinite(r['rms']) else r['rms']/S0
zA=np.array(json.load(open('fwd/eta_vc_ckpt_free.json'))['x'])
if not os.path.exists('fwd/bracket_B.json'):
    bb=[B[n] for n in names]
    r=differential_evolution(cc,bb,seed=3,popsize=8,maxiter=20,polish=False,init='sobol',x0=zA)
    r2=minimize(cc,r.x,method='Nelder-Mead',options=dict(maxiter=400,adaptive=True)); zB=r2.x if r2.fun<r.fun else r.x
    json.dump(dict(x=list(map(float,zB))),open('fwd/bracket_B.json','w'))
zB=np.array(json.load(open('fwd/bracket_B.json'))['x'])
def lin(z): z=z.copy(); z[0]=10**z[0]; z[1]=10**z[1]; z[4]=10**z[4]; return z   # linear-M, linear-(1-G_i), linear-L
def unlin(y): y=y.copy(); y[0]=np.log10(y[0]); y[1]=np.log10(y[1]); y[4]=np.log10(y[4]); return y
rows=[]
for path in ('log','linear'):
    for f in np.linspace(0,1,21):
        z=(1-f)*zA+f*zB if path=='log' else unlin((1-f)*lin(zA)+f*lin(zB))
        r=full(z); d=dict(path=path,f=float(f),z=list(map(float,z)))
        if r is None: d['ok']=False
        else:
            e=r['e'];L=e['L'];lam=L['Hbar']/H70
            d.update(ok=True,chi2=e['chi'],sigma=sig(e['chi'],7),R2=float(1-np.sum((F.Y-F.shifted(e['af'],e['l0'])[F.IMIG])**2)/np.sum((F.Y-F.Y.mean())**2)),
              l0=e['l0'],cancel=r.get('cancel',np.nan),rms_over_floor=r['rms']/1e-5,Tmean=r.get('Tm',np.nan),
              M_kg=e['P']['M']*lam,D=L['D']*lam,vo_c=L['Uo']/F.C,Hpar=L['Hpar']*F.KMSMPC/lam,Hper=L['Hper']*F.KMSMPC/lam,age=e['o']['t']*lam,gk=L['grav_over_kin'])
        rows.append(d); print({k:(round(v,4) if isinstance(v,float) else v) for k,v in d.items() if k!='z'},flush=True)
        json.dump(dict(names=names,zA=list(zA),zB=list(zB),rows=rows),open('fwd/bracket.json','w'),indent=1,default=float)
print('DONE')
