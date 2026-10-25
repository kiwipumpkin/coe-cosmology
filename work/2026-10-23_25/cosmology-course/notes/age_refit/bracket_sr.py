import os,sys,json,numpy as np;os.environ['VMAX']='1';sys.argv=['x','eta']
exec(open('sr_local.py').read().split("if __name__")[0])
from scipy.optimize import differential_evolution,minimize
Q=0.23
def cc(z):
    e=ev(np.array(z),True)
    if e is None: return 1e3
    c=F.cmb_sr_local(e['o'],Q); S=c['S']
    if not np.all(np.isfinite(S)): return 1e3
    return min(score(S,p)[0] for p in np.arange(0,360,3.))/S0
zA=np.array(json.load(open('fwd/eta_vc_ckpt_free.json'))['x'])
if not os.path.exists('fwd/bracket_sr_B.json'):
    r=differential_evolution(cc,[B[n] for n in names],seed=5,popsize=8,maxiter=25,polish=False,init='sobol',x0=zA)
    json.dump(dict(x=list(map(float,r.x)),stage='de'),open('fwd/bracket_sr_B.json','w'))
    r2=minimize(cc,r.x,method='Nelder-Mead',options=dict(maxiter=400,adaptive=True))
    zB=r2.x if r2.fun<r.fun else r.x; json.dump(dict(x=list(map(float,zB)),stage='nm'),open('fwd/bracket_sr_B.json','w'))
zB=np.array(json.load(open('fwd/bracket_sr_B.json'))['x'])
def lin(z): z=z.copy(); z[0]=10**z[0]; z[1]=10**z[1]; z[4]=10**z[4]; return z
def unlin(y): y=y.copy(); y[0]=np.log10(y[0]); y[1]=np.log10(y[1]); y[4]=np.log10(y[4]); return y
rows=[]
for f in np.linspace(0,1,21):
    z=unlin((1-f)*lin(zA)+f*lin(zB)); e=ev(z,True); d=dict(f=float(f),z=list(map(float,z)))
    if e is None: d['ok']=False
    else:
        L=e['L'];lam=L['Hbar']/H70; s=sr_eval(e,Q)
        A=F.shifted(e['af'],e['l0'])[F.IMIG]
        s.pop('ok',None); d.update(ok=True,sr_ok=('cancel_l0' in s),chi2=e['chi'],sigma=sig(e['chi'],7),R2=float(1-np.sum((F.Y-A)**2)/np.sum((F.Y-F.Y.mean())**2)),l0=e['l0'],
          M_kg=e['P']['M']*lam,D=L['D']*lam,vo_c=L['Uo']/F.C,Hpar=L['Hpar']*F.KMSMPC/lam,Hper=L['Hper']*F.KMSMPC/lam,age=e['o']['t']*lam,gk=L['grav_over_kin'],**s)
    rows.append(d); print({k:(round(v,3) if isinstance(v,float) else v) for k,v in d.items() if k not in('z','modes_l0','modes_phifree')},flush=True)
    json.dump(dict(names=names,q=Q,zA=list(zA),zB=list(zB),rows=rows),open('fwd/bracket_sr.json','w'),indent=1,default=float)
print('DONE')
