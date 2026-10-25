import sys,os,json,numpy as np,warnings;warnings.filterwarnings('ignore')
import forward as F
from scipy.optimize import differential_evolution,minimize
from scipy.stats import chi2 as C2,norm
H70=70/F.KMSMPC
def sig(chi,k): d=F.CHI0-chi; return float(norm.isf(C2.sf(d,k))) if d>0 else 0.0
law=sys.argv[1]; VMAX=float(os.environ.get('VMAX','0')); SUF='_vc' if VMAX else ''; os.makedirs('fwd',exist_ok=True)
if law=='jerk':
    def c(z):
        D=10**z[0]; rho=z[1]; Hper=H70/(rho/3+2/3); af=F.jerk_af(D,rho*Hper,Hper,z[2],z[3])
        return F.chi2_l0(af)[0] if np.all(np.isfinite(af)) else 1e3
    r=differential_evolution(c,[(-1,3),(-3,6),(-20,20),(-200,200)],seed=1,popsize=12,maxiter=60,polish=False)
    r2=minimize(c,r.x,method='Nelder-Mead',options=dict(maxiter=2000,adaptive=True)); x=r2.x if r2.fun<r.fun else r.x
    D=10**x[0]; Hper=H70/(x[1]/3+2/3); af=F.jerk_af(D,x[1]*Hper,Hper,x[2],x[3]); chi,l0=F.chi2_l0(af)
    A=F.shifted(af,l0)[F.IMIG]; R2=1-np.sum((F.Y-A)**2)/np.sum((F.Y-F.Y.mean())**2)
    json.dump(dict(D_Gpc=D,Hpar=x[1]*Hper*F.KMSMPC,Hper=Hper*F.KMSMPC,k2=x[2],k3=x[3],l0=l0,chi2=chi,R2=R2,sigma_k5=sig(chi,5),af=af.tolist()),open('fwd/jerk.json','w'),indent=1)
    print('jerk',chi,R2,D,x,l0); sys.exit()
SH=law in('eta','pl','soft'); names=['lM','u','p','lf','lL']+(['s'] if SH else [])
B={'lM':(40,60),'u':(-8,0),'p':(-3,4),'lf':(-4,4),'lL':(-3,5),'s':(0.2,6)}
def P_of(z):
    d=dict(zip(names,z)); return F.params_from(law,d['lM'],1-10**d['u'],d['p'],d['lf'],d['lL'],d.get('s'))
def ev(z,keep=False):
    try:
        P=P_of(z); o=F.integrate(P,keep)
        if o is None: return None
        af,L=F.migkas_forward(o)
        if not np.all(np.isfinite(af)): return None
        if VMAX and L['Uo']>VMAX*F.C: return None
        chi,l0=F.chi2_l0(af); return dict(P=P,o=o,af=af,L=L,chi=chi,l0=l0)
    except Exception as e: return None
def cost(z): e=ev(z); return 1e3 if e is None else e['chi']
def fit(fixu=None,x0=None,de=True):
    if fixu is None:
        bb=[B[n] for n in names]
        r=differential_evolution(cost,bb,seed=2,popsize=10,maxiter=30,polish=False,init='sobol')
        x=r.x; f=r.fun
    else: x=np.array(x0,float); x[1]=fixu; f=cost(x)
    def cc(y):
        z=np.array(y); 
        if fixu is not None: z=np.insert(np.delete(z,1),1,fixu) if False else z
        return cost(z)
    idx=[i for i in range(len(names)) if not(fixu is not None and i==1)]
    def c2(y): z=x.copy(); z[idx]=y; return cost(z)
    r2=minimize(c2,x[idx],method='Nelder-Mead',options=dict(maxiter=800,adaptive=True))
    if r2.fun<f: x[idx]=r2.x
    return x
def report(z,tag):
    e=ev(z,True); P,o,L=e['P'],e['o'],e['L']; lam=L['Hbar']/H70   # rescale so Hbar_today = 70
    c=F.cmb(o); G0=CC.G_FULL; S=c['S']; ok=np.isfinite(S)
    Sf=np.where(ok,S,np.nanmean(S)); Sl=np.interp((CC.ELL-e['l0'])%360,np.arange(361.),np.append(Sf,Sf[0]))  # FIXED: orient theta=ell-l0 (was unshifted)
    res=CC.resid(G0,Sl); rms=CC.rms(res); S0=CC.rms(CC.resid(G0,np.ones(360)))
    fo=CC.fourier(res,6); sdip=CC.fourier(S/np.nanmean(S)-1,4)
    Hdir=lambda a: (L['Hpar']*np.cos(np.radians(a))**2+L['Hper']*np.sin(np.radians(a))**2)*F.KMSMPC/lam
    io=o['io']; k2=L['U2']*L['D']/L['U1']; k3=L['U3']*L['D']**2/L['U1']
    afj=F.jerk_af(L['D'],L['Hpar'],L['Hper'],k2,k3,L['Hbar']); chij,l0j=F.chi2_l0(afj)
    A=F.shifted(e['af'],e['l0'])[F.IMIG]; R2=1-np.sum((F.Y-A)**2)/np.sum((F.Y-F.Y.mean())**2)
    Rrec=float(np.nanmean(c['re']))*lam; dtc,dto=F.dt_ls_formula(L['D']*lam,L['Uo'],Rrec)
    te=c['te']*lam; 
    gk_pm=2*F.GU*P['M']/(L['D']*L['Uo']**2)
    out=dict(tag=tag,law=law,free=names+['l0'],nfree=len(names)+1,z=list(map(float,z)),
      chi2=e['chi'],R2=float(R2),sigma=sig(e['chi'],len(names)+1),l0=e['l0'],
      M_kg=P['M']*lam,start_r_Gpc=1.0*lam,V0_over_c=P['V0']/F.C,p=P['p'],L_Gpc=P['L']*lam,shape=P['s'],f_matter=P['f_matter'],G_i=P['G_i'],
      D_today_Gpc=L['D']*lam,vo_over_c=L['Uo']/F.C,beta_o=float(L['Uo']/np.hypot(L['Uo'],F.C)),
      H_par=L['Hpar']*F.KMSMPC/lam,H_perp=L['Hper']*F.KMSMPC/lam,H_mean=70.0,H_0deg=Hdir(0),H_90deg=Hdir(90),
      age_Gyr=o['t']*lam,q_obs=L['q'],j_obs=L['j'],grav_over_kin_today=L['grav_over_kin'],grav_over_kin_pointmass_today=float(gk_pm),
      jerk_k2=k2,jerk_k3=k3,jerk_chi2_with_forward_coeffs=chij,jerk_l0=l0j,af_jerk_minus_full_rms=float(np.std(afj-e['af'])),
      shell_crossing=o['cross'],n_collapsed_into_CoE=o['ncollapsed'],
      cmb_found=int(ok.sum()),cancel_pct=100*(1-rms/S0),cmb_rms=rms,resid_over_1e5=rms/1e-5,
      Tobs_mean=float(np.nanmean(c['Tobs'])),Tobs_min=float(np.nanmin(c['Tobs'])),Tobs_max=float(np.nanmax(c['Tobs'])),
      cmb_modes={m['m']:m['A'] for m in fo['modes']},S_dipole_modes={m['m']:m['A'] for m in sdip['modes']},
      tLS_0=te[0],tLS_90=te[90],tLS_180=te[180],tLS_range=float(np.nanmax(te)-np.nanmin(te)),
      rLS_0=c['re'][0]*lam,rLS_180=c['re'][180]*lam,Rrec_mean=Rrec,
      dtLS_integrated_CoEframe=float(te[0]-te[180]),dtLS_formula_CoEframe=dtc,dtLS_formula_obsframe=dto,
      dtLS_integrated_obsframe=float(np.sqrt(1+(L['Uo']/F.C)**2)*(te[0]-te[180])-(L['Uo']/F.C)*(c['x'][0]-c['x'][180])*lam/F.C),
      af=e['af'].tolist(),afj=afj.tolist(),Tobs=c['Tobs'].tolist(),tLS=te.tolist(),res=res.tolist())
    return out
import importlib.util as iu
sp=iu.spec_from_file_location('cc','/workspace/cosmology-course/notes/cmb-cancel.py');CC=iu.module_from_spec(sp)
src=open('/workspace/cosmology-course/notes/cmb-cancel.py').read().split('if __name__')[0]; exec(compile(src,'cc','exec'),CC.__dict__)
x=fit() if not os.environ.get('FROMCKPT') else np.array(json.load(open(f'fwd/before/{law}{SUF}_ckpt_free.json'))['x']); json.dump(dict(x=list(map(float,x))),open(f'fwd/{law}{SUF}_ckpt_free.json','w'))
res=dict(free=report(x,'free')); json.dump(res,open(f'fwd/{law}{SUF}.json','w'),indent=1,default=float); print('free',res['free']['chi2'],res['free']['grav_over_kin_today'],flush=True)
res['profile']=[]
for u in ((-1,-2,-3,-4) if VMAX else ()):
    xu=fit(fixu=u,x0=x); r=report(xu,f'u{u}'); res['profile'].append({k:v for k,v in r.items() if k not in('af','afj','Tobs','tLS','res')})
    json.dump(res,open(f'fwd/{law}{SUF}.json','w'),indent=1,default=float); print('u',u,r['chi2'],r['grav_over_kin_today'],flush=True)
print('DONE',law)
