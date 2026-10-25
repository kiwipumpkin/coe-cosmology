import json,glob,os,sys,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline,PchipInterpolator
from numpy.polynomial import legendre as Lg
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;S=0.04103732750534994;Rb=2.010731531164496;SR='../power_spectrum/slope_refit/'
PL=np.loadtxt('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb/planckTT.txt')[:5]
ns=128;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T;psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));l=np.arange(7)
M={}
for f in glob.glob(SR+'mig_*.json')+glob.glob('mig_*.json'):M.update({float(k):v for k,v in json.load(open(f)).items()})
rr=np.array(sorted(M));fc,fr,fl=[PchipInterpolator(rr,[M[r][k] for r in rr]) for k in ('chi2','R2','l0')]
iso=M[min(M,key=lambda r:abs(r-0.9997017))]['chi2'];best=min(M.values(),key=lambda v:v['chi2'])
def shape(d,x):
  z=np.load('sedov_g%.4f.npz'%float(d['GAM']));xi,G=z['xi'],z['G'];xs=1/float(d['Q']);e=np.abs(x)/xs;k=3/(float(d['GAM'])-1)
  g=np.where(e<xi[0],G[0]*(np.maximum(e,1e-300)/xi[0])**k,np.exp(np.interp(e,xi,np.log(G))));return np.where(e<=1,g,float(d['AMB']))*(np.abs(x)<=Rb)
rows=[]
for fn in sorted(glob.glob('runs/*.npz')):
  d=np.load(fn);o=d['o'];o=o[np.argsort(o[:,0])];cs=CubicSpline(o[:,0],np.log(o[:,1]))
  sm=np.exp(-cs(psi));sm=np.asarray(hp.remove_dipole(sm/sm.mean()-1));Dl=l*(l+1)*hp.anafast(sm,lmax=6)*T0**2/(2*np.pi)
  mu=np.linspace(-1,1,2001);f=np.exp(-cs(np.degrees(np.arccos(mu))));f=f/(0.5*np.trapezoid(f,mu))-1
  a=[(2*k+1)/2*np.trapezoid(f*Lg.legval(mu,[0]*k+[1]),mu) for k in range(5)]
  x=np.linspace(1e-4,Rb,4000);s=shape(d,x);xp=float(x[np.argmax(s)]);h=1e-3
  sl=float(((np.log(shape(d,np.array([1+h])))-np.log(shape(d,np.array([1-h]))))/(np.log(1+h)-np.log(1-h)))[0])
  rat=float(d['ratio']);near=min(M,key=lambda r:abs(r-rat));ex=abs(near-rat)<2e-4;inr=rr[0]<=rat<=rr[-1]
  g=lambda k,fun:M[near][k] if ex else (float(fun(rat)) if inr else None)
  rows.append(dict(run=os.path.basename(fn)[:-4],gamma=float(d['GAM']),obs_r_over_Rs=float(d['Q']),x_shock=1/float(d['Q']),ambient=float(d['AMB']),frac=float(d['frac']),x_peak=round(xp,3),rho0_over_peak=float(shape(d,np.array([1e-4]))[0]/s.max()),
    dlnrho_dlnx_at_us=sl,catwise_locked=bool(abs(sl+S)<2e-3),ratio=rat,chi2=g('chi2',fc),R2=g('R2',fr),l0_best=g('l0',fl),migkas_fit='direct' if ex else ('interp(dense ratio profile)' if inr else 'OUT OF RANGE'),
    a1=a[1],a2=a[2],a3=a[3],a4=a[4],D2=float(Dl[2]),D3=float(Dl[3]),D4=float(Dl[4]),D5=float(Dl[5]),D6=float(Dl[6]),z_axis=float(o[0,1]-1),z_anti=float(o[-1,1]-1)))
json.dump(dict(notes='Migkas chi2 from slope_refit dense chi2(ratio) profile (migkas_prof.py model, 71 bins, l0 free, seed 12345); a_l Legendre of dT/T(psi) after monopole; D_l from healpy with dipole removed. Planck 2018 TT low-l (planckTT.txt) D2..D6='+str(PL[:,1].tolist())+'. Density: exact Sedov-Taylor G(r/R_s) from sedov_profile.py; CatWISE slope is dlnrho/dlnx at x=1 (lock -0.041).',
  migkas_ref=dict(best_chi2=best['chi2'],best_l0=best['l0'],ratio1_chi2=iso),rows=rows),open('sedov_results.json','w'),indent=1)
for r in rows:print('%-12s xp=%.2f c0=%.2f sl=%+.3f L=%d rat=%.4f chi2=%s l0=%s a2=%+.2e a3=%+.2e D2=%.3g D3=%.3g'%(r['run'],r['x_peak'],r['rho0_over_peak'],r['dlnrho_dlnx_at_us'],r['catwise_locked'],r['ratio'],
  '%.2f'%r['chi2'] if r['chi2'] else 'NA','%.1f'%r['l0_best'] if r['l0_best'] else 'NA',r['a2'],r['a3'],r['D2'],r['D3']))
fig,ax=plt.subplots(1,3,figsize=(18,5.5));x=np.linspace(0,Rb,1500)
for q,c in zip((0.7,0.85,0.95,0.99),('C0','C1','C2','C3')):
  for gm,ls in (('1.6667','-'),('1.3333','--')):
    fn='runs/s_%s_%s_1e-3_f30.npz'%(gm,q);fn=fn if os.path.exists(fn) else sorted(glob.glob('runs/s_%s_%s_*.npz'%(gm,q)))[0] if glob.glob('runs/s_%s_%s_*.npz'%(gm,q)) else None
    if fn is None: continue
    d=np.load(fn);ax[0].plot(x,shape(d,x),c,ls=ls,lw=1.2,label='γ=%s, us at r/R_s=%s'%(gm[:4],q))
ax[0].axvline(1,color='k',ls=':');ax[0].plot(0,0,'k*',ms=12);ax[0].text(0.02,0.2,'CoE',fontsize=9);ax[0].text(1.01,6.6,'observer',fontsize=9)
ax[0].set_xlim(0,Rb);ax[0].set_xlabel('x (CoE at left, observer x=1)');ax[0].set_ylabel('ρ/ρ_ambient,0 (Sedov G; ambient≈0 shown)');ax[0].legend(fontsize=7);ax[0].set_title('Exact Sedov–Taylor initial density')
for gm,mk in (('1.6667','o'),('1.3333','s')):
  for amb,c in ((1e-3,'C0'),(0.25,'C1'),(1.0,'C2')):
    R=[r for r in rows if abs(r['gamma']-float(gm))<1e-3 and r['ambient']==amb]
    if not R: continue
    for f,al in ((30,1),(10,.6),(3,.3),(0.3,.15)):
      RR=sorted([r for r in R if r['frac']==f],key=lambda r:r['obs_r_over_Rs'])
      if RR:
        ax[1].plot([r['obs_r_over_Rs'] for r in RR],[r['chi2'] if r['chi2'] else np.nan for r in RR],mk+'-',color=c,alpha=al,label='γ=%s amb=%g f=%g'%(gm[:4],amb,f))
        ax[2].semilogy([r['obs_r_over_Rs'] for r in RR],[r['D2'] for r in RR],mk+'-',color=c,alpha=al)
ax[1].axhline(best['chi2'],color='g',ls=':',label='best Migkas %.2f'%best['chi2']);ax[1].axhline(iso,color='gray',ls='--',label='ratio=1 χ² %.2f'%iso);ax[1].set_xlabel('observer position r/R_s');ax[1].set_ylabel('Migkas χ² (71 bins)');ax[1].legend(fontsize=5,ncol=2);ax[1].set_ylim(top=min(ax[1].get_ylim()[1],40))
ax[2].axhline(226,color='r',ls='--',label='Planck 226 μK²');ax[2].set_xlabel('observer position r/R_s');ax[2].set_ylabel('D2 [μK²]');ax[2].legend(fontsize=8)
plt.tight_layout();plt.savefig('sedov.png',dpi=110)
