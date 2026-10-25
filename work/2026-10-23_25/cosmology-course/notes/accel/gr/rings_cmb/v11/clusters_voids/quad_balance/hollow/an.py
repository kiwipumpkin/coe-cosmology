import json,glob,os,sys,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline,PchipInterpolator
from numpy.polynomial import legendre as Lg
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;S=0.04103732750534994;Rb=2.010731531164496;SR='../power_spectrum/slope_refit/'
ns=128;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T;psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));l=np.arange(7)
M={}
for f in glob.glob(SR+'mig_*.json')+glob.glob('mig_*.json'):M.update({float(k):v for k,v in json.load(open(f)).items()})
rr=np.array(sorted(M));fc,fr,fl=[PchipInterpolator(rr,[M[r][k] for r in rr]) for k in ('chi2','R2','l0')]
iso=M[min(M,key=lambda r:abs(r-0.9997017))]['chi2'];best=min(M.values(),key=lambda v:v['chi2'])
def shape(d,x):
  f=str(d['fam'])
  if f=='gsh':return np.exp(-0.5*((x-float(d['XP']))/float(d['W']))**2)*(x<=Rb)
  return x**float(d['p'])*np.clip(1-x/Rb,0,None)**float(d['BB'])
rows=[]
for fn in sorted(glob.glob('runs/*.npz')):
  d=np.load(fn);o=d['o'];o=o[np.argsort(o[:,0])];cs=CubicSpline(o[:,0],np.log(o[:,1]))
  sm=np.exp(-cs(psi));sm=np.asarray(hp.remove_dipole(sm/sm.mean()-1));Dl=l*(l+1)*hp.anafast(sm,lmax=6)*T0**2/(2*np.pi)
  mu=np.linspace(-1,1,2001);f=np.exp(-cs(np.degrees(np.arccos(mu))));f=f/(0.5*np.trapezoid(f,mu))-1
  a=[(2*k+1)/2*np.trapezoid(f*Lg.legval(mu,[0]*k+[1]),mu) for k in range(4)]
  x=np.linspace(1e-4,Rb,4000);s=shape(d,x);xp=float(x[np.argmax(s)]);h=1e-3
  sl=float(((np.log(shape(d,np.array([1+h])))-np.log(shape(d,np.array([1-h]))))/(np.log(1+h)-np.log(1-h)))[0])
  rat=float(d['ratio']);near=min(M,key=lambda r:abs(r-rat));ex=abs(near-rat)<2e-4;inr=rr[0]<=rat<=rr[-1]
  g=lambda k,fun:M[near][k] if ex else (float(fun(rat)) if inr else None)
  rows.append(dict(run=os.path.basename(fn)[:-4],fam=str(d['fam']),a=float(d['p']),b=float(d['BB']),W=float(d['W']),x_peak=round(xp,3),rho0_over_peak=float(shape(d,np.array([1e-4]))[0]/s.max()),
    dlnrho_dlnx_at_us=sl,catwise_locked=bool(abs(sl+S)<2e-3),ratio=rat,chi2=g('chi2',fc),R2=g('R2',fr),l0_best=g('l0',fl),migkas_fit='direct' if ex else ('interp(dense ratio profile)' if inr else 'OUT OF RANGE'),
    a1=a[1],a2=a[2],a3=a[3],D2=float(Dl[2]),D3=float(Dl[3]),D4=float(Dl[4]),z_axis=float(o[0,1]-1),z_anti=float(o[-1,1]-1)))
json.dump(dict(notes='Migkas chi2 from slope_refit dense chi2(ratio) profile (migkas_prof.py model, 71 bins, l0 free, seed 12345); a_l Legendre of dT/T(psi) after monopole; D_l from healpy with dipole removed. Planck D2=226 uK^2.',
  migkas_ref=dict(best_chi2=best['chi2'],best_l0=best['l0'],ratio1_chi2=iso),rows=rows),open('hollow_results.json','w'),indent=1)
for r in rows:print('%-12s xp=%.2f c0=%.2f sl=%+.3f L=%d rat=%.4f chi2=%s l0=%s a2=%+.2e a3=%+.2e D2=%.3g D3=%.3g'%(r['run'],r['x_peak'],r['rho0_over_peak'],r['dlnrho_dlnx_at_us'],r['catwise_locked'],r['ratio'],
  '%.2f'%r['chi2'] if r['chi2'] else 'NA','%.1f'%r['l0_best'] if r['l0_best'] else 'NA',r['a2'],r['a3'],r['D2'],r['D3']))
fig,ax=plt.subplots(1,3,figsize=(18,5.5));x=np.linspace(0,Rb,800)
for fn in sorted(glob.glob('runs/*.npz')):
  d=np.load(fn);s=shape(d,x);ax[0].plot(x,s/s.max(),lw=1,label=os.path.basename(fn)[:-4])
ax[0].axvline(1,color='k',ls=':');ax[0].text(0.01,1.02,'CoE',fontsize=9);ax[0].text(1.01,1.02,'us',fontsize=9);ax[0].set_xlim(0,Rb);ax[0].set_xlabel('x (CoE at left, observer x=1)');ax[0].set_ylabel('ρ/ρ_peak');ax[0].legend(fontsize=5,ncol=2);ax[0].set_title('hollow-core density profiles')
for fam,mk in [('g_','o'),('b_','s'),('bL_','^'),('gL_','v')]:
  R=[r for r in rows if r['run'].startswith(fam)]
  ax[1].plot([r['x_peak'] for r in R],[r['chi2'] if r['chi2'] else np.nan for r in R],mk,label=fam.rstrip('_'))
  ax[2].semilogy([r['x_peak'] for r in R],[r['D2'] for r in R],mk,label=fam.rstrip('_'))
ax[1].axhline(best['chi2'],color='g',ls=':',label='best Migkas %.2f'%best['chi2']);ax[1].axhline(iso,color='gray',ls='--',label='ratio=1 χ² %.2f'%iso);ax[1].set_xlabel('x_peak');ax[1].set_ylabel('Migkas χ² (71 bins)');ax[1].legend(fontsize=8);ax[1].set_ylim(top=min(ax[1].get_ylim()[1],60))
ax[2].axhline(226,color='r',ls='--',label='Planck 226 μK²');ax[2].set_xlabel('x_peak');ax[2].set_ylabel('D2 [μK²]');ax[2].legend(fontsize=8)
plt.tight_layout();plt.savefig('hollow.png',dpi=120)
