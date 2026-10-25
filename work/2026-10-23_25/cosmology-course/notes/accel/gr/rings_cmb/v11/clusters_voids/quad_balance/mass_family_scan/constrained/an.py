import json,glob,os,sys,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline,PchipInterpolator
from numpy.polynomial.legendre import Legendre
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;Q='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/'
ns=128;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T;psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));l=np.arange(7);mu=np.linspace(-1,1,4001)
M={float(k):v for k,v in json.load(open(Q+'power_spectrum/slope_refit/migkas_vs_ratio_dense.json')).items()}
for f in glob.glob('mig_*.json'):M.update({float(k):v for k,v in json.load(open(f)).items()})
rr=np.array(sorted(M));fc=PchipInterpolator(rr,[M[r]['chi2'] for r in rr]);fr=PchipInterpolator(rr,[M[r]['R2'] for r in rr]);fl=PchipInterpolator(rr,[M[r]['l0'] for r in rr])
rows=[]
fs=sorted(glob.glob('runs/*.npz'))+[Q+'mass_family_scan/runs_sl/slc_0_0.npz']+sorted(glob.glob(Q+'mass_family_scan/runs_sl2/*.npz'))
for fn in fs:
  d=np.load(fn);o=d['o'];o=o[np.argsort(o[:,0])];cs=CubicSpline(o[:,0],np.log(o[:,1]))
  sm=np.exp(-cs(psi));sm=np.asarray(hp.remove_dipole(sm/sm.mean()-1));Dl=l*(l+1)*hp.anafast(sm,lmax=6)*T0**2/(2*np.pi)
  y=-cs(np.degrees(np.arccos(mu)));y=y-0.5*np.trapezoid(y,mu)   # dT/T = -ln(1+z) about monopole
  a=[(2*L+1)/2*np.trapezoid(y*Legendre.basis(L)(mu),mu) for L in (1,2,3)]
  rat=float(d['ratio']);near=min(M,key=lambda r:abs(r-rat));ex=abs(near-rat)<2e-4;inr=rr[0]<=rat<=rr[-1]
  rows.append(dict(run=os.path.relpath(fn,Q+'mass_family_scan'),fam=str(d['fam']),p=float(d['p']),W=float(d['W']) if 'W' in d.files else None,SL=float(d['sl']),
    ratio_local=rat,a1=a[0],a2=a[1],a3=a[2],int_stretch_cmb=1.5*abs(a[1]),int_stretch_ratio=float(np.exp(1.5*abs(a[1]))),passes_planck=bool(abs(a[1])<=8.0e-6),D2=float(Dl[2]),D3=float(Dl[3]),
    chi2=float(M[near]['chi2'] if ex else fc(rat)) if inr else None,R2=float(M[near]['R2'] if ex else fr(rat)) if inr else None,l0=float(M[near]['l0'] if ex else fl(rat)) if inr else None,
    migkas='direct' if ex else ('interp dense chi2(ratio) profile' if inr else 'out of range')))
json.dump(rows,open('constrained_results.json','w'),indent=1)
for r in rows:print('%-34s rat=%.5f a2=%+.2e int=%.2e D2=%.3g chi2=%s l0=%s %s'%(r['run'],r['ratio_local'],r['a2'],r['int_stretch_cmb'],r['D2'],r['chi2'] and '%.2f'%r['chi2'],r['l0'] and '%.1f'%r['l0'],r['migkas']))
fig,AX=plt.subplots(1,2,figsize=(11,4.2));ax=AX[0]
for r in rows:
  c='tab:blue' if r['fam']=='slc' else ('tab:red' if r['fam']=='shell' else ('tab:purple' if r['fam']=='combo' else 'tab:green'));ax.loglog(abs(r['ratio_local']-1)+1e-6,r['int_stretch_cmb'],'o',color=c,ms=4)
ax.axhline(1.2e-5,color='k',ls='--',lw=.8,label='Planck bound 1.2e-5');ax.axvspan(0.03,0.13,color='orange',alpha=.15,label='Migkas χ²≲9 band (ratio 1.03–1.13)')
for c,n in [('tab:blue','slc (free slope+curv)'),('tab:red','shell (Gaussian at x=1)'),('tab:green','hat (local gradient at x=1)'),('tab:purple','combo (thin shell + outward-rising slope)')]:ax.plot([],[],'o',color=c,label=n)
ax.set_xlabel('|local stretch ratio − 1| at observer');ax.set_ylabel('|∫(H_r−H_⊥)dt| on CMB rays = 1.5|a2|');ax.legend(fontsize=7);ax.grid(alpha=.3,which='both')
b=AX[1];x=np.linspace(0.02,1.5,1500)
for lab,f,c in [('uniform (ratio 1.000, D2=101)',lambda x:np.ones_like(x),'gray'),('slc s=0.041 p=0.035 (ratio 1.063, D2=1.6e8)',lambda x:x**-0.041037*np.exp(-0.035*(x-1)**2),'tab:blue'),('combo shell p=-0.026 W=0.03 + x^+0.0386 (ratio 1.064, D2=0.44)',lambda x:x**0.0386*np.exp(-0.026*np.exp(-0.5*((x-1)/0.03)**2)),'tab:red')]:
  b.plot(x,f(x)/f(1.0+1e-9),color=c,label=lab)
b.axvline(1,color='k',ls=':');b.text(1.01,b.get_ylim()[0],' us',fontsize=8);b.text(0.03,b.get_ylim()[1],'CoE',fontsize=8,va='top')
b.set_xlabel('distance from CoE x = r/D (CoE at left)');b.set_ylabel('ρ / ρ(us) today');b.legend(fontsize=7);b.grid(alpha=.3)
plt.tight_layout();plt.savefig('constrained.png',dpi=100)
