import json,glob,os,numpy as np,healpy as hp,sys,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline,PchipInterpolator
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;Q='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/';LOCK=0.04103732750534994
ns=128;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T;psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));l=np.arange(7)
M={}
for f in glob.glob('mig_*.json'):M.update({float(k):v for k,v in json.load(open(f)).items()})
rr=np.array(sorted(M));c2=np.array([M[r]['chi2'] for r in rr]);r2=np.array([M[r]['R2'] for r in rr]);l0=np.array([M[r]['l0'] for r in rr])
iso=M[min(M,key=lambda r:abs(r-0.9997017))]['chi2']  # ratio=1 (no stretch anisotropy) reference
json.dump({str(r):M[r] for r in rr},open('migkas_vs_ratio_dense.json','w'),indent=1)
fc=PchipInterpolator(rr,c2);fr=PchipInterpolator(rr,r2);fl=PchipInterpolator(rr,l0)
rows=[]
for fn in sorted(glob.glob(Q+'runs_sl/*.npz')+glob.glob(Q+'runs_sl2/*.npz')):
  d=np.load(fn);o=d['o'];o=o[np.argsort(o[:,0])];cs=CubicSpline(o[:,0],np.log(o[:,1]))
  sm=np.exp(-cs(psi));sm=np.asarray(hp.remove_dipole(sm/sm.mean()-1));Dl=l*(l+1)*hp.anafast(sm,lmax=6)*T0**2/(2*np.pi)
  rat=float(d['ratio']);sl=float(d['sl']);near=min(M,key=lambda r:abs(r-rat));ex=abs(near-rat)<2e-4
  ch=M[near]['chi2'] if ex else float(fc(rat))
  rows.append(dict(s=-sl,catwise_pct=3*sl/LOCK,p=float(d['p']),ratio=rat,l0_best=M[near]['l0'] if ex else float(fl(rat)),chi2=ch,R2=M[near]['R2'] if ex else float(fr(rat)),
    sigma_vs_ratio1=float(np.sqrt(max(iso-ch,0))),dchi2_vs_current=ch-M[1.051]['chi2'],D2=float(Dl[2]),D3=float(Dl[3]),D4=float(Dl[4]),migkas_fit='direct' if ex else 'interp(dense ratio profile)'))
rows.sort(key=lambda r:(r['s'],r['p']));json.dump(rows,open('slope_refit_rows.json','w'),indent=1)
i=np.argmin(c2);print('best ratio',rr[i],c2[i],r2[i],l0[i],'iso',iso,'cur',M[1.051])
for r in rows:print('s=%+.3f(%4.1f%%) p=%.3f rat=%.4f chi2=%.2f R2=%.3f l0=%.1f sig=%.2f D2=%.3g D3=%.3g'%(r['s'],r['catwise_pct'],r['p'],r['ratio'],r['chi2'],r['R2'],r['l0_best'],r['sigma_vs_ratio1'],r['D2'],r['D3']))
fig,ax=plt.subplots(1,2,figsize=(14,5.5))
a=ax[0];a.plot(rr,c2,'k.-',label='Migkas best χ² (ℓ0 free)');a.axhline(M[1.051]['chi2'],color='g',ls=':');a.axvline(1.051,color='g',ls=':',label='current soft p=2 (ratio 1.051, χ²=%.2f)'%M[1.051]['chi2'])
a.set_xlabel('transverse/radial stretching ratio at observer');a.set_ylabel('Migkas χ² (71 bins)');b=a.twinx()
R=np.array([[r['ratio'],r['D2']] for r in rows]);b.semilogy(R[:,0],R[:,1],'r^',ms=5,label='D2 of free-slope runs');b.axhline(226,color='r',ls='--',lw=.8);b.axhline(1.08e8,color='r',ls=':',lw=.8);b.set_ylabel('D2 [μK²]',color='r')
a.legend(loc='upper left',fontsize=8);b.legend(loc='upper right',fontsize=8);a.set_title('Migkas χ² vs ratio; D2 (red) on twin axis')
a=ax[1];b=a.twinx()
for p,cc in zip(sorted(set(r['p'] for r in rows)),plt.cm.viridis(np.linspace(0,1,len(set(r['p'] for r in rows))))):
  S=sorted([r for r in rows if r['p']==p],key=lambda r:r['s'])
  a.plot([r['s'] for r in S],[r['chi2'] for r in S],'o-',color=cc,label='p=%g'%p);b.semilogy([r['s'] for r in S],[r['D2'] for r in S],'x--',color=cc,alpha=.6)
a.axhline(M[1.051]['chi2'],color='g',ls=':');a.axvline(-LOCK,color='gray',ls=':');a.text(-LOCK,a.get_ylim()[1]*0.95,' CatWISE 3% lock',fontsize=8)
a.set_xlabel('near-observer log slope s (ln ρ = -s... ; shown as d ln ρ/d ln x)');a.set_ylabel('Migkas χ² (solid)');b.set_ylabel('D2 μK² (dashed x)');b.axhline(226,color='r',ls='--',lw=.8)
a.legend(fontsize=7,loc='lower left');a.set_ylim(top=min(a.get_ylim()[1],60))
plt.tight_layout();plt.savefig('migkas_refit_vs_slope.png',dpi=130)
