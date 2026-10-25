import sys,json,numpy as np,healpy as hp;sys.path.insert(0,'../rotfit')
import fit
from scipy.interpolate import CubicSpline
T0=2.7255e6;PD=3362/T0;nside=64;V=np.array(hp.pix2vec(nside,np.arange(hp.nside2npix(nside)))).T
o=np.load('../shrink_rec/b_1100_k0.npz')['o'];S=np.load('s.npy')
aC=fit.vec(0,15);psi=np.degrees(np.arccos(np.clip(V@aC,-1,1)))
lnz=CubicSpline(o[:,0],np.log(o[:,1]))(psi);s_pix=np.interp(psi,S[:,0],S[:,1])
smooth=np.exp(-lnz);smooth=smooth/smooth.mean()-1          # dT/T of smooth model
def dipvec(m):return 3*(V*m[:,None]).mean(0)
def hi_rms(m):a=hp.map2alm(m,lmax=3*nside-1,iter=1);fl=np.ones(3*nside);fl[:2]=0;return float(np.std(hp.alm2map(hp.almxfl(a,fl),nside)))
Dm=dipvec(smooth);out=dict(smooth_xPlanck=float(np.linalg.norm(Dm)/PD),s_mean=float(s_pix.mean()),s_min=float(s_pix.min()),s_max=float(s_pix.max()))
def field(rng,fwhm_deg,drms):
    w=rng.normal(size=len(V));m=hp.smoothing(w,fwhm=np.radians(fwhm_deg),iter=0);m-=m.mean();return m/np.std(m)*drms
def total(delta):   # threshold rule: emission when rho(1+delta)=rho_c  -> ln(1+z) shifts by -s*delta
    T=np.exp(-(lnz-s_pix*delta));return T/T.mean()-1
rng=np.random.default_rng(7);res={};N=int(sys.argv[1]) if len(sys.argv)>1 else 300
for fw in (1.0,5.0):
  for drms in (1e-5,1e-4,1e-3,1e-2):
    d=[];h=[];hot=[]
    for i in range(N):
        m=total(field(rng,fw,drms));v=dipvec(m);d.append(np.linalg.norm(v)/PD);hot.append(fit.ang(v,-aC))
        if i<20:h.append(hi_rms(m))
    d=np.array(d);res[f'fwhm{fw}_rms{drms}']=dict(dipole_median=float(np.median(d)),p05=float(np.percentile(d,5)),p95=float(np.percentile(d,95)),min=float(d.min()),hot_offset_from_CoE_median=float(np.median(hot)),small_scale_rms_dTT=float(np.mean(h)),small_scale_rms_uK=float(np.mean(h)*T0),frac_below_1x=float(np.mean(d<=1)))
    np.save(f'd_fwhm{fw}_rms{drms}.npy',d);print(fw,drms,res[f'fwhm{fw}_rms{drms}'],flush=True)
out['random']=res;out['smooth_hi_rms_dTT']=hi_rms(smooth)
# arranged: add a net density dipole (hemisphere imbalance of clusters vs voids) that cancels the excess
json.dump(dict(out,random=res),open("cv_results.json","w"),indent=1)
need=(1-1/out['smooth_xPlanck'])  # fraction of model dipole to remove
# density pattern A*cos(angle from CoE): its dT/T dipole is s*A (cos has dipole amplitude A). Solve numerically
cosC=V@(-aC);A=lambda a:np.linalg.norm(dipvec(total(a*cosC)))/PD
from scipy.optimize import brentq
a1=brentq(lambda a:A(a)-1,-0.05,0);a0=brentq(lambda a:A(a)-0.0,-0.05,0) if A(-0.05)<0 else None
from scipy.optimize import minimize_scalar
mn=minimize_scalar(A,bounds=(-0.05,0),method='bounded')
out['arranged']=dict(density_dipole_for_1x=float(a1),hemisphere_contrast_for_1x=float(2*abs(a1)),density_dipole_for_min=float(mn.x),min_xPlanck=float(mn.fun))
for fw,drms in((1.0,1e-3),(5.0,1e-3),(1.0,1e-2)):
    dd=[];
    for i in range(100):
        m=total(field(rng,fw,drms)+a1*cosC);dd.append(np.linalg.norm(dipvec(m))/PD)
    out['arranged'][f'with_random_fwhm{fw}_rms{drms}']=dict(median=float(np.median(dd)),p05=float(np.percentile(dd,5)),p95=float(np.percentile(dd,95)),small_scale_rms_uK=float(hi_rms(total(field(rng,fw,drms)+a1*cosC))*T0))
print(json.dumps(out,indent=1));json.dump(out,open('cv_results.json','w'),indent=1)
# plots
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
r1=np.random.default_rng(11);m=total(field(r1,1.0,1e-3));v=dipvec(m)
fig=plt.figure(figsize=(16,10));hp.mollview(m*T0,fig=fig.number,cmap='RdBu_r',unit='ΔT (μK)',title=f'one random realisation: smooth model + clusters/voids (density rms 0.1%, 1° scale)\ndipole {np.linalg.norm(v)/PD:.2f}× Planck, hot end {fit.ang(v,-aC):.1f}° from CoE; small-scale rms {hi_rms(m)*T0:.0f} μK')
hp.graticule(dpar=30,dmer=60,alpha=.3)
for (l,b),c,mk,lab in (((264.02,48.25),'r','^','Planck dipole hot'),((84.02,-48.25),'b','v','Planck dipole cold'),((180,-15),'k','*','CoE'),((0,15),'k','X','anti-CoE')):hp.projscatter(l,b,lonlat=True,c=c,marker=mk,s=300,edgecolors='w',label=lab)
plt.gca().title.set_fontsize(15);fig.legend(loc='upper center',ncol=4,fontsize=14,frameon=False,bbox_to_anchor=(.5,1.04));fig.savefig('map_realisation.png',dpi=80,bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(16,9))
for k,c,lab in (('fwhm1.0_rms0.001','C0','random, 1° scale, density rms 0.1%'),('fwhm1.0_rms0.01','C1','random, 1° scale, density rms 1%'),('fwhm5.0_rms0.01','C3','random, 5° scale, density rms 1%')):
    d=np.load(f'd_{k}.npy');ax.hist(d,bins=60,histtype='step',lw=3,color=c,label=f"{lab}: median {np.median(d):.3f}×, 5–95% {np.percentile(d,5):.3f}–{np.percentile(d,95):.3f}×")
ax.axvline(1,color='k',ls='--',lw=2,label='1× Planck');ax.axvline(out['smooth_xPlanck'],color='m',lw=2,label=f"smooth model {out['smooth_xPlanck']:.2f}×")
ax.set_xlabel('combined model dipole (× Planck), no 370 km/s motion',fontsize=16);ax.set_ylabel(f'number of realisations (of {N})',fontsize=16);ax.tick_params(labelsize=14);ax.grid(alpha=.3)
ax.legend(loc='lower center',bbox_to_anchor=(.5,1.0),fontsize=13,frameon=False);fig.savefig('dipole_hist.png',dpi=90,bbox_inches='tight')
