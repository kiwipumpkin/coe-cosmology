import json,sys,glob,os,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from numpy.polynomial.legendre import Legendre
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;Q='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/'
ns=128;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T
psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));LMAX=6;l=np.arange(LMAX+1);mu=np.linspace(-1,1,4001)
fs=glob.glob(Q+'runs_shape/*.npz')+glob.glob(Q+'runs_amp/*.npz')+[Q+'runs/soft_2_30.npz',Q+'runs/soft_0.3_30.npz'];rows=[]
for fn in fs:
  d=np.load(fn);o=d['o'];o=o[np.argsort(o[:,0])];cs=CubicSpline(o[:,0],np.log(o[:,1]))
  sm=np.exp(-cs(psi));sm=np.asarray(hp.remove_dipole(sm/sm.mean()-1));Dl=l*(l+1)*hp.anafast(sm,lmax=LMAX)*T0**2/(2*np.pi)
  A=np.exp(-cs(np.degrees(np.arccos(mu))));dl=A/np.trapezoid(A,mu)*2-1
  a=[(2*L+1)/2*np.trapezoid(dl*Legendre.basis(L)(mu),mu) for L in (1,2)]
  rows.append(dict(shape=os.path.basename(fn)[:-4],AMP=float(d['AMP']) if 'AMP' in d else 1.0,fam=str(d['fam']),p=float(d['p']),M=float(d['frac']),ratio=float(d['ratio']),D2=float(Dl[2]),D2_planck=float(Dl[2]/226.1),a1=float(a[0]),a2=float(a[1]),D3=float(Dl[3]),D4=float(Dl[4]),z0=float(o[0,1]),z90=float(np.interp(90,o[:,0],o[:,1])),z180=float(o[-1,1]),dt_emit_range=float(o[:,2].max()-o[:,2].min()),obs_logslope=float(-float(d['p']) if str(d['fam'])=='slope' else -0.04103732750534994)))
rows.sort(key=lambda r:r['D2']);[print('%-14s ratio=%.5f D2=%.4e x%.4g a2=%.4e a1=%.4e D3=%.3g D4=%.3g dt=%.2e z=%.2f/%.2f/%.2f amp=%g'%(r['shape'],r['ratio'],r['D2'],r['D2_planck'],r['a2'],r['a1'],r['D3'],r['D4'],r['dt_emit_range'],r['z0'],r['z90'],r['z180'],r['AMP'])) for r in rows]
json.dump(rows,open('density_shape_quad.json','w'),indent=1)
fig,ax=plt.subplots(figsize=(9,6));mk={'soft':'o','eta':'s','exp':'^','pow':'v','curv':'D','core':'P','outer':'*','slope':'X'}
for r in rows: ax.loglog(abs(r['ratio']-1),r['D2'],mk.get(r['fam'],'x'),ms=8,c='C%d'%list(mk).index(r['fam']) if r['fam'] in mk else 'k');ax.annotate(r['shape'],(abs(r['ratio']-1),r['D2']),fontsize=7)
ax.axhline(226.1,c='gold',lw=2,label='Planck D_2=226 μK²');ax.set_xlabel('|ratio − 1|  (transverse/radial stretching anisotropy at observer)');ax.set_ylabel('D_2 [μK²] (dipole removed)');ax.legend();ax.set_title('Smooth-sky quadrupole vs density shape (M=30, D=42)')
plt.savefig('../density_shape_quad.png',dpi=110,bbox_inches='tight')
