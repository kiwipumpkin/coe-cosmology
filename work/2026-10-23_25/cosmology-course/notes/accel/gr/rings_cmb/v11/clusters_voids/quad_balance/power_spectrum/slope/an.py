import json,sys,glob,os,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from numpy.polynomial.legendre import Legendre
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;Q='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/'
ns=128;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T
psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));LMAX=6;l=np.arange(LMAX+1);mu=np.linspace(-1,1,4001);LOCK=0.04103732750534994
rows=[]
for fn in sorted(glob.glob(Q+'runs_sl/*.npz')):
  d=np.load(fn);o=d['o'];o=o[np.argsort(o[:,0])];cs=CubicSpline(o[:,0],np.log(o[:,1]))
  sm=np.exp(-cs(psi));sm=np.asarray(hp.remove_dipole(sm/sm.mean()-1));Dl=l*(l+1)*hp.anafast(sm,lmax=LMAX)*T0**2/(2*np.pi)
  A=np.exp(-cs(np.degrees(np.arccos(mu))));dl=A/np.trapezoid(A,mu)*2-1
  a=[(2*L+1)/2*np.trapezoid(dl*Legendre.basis(L)(mu),mu) for L in (1,2,3)]
  sl=float(d['sl']);rows.append(dict(run=os.path.basename(fn)[:-4],near_logslope=-sl,catwise_pct_scaled=3.0*sl/LOCK,p=float(d['p']),ratio=float(d['ratio']),a1=a[0],a2=a[1],a3=a[2],D2=float(Dl[2]),D3=float(Dl[3]),D4=float(Dl[4]),z0=float(o[0,1]),z90=float(np.interp(90,o[:,0],o[:,1])),z180=float(o[-1,1])))
rows.sort(key=lambda r:(r['near_logslope'],r['p']))
for r in rows: print('slope=%+.4f (%.2f%%) p=%.2f ratio=%.5f a2=%+.3e a3=%+.3e D2=%.3g D3=%.3g'%(r['near_logslope'],r['catwise_pct_scaled'],r['p'],r['ratio'],r['a2'],r['a3'],r['D2'],r['D3']))
json.dump(rows,open('slope/slope_release_rows.json','w'),indent=1)
