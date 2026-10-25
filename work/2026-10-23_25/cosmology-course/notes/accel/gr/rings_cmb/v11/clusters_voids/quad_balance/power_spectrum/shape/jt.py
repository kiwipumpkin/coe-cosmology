import json,sys,glob,os,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from numpy.polynomial.legendre import Legendre
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;Q='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/'
ns=128;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T
psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));LMAX=8;l=np.arange(LMAX+1);mu=np.linspace(-1,1,4001)
fs=sorted(glob.glob(Q+'runs_cub/*.npz'))+glob.glob(Q+'runs_T/*.npz')+[Q+'runs_shape/curv_0.14.npz',Q+'runs/soft_2_30.npz',Q+'runs_rec/f0.npz'];rows=[]
for fn in fs:
  d=np.load(fn);o=d['o'];o=o[np.argsort(o[:,0])];cs=CubicSpline(o[:,0],np.log(o[:,1]))
  sm=np.exp(-cs(psi));sm=np.asarray(hp.remove_dipole(sm/sm.mean()-1));Dl=l*(l+1)*hp.anafast(sm,lmax=LMAX)*T0**2/(2*np.pi)
  A=np.exp(-cs(np.degrees(np.arccos(mu))));dl=A/np.trapezoid(A,mu)*2-1
  a=[(2*L+1)/2*np.trapezoid(dl*Legendre.basis(L)(mu),mu) for L in (1,2,3,4)]
  r=dict(run=os.path.relpath(fn,Q)[:-4],p=float(d['p']),q=float(d['q']) if 'q' in d else 0.0,ratio=float(d['ratio']),a1=a[0],a2=a[1],a3=a[2],a4=a[3],**{f'D{L}':float(Dl[L]) for L in range(2,9)},
    z0=float(o[0,1]),z90=float(np.interp(90,o[:,0],o[:,1])),z180=float(o[-1,1]),dt_emit_range=float(o[:,2].max()-o[:,2].min()))
  if 'obs_slope' in d: r['obs_slope']=float(d['obs_slope']);r['iter_hist']=[float(h) for h in d['hist']]
  rows.append(r)
for r in rows: print('%-22s p=%.3f q=%5.2f ratio=%.5f a2=%+.3e a3=%+.3e D2=%.3g D3=%.3g D4=%.3g D5=%.3g D6=%.3g z=%.2f/%.2f/%.2f dt=%.2e'%(r['run'],r['p'],r['q'],r['ratio'],r['a2'],r['a3'],r['D2'],r['D3'],r['D4'],r['D5'],r['D6'],r['z0'],r['z90'],r['z180'],r['dt_emit_range']))
json.dump(rows,open('joint_shape_rows.json','w'),indent=1)
